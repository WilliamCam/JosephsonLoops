using Symbolics
using SymbolicUtils
using NonlinearSolve
struct HarmonicSystem
    complete_sys::ModelingToolkit.AbstractSystem
    ω_var::Num
    t_var::Num          # independent variable from the original ODE system
    N::Int
    variable_map::Dict{Tuple{String, Int, Symbol}, Num}
    observed_equations::Dict{Any, Any}
end

struct HarmonicProblem
    harmonic_system::HarmonicSystem
    params::Dict
end

struct HarmonicSweepProblem
    prob::HarmonicProblem
    sweep_var::Num
    sweep_vals::AbstractVector
    
    function HarmonicSweepProblem(prob::HarmonicProblem, sweep_var::Num, sweep_vals::AbstractVector)
        new_params = copy(prob.params)
        if haskey(new_params, sweep_var)
            delete!(new_params, sweep_var)
        end
        sweep_var_uw = Symbolics.unwrap(sweep_var)
        if haskey(new_params, sweep_var_uw)
            delete!(new_params, sweep_var_uw)
        end
        new_prob = HarmonicProblem(prob.harmonic_system, new_params)
        new(new_prob, sweep_var, sweep_vals)
    end
end

HarmonicSweepProblem(prob::HarmonicProblem, sweep_var::Nothing, sweep_vals) = prob
HarmonicSweepProblem(prob::HarmonicProblem) = prob

struct HarmonicResult
    sweep_var::Num
    sweep_vals::AbstractVector
    results::Dict{Num, Vector{Float64}}
    observed_results::Dict{Any, Vector{Float64}}
end

function var_is_in(vars::Vector, target_var::SymbolicUtils.BasicSymbolic{Real})
    ret = false
    for var in vars
        if isequal(var, target_var)
            ret = true
            break
        end
    end
return ret
end

function var_is_in(vars::Vector, target_var::Num)
    ret = false
    for var in vars
        if isequal(var, target_var)
            ret = true
            break
        end
    end
return ret
end

function var_index(vars::Vector, target_var::SymbolicUtils.BasicSymbolic{Real})
    return findfirst(x->isequal(x, target_var),vars)
end

function get_HB_scattering_matrix(model::System,i::Char,j::Char)
    port_i_sym = Symbol('P'*i)
    port_j_sym = Symbol('P'*j)
    @assert (hasproperty(model,port_i_sym)) "Error: Port $i not defined in variable map"
    @assert (hasproperty(model,port_j_sym)) "Error: Port $j not defined in variable map"
    N_ports = maximum([parse(Int,port_index) for port_index in [i,j]])
    @assert (N_ports < 3) "Maximum of 2 port networks supported"
    a = zeros(Num,1,N_ports)
    b = zeros(Num,1,N_ports)
    #TODO: assert port reference impedances are equal
    for k in N_ports
        port_k_sym = Symbol('P'*string(k))
        port_k = getproperty(model, port_k_sym)
        a[k] = 0.5/sqrt(port_k.Rsrc.R)*(port_k.dθ*Φ₀/2π+port_k.Rsrc.R*port_k.i)
        b[k] = 0.5/sqrt(port_k.Rsrc.R)*(port_k.dθ*Φ₀/2π-port_k.Rsrc.R*port_k.i)
    end
    return b / a
end

function only_derivatives(expr, var, tvar)
    # The maximum order of derivative to check against.
    # We must substitute all known derivative terms of `var`.
    sub_rules = Dict(
        Differential(tvar)(var) => 0,
        Differential(tvar)(Differential(tvar)(var)) => 0,
    )
    expr_sub = substitute(expr, sub_rules)
    remaining_vars = get_variables(expr_sub)
    return !var_is_in(remaining_vars, var)
end

function get_derivatives(X, t)
    D = Differential(t)
    dXdt = Symbolics.expand_derivatives(D(X))
    d2Xdt2 = Symbolics.expand_derivatives(D(dXdt))
    return dXdt, d2Xdt2
end

function get_full_equations(model::ModelingToolkit.System, tvar::Num)
    eqs = full_equations(model)
    states = unknowns(model)

    diff2vars = Vector{Num}()
    diffvars = Vector{Num}()
    remove_idxs = Int[]
    for (i, eq) in enumerate(eqs)
        vars = get_variables(eq.rhs)
        if length(vars) == 1 && var_is_in(states, vars[1])
            push!(diff2vars, vars[1])
            push!(diffvars, get_variables(eq.lhs)[1])
            push!(remove_idxs, i)
        end
    end

    for i in reverse(remove_idxs)
        deleteat!(eqs, i)
    end

    for (i,var) in enumerate(diffvars)
        eqs = substitute(eqs, Dict(diff2vars[i]=>Differential(tvar)(diffvars[i])))
    end
    remove_idxs = Int[]
    for (i,var) in enumerate(states)
        if var_is_in(diff2vars, var)
            push!(remove_idxs, i)   
        end
    end
    for i in reverse(remove_idxs)
        deleteat!(states, i)
    end
    return eqs, states
end

function is_term(set, target_term)
    if typeof(set) == Equation
        vars = get_variables(set)
    elseif typeof(set) == SymbolicUtils.BasicSymbolic{Real}
        vars = get_variables(set)
    elseif typeof(set) == Num
        vars = get_variables(set)
    else
        vars = set
    end
    ret = false
    for term in vars
        if isequal(term, target_term)
            ret = true
            break
        else
            ret = false
        end
    end
    return ret
end


function build_jacobians(rotated_system, vars, dvars)
    #TODO check ordering
    _jac = Symbolics.jacobian(rotated_system, vars)
    jac_0 = Num.(simplify(substitute(_jac, Dict(dvars .=> 0))))
    jac_1 = Symbolics.jacobian(rotated_system, dvars)
    return jac_0, jac_1
end
function rotate_to_harmonic_frame(N, Nt, harmonic_system)
    Γ = Matrix{Num}(undef, 2N + 1, Nt)
    # Currently only works for one harmonic ansatz. e.g. M=1
    for j in 1:Nt
        # 1. Place DC at the first index [1, j]
        Γ[1, j] = Num(1//Nt)
        # 2. Place Cosines at indices [2 to N+1]
        for n in 1:N
            phase = n * (j - 1) * (2π / Nt)
            Γ[n + 1, j] = Num((2//Nt) * cos(phase))
        end
        # 3. Place Sines at indices [N+2 to 2N+1]
        for n in 1:N
            phase = n * (j - 1) * (2π / Nt)
            Γ[N + 1 + n, j] = Num((2//Nt) * sin(phase))
        end
    end
    
    # Perform the rotation/projection
    rotated_system = Γ * [equation.lhs for equation in harmonic_system]
    
    return simplify.(rotated_system)
end

#  Shared Helpers
"""Strip `(t)` suffix and whitespace from a symbolic variable name."""
clean_name(x) = replace(string(x), "(t)" => "", " " => "")

"""Find a matching key in a variable_map by exact or fuzzy name match."""
function find_varmap_key(variable_map, vname, order, component)
    haskey(variable_map, (vname, order, component)) && return (vname, order, component)
    for k in keys(variable_map)
        k[2] == order && k[3] == component && occursin(vname, k[1]) && return k
    end
    return nothing
end

"""Unified symbolic lookup across one or more Dicts. Tries hash, isequal, then string match."""
function _find_sym(sym::Num, dicts...)
    sym_uw = Symbolics.unwrap(sym)
    sym_str = clean_name(sym)
    for d in dicts
        haskey(d, sym) && return d[sym]
        for (k, v) in d
            (isequal(k, sym) || isequal(Symbolics.unwrap(k), sym_uw)) && return v
        end
        for (k, v) in d
            clean_name(k) == sym_str && return v
        end
    end
    return nothing
end


"""
    HarmonicSystem(sys, ωvar; tearing=true, N=1) -> HarmonicSystem

Construct a `HarmonicSystem` by applying the harmonic balance method to a
time-domain ODE system.

Each state variable `x(t)` is expanded into a truncated Fourier series
    x(t) = A₀ + Σₙ₌₁ᴺ [Aₙ cos(nωt) + Bₙ sin(nωt)]
and the ODE residual is projected onto each harmonic basis function, converting
the differential equations into a set of nonlinear algebraic equations in the
Fourier coefficients {A₀, Aₙ, Bₙ}.

# Arguments
- `sys`: A `ModelingToolkit.ODESystem` (or compatible system) containing the
  time-domain differential equations.
- `ωvar::Num`: Symbolic variable representing the fundamental angular frequency ω.
  Must already be declared as a parameter in `sys`.

# Keywords
- `N::Int=1`: Harmonic truncation order. `N=1` retains only the fundamental; larger
  values improve accuracy at the cost of a larger nonlinear system.
- `tearing::Bool=true`: Whether to apply structural simplification (`mtkcompile`) to
  the resulting `NonlinearSystem`. Tearing reduces the number of active unknowns and
  typically speeds up the solve.

# Returns
A `HarmonicSystem` containing the compiled nonlinear system, frequency metadata,
the variable map, and the original observed equations.

# Details
- The independent variable (time) is inferred automatically via
  `ModelingToolkit.get_iv(sys)`.
- If the expanded system is over-determined (more equations than unknowns — which
  can occur when symmetry constraints force certain coefficients to zero), trailing
  equations are dropped with a warning. The ordering is determined by
  `harmonic_equation`, so inspect results carefully in this case.
- Observed equations from `sys` are stored verbatim in `HarmonicSystem.observed_equations`
  so that non-state quantities (e.g. branch currents, voltages across elements) can be
  reconstructed from the HB solution without re-solving the system.
"""
function HarmonicSystem(sys, ωvar::Num; tearing::Bool=true, N::Int=1)
    tvar = Num(ModelingToolkit.get_iv(sys))
    eqs, states = get_full_equations(sys, tvar)

    nonlinear_sys, _, variable_map = harmonic_equation(eqs, states, tvar, ωvar, N)
    
    sys_eqs = equations(nonlinear_sys)
    sys_vars = unknowns(nonlinear_sys)
    
    if length(sys_eqs) > length(sys_vars)
        n_drop = length(sys_eqs) - length(sys_vars)
        @warn "System is overdetermined: $(length(sys_eqs)) equations for $(length(sys_vars)) variables. " *
              "Dropping the last equation(s). Caution: This behavior depends on variable order."
        sys_eqs = sys_eqs[1:end-n_drop]
    end
        
    @named nonlinear_sys = NonlinearSystem(sys_eqs, sys_vars, parameters(sys))
    complete_sys = tearing ? mtkcompile(nonlinear_sys; fully_determined=false) : nonlinear_sys
    
    # Store original observed equations to enable reconstruction of non-state variables (eg. C1.i)
    observed_eqs = Dict{Any, Any}(eq.lhs => eq.rhs for eq in observed(sys))
  
    return HarmonicSystem(complete_sys, ωvar, tvar, N, variable_map, observed_eqs)
end

"""
    HarmonicProblem(sys, ωvar, params; tearing=true, N=1) -> HarmonicProblem

Convenience constructor that builds a `HarmonicSystem` from a time-domain ODE system
and immediately pairs it with parameter values.

# Arguments
- `sys::ModelingToolkit.AbstractSystem`: The time-domain ODE system to expand.
- `ωvar::Num`: Symbolic variable for the fundamental angular frequency ω.
- `params`: Parameter map from symbolic variables to numeric values (e.g. `Dict(ω => 1.0, R => 50.0)`).

# Keywords
- `N::Int=1`: Harmonic truncation order passed to `HarmonicSystem`.
- `tearing::Bool=true`: Whether to structurally simplify the nonlinear system via `mtkcompile`.

# Returns
A `HarmonicProblem` ready to be passed to `solve` or wrapped in a `HarmonicSweepProblem`.
"""
function HarmonicProblem(sys::ModelingToolkit.AbstractSystem, ωvar::Num, params; tearing::Bool=true, N::Int=1)
    harmonic_sys = HarmonicSystem(sys, ωvar; tearing=tearing, N=N)
    return HarmonicProblem(harmonic_sys, params)
end

"""
    solve(prob::HarmonicProblem; kwargs...) -> NonlinearSolution

Solve the harmonic balance system at a single parameter point.

Constructs a `NonlinearProblem` from the compiled HB system using a zero initial guess
for all unknowns, then calls `ModelingToolkit.solve`. All keyword arguments are
forwarded to the underlying solver (e.g. `alg`, `abstol`, `reltol`).

# Returns
A `NonlinearSolution` as returned by `ModelingToolkit.solve`. Harmonic coefficients
can be read directly from `sol[var]` using the symbolic unknowns of
`prob.harmonic_system.complete_sys`.
"""
function solve(prob::HarmonicProblem; kwargs...)
    hsys = prob.harmonic_system
    sys = hsys.complete_sys
    u0_guess = fill(0.0, length(unknowns(sys)))
    combined_args = merge(Dict(unknowns(sys) .=> u0_guess), prob.params)
    nl_prob = NonlinearProblem(sys, combined_args)
    return ModelingToolkit.solve(nl_prob; kwargs...)
end

"""
    solve(sweepprob::HarmonicSweepProblem; kwargs...) -> HarmonicResult

Solve the harmonic balance system across a parameter sweep using continuation.

Iterates over each value in `sweepprob.sweep_vals`, updating `sweep_var` at each step
and using the previous solution as the initial guess for the next (`u0` continuation).
All keyword arguments are forwarded to `ModelingToolkit.solve` at each point.

# Returns
A `HarmonicResult` containing:
- `sweep_var`: The swept symbolic parameter.
- `sweep_vals`: The vector of sweep values as collected.
- `results`: `Dict{Num, Vector{Float64}}` mapping each HB unknown to its values across
  the sweep.
- `observed_results`: `Dict` for pre-computed observed quantities (populated downstream
  by `get_harmonic` / `fetch_harmonic_coeff` on demand).

# Details
- Continuation means the solver is warm-started from the previous sweep point's solution.
  This improves convergence along smooth solution branches but may cause the solver to
  track an unstable branch through a bifurcation rather than jump to the stable one.
- Progress is printed to stdout for each sweep point.
"""
function solve(sweepprob::HarmonicSweepProblem; kwargs...)
    prob = sweepprob.prob
    hsys = prob.harmonic_system
    sweep_var = sweepprob.sweep_var
    sweep_vals = sweepprob.sweep_vals
    sys = hsys.complete_sys

    # Setup Parameters
    current_params = copy(prob.params)
    system_unknowns = unknowns(sys)
    
    # Initial guess
    u0_guess = fill(0.0, length(system_unknowns))

    current_params[sweep_var] = first(sweep_vals)

    combined_args = merge(Dict(system_unknowns .=> u0_guess), current_params)
    nl_prob = NonlinearProblem(sys, combined_args)
    
    # Results containers for unknowns
    results = Dict{Num, Vector{Float64}}()
    for v in system_unknowns
        results[v] = Float64[]
        sizehint!(results[v], length(sweep_vals))
    end
    
    # Results containers for observed symbols and extra varmap symbols
    # Key by Num symbol (eq.lhs) for exact isequal matching
    observed_results = Dict{Any, Vector{Float64}}()

    # Also register extra symbols pre-identified in HarmonicSystem
   
    println("Sweeping $(sweep_var) over $(length(sweep_vals)) points...")
    # Initialize last_u with the numeric initial guess
    last_u = nl_prob.u0
  
    for  val in sweep_vals     
        # Continuation: use previous solution (last_u)
        nl_prob = remake(nl_prob; u0 = last_u, p = [sweep_var => val])
        sol = ModelingToolkit.solve(nl_prob, kwargs...)
        last_u = sol.u
        for (i, v) in enumerate(system_unknowns)
            push!(results[v], sol.u[i])
        end
    end
    return HarmonicResult(sweep_var, collect(sweep_vals), results, observed_results)
end


"""
    fetch_harmonic_coeff(sym, h_prob, sweep_res) -> Vector{Float64}

Retrieve the sweep-resolved numeric values for a single HB coefficient symbol `sym`.

# Arguments
- `sym::Num`: A symbolic variable representing one harmonic coefficient (e.g. the cosine
  or sine amplitude of a particular state at a particular harmonic order), as stored in
  `h_prob.harmonic_system.variable_map`.
- `h_prob::HarmonicProblem`: The harmonic balance problem, used to access the system
  definition and parameter values when a fallback evaluation is required.
- `sweep_res::HarmonicResult`: The solver output containing `results` (primary HB
  unknowns) and `observed_results` (pre-computed observed quantities).

# Returns
A `Vector{Float64}` of length `length(sweep_res.sweep_vals)`, giving the value of `sym`
at each sweep point.

# Lookup strategy
1. **Direct lookup** — searches `sweep_res.results` then `sweep_res.observed_results`
   using hash equality, `isequal`, and finally string-name matching via `_find_sym`.
2. **Symbolic fallback** — if the symbol is not found in the cached results, locates the
   corresponding observed equation in `ModelingToolkit.observed(sys)`, recursively
   substitutes any observed-in-observed dependencies (up to 20 levels), and evaluates
   the resulting expression numerically at each sweep point by substituting all HB
   unknowns and sweep parameters.

# Errors
Throws if `sym` cannot be located in either the results or the observed equations, or if
the fallback expression cannot be fully reduced to a numeric value.
"""
function fetch_harmonic_coeff(sym::Num, h_prob::HarmonicProblem, sweep_res::HarmonicResult)
    result = _find_sym(sym, sweep_res.results, sweep_res.observed_results)
    result !== nothing && return result

    # Last resort: manual symbolic substitution from observed equations
    sys = h_prob.harmonic_system.complete_sys
    sym_str = clean_name(sym)
    target_eq = nothing
    for eq in observed(sys)
        if isequal(eq.lhs, sym) || clean_name(eq.lhs) == sym_str
            target_eq = eq
            break
        end
    end
    target_eq === nothing && error("Harmonic coefficient $sym not found in results or observed equations")

    expr = Symbolics.unwrap(target_eq.rhs)
    obs_dict = Dict{Any, Any}(Symbolics.unwrap(eq.lhs) => Symbolics.unwrap(eq.rhs) for eq in observed(sys))
    states_set = Set(Symbolics.unwrap(v) for v in unknowns(sys))
    params_set = Set(Symbolics.unwrap(p) for p in parameters(sys))
    tvar_uw = Symbolics.unwrap(h_prob.harmonic_system.t_var)

    for _ in 1:20
        vars_in_expr = Symbolics.get_variables(expr)
        unknown_vars = filter(v -> !in(v, states_set) && !in(v, params_set) && !isequal(v, tvar_uw), vars_in_expr)
        isempty(unknown_vars) && break
        
        subs_dict_expr = Dict()
        for u in unknown_vars
            # Handle derivatives like Differential(t)(J1₊out₊Φ(t))
            if SymbolicUtils.iscall(u) && SymbolicUtils.operation(u) isa Differential
                base_var = SymbolicUtils.arguments(u)[1]
                if haskey(obs_dict, base_var)
                    subs_dict_expr[u] = Symbolics.expand_derivatives(Differential(tvar_uw)(obs_dict[base_var]))
                end
            elseif haskey(obs_dict, u)
                subs_dict_expr[u] = obs_dict[u]
            end
        end
        
        isempty(subs_dict_expr) && break
        expr = Symbolics.substitute(expr, subs_dict_expr)
        expr = Symbolics.expand_derivatives(expr)
        expr = Symbolics.simplify(expr)
    end

    n_points = length(sweep_res.sweep_vals)
    out_vals = zeros(Float64, n_points)
    
    param_syms = parameters(sys)
    p_vals_base = [h_prob.params[p] for p in param_syms]
    sweep_idx = findfirst(isequal(sweep_res.sweep_var), param_syms)
    states_vec = unknowns(sys)
    
    for i in 1:n_points
        p_dict = Dict{Any, Any}(Symbolics.unwrap(param_syms[j]) => p_vals_base[j] for j in eachindex(param_syms))
        sweep_idx !== nothing && (p_dict[Symbolics.unwrap(sweep_res.sweep_var)] = sweep_res.sweep_vals[i])
        u_dict = Dict{Any, Any}(Symbolics.unwrap(v) => sweep_res.results[v][i] for v in states_vec)

        subs_dict = merge(p_dict, u_dict)
        val = Symbolics.substitute(expr, subs_dict)
        if val isa SymbolicUtils.BasicSymbolic || val isa Num
            val = Symbolics.substitute(val, subs_dict)
        end
        out_vals[i] = Float64(Symbolics.unwrap(val))
    end
    return out_vals
end


"""
    get_harmonic(h_prob, sweep_res, var_name, order; output_type=:complex) -> Vector

Return the harmonic phasor (or a derived quantity) for a named variable at a given
harmonic order across the full parameter sweep.

# Arguments
- `h_prob::HarmonicProblem`: The harmonic balance problem.
- `sweep_res::HarmonicResult`: Sweep results returned by the HB solver.
- `var_name::String`: Name of the variable of interest (namespace prefixes are stripped
  automatically by `clean_name`).
- `order::Int`: Harmonic order to extract. `0` returns the DC component; `n ≥ 1` returns
  the n-th harmonic.
- `output_type::Symbol` (keyword, default `:complex`): Controls the form of the output.
  - `:complex`   — phasor  `Aₙ - im·Bₙ`  (complex amplitudes, `Vector{ComplexF64}`)
  - `:magnitude` — `|Aₙ - im·Bₙ|`  (`Vector{Float64}`)
  - `:phase`     — `∠(Aₙ - im·Bₙ)` in radians  (`Vector{Float64}`)

# Returns
A vector of length `length(sweep_res.sweep_vals)` containing the requested quantity at
each sweep point.

# Details
- For `order == 0` the DC value is the cosine coefficient A₀ (sine is identically zero
  by convention).
- For `order ≥ 1` the phasor is constructed as `Aₙ - im·Bₙ`, consistent with the
  convention `x(t) = Re[(Aₙ - im·Bₙ) e^{inωt}]`.
- Coefficients are fetched via `fetch_harmonic_coeff` when the variable appears directly
  in the HB variable map, or via `reconstruct_from_observed` when the variable is only
  defined through observed (algebraic) equations.
"""
function get_harmonic(h_prob::HarmonicProblem, sweep_res::HarmonicResult, var_name::String, order::Int; output_type::Symbol=:complex)

    var_clean = clean_name(var_name)
    vmap = h_prob.harmonic_system.variable_map
    
    if order == 0
        key = find_varmap_key(vmap, var_clean, 0, :Cos)
        if key !== nothing
            phasor = complex.(fetch_harmonic_coeff(vmap[key], h_prob, sweep_res))
        else
            phasor = complex.(reconstruct_from_observed(h_prob, sweep_res, var_clean, 0, :Cos))
        end
    else
        key_cos = find_varmap_key(vmap, var_clean, order, :Cos)
        key_sin = find_varmap_key(vmap, var_clean, order, :Sin)
        
        if key_cos !== nothing && key_sin !== nothing
            A_vec = fetch_harmonic_coeff(vmap[key_cos], h_prob, sweep_res)
            B_vec = fetch_harmonic_coeff(vmap[key_sin], h_prob, sweep_res)
            phasor = @. (A_vec - im*B_vec)
        else
            A_vec = reconstruct_from_observed(h_prob, sweep_res, var_clean, order, :Cos)
            B_vec = reconstruct_from_observed(h_prob, sweep_res, var_clean, order, :Sin)
            phasor = @. (A_vec - im*B_vec)
        end
    end

    if output_type == :magnitude
        return abs.(phasor)
    elseif output_type == :phase
        return angle.(phasor)
    else
        return phasor
    end
end

"""
    reconstruct_from_observed(h_prob, sweep_res, var_name, order, component) -> Vector{Float64}

Extract a harmonic coefficient for an **observed** (algebraically-defined) variable
across a parameter sweep by symbolically reconstructing the variable's time-domain
expression and then projecting it onto the requested harmonic.

# Arguments
- `h_prob::HarmonicProblem`: The harmonic balance problem containing the system definition,
  variable map, and sweep parameters.
- `sweep_res::HarmonicResult`: The solution returned by the HB solver, containing the
  sweep variable, sweep values, and result arrays for each HB unknown.
- `var_name::String`: Name of the observed variable to reconstruct (matched after
  cleaning to remove namespace prefixes).
- `order::Int`: Harmonic order to extract. `0` gives the DC (constant) component;
  `n ≥ 1` gives the n-th harmonic.
- `component::Symbol`: Which quadrature to extract at `order ≥ 1`.
  - `:Cos` — cosine coefficient  Aₙ
  - `:Sin` — sine coefficient    Bₙ

# Returns
A `Vector{Float64}` of length `length(sweep_res.sweep_vals)`, where each element is
the requested harmonic coefficient of `var_name` evaluated at the corresponding sweep
point.

# Method
1. Collects all observed equations (both from the `HarmonicSystem` and from
   `ModelingToolkit.observed(sys)`) into a substitution dictionary.
2. Locates `var_name` in the dictionary and recursively flattens any
   observed-in-observed dependencies (up to 20 levels deep) so the expression is
   written purely in terms of HB unknowns, parameters, and the time variable `t`.
3. For each sweep point, builds the full time-domain harmonic ansatz
       x(t) = A₀ + Σₙ Aₙ cos(nωt) + Bₙ sin(nωt)
   for every original ODE state variable still present in the expression, using the
   numerical HB results, then substitutes all parameters and HB unknowns.
4. If `t` remains after substitution, the harmonic coefficient is extracted by
   evaluating the resulting function at specific phase points:
   - DC  (`order == 0`):  `(f(0) + f(T/2)) / 2`
   - Cos (`order ≥ 1`):   `(f(0) − f(T/2)) / 2`
   - Sin (`order ≥ 1`):   `(f(T/4) − f(3T/4)) / 2`
   where `T = 2π/ω` is the fundamental period.

# Errors
- Throws an error if `var_name` is not found in the observed equations.
- Throws an error if the expression cannot be reduced to a numeric value at a sweep
  point (i.e. unresolved symbolic variables remain after all substitutions).
"""
function reconstruct_from_observed(h_prob::HarmonicProblem, sweep_res::HarmonicResult, var_name::String, order::Int, component::Symbol)
    sys = h_prob.harmonic_system.complete_sys
    t_sym = h_prob.harmonic_system.t_var
    tvar_uw = Symbolics.unwrap(t_sym)

    obs_dict = Dict{Any, Any}()
    for (k, v) in h_prob.harmonic_system.observed_equations
        obs_dict[Symbolics.unwrap(k)] = Symbolics.unwrap(v)
    end
    for eq in observed(sys)
        obs_dict[Symbolics.unwrap(eq.lhs)] = Symbolics.unwrap(eq.rhs)
    end

    # Find target variable
    var_name_clean = clean_name(var_name)
    target_sym = nothing
    for k in keys(obs_dict)
        k_str = clean_name(k)
        if k_str == var_name_clean || occursin(var_name_clean, k_str)
            target_sym = k
            break
        end
    end
    target_sym === nothing && error("Variable $var_name not found in observed equations. Available: $([clean_name(k) for k in keys(obs_dict)])")

    # Start with the RHS of the observed equation
    expr = obs_dict[target_sym]
    states_set = Set(Symbolics.unwrap(v) for v in unknowns(sys))
    params_set = Set(Symbolics.unwrap(p) for p in parameters(sys))

    # Flatten observed-in-observed dependencies 
    for _ in 1:20
        vars_in_expr = Symbolics.get_variables(expr)
        unknown_vars = filter(v -> !in(v, states_set) && !in(v, params_set) && !isequal(v, tvar_uw), vars_in_expr)
        isempty(unknown_vars) && break
        subs = Dict()
        for u in unknown_vars
            if SymbolicUtils.iscall(u) && SymbolicUtils.operation(u) isa Differential
                base_var = SymbolicUtils.arguments(u)[1]
                haskey(obs_dict, base_var) && (subs[u] = Symbolics.expand_derivatives(Differential(tvar_uw)(obs_dict[base_var])))
            elseif haskey(obs_dict, u)
                subs[u] = obs_dict[u]
            end
        end
        isempty(subs) && break
        expr = Symbolics.substitute(expr, subs)
        expr = Symbolics.expand_derivatives(expr)
    end
    # Identify original ODE state variables still present in expr 
    orig_state_names = Set(k[1] for k in keys(h_prob.harmonic_system.variable_map))
    orig_state_map = Dict{Any, Tuple{String, Int}}()  # symbolic_var => (state_name, deriv_order)
    for v in Symbolics.get_variables(expr)
        in(v, states_set) && continue
        in(v, params_set) && continue
        isequal(v, tvar_uw) && continue
        v_str = clean_name(v)
        found = false
        for name in orig_state_names
            if v_str == name
                orig_state_map[v] = (name, 0); found = true; break
            elseif v_str == name * "ˍt" || endswith(v_str, "₊" * name * "ˍt")
                orig_state_map[v] = (name, 1); found = true; break
            elseif v_str == name * "ˍtt" || endswith(v_str, "₊" * name * "ˍtt")
                orig_state_map[v] = (name, 2); found = true; break
            end
            if !found && SymbolicUtils.iscall(v) && SymbolicUtils.operation(v) isa Differential
                base = SymbolicUtils.arguments(v)[1]
                if clean_name(base) == name
                    orig_state_map[v] = (name, 1); break
                end
            end
        end
    end

    vmap = h_prob.harmonic_system.variable_map
    param_vec = collect(parameters(sys))
    ωvar_sym = h_prob.harmonic_system.ω_var
    N = h_prob.harmonic_system.N
    base_params = h_prob.params
    sweep_idx = findfirst(isequal(sweep_res.sweep_var), param_vec)

    n_points = length(sweep_res.sweep_vals)
    out_vals = zeros(Float64, n_points)

    for i in 1:n_points
        p_dict = Dict{Any, Any}(Symbolics.unwrap(p) => base_params[p] for p in param_vec)
        sweep_idx !== nothing && (p_dict[Symbolics.unwrap(sweep_res.sweep_var)] = sweep_res.sweep_vals[i])

        ω_val = Float64(Symbolics.unwrap(Symbolics.substitute(Symbolics.unwrap(ωvar_sym), p_dict)))
        T_period = 2π / ω_val

        subs = Dict{Any, Any}(p_dict)

        # HB system unknowns → their numeric values from results
        for state_j in collect(unknowns(sys))
            state_j_uw = Symbolics.unwrap(state_j)
            haskey(sweep_res.results, state_j) && (subs[state_j_uw] = sweep_res.results[state_j][i])
        end

        # Original ODE state variables 
        for (orig_var, (orig_name, deriv_order)) in orig_state_map
            x_t   = Num(0.0)
            dx_t  = Num(0.0)
            ddx_t = Num(0.0)

            # Look up a harmonic coefficient from results; default to 0.0 if the
            # symbol was eliminated during compilation (e.g. DC forced to zero).
            lookup = sym -> begin
                haskey(sweep_res.results, sym) && return sweep_res.results[sym][i]
                sym_str = clean_name(sym)
                for (k, v) in sweep_res.results
                    (isequal(k, sym) || clean_name(k) == sym_str) && return v[i]
                end
                return 0.0
            end

            key_dc = find_varmap_key(vmap, orig_name, 0, :Cos)
            key_dc !== nothing && (x_t += lookup(vmap[key_dc]))

            for n in 1:N
                key_c = find_varmap_key(vmap, orig_name, n, :Cos)
                key_s = find_varmap_key(vmap, orig_name, n, :Sin)
                An = key_c !== nothing ? lookup(vmap[key_c]) : 0.0
                Bn = key_s !== nothing ? lookup(vmap[key_s]) : 0.0
                nω = n * ω_val
                x_t   +=  An * cos(nω * t_sym) + Bn * sin(nω * t_sym)
                dx_t  += -An * nω * sin(nω * t_sym) + Bn * nω * cos(nω * t_sym)
                ddx_t += -An * nω^2 * cos(nω * t_sym) - Bn * nω^2 * sin(nω * t_sym)
            end

            subs[orig_var] = Symbolics.unwrap(
                deriv_order == 0 ? x_t : deriv_order == 1 ? dx_t : ddx_t)
        end

        # Substitute parameters, HB unknowns, and full state ansatz
        val_t = Symbolics.substitute(expr, subs)
        val_t = Symbolics.expand_derivatives(val_t)

        # If t remains  extract the  harmonic coefficient by evaluating at two time points first-harmonic expressions.
        remaining = Symbolics.get_variables(val_t)
        if any(isequal(v, tvar_uw) for v in remaining)
            eval_at = tv -> Float64(Symbolics.unwrap(
                Symbolics.substitute(Symbolics.unwrap(val_t), Dict(tvar_uw => tv))))
            if order == 0
                out_vals[i] = (eval_at(0.0) + eval_at(T_period / 2)) / 2
            elseif component == :Cos
                out_vals[i] = (eval_at(0.0) - eval_at(T_period / 2)) / 2
            else  # :Sin
                out_vals[i] = (eval_at(T_period / 4) - eval_at(3 * T_period / 4)) / 2
            end
        else
            try
                out_vals[i] = Float64(Symbolics.unwrap(val_t))
            catch
                vars = Symbolics.get_variables(val_t)
                error("Cannot convert to Float64. Unresolved vars: $vars in expr: $val_t")
            end
        end
    end

    return out_vals
end
    