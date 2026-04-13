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
    sweep_var::Union{Num, Nothing}
    sweep_vals::Union{AbstractVector, Nothing}
    u0::Vector{Float64}
end


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
    HarmonicProblem(harmonic_sys, params; sweep_var=nothing, sweep_vals=nothing) -> HarmonicProblem
    HarmonicProblem(sys, ωvar, params; tearing=true, N=1, sweep_var=nothing, sweep_vals=nothing) -> HarmonicProblem

Construct a `HarmonicProblem` from either a pre-built `HarmonicSystem` or directly
from a time-domain ODE system (which will be expanded into a `HarmonicSystem` internally).

Prefer passing a pre-built `HarmonicSystem` when solving multiple times with different
parameters, so the (potentially expensive) harmonic expansion is not repeated.

# Arguments
- `harmonic_sys::HarmonicSystem`: A pre-built harmonic system (first method).
- `sys::ModelingToolkit.AbstractSystem`: The time-domain ODE system to expand (second method).
- `ωvar::Num`: Symbolic variable for the fundamental angular frequency ω (second method only).
- `params`: Parameter map from symbolic variables to numeric values (e.g. `Dict(ω => 1.0, R => 50.0)`).

# Keywords
- `sweep_var::Union{Num,Nothing}=nothing`: Parameter to sweep over.
- `sweep_vals::Union{AbstractVector,Nothing}=nothing`: Values for the sweep parameter.
- `N::Int=1`: Harmonic truncation order — second method only, passed to `HarmonicSystem`.
- `tearing::Bool=true`: Whether to structurally simplify via `mtkcompile` — second method only.

# Returns
A `HarmonicProblem` ready to be passed to `solve`.
"""
function HarmonicProblem(harmonic_sys::HarmonicSystem, params; sweep_var::Union{Num, Nothing}=nothing, sweep_vals::Union{AbstractVector, Nothing}=nothing)
    u0 = fill(0.0, length(unknowns(harmonic_sys.complete_sys)))
    return HarmonicProblem(harmonic_sys, params, sweep_var, sweep_vals, u0)
end
"""
    solve(prob::HarmonicProblem; kwargs...) -> NonlinearSolution | HarmonicResult

Solve the harmonic balance system. If `sweep_var` and `sweep_vals` are set on the
problem, performs a parameter sweep using continuation (warm-starting each point from
the previous solution). Otherwise, solves at a single parameter point.

All keyword arguments are forwarded to `ModelingToolkit.solve` (e.g. `alg`, `abstol`).

# Returns
- **No sweep**: A `NonlinearSolution` whose harmonic coefficients can be read via `sol[var]`.
- **With sweep**: A `HarmonicResult` containing the swept parameter, sweep values, and
  result arrays for each HB unknown.
"""
function solve(prob::HarmonicProblem; kwargs...)
    hsys = prob.harmonic_system
    sys = hsys.complete_sys
    system_unknowns = unknowns(sys)

    # no sweep
    if prob.sweep_var === nothing || prob.sweep_vals === nothing
        combined_args = merge(Dict(system_unknowns .=> prob.u0), prob.params)
        nl_prob = NonlinearProblem(sys, combined_args)
        return ModelingToolkit.solve(nl_prob; kwargs...)
    end

    #Sweep solve 
    sweep_var = prob.sweep_var
    sweep_vals = prob.sweep_vals

    # Remove sweep variable from fixed params to avoid conflicts
    current_params = copy(prob.params)
    if haskey(current_params, sweep_var)
        delete!(current_params, sweep_var)
    end
    sweep_var_uw = Symbolics.unwrap(sweep_var)
    if haskey(current_params, sweep_var_uw)
        delete!(current_params, sweep_var_uw)
    end

    current_params[sweep_var] = first(sweep_vals)
    combined_args = merge(Dict(system_unknowns .=> prob.u0), current_params)
    nl_prob = NonlinearProblem(sys, combined_args)

    # Results containers
    results = Dict{Num, Vector{Float64}}()
    for v in system_unknowns
        results[v] = Float64[]
        sizehint!(results[v], length(sweep_vals))
    end
    observed_results = Dict{Any, Vector{Float64}}()

    println("Sweeping $(sweep_var) over $(length(sweep_vals)) points...")
    last_u = nl_prob.u0

    for val in sweep_vals
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
    _get_coeff_expr(sym, h_prob) -> Num

Return a symbolic expression for a single HB coefficient symbol `sym`.

If `sym` is a direct unknown of the compiled system it is returned as-is (it will be
substituted numerically later by `get_harmonic`).  If it was eliminated by `mtkcompile`,
its expression is recovered from the compiled system's observed equations.
"""
function _get_coeff_expr(sym::Num, h_prob::HarmonicProblem)
    sys = h_prob.harmonic_system.complete_sys
    # If sym is a system unknown, return it directly — it will be substituted numerically
    for u in unknowns(sys)
        isequal(u, sym) && return sym
    end
    # isequal can fail for array-indexed vars (different underlying array objects after
    # mtkcompile).  Retry with string comparison, returning the COMPILED system's symbol
    # so that build_function can map it to the correct input position.
    sym_str = clean_name(sym)
    for u in unknowns(sys)
        clean_name(u) == sym_str && return Num(Symbolics.unwrap(u))
    end
    # sym was eliminated by tearing — find it in observed equations
    for eq in observed(sys)
        isequal(eq.lhs, sym) && return Num(Symbolics.unwrap(eq.rhs))
    end
    for eq in observed(sys)
        clean_name(eq.lhs) == sym_str && return Num(Symbolics.unwrap(eq.rhs))
    end
    # Not found anywhere — return the symbol itself (may be zero-forced by structure)
    @warn "_get_coeff_expr: could not resolve $(clean_name(sym)) — returning as-is. unknowns=$(clean_name.(unknowns(sys))), observed_lhs=$(clean_name.(eq.lhs for eq in observed(sys)))"
    return sym
end



"""
    get_harmonic(h_prob, var_name, order) -> Function

Build a compiled harmonic extractor for a named variable at a given harmonic order.

Returns a **function** `f(sweep_res::HarmonicResult) -> Vector{ComplexF64}` that
evaluates the phasor at every sweep point using fast numeric evaluation (no symbolic
work at call time).

The symbolic expression is built and compiled via `Symbolics.build_function` once
when `get_harmonic` is called, so the returned function can be reused cheaply on
different sweep results from the same `HarmonicProblem`.

# Usage
```julia
# Compile once (does symbolic work + build_function)
extract_I = get_harmonic(h_prob, "C1₊i", 1)

# Evaluate on any sweep result (fast numeric only)
phasor = extract_I(sweep_res)
mag    = abs.(extract_I(sweep_res))
phase  = angle.(extract_I(sweep_res))
```

# Arguments
- `h_prob::HarmonicProblem`: The harmonic balance problem.
- `var_name::String`: Name of the variable of interest.
- `order::Int`: Harmonic order to extract. `0` = DC, `n ≥ 1` = n-th harmonic.
"""
function get_harmonic(h_prob::HarmonicProblem, var_name::String, order::Int)
    var_clean = clean_name(var_name)
    vmap = h_prob.harmonic_system.variable_map
    sys  = h_prob.harmonic_system.complete_sys
    #Build symbolic phasor expression using raw vmap symbols 
    if order == 0
        key = find_varmap_key(vmap, var_clean, 0, :Cos)
        phasor_re = key !== nothing ? vmap[key] : reconstruct_from_observed(h_prob, var_clean, 0, :Cos)
        phasor_im = Num(0.0)
    else
        key_cos = find_varmap_key(vmap, var_clean, order, :Cos)
        key_sin = find_varmap_key(vmap, var_clean, order, :Sin)
        if key_cos !== nothing && key_sin !== nothing
            A_expr = vmap[key_cos]
            B_expr = vmap[key_sin]
        else
            A_expr = reconstruct_from_observed(h_prob, var_clean, order, :Cos)
            B_expr = reconstruct_from_observed(h_prob, var_clean, order, :Sin)
        end
        phasor_re = A_expr
        phasor_im = -B_expr
    end

    all_vmap_syms = collect(values(vmap))
    param_vec     = collect(parameters(sys))
    all_vars_uw   = vcat(Symbolics.unwrap.(all_vmap_syms), Symbolics.unwrap.(param_vec))

    #Compile build_function 
    f_re = Symbolics.build_function(Symbolics.unwrap(phasor_re), all_vars_uw; expression=Val{false})
    f_im = Symbolics.build_function(Symbolics.unwrap(phasor_im), all_vars_uw; expression=Val{false})

    # Return closure for fast numeric evaluation ─
    n_vmap   = length(all_vmap_syms)
    n_params = length(param_vec)

    # Pre-compute fixed parameter values
    param_vals = zeros(Float64, n_params)
    for (j, p) in enumerate(param_vec)
        val = _find_sym(Num(p), h_prob.params)
        val !== nothing && (param_vals[j] = Float64(val))
    end

    return function(sweep_res::HarmonicResult)
        n_points = length(sweep_res.sweep_vals)
        input_vec = zeros(Float64, n_vmap + n_params)
        input_vec[n_vmap+1:end] .= param_vals
        vmap_vecs = Vector{Union{Nothing, Vector{Float64}}}(undef, n_vmap)
        for (j, sym) in enumerate(all_vmap_syms)
            vmap_vecs[j] = _find_sym(sym, sweep_res.results)
        end
        # Locate sweep parameter index in param section
        sweep_j = nothing
        sweep_str = clean_name(sweep_res.sweep_var)
        for (j, p) in enumerate(param_vec)
            if clean_name(p) == sweep_str
                sweep_j = j; break
            end
        end
        phasor = Vector{ComplexF64}(undef, n_points)
        for i in 1:n_points
            for j in 1:n_vmap
                input_vec[j] = vmap_vecs[j] !== nothing ? vmap_vecs[j][i] : 0.0
            end
            sweep_j !== nothing && (input_vec[n_vmap + sweep_j] = sweep_res.sweep_vals[i])
            phasor[i] = complex(f_re(input_vec), f_im(input_vec))
        end
        return phasor
    end
end

function get_harmonic_Will(h_prob::HarmonicProblem, var_name::String, order::Int)
    var_clean = clean_name(var_name) #drops (t)
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
            A = vmap[key_cos] # Num
            B = vmap[key_sin] # Num
            phasor_expr = (A - 1im * B) # Equation / Expression
            phasor = build_function(phasor_expr, [A, B], expression=Val{false}) # Function
        else
            A_vec = reconstruct_from_observed(h_prob, sweep_res, var_clean, order, :Cos)
            B_vec = reconstruct_from_observed(h_prob, sweep_res, var_clean, order, :Sin)
            phasor = @. (A_vec - im*B_vec)
        end
    end
    return phasor
end

"""
    reconstruct_from_observed(h_prob, var_name, order, component) -> Num

Return a **symbolic expression** for the requested harmonic coefficient of an
observed (algebraically-defined) variable.

The expression is written purely in terms of HB system unknowns and parameters,
ready to be compiled via `Symbolics.build_function` by the caller (`get_harmonic`).

# Arguments
- `h_prob::HarmonicProblem`: The harmonic balance problem.
- `var_name::String`: Name of the observed variable to reconstruct.
- `order::Int`: Harmonic order (`0` = DC, `n ≥ 1` = n-th harmonic).
- `component::Symbol`: `:Cos` or `:Sin`.

# Returns
A `Num` symbolic expression.
"""
function reconstruct_from_observed(h_prob::HarmonicProblem, var_name::String, order::Int, component::Symbol)
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
    orig_state_map = Dict{Any, Tuple{String, Int}}()
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
    ωvar_sym = h_prob.harmonic_system.ω_var
    N = h_prob.harmonic_system.N

    # Build symbolic ansatz substitution 
    # Use _get_coeff_expr to pre-resolve any vmap symbols that were eliminatedby mtkcompile 
    ansatz_subs = Dict{Any, Any}()
    for (orig_var, (orig_name, deriv_order)) in orig_state_map
        x_t = dx_t = ddx_t = Num(0.0)

        key_dc = find_varmap_key(vmap, orig_name, 0, :Cos)
        key_dc !== nothing && (x_t += _get_coeff_expr(vmap[key_dc], h_prob))

        for n in 1:N
            key_c = find_varmap_key(vmap, orig_name, n, :Cos)
            key_s = find_varmap_key(vmap, orig_name, n, :Sin)
            An = key_c !== nothing ? _get_coeff_expr(vmap[key_c], h_prob) : Num(0.0)
            Bn = key_s !== nothing ? _get_coeff_expr(vmap[key_s], h_prob) : Num(0.0)
            nω = n * ωvar_sym
            x_t   +=  An * cos(nω * t_sym) + Bn * sin(nω * t_sym)
            dx_t  += -An * nω * sin(nω * t_sym) + Bn * nω * cos(nω * t_sym)
            ddx_t += -An * nω^2 * cos(nω * t_sym) - Bn * nω^2 * sin(nω * t_sym)
        end

        ansatz_subs[orig_var] = Symbolics.unwrap(
            deriv_order == 0 ? x_t : deriv_order == 1 ? dx_t : ddx_t)
    end

    # Substitute symbolic ansatz into expression
    val_t = Symbolics.substitute(expr, ansatz_subs)
    val_t = Symbolics.expand_derivatives(val_t)

    #  Isolate the target harmonic component 
    remaining = Symbolics.get_variables(val_t)
    if any(isequal(v, tvar_uw) for v in remaining)
        val_expanded = Symbolics.expand(val_t)
        isolate = Dict{Any, Any}()
        for k in 1:(2*N+1)
            cos_k = Symbolics.unwrap(cos(Num(k) * ωvar_sym * t_sym))
            sin_k = Symbolics.unwrap(sin(Num(k) * ωvar_sym * t_sym))
            if k == order
                isolate[cos_k] = (component == :Cos ? 1.0 : 0.0)
                isolate[sin_k] = (component == :Sin ? 1.0 : 0.0)
            else
                isolate[cos_k] = 0.0
                isolate[sin_k] = 0.0
            end
        end
        final_expr = Symbolics.substitute(Symbolics.unwrap(val_expanded), isolate)
    else
        final_expr = val_t
    end

    return Num(final_expr)
end
  