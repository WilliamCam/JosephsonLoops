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
    J0::Matrix{Num}     # empty unless built with linear=true
    J1::Matrix{Num}     # empty unless built with linear=true
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
end

struct LinearizedProblem
    harmonic_system::HarmonicSystem
    params::Dict
    ωp::Float64
    U₀::Dict{Num, Float64}      # working point: every harmonic coefficient that appears in J0/J1
    Ω_vals::AbstractVector
    perturb::Vector{Float64}
end

struct LinearizedResult
    sweep_var::Num
    sweep_vals::AbstractVector
    results::Dict{Num, Vector{ComplexF64}}
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

function rotate_to_harmonic_frame(M, N, Nt, harmonic_system)
    # M: Number of state variables (dimension of state vector)
    block_rows = 2N + 1
    # Total Matrix Dimensions
    total_rows = M * block_rows
    total_cols = M * Nt
    # Initialize the large block-diagonal matrix
    Γ_total = zeros(Num, total_rows, total_cols)
    # Generate the single-variable operator (the logic you already have)
    Γ_single = Matrix{Num}(undef, block_rows, Nt)
    for j in 1:Nt
        # DC Term
        Γ_single[1, j] = Num(1//Nt)
        # AC Terms
        for n in 1:N
            phase = n * (j - 1) * (2π / Nt)
            Γ_single[n + 1, j] = Num((2//Nt) * cos(phase))      # Cosine
            Γ_single[N + 1 + n, j] = Num((2//Nt) * sin(phase))  # Sine
        end
    end
    # Place the single-variable operator into the block diagonal
    for d in 1:M
        row_range = (d-1)*block_rows + 1 : d*block_rows
        col_range = (d-1)*Nt + 1 : d*Nt
        Γ_total[row_range, col_range] .= Γ_single
    end
    rotated_system = Γ_total * [equation.lhs for equation in harmonic_system]    
    return simplify.(rotated_system)
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
- `linear::Bool=false`: If `true`, also populate `J0`, `J1` (Kosata 2022, eq. 5.12)
  via `harmonic_equation(...; jac=true)`. Otherwise those fields are left empty.

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
function HarmonicSystem(sys, ωvar::Num; tearing::Bool=true, N::Int=1, linear::Bool=false)
    tvar = Num(ModelingToolkit.get_iv(sys))
    eqs, states = get_full_equations(sys, tvar)

    # harmonic_equation's M==1 branch wraps `eqs = [eqs]`, expecting a scalar input.
    # Unwrap here so we can pass Vectors uniformly from get_full_equations.
    # `Num(states[1])` matches the prototype's input type, so `states[k]` indexing works inside.
    eqs_arg    = length(states) == 1 ? eqs[1]         : eqs
    states_arg = length(states) == 1 ? Num(states[1]) : states

    if linear
        nonlinear_sys, _, variable_map, (J0, J1) = harmonic_equation(eqs_arg, states_arg, tvar, ωvar, N; jac=true)
    else
        nonlinear_sys, _, variable_map = harmonic_equation(eqs_arg, states_arg, tvar, ωvar, N)
        J0 = Matrix{Num}(undef, 0, 0)
        J1 = Matrix{Num}(undef, 0, 0)
    end
    
    sys_eqs = equations(nonlinear_sys)
    sys_vars = unknowns(nonlinear_sys)
    
    if length(sys_eqs) > length(sys_vars)
        n_drop = length(sys_eqs) - length(sys_vars)
        @warn "System is overdetermined: $(length(sys_eqs)) equations for $(length(sys_vars)) variables. " *
              "Dropping the last equation(s). Caution: This behavior depends on variable order."
        sys_eqs = sys_eqs[1:end-n_drop]
    end
        
    # When skipping mtkcompile, MTK 10.x requires equations in `0 ~ residual` form;
    # harmonic_equation produces `residual ~ 0`, so flip when we won't tear.
    sys_eqs_built = tearing ? sys_eqs : [0 ~ eq.lhs - eq.rhs for eq in sys_eqs]
    @named nonlinear_sys = NonlinearSystem(sys_eqs_built, sys_vars, parameters(sys))
    complete_sys = tearing ? mtkcompile(nonlinear_sys) : complete(nonlinear_sys)
    
    # Store original observed equations to enable reconstruction of non-state variables (eg. C1.i)
    observed_eqs = Dict{Any, Any}(eq.lhs => eq.rhs for eq in observed(sys))
    
    return HarmonicSystem(complete_sys, ωvar, tvar, N, variable_map, observed_eqs, J0, J1)
end

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
    return HarmonicResult(sweep_var, collect(sweep_vals), results)
end

"""
    solve(prob::LinearizedProblem) -> LinearizedResult

Solve the linearized harmonic-balance system around the working point. At each `Ω` in
`prob.Ω_vals`, form `mat = J0 - i(Ω-ωp) J1` (Kosata 2022, eq. 5.12) with `J0, J1`
substituted at the working point, solve `mat \\ perturb`, and store the complex
response of every harmonic coefficient.
"""
function solve(prob::LinearizedProblem)
    h_sys = prob.harmonic_system

    # Numeric substitution: working point + parameters + ω => ωp
    sub_rules = merge(prob.params, prob.U₀, Dict(h_sys.ω_var => prob.ωp))
    J0_num = ComplexF64.(Symbolics.value.(simplify.(substitute(h_sys.J0, sub_rules))))
    J1_num = ComplexF64.(Symbolics.value.(simplify.(substitute(h_sys.J1, sub_rules))))

    perturb_c = ComplexF64.(prob.perturb)
    var_syms = _ordered_harmonic_vars(h_sys)

    results = Dict{Num, Vector{ComplexF64}}()
    for sym in var_syms
        results[sym] = Vector{ComplexF64}(undef, length(prob.Ω_vals))
    end

    for (i, Ω) in enumerate(prob.Ω_vals)
        mat = J0_num - 1im * (Ω - prob.ωp) * J1_num
        # `pinv` instead of `\`: J0/J1 carry gauge-redundant zero columns (e.g. capacitor
        # DC phase), making `mat` rank-deficient. Min-norm solution gives 0 along those
        # gauge directions, which is the physically correct choice.
        resp = pinv(mat) * perturb_c
        for (j, sym) in enumerate(var_syms)
            results[sym][i] = resp[j]
        end
    end

    @variables Ω_signal
    return LinearizedResult(Ω_signal, collect(prob.Ω_vals), results)
end

# Index of (var_name, order, component) in the J0/J1 column ordering.
# State k uses letter pair (coeff_labels[2k-1], coeff_labels[2k]) — so 'A','C','E','G' map to states 1..4.
# Within a state, ordering is [DC, cos₁, sin₁, cos₂, sin₂, …].
function _harmonic_var_index(h_sys::HarmonicSystem, var_name::String, order::Int, component::Symbol)
    haskey(h_sys.variable_map, (var_name, 0, :Cos)) ||
        error("Variable $var_name not found in variable_map")
    dc_sym = h_sys.variable_map[(var_name, 0, :Cos)]
    letter = string(Symbolics.unwrap(dc_sym))[1]
    state_idx = (Int(letter) - Int('A')) ÷ 2 + 1
    offset = order == 0 ? 0 : (component == :Cos ? 2*order - 1 : 2*order)
    return (state_idx - 1) * (2*h_sys.N + 1) + offset + 1
end

# Enumerate all harmonic-coefficient symbols in J0/J1 column order.
# State ordering is recovered from the letter assignment in variable_map values.
function _ordered_harmonic_vars(h_sys::HarmonicSystem)
    state_names = unique([k[1] for k in keys(h_sys.variable_map)])
    sort!(state_names; by = name -> begin
        sym = h_sys.variable_map[(name, 0, :Cos)]
        Int(string(Symbolics.unwrap(sym))[1])
    end)
    var_syms = Num[]
    for name in state_names
        push!(var_syms, h_sys.variable_map[(name, 0, :Cos)])
        for n in 1:h_sys.N
            push!(var_syms, h_sys.variable_map[(name, n, :Cos)])
            push!(var_syms, h_sys.variable_map[(name, n, :Sin)])
        end
    end
    return var_syms
end

"""
    LinearizedProblem(h_sys, params, ωp, U₀, Ω_vals;
                      perturb_var, perturb_order=1, perturb_component=:Cos, perturb_amplitude=1.0)
        -> LinearizedProblem

Construct a linear-analysis problem around the working point `U₀` at pump frequency `ωp`.

# Arguments
- `h_sys::HarmonicSystem`: must be built with `linear=true` so `J0, J1` are populated.
- `params::Dict`: parameter map. Must include the pump-frequency parameter set to `ωp`.
- `ωp::Real`: pump frequency.
- `U₀::Dict{Num, Float64}`: working-point value for every harmonic coefficient appearing in J0/J1.
- `Ω_vals::AbstractVector`: small-signal frequency sweep.

# Keywords (perturb spec)
- `perturb_var::String`: name of the physical variable to perturb (e.g. `"P1₊i"`).
- `perturb_order::Int=1`: harmonic order of the perturbation.
- `perturb_component::Symbol=:Cos`: `:Cos` or `:Sin`.
- `perturb_amplitude::Real=1.0`: kick magnitude.
"""
function LinearizedProblem(h_sys::HarmonicSystem, params::Dict, ωp::Real,
                           U₀::Dict{Num, Float64}, Ω_vals::AbstractVector;
                           perturb_var::String, perturb_order::Int=1,
                           perturb_component::Symbol=:Cos, perturb_amplitude::Real=1.0)
    isempty(h_sys.J0) && error("HarmonicSystem must be built with linear=true")
    n_vars = size(h_sys.J0, 1)
    perturb = zeros(Float64, n_vars)
    idx = _harmonic_var_index(h_sys, perturb_var, perturb_order, perturb_component)
    perturb[idx] = Float64(perturb_amplitude)
    return LinearizedProblem(h_sys, params, Float64(ωp), U₀, Ω_vals, perturb)
end

