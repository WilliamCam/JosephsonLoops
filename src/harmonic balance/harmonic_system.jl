using Symbolics
using SymbolicUtils
using NonlinearSolve

struct HarmonicSystem
    system::ModelingToolkit.AbstractSystem
    time_domain_system::ModelingToolkit.System
    ω::Num
    N::Int
    variable_map::Dict{Tuple{String, Int, Symbol}, Num}
    jacobian::Union{Tuple{Matrix{Num}, Matrix{Num}},Nothing}
end

struct HarmonicProblem
    harmonic_system::HarmonicSystem
    parameters::Dict
    swept_parameters::Union{Dict{Num, Vector{Float64}},Nothing}
    U₀::Vector{Float64}
    result::Dict{Num, Matrix{Float64}}
end

struct HarmonicResult
    result::Dict{Num, Vector{Float64}}
    #TODO: add retcodes
end

struct LinearizedProblem
    harmonic_system::HarmonicSystem
    params::Dict
    ωp::Float64
    U₀::Dict{Num, Float64}   
    Ω_vals::AbstractVector
    perturb::Vector{Float64}
end

struct LinearizedResult
    sweep_var::Num
    sweep_vals::AbstractVector
    results::Dict{Num, Vector{ComplexF64}}
end

function solve!(prob::HarmonicProblem; continuation::Bool = true, kwargs...)
    harmonic_system = prob.harmonic_system
    system = harmonic_system.system
    system_unknowns = unknowns(system)
    result = prob.result
    
    #Initialise NonlinearProblem
    system_parameters = merge(Dict(system_unknowns .=> prob.U₀), prob.parameters)
    nonlinear_prob = NonlinearProblem(system, system_parameters)

    #Solve if no parameter sweep
    if prob.swept_parameters === nothing
        return ModelingToolkit.solve(nonlinear_prob; kwargs...)
    end

    #Swept parameters
    for sweep in prob.swept_parameters
        sweep_values = sweep.second
        sweep_parameter = sweep.first

        println("Sweeping $(sweep_parameter) over $(length(sweep_values)) points...")
        output_matrix = result[sweep_parameter]

        working_prob = remake(nonlinear_prob; u0 = prob.U₀, p = [sweep_parameter => first(sweep_values)])
        for (column_index,_val) in enumerate(sweep_values)

            working_prob.ps[sweep_parameter] = _val

            sol = ModelingToolkit.solve(working_prob, kwargs...)

            if continuation
                working_prob.u0 .= sol.u
            end
            output_matrix[:, column_index] .= sol.u  
        end
    end

    return result
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
    J₀ = substitute(h_sys.J0, sub_rules)
    J₁ = substitute(h_sys.J1, sub_rules)

    perturb_c = (prob.perturb)
    var_syms = _ordered_harmonic_vars(h_sys)

    results = Dict{Num, Vector{ComplexF64}}()
    for sym in var_syms
        results[sym] = Vector{ComplexF64}(undef, length(prob.Ω_vals))
    end

    for (i, Ω) in enumerate(prob.Ω_vals)
        mat = J₀ - 1im * (Ω - prob.ωp) * J₁
        #resp = pinv(mat) * perturb_c #Potnetially faster ? approximate method
        resp = mat \ perturb_c
        for (j, sym) in enumerate(var_syms)
            results[sym][i] = resp[j]
        end
    end

    @variables Ω_signal
    return LinearizedResult(Ω_signal, collect(prob.Ω_vals), results)
end

function HarmonicProblem(harmonic_system::HarmonicSystem, parameters::Dict; 
        swept_parameters::Union{Dict{Num, Vector{Float64}}, Nothing}=nothing, 
        U₀::Union{Vector{Float64},Nothing} = nothing, 
        linearise_problem::Bool = false
    )
    system = harmonic_system.system
    system_unknowns = unknowns(system)
    _Nvars = length(system_unknowns)

    #Preallocate results object
    keys_list = collect(keys(swept_parameters))
    _get_result_size(parameter::Num) = (_Nvars, length(swept_parameters[parameter])) 
    results = Dict{Num, Matrix{Float64}}(key => Matrix{Float64}(undef, _get_result_size(key)...) for key in keys_list)

    #Determine inital condition state vector
     if isnothing(U₀)  
        U₀ = fill(0.0, length(unknowns(harmonic_system.system)))
     end

     if linearise_problem
        #TODO: Construct working point for linearised system
        nothing
     end

    return HarmonicProblem(harmonic_system, parameters, swept_parameters, U₀, results)
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
function HarmonicSystem(sys, ωvar::Num, N::Int; tearing::Bool=true, determine_jacobian::Bool=false)
    tvar = Num(ModelingToolkit.get_iv(sys))
    eqs, states = get_full_equations(sys, tvar)

    # harmonic_equation's M==1 branch wraps `eqs = [eqs]`, expecting a scalar input.
    # Unwrap here so we can pass Vectors uniformly from get_full_equations.
    # `Num(states[1])` matches the prototype's input type, so `states[k]` indexing works inside.
    eqs_arg    = length(states) == 1 ? eqs[1]         : eqs
    states_arg = length(states) == 1 ? Num(states[1]) : states

    if determine_jacobian
        nonlinear_sys, _, variable_map, jac = harmonic_equation(eqs_arg, states_arg, tvar, ωvar, N; jac=true)
    else
        nonlinear_sys, _, variable_map = harmonic_equation(eqs_arg, states_arg, tvar, ωvar, N)
        jac = nothing
    end
    
    sys_eqs, sys_vars = equations(nonlinear_sys), unknowns(nonlinear_sys)
    
    if length(sys_eqs) > length(sys_vars)
        n_drop = length(sys_eqs) - length(sys_vars)
        @warn "Harmonic system is overdetermined: $(length(sys_eqs)) equations for $(length(sys_vars)) variables. " *
              "Dropping the last equation(s). Caution: This behavior depends on variable order."
        sys_eqs = sys_eqs[1:end-n_drop]
    end
        
    sys_eqs_built = tearing ? sys_eqs : [0 ~ eq.lhs - eq.rhs for eq in sys_eqs]

    @named nonlinear_sys = NonlinearSystem(sys_eqs_built, sys_vars, parameters(sys))
    complete_sys = tearing ? mtkcompile(nonlinear_sys) : complete(nonlinear_sys)
    
    return HarmonicSystem(complete_sys, sys, ωvar, N, variable_map, jac)
end