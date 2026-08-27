using Symbolics
using SymbolicUtils
using NonlinearSolve
using BenchmarkTools
using StaticArrays

struct FourierBasis
    dc_coeff::Num
    d_dc_coeff::Num
    cos_coeffs::Symbolics.Arr{Num, 1}
    sin_coeffs::Symbolics.Arr{Num, 1}

    # Symbolic first derivative variable e.g (d/dt A₁) needed for linearised jacobian J₁
    d_cos_coeffs::Symbolics.Arr{Num, 1}
    d_sin_coeffs::Symbolics.Arr{Num, 1}

    # Symbolic second derivative variables (d²/dt² A₁) needed for the δ²·J₂ term that
    # makes the linearised response exact in detuning (previously first-order only)
    d2_dc_coeff::Num
    d2_cos_coeffs::Symbolics.Arr{Num, 1}
    d2_sin_coeffs::Symbolics.Arr{Num, 1}

    fourier_indicies::Vector{Tuple{Int,Int}}
    coeff_map::Dict{Tuple{Tuple{Int,Int},Symbol}, Num}
end

struct HarmonicSystem
    system::ModelingToolkit.AbstractSystem
    time_domain_system::ModelingToolkit.System
    ω::Tuple{Num,Num} # (ω1, ω2) in multi-tone systems
    N::Int
    harmonic_ansatz::Vector{Num}
    variable_map::Dict{Num, FourierBasis}
    jacobian::Union{Tuple{Matrix{Num}, Matrix{Num}}, Tuple{Matrix{Num}, Matrix{Num}, Matrix{Num}}, Nothing}
end

struct HarmonicResult
    dependent_parameters::Union{Num,Vector{Num}}
    solution::AbstractArray
    #TODO: add retcodes
end

struct HarmonicProblem
    harmonic_system::HarmonicSystem
    problem::NonlinearSolve.NonlinearProblem
    parameters::Dict
    parameter_sweep::Union{Vector{Pair{Num, Vector{Float64}}}, Nothing}
    U₀::Vector{Float64}
    result::HarmonicResult
end


function solve!(harmonic_problem::HarmonicProblem; continuation::Bool = true, kwargs...)
    result = harmonic_problem.result.solution
    nonlinear_prob = harmonic_problem.problem
    sweep_space = harmonic_problem.parameter_sweep
    continuation_axis = 1

 #continuation, not sure if needed ask will 
    working_prob = remake(nonlinear_prob; u0 = copy(harmonic_problem.U₀))

    if isnothing(sweep_space)
        sol = ModelingToolkit.solve(working_prob, kwargs...)
        result[:] .= sol.u
    else
        param_keys = [sweep_space[d].first for d in eachindex(sweep_space)]
        update_parameters! = ModelingToolkit.setp(working_prob, param_keys)
        param_space = CartesianIndices(axes(result)[2:end]) 
        for cart_index in param_space
            current_params = [sweep_space[axis].second[cart_index[axis]] for axis in eachindex(sweep_space)]

            if cart_index[continuation_axis] == 1 && cart_index != first(param_space)
                _cont_step = false
                print("Paused continuation")
            else
                _cont_step = continuation
            end
  
            update_parameters!(working_prob, current_params)
            #TODO: This should be made into a parallel Ensemble sweep using DiffEq.jl
            sol = ModelingToolkit.solve(working_prob, kwargs...)
            for var_idx in axes(result, 1)
                full_idx = CartesianIndex(var_idx, cart_index.I...)
                result[full_idx] = sol.u[var_idx]
            end
            _cont_step ? working_prob.u0 .= sol.u : working_prob.u0 .= harmonic_problem.U₀
        end
    end

    return result
end

function HarmonicProblem(harmonic_system::HarmonicSystem, parameters::Dict; 
        parameter_sweep::Union{Vector{Pair{Num, Vector{Float64}}}, Nothing}=nothing, 
        U₀::Union{Vector{Float64},Nothing} = nothing, 
        kwargs...
    )
    system = harmonic_system.system
    system_unknowns = unknowns(system)
    _Nvars = length(system_unknowns)
    results_size = [_Nvars]
    dep_params = Num[]
    #Preallocate results object
    if !isnothing(parameter_sweep)
        for sweep in parameter_sweep
            sweep_param, vals = sweep
            results_size = vcat(results_size, length(vals))
            append!(dep_params, sweep_param)
        end
    end
    result_arr = Array{ComplexF64}(undef,results_size...) 
    output = HarmonicResult(dep_params, result_arr)

    #Determine inital condition state vector
    if isnothing(U₀)  
        U₀ = fill(0.0, length(unknowns(system)))
    end

    #Initialise NonlinearProblem
    system_parameters = merge(Dict(system_unknowns .=> U₀), parameters)
    nonlinear_prob = NonlinearProblem(system, system_parameters)

    return HarmonicProblem(harmonic_system, nonlinear_prob, parameters, parameter_sweep, U₀, output)
end

function HarmonicSystem(sys, ω::Union{Num,Tuple{Num,Num}}, N::Int; tearing::Bool=true, determine_jacobian::Bool=false,
        intermod_order::Int=0, tones::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
        commensurate_tol::Real=1e-6, max_denominator::Int=1000, oversample::Int=2)
    # `typeof(ω) !== Tuple` was always true (Tuple{Num,Num} !== the UnionAll Tuple),
    # double-wrapping real two-tone inputs.
    if !(ω isa Tuple)
        ω = (ω, Num(0))
    end

    # Two-tone grid selection, all in the backend (no integers in the API):
    #  * `tones` given (the pump frequencies as plain numbers, any consistent units —
    #    only the ratio is used) and the ratio is within commensurate_tol of a small
    #    rational: 1D commensurate grid (ω1 = p·ω0, ω2 = q·ω0), tone 2 snapped.
    #  * otherwise (no tones, or a nearly incommensurate ratio): 2D torus grid — both
    #    tones stay fully symbolic, any ratio solvable, no snapping.
    commensurate = nothing
    if !isequal(ω[2], Num(0))
        if tones === nothing
            @info "two-tone HB: no tones given — using the 2D torus grid (any tone ratio, both tones symbolic)"
        else
            rat = rationalize(Int, float(tones[2] / tones[1]), tol = commensurate_tol)
            p, q = denominator(rat), numerator(rat)
            if max(p, q) > max_denominator
                @info "tone ratio $(tones[2]/tones[1]) has no rational form with integers ≤ $max_denominator within tol = $commensurate_tol — using the 2D torus grid"
            else
                snapped = (q / p) * tones[1]
                @info "commensurate two-tone grid: ω1:ω2 = $p:$q, tone 2 snapped to $(snapped) (relative shift $(abs(snapped / tones[2] - 1)))"
                commensurate = (p, q)
            end
        end
    end

    #use of jacobian in linearisation requires all system variable to be present in final model
    determine_jacobian && (tearing = false)

    tvar = Num(ModelingToolkit.get_iv(sys))
    eqs, states, _, _ = get_full_equations(sys)

    eqs_arg    = length(states) == 1 ? eqs[1]         : eqs
    states_arg = length(states) == 1 ? Num(states[1]) : states
    nonlinear_sys, X, variable_map, jac = harmonic_equation(eqs_arg, Num.(states_arg), tvar, ω, N;
        jac=determine_jacobian, intermod_order=intermod_order, commensurate=commensurate,
        oversample=oversample)
    
    sys_eqs, sys_vars = equations(nonlinear_sys), unknowns(nonlinear_sys)
    
    #TODO: Remove the correct equation
    if length(sys_eqs) > length(sys_vars)
        n_drop = length(sys_eqs) - length(sys_vars)
        @warn "Harmonic system is overdetermined: $(length(sys_eqs)) equations for $(length(sys_vars)) variables. " *
              "Dropping the last equation(s). Caution: This behavior depends on variable order."
        sys_eqs = sys_eqs[1:end-n_drop]
    end
        
    built_equations = tearing ? sys_eqs : [0 ~ eq.lhs - eq.rhs for eq in sys_eqs]

    @named nonlinear_sys = NonlinearSystem(built_equations, sys_vars, parameters(sys))
    complete_sys = tearing ? mtkcompile(nonlinear_sys) : complete(nonlinear_sys)

    return HarmonicSystem(complete_sys, sys, ω, N, X, variable_map, jac)
end