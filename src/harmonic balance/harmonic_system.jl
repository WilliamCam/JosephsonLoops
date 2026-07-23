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

    fourier_indicies::Vector{Tuple{Int,Int}}
    coeff_map::Dict{Tuple{Tuple{Int,Int},Symbol}, Num}
end

struct HarmonicSystem
    system::ModelingToolkit.AbstractSystem
    time_domain_system::ModelingToolkit.System
    ω::Tuple{Num,Num} # (ωp, ωs)
    N::Int
    harmonic_ansatz::Vector{Num}
    variable_map::Dict{Num, FourierBasis}
    jacobian::Union{Tuple{Matrix{Num}, Matrix{Num}},Nothing}
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


struct LinearisedProblem
    harmonic_system::HarmonicSystem
    jacobian::Tuple{Matrix{Num},Matrix{Num}}
    Ωs::Union{Float64,Vector{Float64}}
    Ωp::Union{Float64,Vector{Float64}}
    parameters::Dict
    parameter_sweep::Union{Dict{Num, Vector{Float64}},Nothing}
    δU::Union{ComplexF64,Vector{ComplexF64}}
    result::HarmonicResult
end

function solve!(harmonic_problem::HarmonicProblem; continuation::Bool = true, kwargs...)
    result = harmonic_problem.result.solution
    nonlinear_prob = harmonic_problem.problem
    sweep_space = harmonic_problem.parameter_sweep
    continuation_axis = 1

    working_prob = remake(nonlinear_prob; u0 = harmonic_problem.U₀)

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

function solve!(linear_problem::LinearisedProblem)
    ω, ω_values = linear_problem.ω_sweep
    result = linear_problem.result

    if isnothing(linear_problem.parameter_sweep)
        println("Performing sweep $(ω) over $(length(ω_values)) points...")
        output_array = result.solution[ω]
        K = size(linear_problem.jacobian[1],1)
        numeric_substitution = merge(linear_problem.parameters, linear_problem.working_point, Dict(ω => linear_problem.ω_pump))
        # Float64 conversion errors loudly if any symbol survived the substitution
        # (e.g. a missing working-point entry) instead of silently solving symbolically.
        J₀ = Float64.(Symbolics.value.(substitute(linear_problem.jacobian[1], numeric_substitution)))
        J₁ = Float64.(Symbolics.value.(substitute(linear_problem.jacobian[2], numeric_substitution)))

        U_small_signal = SVector{K}(linear_problem.U_perturbation)

        for (column_index, Ω) in enumerate(ω_values)
            mat = J₀ - 1im * (Ω - linear_problem.ω_pump) * J₁
            output_array[:, column_index] .= _linear_solve(mat, U_small_signal, column_index == 1)
        end
    else
        for param_sweep in linear_problem.parameter_sweep
            sweep_parameter = param_sweep.first
            sweep_values = param_sweep.second
            println("Performing 2D Sweep $(sweep_parameter) over $(length(ω_values)*length(sweep_values)) points...")
            output_array = result.solution[sweep_parameter]

            U_small_signal = (linear_problem.U_perturbation)

            for (parameter_index,_val) in enumerate(sweep_values)
                numeric_substitution = merge(linear_problem.parameters, linear_problem.working_point,
                    Dict(ω => linear_problem.ω_pump, sweep_parameter => _val))
                J₀ = Float64.(Symbolics.value.(substitute(linear_problem.jacobian[1], numeric_substitution)))
                J₁ = Float64.(Symbolics.value.(substitute(linear_problem.jacobian[2], numeric_substitution)))

                for (column_index, Ω) in enumerate(ω_values)
                    mat = J₀ - 1im * (Ω - linear_problem.ω_pump) * J₁
                    output_array[:, parameter_index, column_index] .= _linear_solve(mat, U_small_signal, column_index == 1)
                end
            end
        end

    end
    return result
end


function _linear_solve(mat::AbstractMatrix, U::AbstractVector, warn_once::Bool)
    F = LinearAlgebra.lu(mat; check = false)
    LinearAlgebra.issuccess(F) && return F \ U
    warn_once && @warn "Linearised jacobian is singular (free coefficient / gauge freedom); using the minimum-norm solution."
    return LinearAlgebra.pinv(mat) * U
end

function _nl_solve_method!(prealloc_array::Array{ComplexF64}, problem::NonlinearSolve.NonlinearProblem, sweep_parameter::Num, sweep_values::Union{Float64, Vector{Float64}}; 
    continuation::Bool=true, parameter_index::Union{Int, Vector{Int}}, kwargs...
        )
    for (column_index, val) in enumerate(sweep_values)

        problem.ps[sweep_parameter] = val

        #TODO: Benchmark allocations for performance
        sol = ModelingToolkit.solve(problem, kwargs...)

        if continuation 
            problem.u0 .= sol.u
        end
        if parameter_index != 0
            prealloc_array[:, parameter_index..., column_index] .= sol.u
        else
            prealloc_array[:, column_index] .= sol.u
        end
    end
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

function LinearisedProblem(harmonic_system::HarmonicSystem, parameters::Dict,
    δU::AbstractVector; 
    parameter_sweep::Union{Dict{Num, Vector{Float64}}, Nothing}=nothing, 
    U₀::Union{Vector{Float64},Nothing} = nothing,
    kwargs...)
    @assert !isnothing(harmonic_system.jacobian) "Linearised problem requires a harmonic system with determined jacobian"

    system = harmonic_system.system
    ωp, ωs = harmonic_system.ω
    system_unknowns = unknowns(system)
    _Nvars = size(harmonic_system.jacobian[1], 1)
    results_size = [_Nvars]
    #Preallocate results object
    key_map = Dict{Num,Int}()
    if !isnothing(parameter_sweep)
        if has_key(parameter_sweep, ωs)
            Ωs = pop!(parameter_sweep, ωs)
            results_size = vcat(results_size, length(Ωs))
        end
        if has_key(parameter_sweep, ωp)
            Ωp = pop!(parameter_sweep, ωp)
            results_size = vcat(results_size, length(Ωp))
        end

        keys_list = collect(keys(parameter_sweep))
        for (sweep_key, index) in enumerate(keys_list)
            results_size = vcat(result_size, length(parameter_sweep[sweep_key]))
            key_map[sweep_key] = index
        end
    else
        Ωp = parameters[ωp]
        Ωs = parameters[ωs]
    end
    result_arr = Array{ComplexF64}(undef,results_size...)
    output = HarmonicResult(key_map, result_arr)

    #Determine inital condition state vector
    if isnothing(U₀)  
        U₀ = fill(0.0, length(unknowns(system)))
    end
    return LinearisedProblem(harmonic_system, harmonic_system.jacobian, Ωp, Ωs, parameters, parameter_sweep, ComplexF64.(δU), output)
end

#TODO:: LinearisedProbelm(::HarmonicProblem)

function HarmonicSystem(sys, ω::Union{Num,Tuple{Num,Num}}, N::Int; tearing::Bool=true, determine_jacobian::Bool=false)
    if typeof(ω) !== Tuple
        ω = (ω, Num(0))
    end
    determine_jacobian && (tearing = false)

    tvar = Num(ModelingToolkit.get_iv(sys))
    eqs, states, _, _ = get_full_equations(sys)

    eqs_arg    = length(states) == 1 ? eqs[1]         : eqs
    states_arg = length(states) == 1 ? Num(states[1]) : states
    nonlinear_sys, X, variable_map, jac = harmonic_equation(eqs_arg, Num.(states_arg), tvar, ω, N; jac=determine_jacobian) 
    
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