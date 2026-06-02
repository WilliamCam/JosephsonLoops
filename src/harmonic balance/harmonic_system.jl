using Symbolics
using SymbolicUtils
using NonlinearSolve
using BenchmarkTools

struct HarmonicSystem
    system::ModelingToolkit.AbstractSystem
    time_domain_system::ModelingToolkit.System
    ω::Num
    N::Int
    harmonic_ansatz::Any
    variable_map::Dict{Tuple{String, Int, Symbol}, Num}
    jacobian::Union{Tuple{Matrix{Num}, Matrix{Num}},Nothing}
end

struct HarmonicResult
    solution::Dict{Num, AbstractArray}
    #TODO: add retcodes
end

struct HarmonicProblem
    problem::NonlinearSolve.NonlinearProblem
    ω_sweep::Tuple{Num, Union{Float64,Vector{Float64}}}
    parameters::Dict
    parameter_sweep::Union{Dict{Num, Vector{Float64}},Nothing}
    U₀::Vector{Float64}
    result::HarmonicResult
end

struct LinearisedProblem
    jacobian::Tuple{Matrix{Num},Matrix{Num}}
    ω_sweep::Tuple{Num, Union{Float64,Vector{Float64}}}
    parameters::Dict
    parameter_sweep::Union{Dict{Num, Vector{Float64}},Nothing}
    working_point::Dict{Num, Float64}
    ω_pump::Float64
    U_perturbation::Vector{Float64}
    result::HarmonicResult
end

function solve!(harmonic_problem::HarmonicProblem; continuation::Bool = true, kwargs...)
    result = harmonic_problem.result
    nonlinear_prob = harmonic_problem.problem
    ω, ω_values = harmonic_problem.ω_sweep

    if isnothing(harmonic_problem.parameter_sweep)
        println("Performing sweep $(ω) over $(length(ω_values)) points...")
        output_array = result.solution[ω]
        print(typeof(output_array))
        working_prob = remake(nonlinear_prob; u0 = harmonic_problem.U₀, p = [ω => first(ω_values)])
    
        _nl_solve_method!(output_array, working_prob, ω, ω_values, continuation=continuation)
    else
        for param_sweep in harmonic_problem.parameter_sweep
            sweep_parameter = param_sweep.first
            sweep_values = param_sweep.second

            println("Performing 2D Sweep $(sweep_parameter) over $(length(ω_values)*length(sweep_values)) points...")

            output_array = result.solution[sweep_parameter]
            #output array indexed as [state_variable, parameter_value, ω_value]

            working_prob = remake(nonlinear_prob; u0 = harmonic_problem.U₀, p = [sweep_parameter => first(sweep_values)])
            #TODO: This should be made into a parallel Ensemble sweep using DiffEq.jl
            for (parameter_index,_val) in enumerate(sweep_values)
                working_prob.ps[sweep_parameter] = _val
                _nl_solve_method!(output_array, working_prob, ω, ω_values, continuation=continuation, parameter_index = parameter_index)
            end
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
        numeric_substitution = merge(linear_problem.parameters, linear_problem.working_point, Dict(ω => linear_problem.ω_pump))
        J₀ = simplify(substitute(linear_problem.jacobian[1], numeric_substitution))
        J₁ = simplify(substitute(linear_problem.jacobian[2], numeric_substitution))

        U_small_signal = (linear_problem.U_perturbation)

        for (column_index, Ω) in enumerate(ω_values)
            mat = simplify(J₀ - 1im * (Ω - linear_problem.ω_pump) * J₁)
            resp = mat \ U_small_signal
            output_array[:, column_index] .= resp
        end
    else
        for param_sweep in linear_problem.parameter_sweep
            sweep_parameter = param_sweep.first
            sweep_values = param_sweep.second
            println("Performing 2D Sweep $(sweep_parameter) over $(length(ω_values)*length(sweep_values)) points...")
            ps = copy(linear_problem.parameter_sweep)
            output_array = result.solution[sweep_parameter]

            for (parameter_index,_val) in enumerate(sweep_values)
                ps[sweep_parameter] = _val
                numeric_substitution = merge(ps, linear_problem.U₀, Dict(ω => linear_problem.ω_pump))
                J₀ = substitute(linear_problem.jacobian[1], numeric_substitution)
                J₁ = substitute(linear_problem.jacobian[2], numeric_substitution)

                U_small_signal = (linear_problem.U_perturbation)

                for (column_index, Ω) in enumerate(ω_values)
                    mat = J₀ - 1im * (Ω - linear_problem.ω_pump) * J₁
                    resp = mat \ perturb_c
                    output_array[:, parameter_index, column_index] .= resp
                end
            end
        end

    end
end

function _nl_solve_method!(prealloc_array::Array{Float64}, problem::NonlinearSolve.NonlinearProblem, ω_variable::Num, ω_sweep_values::Union{Float64, Vector{Float64}}; 
    continuation::Bool=true, parameter_index::Int = 0, kwargs...
        )
    for (column_index, frequency) in enumerate(ω_sweep_values)

        problem.ps[ω_variable] = frequency

        #TODO: Benchmark allocations for performance
        sol = ModelingToolkit.solve(problem, kwargs...)

        if continuation 
            problem.u0 .= sol.u
        end
        if parameter_index != 0
            prealloc_array[:, parameter_index, column_index] .= sol.u
        else
            prealloc_array[:, column_index] .= sol.u
        end
    end
end

function HarmonicProblem(harmonic_system::HarmonicSystem, ω_values::Union{Float64, Vector{Float64}}, parameters::Dict; 
        parameter_sweep::Union{Dict{Num, Vector{Float64}}, Nothing}=nothing, 
        U₀::Union{Vector{Float64},Nothing} = nothing, 
        linear_response::Union{Tuple{Float64,Vector{Float64}},Nothing} = nothing,
        kwargs...
    )
    system = harmonic_system.system
    ωvar = harmonic_system.ω
    ω_sweep = (ωvar, ω_values)
    system_unknowns = unknowns(system)
    _Nvars = isnothing(linear_response) ? length(system_unknowns) : size(harmonic_system.jacobian[1], 1)

    #Preallocate results object
    if isnothing(parameter_sweep)
        result_dict = Dict{Num, Array{Float64}}(ωvar => Array{Float64}(undef, _Nvars, length(ω_values)))
    else
        keys_list = collect(keys(parameter_sweep))
        _get_result_size(parameter::Num) = (_Nvars, length(parameter_sweep[parameter]), length(ω_values)) 
        result_dict = Dict{Num, Array{Float64}}(key => Array{Float64}(undef, _get_result_size(key)...) for key in keys_list)
    end
    output = HarmonicResult(result_dict)

    #Determine inital condition state vector
     if isnothing(U₀)  
        U₀ = fill(0.0, length(unknowns(system)))
     end

    #Initialise NonlinearProblem
    system_parameters = merge(Dict(system_unknowns .=> U₀), parameters)
    nonlinear_prob = NonlinearProblem(system, system_parameters)
    if isnothing(linear_response)
        return HarmonicProblem(nonlinear_prob, ω_sweep, parameters, parameter_sweep, U₀, output)
    else
        #TODO: Assert jacobian is generated for linear response
        pump_frequency, perturbation = linear_response
        working_prob = remake(nonlinear_prob; u0 = U₀, p = [harmonic_system.ω => pump_frequency])
        sol = ModelingToolkit.solve(working_prob, kwargs...)
        working_point = Dict(
            key => begin
                try
                    sol[key]
                catch e
                    0.0
                end
            end 
            for key in values(harmonic_system.variable_map)
        )
        return LinearisedProblem(harmonic_system.jacobian, ω_sweep, parameters, parameter_sweep, working_point, pump_frequency, perturbation, output)
     end
end

function HarmonicSystem(sys, ωvar::Num, N::Int; tearing::Bool=true, determine_jacobian::Bool=false)
    tvar = Num(ModelingToolkit.get_iv(sys))
    eqs, states = get_full_equations(sys, tvar)

    # harmonic_equation's M==1 branch wraps `eqs = [eqs]`, expecting a scalar input.
    # Unwrap here so we can pass Vectors uniformly from get_full_equations.
    # `Num(states[1])` matches the prototype's input type, so `states[k]` indexing works inside.
    eqs_arg    = length(states) == 1 ? eqs[1]         : eqs
    states_arg = length(states) == 1 ? Num(states[1]) : states
    if determine_jacobian
        nonlinear_sys, X, variable_map, jac = harmonic_equation(eqs_arg, states_arg, tvar, ωvar, N; jac=true)
    else
        nonlinear_sys, X, variable_map = harmonic_equation(eqs_arg, states_arg, tvar, ωvar, N)
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

    return HarmonicSystem(complete_sys, sys, ωvar, N, X, variable_map, jac)
end