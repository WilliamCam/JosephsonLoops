struct LinearisedProblem
    harmonic_system::HarmonicSystem
    problem::NonlinearSolve.NonlinearProblem
    parameters::Dict
    parameter_sweep::Union{Vector{Pair{Num, Vector{Float64}}}, Nothing}
    Ωs::Vector{Float64}
    δU::Union{ComplexF64,Vector{ComplexF64}}
    result::HarmonicResult
end

function solve!(linear_problem::LinearisedProblem; kwargs...)
    result = linear_problem.result.solution
    nonlinear_prob = linear_problem.problem
    J₀, J₁ = linear_problem.harmonic_system.jacobian
    sweep_space = linear_problem.parameter_sweep

    Ωp = linear_problem.harmonic_system.ω
    single_tone = Ωp[2] == Num(0)

    Ωs = linear_problem.Ωs
    δU = linear_problem.δU

    if isnothing(sweep_space)
        K = size(J₀,1)
        sol = ModelingToolkit.solve(nonlinear_prob, kwargs...)

        working_point = Dict([var=>sol[var] for var in unknowns(linear_problem.harmonic_system.system)])
        numeric_substitution = merge(linear_problem.parameters, working_point)

        J₀ = Float64.(Symbolics.value.(substitute(J₀, numeric_substitution)))
        J₁ = Float64.(Symbolics.value.(substitute(J₁, numeric_substitution)))

        U_small_signal = SVector{K}(δU)
        if single_tone
            Ωp_value = linear_problem.parameters[Ωp[1]]
            for (column_index, Ω) in enumerate(Ωs)
                mat = J₀ - 1im * (Ω - Ωp_value) * J₁
                result[:, column_index] .= _linear_solve(mat, U_small_signal, column_index == 1)
            end
        else
            for (pump_index, Ω_pump) in enumerate(Ωp)
                Ωp_value = linear_problem.parameters[Ω_pump]
                for (column_index, Ω) in enumerate(Ωs)
                    mat = J₀ - 1im * (Ω - Ωp_value) * J₁
                    output_array[:, pump_index, column_index] .= _linear_solve(mat, U_small_signal, column_index == 1)
                end
            end
        end
    else
        #parameter sweep

    end
    return result
end

function _linear_solve(mat::AbstractMatrix, U::AbstractVector, warn_once::Bool)
    F = LinearAlgebra.lu(mat; check = false)
    LinearAlgebra.issuccess(F) && return F \ U
    warn_once && @warn "Linearised jacobian is singular (free coefficient / gauge freedom); using the minimum-norm solution."
    return LinearAlgebra.pinv(mat) * U
end

function LinearisedProblem(harmonic_system::HarmonicSystem, parameters::Dict,
    δU::Vector{ComplexF64}, Ωs::Vector{Float64}; 
    parameter_sweep::Union{Vector{Pair{Num, Vector{Float64}}}, Nothing}=nothing, 
    U₀::Union{Vector{Float64},Nothing} = nothing,
    kwargs...)

    @assert !isnothing(harmonic_system.jacobian) "Linearised problem requires a harmonic system with determined jacobian"

    system = harmonic_system.system
    system_unknowns = unknowns(system)
    _Nvars = size(harmonic_system.jacobian[1], 1)

    # Will linearise around each pump tone if there is more than one
    ω1, ω2 = harmonic_system.ω
    ω2 == Num(0) ? results_size = [_Nvars, length(Ωs)] : results_size = [_Nvars, 2, length(Ωs)]  
    dep_params = Num[]
    #Preallocate results object
    #TODO: Support for linear parameter sweeps i.e. capacitance/resistance -> doesnt require re-solving NL system working point!
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
    working_prob = remake(nonlinear_prob; u0 = U₀)

    return LinearisedProblem(harmonic_system, working_prob, parameters, parameter_sweep, Ωs, (δU), output)
end

#TODO:: LinearisedProbelm(::HarmonicProblem)

#Creates system pertubation response in harmonic frame
function perturbation_response(h_sys::HarmonicSystem, source_param::Num, parameters::Dict; amplitude::Float64 = 1.0)
    t = Num(ModelingToolkit.get_iv(h_sys.time_domain_system))
    eqs, _, _, _ = get_full_equations(h_sys.time_domain_system)
    ω1, ω2 = h_sys.ω
    fourier_indicies = h_sys.variable_map[unknowns(h_sys.time_domain_system)[1]].fourier_indicies
    Nt = 2*length(fourier_indicies)-1

    U = zeros(ComplexF64, length(eqs) * Nt)

    fixed_params = copy(parameters)

    delete!(fixed_params, ω1)
    delete!(fixed_params, ω2)

    system_response_symbolic = Symbolics.jacobian([eq.lhs - eq.rhs for eq in eqs], [source_param])
    system_response = Symbolics.substitute(system_response_symbolic, fixed_params)
    @assert any(x -> !isequal(x, Num(0)), system_response) "Parameter $source_param not found in time domain system"
    zero_subs = Dict{Num, Float64}()
    for n in eachindex(fourier_indicies)
        pump1_index, pump2_index = fourier_indicies[n]
        zero_subs[cos((pump1_index * ω1 + pump2_index * ω2)*t)] = 0.0
        zero_subs[sin((pump1_index * ω1 + pump2_index * ω2)*t)] = 0.0
    end
  
    for (k, eq) in enumerate(system_response)
        isequal(Symbolics.simplify(eq), Num(0)) && continue
        base = (k - 1) * Nt
        U[base + 1] -= amplitude * Symbolics.substitute(eq, zero_subs)
        print(Symbolics.substitute(eq, zero_subs)) 
        for n in eachindex(fourier_indicies)
            if fourier_indicies[n] == (0,0)
                continue
            end
            pump1_index, pump2_index = fourier_indicies[n]
            I = Symbolics.coeff(eq, cos((pump1_index * ω1 + pump2_index * ω2)*t))
            Q = Symbolics.coeff(eq, sin((pump1_index * ω1 + pump2_index * ω2)*t)) 
            if isequal(I, 0.0) && isequal(Q, 0.0)
                continue
            end
            U[base + 2n + 1] -= amplitude * (I - im * Q)
            U[base + 2n] -= amplitude * (Q + im * I)
        end
    end
    return U
end

