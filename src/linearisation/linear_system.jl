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
    J₀, J₁, J₂ = linear_problem.harmonic_system.jacobian
    sweep_space = linear_problem.parameter_sweep

    Ωp = linear_problem.harmonic_system.ω

    Ωs = linear_problem.Ωs
    δU = linear_problem.δU

    if isnothing(sweep_space)
        K = size(J₀,1)
        #TODO: ????
        # gauge-free DC coefficients can make the Newton jacobian singular; fall back to LM
        sol = try
            ModelingToolkit.solve(nonlinear_prob; kwargs...)
        catch
            ModelingToolkit.solve(nonlinear_prob, NonlinearSolve.LevenbergMarquardt())
        end

        working_point = Dict([var=>sol[var] for var in unknowns(linear_problem.harmonic_system.system)])
        numeric_substitution = merge(linear_problem.parameters, working_point)

        J₀ = Float64.(Symbolics.value.(substitute(J₀, numeric_substitution)))
        J₁ = Float64.(Symbolics.value.(substitute(J₁, numeric_substitution)))
        J₂ = Float64.(Symbolics.value.(substitute(J₂, numeric_substitution)))

        U_small_signal = SVector{K}(δU)
        # linearised around the single declared pump (tone 1): δ = Ω − ω1; the δ²J₂ term
        # makes the response exact in δ for the truncated basis
        Ωp_value = linear_problem.parameters[Ωp[1]]
        #TODO: optional argument to only compute taylor expansion to second order (for speed)
        for (column_index, Ω) in enumerate(Ωs)
            δ = Ω - Ωp_value
            mat = J₀ - 1im * δ * J₁ - δ^2 * J₂
            result[:, column_index] .= _linear_solve(mat, U_small_signal, column_index == 1)
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

    # Linearised around one declared pump (tone 1), so the response is always [vars, Ω]
    results_size = [_Nvars, length(Ωs)]
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

function get_solution(lin_prob::LinearisedProblem, expression::Num, order::Union{Int,Tuple{Int,Int}} = 1)
    h_sys = lin_prob.harmonic_system
    expr_vars = Num.(get_variables(expression))
    if length(expr_vars) > 1
        sub_dict = Dict{Num,Num}()
        for var in expr_vars
            #TODO assert var is in parameters or observed or states
            if var_is_in(parameters(h_sys.system), var)
                continue
            end
            harmonic_var = expression_for_var(h_sys, var, order)
            sub_dict[var] = harmonic_var
        end
        harmonic_expr = substitute(expression, sub_dict)
    else
        harmonic_expr = expression_for_var(h_sys, expression, order)
    end
    harmonic_expr = simplify(expand(harmonic_expr))
    sol = lin_prob.result.solution
    output_arr = Array{ComplexF64}(undef, size(sol)[2:end]...)
    apply_harmonic_expression!(output_arr, lin_prob, harmonic_expr)
    return output_arr
end

# Linearised counterpart: the response array rows follow the jacobian's `vars` ordering
# (see harmonic_equation), not unknowns(system), and the responses are complex, so the
# same phasor expression is evaluated with complex inputs. The signal frequency enters
# through the tone symbols (the differentiated ansatz carries them), so they are compiled
# as inputs and fed per column. For a plain state this reduces to
# resp[row(cos)] + im*resp[row(sin)].
function apply_harmonic_expression!(output::AbstractArray, lin_prob::LinearisedProblem, expression::Complex{Num})
    h_sys = lin_prob.harmonic_system
    result = lin_prob.result.solution
    ω1, _ = h_sys.ω

    _, states, _, _ = get_full_equations(h_sys.time_domain_system)
    input_syms = Num[]
    for state in states
        basis = h_sys.variable_map[state]
        push!(input_syms, basis.dc_coeff)
        for n in 1:length(basis.fourier_indicies)-1
            push!(input_syms, basis.cos_coeffs[n], basis.sin_coeffs[n])
        end
    end
    @assert length(input_syms) == size(result, 1) "Jacobian variable ordering does not match linear response rows"

    # Linearisation is always around a single declared pump (any other tone of a
    # multi-tone system stays at its parameter value), so the response is [vars, Ω] and
    # the signal frequency is fed through the pump symbol per column.
    @assert ndims(result) == 2 "Linearised response must be [vars, Ω] — declare a single pump when solving"
    output_func = compile_expression(expression, h_sys.system, lin_prob.parameters;
        sweep_parameters = [ω1], input_syms = Symbolics.unwrap.(input_syms))

    state_len = length(input_syms)
    input_vec = Vector{ComplexF64}(undef, state_len + 1)
    @views for (column_index, Ω) in enumerate(lin_prob.Ωs)
        input_vec[1:state_len] .= result[:, column_index]
        input_vec[end] = Ω
        output[column_index] = output_func(input_vec)
    end
end

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
            # fourier_indicies[n] fills coefficient slot n-1 (entry 1 is DC), whose
            # projected rows are base + 2(n-1) (cos) and base + 2(n-1)+1 (sin).
            U[base + 2n - 1] -= amplitude * (I - im * Q)
            U[base + 2n - 2] -= amplitude * (Q + im * I)
        end
    end
    return U
end

function apply_harmonic_expression!(output::AbstractArray, h_prob::LinearProblem, expression::Complex{Num})
    system = h_prob.harmonic_system.system
    result = h_prob.result.solution
    sweep = h_prob.parameter_sweep
    sweep_params = h_prob.result.dependent_parameters
    if !isnothing(sweep)
        output_func = compile_expression(expression, system, h_prob.parameters, sweep_parameters = sweep_params)
        input_vec = similar(result, size(result, 1) + length(sweep))
        state_len = size(result, 1)

        @views for (linear_idx, cart_index) in enumerate(CartesianIndices(axes(result)[2:end]))

            input_vec[1:state_len] .= result[:, cart_index]

            for (idx, param_sweep) in enumerate(sweep)
                p_vals = last(param_sweep)
                input_vec[state_len + idx] = p_vals[linear_idx]
            end

            output[cart_index] = output_func(input_vec)
        end
    else
        print("Non sweep output")
    end
end