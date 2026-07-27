
function substitute_to_fixpoint(expr, subs; max_iters = 20, transform = identity)
    for _ in 1:max_iters
        next = transform(Symbolics.substitute(expr, subs))
        isequal(next, expr) && break
        expr = next
    end
    return expr
end

function get_solution(h_prob::HarmonicProblem, expression::Num, order::Union{Int,Tuple{Int,Int}} = 1;
    time_domain_solution::Bool = false)
    expr_vars = Num.(get_variables(expression))
    if time_domain_solution
        order = nothing
    end
    if length(expr_vars) > 1
        sub_dict = Dict{Num,Num}()
        for var in expr_vars
            #TODO assert var is in parameters or observed or states
            if var_is_in(parameters(h_prob.harmonic_system.system), var)
                continue
            end
            harmonic_var = expression_for_var(h_prob.harmonic_system, var, order)
            sub_dict[var] = harmonic_var
        end
        harmonic_expr = substitute(expression, sub_dict)
    else
        harmonic_expr = expression_for_var(h_prob.harmonic_system, expression, order)
    end
    harmonic_expr = simplify(expand(harmonic_expr))
    sol = h_prob.result.solution
    output_arr = Array{ComplexF64}(undef, size(sol)[2:end]...)
    apply_harmonic_expression!(output_arr, h_prob, harmonic_expr)
    return output_arr
end


#  Harmonic expressions
"""
Retrieves the expression for target harmonic specified to some order by either pulling directly from
the state variables via h_sys.variable_map, or by using reconstruct_from_observed() if the the target variable 'var'
was removed via system tearing in MTK.structural_simplify().

Order is a tuple that specifies the fourier indicies of the harmonic ansatz in the format (ωp,ωs). i.e. (2,-1) will return 
the expression for the coefficient of the (2*ωp - ωs)*t term in the solution for var

If no order specified returns the full expression (i.e. all fourier components of expanded solution)
"""
function expression_for_var(h_sys::HarmonicSystem, var::Num, order::Union{Int,Tuple{Int,Int},Nothing})
    if typeof(order)==Int
        order = (order, 0)
    end
    _, states, _, _ = get_full_equations(h_sys.time_domain_system)
   
    if var_is_in(states, var)
        state_index = var_index(states, Symbolics.unwrap(var))
        expr = h_sys.harmonic_ansatz[state_index]
    else
        expr = reconstruct_from_observed(h_sys, var)
    end
    if order !== nothing
        pump_index, signal_index = order
        ωp, ωs = h_sys.ω
        I = Symbolics.coeff(expr, cos((pump_index * ωp + signal_index * ωs)*t))
        Q = Symbolics.coeff(expr, sin((pump_index * ωp + signal_index * ωs)*t)) 
        return Num.(I + 1im*Q)
    else
        return Num.(expr)
    end
end


function reconstruct_from_observed(h_sys::HarmonicSystem, observed_var::Num)
    hs = h_sys
    #TODO: Are the observed equations of the harmonic system (i.e. observed(hs.system)) necessary ?
    observed_map = Dict{Any, Any}(Symbolics.unwrap(eq.lhs) => Symbolics.unwrap(eq.rhs)
                                  for sys in (hs.time_domain_system, hs.system) for eq in observed(sys))
    var_symbol = Symbolics.unwrap(observed_var)
    target_rhs = get(observed_map, var_symbol, nothing)
    target_rhs === nothing && error("Variable $var not found in observed equations")

    flat_subs = Dict(k => v for (k, v) in observed_map if !(Symbolics.symtype(k) <: AbstractArray))
    expr = substitute_to_fixpoint(target_rhs, flat_subs; transform = Symbolics.expand_derivatives)
    (_, states, diffvars, diff2vars) = get_full_equations(hs.time_domain_system)
    ansatz_subs = Dict{Any, Any}()
    for (k, s) in enumerate(states)
        ansatz_subs[Symbolics.unwrap(s)] = Symbolics.unwrap(hs.harmonic_ansatz[k])
    end
    for (i, dv) in enumerate(diff2vars)
        k = var_index(states, Symbolics.unwrap(diffvars[i]))
        ansatz_subs[Symbolics.unwrap(dv)] = Symbolics.unwrap(Symbolics.expand_derivatives(D(hs.harmonic_ansatz[k]))) 
    end

    reconstructed_expr = Num(Symbolics.expand(Symbolics.substitute(expr, ansatz_subs)))
    return simplify(reconstructed_expr)
end


function apply_harmonic_expression!(output::AbstractArray, h_prob::HarmonicProblem, expression::Complex{Num})
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

# Linearised counterpart: the response array rows follow the jacobian's `vars` ordering,
# and the responses are complex, so evaluate the same expressions with complex inputs.
# For a plain state this reduces to resp[row(cos)] + im*resp[row(sin)].

function compile_expression(expr::Complex{Num}, system::ModelingToolkit.System, parameter_values::Dict{Num,Float64};
        sweep_parameters::Vector{Num} = [Num(nothing)]
    )
    input_syms = Symbolics.unwrap.(unknowns(system))
    observed_subs = Dict(eq.lhs => eq.rhs for eq in observed(system)
                         if !(Symbolics.symtype(Symbolics.unwrap(eq.lhs)) <: AbstractArray))
    fixed_params = Dict(Num(p) => Float64(parameter_values[Num(p)]) for p in parameters(system)
                        if !var_is_in(sweep_parameters, Num(p)) && haskey(parameter_values, Num(p)))
    inputs = vcat(input_syms, Symbolics.unwrap.(sweep_parameters)...)
    compiled_function = Symbolics.build_function(
        Symbolics.unwrap(Symbolics.substitute(substitute_to_fixpoint(expr, observed_subs), fixed_params)),
        inputs; expression = Val{false})
    return compiled_function
end
