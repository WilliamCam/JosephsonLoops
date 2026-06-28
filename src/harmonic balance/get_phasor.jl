# Extracts the numeric phasor Aₙ + iBₙ from a harmonic-balance result for some order in ω.
using ModelingToolkit

# Repeatedly substitute `subs` into `expr` until it stops changing. Observed equations can
# reference other observed equations, so a single pass may not fully flatten the result.
# `transform` is applied after each substitution (e.g. expand_derivatives).
function substitute_to_fixpoint(expr, subs; max_iters = 20, transform = identity)
    for _ in 1:max_iters
        next = transform(Symbolics.substitute(expr, subs))
        isequal(next, expr) && break
        expr = next
    end
    return expr
end

function get_output(h_prob::HarmonicProblem, result::HarmonicResult, var::Num, order::Int = 1)
    expression = get_harmonic_expression(h_prob.harmonic_system, var, order)
    return apply_harmonic_expression(h_prob, result, expression)
end

# Small-signal response phasor of `var` (the variable symbol) at the requested order. The harmonic system
# is passed in (rather than read off the problem) so the harmonic_system field can later
# be dropped from the problem structs. Restricted to state variables: an observed
# quantity's expression would have to be differentiated at the working point (δf = ∇f·δc),
# not evaluated at the perturbation — TODO.
function get_output(h_sys::HarmonicSystem, lin_prob::LinearisedProblem, result::HarmonicResult, var::Num, order::Int = 1)
    haskey(h_sys.variable_map, (Symbolics.unwrap(var), max(order, 1), :Cos)) ||
        error("Linearised get_output supports state variables only; $var is not a state.")
    expression = get_harmonic_expression(h_sys, var, order)
    return apply_harmonic_expression(h_sys, lin_prob, result, expression)
end

#  Harmonic expressions

function get_harmonic_expression(h_sys::HarmonicSystem, var_name::String, order::Int)
    variable_map = h_sys.variable_map
    var_symbol = Symbolics.unwrap(var)  
    if order == 0
        phasor_re = get(variable_map, (var_symbol, 0, :Cos), nothing)
        phasor_re === nothing && (phasor_re = reconstruct_from_observed(h_sys, var, 0, :Cos))
        return phasor_re, Num(0.0)
    end
    phasor_re = get(variable_map, (var_symbol, order, :Cos), nothing)
    phasor_im = get(variable_map, (var_symbol, order, :Sin), nothing)
    phasor_re === nothing && (phasor_re = reconstruct_from_observed(h_sys, var, order, :Cos))
    phasor_im === nothing && (phasor_im = reconstruct_from_observed(h_sys, var, order, :Sin))
    return phasor_re, phasor_im
end

get_harmonic_expression(h_prob::HarmonicProblem, var::Num, order::Int) =
    get_harmonic_expression(h_prob.harmonic_system, var, order)

# Compile the symbolic phasor and evaluate it at every ω point. The expression is reduced
# to the system unknowns and ω, then mapped over the solved `[unknown × ω point]` array.
function apply_harmonic_expression(h_prob::HarmonicProblem, result::HarmonicResult, expression::Tuple{Num,Num})
    system = h_prob.harmonic_system.system
    isnothing(h_prob.parameter_sweep) || error("get_output supports ω-only sweeps")
    ω, ω_values = h_prob.ω_sweep
    ω_vec = ω_values isa Number ? [Float64(ω_values)] : ω_values

    #TODO: functionality for parameter sweeps
    solution = result.solution[ω]
    input_syms = Symbolics.unwrap.(unknowns(system))

    #A bit unsure if this is computationally efficient, need to ask will - ai suggested
    #if its purely a state then you just grab the cos/sin rows from sol array, otherwise compile
    cos_row = findfirst(x -> isequal(x, Symbolics.unwrap(expression[1])), input_syms)
    sin_row = findfirst(x -> isequal(x, Symbolics.unwrap(expression[2])), input_syms)
    if cos_row !== nothing && sin_row !== nothing
        return complex.(real.(solution[cos_row, :]), real.(solution[sin_row, :]))
    end

    f_re, f_im = compile_phasor(expression, input_syms, h_prob.parameters, system, ω)

    # Harmonic-balance solutions are real; evaluate both coefficients at the solved state.
    return map(axes(solution, 2)) do i
        input_vec = real.([solution[:, i]; ω_vec[i]])
        complex(f_re(input_vec), f_im(input_vec))
    end
end

# Linearised counterpart: the response array rows follow the jacobian's `vars` ordering,
# and the responses are complex, so evaluate the same expressions with complex inputs.
# For a plain state this reduces to resp[row(cos)] + im*resp[row(sin)].
function apply_harmonic_expression(h_sys::HarmonicSystem, lin_prob::LinearisedProblem, result::HarmonicResult, expression::Tuple{Num,Num})
    isnothing(lin_prob.parameter_sweep) || error("get_output supports ω-only sweeps")
    ω, ω_values = lin_prob.ω_sweep
    ω_vec = ω_values isa Number ? [Float64(ω_values)] : ω_values

    solution = result.solution[ω]
    states = h_sys.full_eqs.states   # cache
    vmap = h_sys.variable_map
    vars = Num[]
    for s in states
        state_symbol = Symbolics.unwrap(s)
        push!(vars, vmap[(state_symbol, 0, :Cos)])
        for n in 1:h_sys.N
            push!(vars, vmap[(state_symbol, n, :Cos)])
            push!(vars, vmap[(state_symbol, n, :Sin)])
        end
    end

    input_syms = Symbolics.unwrap.(vars)
    cos_row = findfirst(x -> isequal(x, Symbolics.unwrap(expression[1])), input_syms)
    sin_row = findfirst(x -> isequal(x, Symbolics.unwrap(expression[2])), input_syms)
    if cos_row !== nothing && sin_row !== nothing
        return solution[cos_row, :] .+ im .* solution[sin_row, :]
    end

    f_re, f_im = compile_phasor(expression, input_syms,
                                lin_prob.parameters, h_sys.system, ω)

    return map(axes(solution, 2)) do i
        input_vec = [solution[:, i]; complex(ω_vec[i])]
        f_re(input_vec) + im * f_im(input_vec)
    end
end

# Shared compile step: reduce both coefficient expressions to `input_syms` plus ω and
# build callable functions.
#   observed_subs — observed equations, skipping whole-array reconstructions
#     (mtkcompile's `E => change_origin(...)`) that would become uncompilable; and
#   fixed_params  — every parameter value in `parameter_values` except the swept ω.
function compile_phasor(expression::Tuple{Num,Num}, input_syms, parameter_values, system, ω)
    observed_subs = Dict(eq.lhs => eq.rhs for eq in observed(system)
                         if !(Symbolics.symtype(Symbolics.unwrap(eq.lhs)) <: AbstractArray))
    fixed_params = Dict(Num(p) => Float64(parameter_values[Num(p)]) for p in parameters(system)
                        if !isequal(Num(p), ω) && haskey(parameter_values, Num(p)))

    inputs = vcat(input_syms, Symbolics.unwrap(ω))
    compile(expr) = Symbolics.build_function(
        Symbolics.unwrap(Symbolics.substitute(substitute_to_fixpoint(expr, observed_subs), fixed_params)),
        inputs; expression = Val{false})
    return compile(expression[1]), compile(expression[2])
end

#  Observed-variable reconstruction

# Rebuild the harmonic coefficient of an observed (non-state) variable: flatten its
# defining equation down to the original states, substitute each state's harmonic ansatz,
# then read off the requested cos/sin coefficient.
function reconstruct_from_observed(h_sys::HarmonicSystem, var::Num, order::Int, component::Symbol)
    hs = h_sys
    t  = ModelingToolkit.get_iv(hs.time_domain_system)
    observed_map = Dict{Any, Any}(Symbolics.unwrap(eq.lhs) => Symbolics.unwrap(eq.rhs)
                                  for sys in (hs.time_domain_system, hs.system) for eq in observed(sys))

    var_symbol = Symbolics.unwrap(var)
    target_rhs = get(observed_map, var_symbol, nothing)
    target_rhs === nothing && error("Variable $var not found in observed equations")

    flat_subs = Dict(k => v for (k, v) in observed_map if !(Symbolics.symtype(k) <: AbstractArray))
    expr = substitute_to_fixpoint(target_rhs, flat_subs; transform = Symbolics.expand_derivatives)
    (; states, diffvars, diff2vars) = hs.full_eqs   # cached 
    ansatz_subs = Dict{Any, Any}()
    for (k, s) in enumerate(states)
        ansatz_subs[Symbolics.unwrap(s)] = Symbolics.unwrap(hs.harmonic_ansatz[k])
    end
    for (i, dv) in enumerate(diff2vars)
        k = var_index(states, Symbolics.unwrap(diffvars[i]))
        ansatz_subs[Symbolics.unwrap(dv)] = Symbolics.unwrap(hs.harmonic_ansatz_dt[k])   # cached 
    end

    substituted = Num(Symbolics.expand(Symbolics.substitute(expr, ansatz_subs)))

    if order == 0
        return Num(Symbolics.substitute(Symbolics.unwrap(substituted), zero_harmonics(hs.N, hs.ω, t)))
    end

    basis = component == :Cos ? cos(Num(order) * hs.ω * t) : sin(Num(order) * hs.ω * t)
    return Num(Symbolics.coeff(substituted, basis))
end

reconstruct_from_observed(h_prob::HarmonicProblem, var::Num, order::Int, component::Symbol) =
    reconstruct_from_observed(h_prob.harmonic_system, var, order, component)
