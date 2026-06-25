# Extracts the numeric phasor Aₙ + iBₙ from a harmonic-balance result for some order in ω.
using ModelingToolkit

# Drop the `(t)` time argument from a symbol's printed name, e.g. "C1₊i(t)" -> "C1₊i".
strip_t(x) = replace(string(x), "(t)" => "")

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

function get_output(h_prob::HarmonicProblem, result::HarmonicResult, var_name::String, order::Int = 1)
    expression = get_harmonic_expression(h_prob.harmonic_system, var_name, order)
    return apply_harmonic_expression(h_prob, result, expression)
end

# Small-signal response phasor of `var_name` at the requested order. The harmonic system
# is passed in (rather than read off the problem) so the harmonic_system field can later
# be dropped from the problem structs. Restricted to state variables: an observed
# quantity's expression would have to be differentiated at the working point (δf = ∇f·δc),
# not evaluated at the perturbation — TODO.
function get_output(h_sys::HarmonicSystem, lin_prob::LinearisedProblem, result::HarmonicResult, var_name::String, order::Int = 1)
    haskey(h_sys.variable_map, (var_name, max(order, 1), :Cos)) ||
        error("Linearised get_output supports state variables only; $var_name is not a state. " *
              "States: $(sort!(unique(k[1] for k in keys(h_sys.variable_map))))")
    expression = get_harmonic_expression(h_sys, var_name, order)
    return apply_harmonic_expression(h_sys, lin_prob, result, expression)
end

#  Harmonic expressions

# Symbolic (cos, sin) Fourier coefficients of `var_name` at the requested order. A state's
# coefficients live in the variable_map; anything else (an observed quantity such as a
# branch current) is rebuilt from the observed equations.
function get_harmonic_expression(h_sys::HarmonicSystem, var_name::String, order::Int)
    variable_map = h_sys.variable_map
    if order == 0
        phasor_re = get(variable_map, (var_name, 0, :Cos), nothing)
        phasor_re === nothing && (phasor_re = reconstruct_from_observed(h_sys, var_name, 0, :Cos))
        return phasor_re, Num(0.0)
    end
    phasor_re = get(variable_map, (var_name, order, :Cos), nothing)
    phasor_im = get(variable_map, (var_name, order, :Sin), nothing)
    phasor_re === nothing && (phasor_re = reconstruct_from_observed(h_sys, var_name, order, :Cos))
    phasor_im === nothing && (phasor_im = reconstruct_from_observed(h_sys, var_name, order, :Sin))
    return phasor_re, phasor_im
end

get_harmonic_expression(h_prob::HarmonicProblem, var_name::String, order::Int) =
    get_harmonic_expression(h_prob.harmonic_system, var_name, order)

# Compile the symbolic phasor and evaluate it at every ω point. The expression is reduced
# to the system unknowns and ω, then mapped over the solved `[unknown × ω point]` array.
function apply_harmonic_expression(h_prob::HarmonicProblem, result::HarmonicResult, expression::Tuple{Num,Num})
    system = h_prob.harmonic_system.system
    isnothing(h_prob.parameter_sweep) || error("get_output supports ω-only sweeps")
    ω, ω_values = h_prob.ω_sweep
    ω_vec = ω_values isa Number ? [Float64(ω_values)] : ω_values

    #TODO: functionality for parameter sweeps
    solution = result.solution[ω]

    f_re, f_im = compile_phasor(expression, Symbolics.unwrap.(unknowns(system)),
                                h_prob.parameters, system, ω)

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

    # Rebuild the response-row ordering (= jacobian columns): per state in get_full_equations
    # order, the coefficients [DC, cos₁, sin₁, cos₂, sin₂, ...] pulled from variable_map. This
    # matches how harmonic_equation assembled `vars`, so the compiled phasor reads each
    # coefficient from the correct response row.
    tvar = Num(ModelingToolkit.get_iv(h_sys.time_domain_system))
    _, states = get_full_equations(h_sys.time_domain_system, tvar)
    vmap = h_sys.variable_map
    vars = Num[]
    for s in states
        name = strip_t(s)
        push!(vars, vmap[(name, 0, :Cos)])
        for n in 1:h_sys.N
            push!(vars, vmap[(name, n, :Cos)])
            push!(vars, vmap[(name, n, :Sin)])
        end
    end

    f_re, f_im = compile_phasor(expression, Symbolics.unwrap.(vars),
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
function reconstruct_from_observed(h_sys::HarmonicSystem, var_name::String, order::Int, component::Symbol)
    hs = h_sys
    t  = ModelingToolkit.get_iv(hs.time_domain_system)

    # Every observed equation (time-domain and harmonic), keyed by its unwrapped lhs. The
    # harmonic system is listed last so its definitions win on any shared key.
    observed_map = Dict{Any, Any}(Symbolics.unwrap(eq.lhs) => Symbolics.unwrap(eq.rhs)
                                  for sys in (hs.time_domain_system, hs.system) for eq in observed(sys))

    # Locate the observable's defining equation by exact (t-stripped) name.
    target = nothing
    for key in keys(observed_map)
        if strip_t(key) == var_name
            target = key
            break
        end
    end
    target === nothing && error("Variable $var_name not found in observed equations")

    # Flatten observed-in-observed down to states/params. Skip whole-array reconstructions.
    flat_subs = Dict(k => v for (k, v) in observed_map if !(Symbolics.symtype(k) <: AbstractArray))
    expr = substitute_to_fixpoint(observed_map[target], flat_subs; transform = Symbolics.expand_derivatives)

    # Substitute each original state by its harmonic ansatz (and derivatives). MTK lowers
    # D(state) either to a renamed `stateˍt(t)` symbol or a Differential wrapper.
    ansatz_subs = Dict{Any, Any}()
    for state_name in unique(key[1] for key in keys(hs.variable_map))
        A_coeffs = [hs.variable_map[(state_name, n, :Cos)] for n in 0:hs.N]
        B_coeffs = [hs.variable_map[(state_name, n, :Sin)] for n in 1:hs.N]
        X = harmonic_solution(hs.N, t, hs.ω, A_coeffs, B_coeffs)
        dX, d2X = get_derivatives(X, t)
        for v in Symbolics.get_variables(expr)
            name = strip_t(v)
            if name == state_name
                ansatz_subs[v] = Symbolics.unwrap(X)
            elseif name == state_name * "ˍt" || endswith(name, "₊" * state_name * "ˍt")
                ansatz_subs[v] = Symbolics.unwrap(dX)
            elseif name == state_name * "ˍtt" || endswith(name, "₊" * state_name * "ˍtt")
                ansatz_subs[v] = Symbolics.unwrap(d2X)
            elseif SymbolicUtils.iscall(v) && SymbolicUtils.operation(v) isa Differential &&
                   strip_t(SymbolicUtils.arguments(v)[1]) == state_name
                ansatz_subs[v] = Symbolics.unwrap(dX)
            end
        end
    end

    substituted = Num(Symbolics.expand(Symbolics.expand_derivatives(Symbolics.substitute(expr, ansatz_subs))))

    # Read off the requested harmonic. Symbolics.coeff returns the basis-function
    # coefficient; the DC term is what remains once every cosₖ/sinₖ is zeroed.
    if order == 0
        zero_harmonics = Dict{Any, Any}()
        for k in 1:(2*hs.N + 1)
            zero_harmonics[Symbolics.unwrap(cos(Num(k) * hs.ω * t))] = 0.0
            zero_harmonics[Symbolics.unwrap(sin(Num(k) * hs.ω * t))] = 0.0
        end
        return Num(Symbolics.substitute(Symbolics.unwrap(substituted), zero_harmonics))
    end

    basis = component == :Cos ? cos(Num(order) * hs.ω * t) : sin(Num(order) * hs.ω * t)
    return Num(Symbolics.coeff(substituted, basis))
end

reconstruct_from_observed(h_prob::HarmonicProblem, var_name::String, order::Int, component::Symbol) =
    reconstruct_from_observed(h_prob.harmonic_system, var_name, order, component)
