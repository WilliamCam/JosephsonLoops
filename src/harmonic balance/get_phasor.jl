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

#  Linearised-problem indexing
#
# The linearised solve works in the coefficient ordering of the jacobians, NOT in
# unknowns(system) order: build_jacobians differentiates with respect to `vars`, which
# harmonic_equation assembles state by state as [DC, cos₁, sin₁, cos₂, sin₂, ...]. Rows of
# the rotated equations carry the same interleaved ordering (rotate_to_harmonic_frame), so
# the maps below address both the perturbation vector and the response array.

# Position k of each state in the assembly order of harmonic_equation. The k-th entry of
# harmonic_ansatz is the k-th state's ansatz, so a state's position is wherever its DC
# coefficient symbol appears.
function linearised_state_order(h_sys::HarmonicSystem)
    ansatz = h_sys.harmonic_ansatz
    order = Dict{String, Int}()
    for name in unique(key[1] for key in keys(h_sys.variable_map))
        dc = Symbolics.unwrap(h_sys.variable_map[(name, 0, :Cos)])
        k = findfirst(X -> any(v -> isequal(v, dc), Symbolics.get_variables(X)), ansatz)
        k === nothing && error("State $name not found in harmonic ansatz")
        order[name] = k
    end
    return order
end

# (state name, order, :Cos/:Sin) -> row in the linearised perturbation/response vectors.
function linearised_row_map(h_sys::HarmonicSystem)
    block = 2 * h_sys.N + 1
    rows = Dict{Tuple{String, Int, Symbol}, Int}()
    for (name, k) in linearised_state_order(h_sys)
        base = (k - 1) * block
        rows[(name, 0, :Cos)] = base + 1
        for n in 1:h_sys.N
            rows[(name, n, :Cos)] = base + 2n
            rows[(name, n, :Sin)] = base + 2n + 1
        end
    end
    return rows
end

# The coefficient symbols in linearised row order (the `vars` vector the jacobians were
# built against).
function linearised_vars(h_sys::HarmonicSystem)
    rows = linearised_row_map(h_sys)
    vars = Vector{Num}(undef, length(rows))
    for (key, idx) in rows
        vars[idx] = h_sys.variable_map[key]
    end
    return vars
end

# Small-signal injection vector: a kick of `amplitude` on the row addressed by
# (var_name, order, component). NB: rows are EQUATION-block rows — block k holds the k-th
# time-domain equation, which for multi-state systems is not necessarily var_name's own
# equation (states and equations are ordered independently by get_full_equations). For a
# physical source perturbation prefer source_perturbation_vector, which locates the right
# equation by the source parameter itself.
function perturbation_vector(h_sys::HarmonicSystem, var_name::String, order::Int = 1,
                             component::Symbol = :Sin; amplitude::Float64 = 1.0)
    rows = linearised_row_map(h_sys)
    haskey(rows, (var_name, order, component)) ||
        error("($var_name, $order, $component) is not a state harmonic. States: $(sort!(unique(k[1] for k in keys(rows))))")
    U = zeros(Float64, length(rows))
    U[rows[(var_name, order, component)]] = amplitude
    return U
end

# Injection vector for perturbing a sinusoidal source amplitude (e.g. P1₊Isrc₊I):
# U = -∂F/∂param · amplitude, where F are the rotated harmonic equations. The forcing
# profile ∂(residual)/∂param is read from each time-domain equation containing the
# parameter and projected onto the harmonic basis, so the equation rows, quadratures,
# signs AND equation scalings (e.g. the 2π/(Φ₀C) of a junction equation) are all
# determined automatically. `parameters` supplies numeric values for parameters that
# multiply the source term inside an equation. Perturbing a current source by δI then
# gives the small-signal response to an extra δI·(source waveform at Ω) drive, in
# source units.
function source_perturbation_vector(h_sys::HarmonicSystem, source_param::Num, parameters::Dict; amplitude::Float64 = 1.0)
    t = Num(ModelingToolkit.get_iv(h_sys.time_domain_system))
    eqs, _states = get_full_equations(h_sys.time_domain_system, t)
    N, ω = h_sys.N, h_sys.ω
    block = 2N + 1
    U = zeros(Float64, length(eqs) * block)

    # numeric parameter values, never the swept ω (it is part of the harmonic basis)
    fixed_params = Dict(k => v for (k, v) in parameters if !isequal(Num(k), ω))
    as_number(x) = Float64(Symbolics.value(Symbolics.simplify(x)))
    zero_harmonics = Dict{Any, Any}()
    for n in 1:(2N + 1)
        zero_harmonics[Symbolics.unwrap(cos(Num(n) * ω * t))] = 0.0
        zero_harmonics[Symbolics.unwrap(sin(Num(n) * ω * t))] = 0.0
    end

    found = false
    for (k, eq) in enumerate(eqs)
        profile = Symbolics.derivative(eq.lhs - eq.rhs, source_param)
        profile = Symbolics.expand(Symbolics.substitute(profile, fixed_params))
        isequal(Symbolics.simplify(profile), Num(0)) && continue
        found = true
        base = (k - 1) * block
        U[base + 1] -= amplitude * as_number(Symbolics.substitute(profile, zero_harmonics))
        for n in 1:N
            U[base + 2n]     -= amplitude * as_number(Symbolics.coeff(profile, cos(Num(n) * ω * t)))
            U[base + 2n + 1] -= amplitude * as_number(Symbolics.coeff(profile, sin(Num(n) * ω * t)))
        end
    end
    found || error("Parameter $source_param does not appear in any time-domain equation")
    return U
end

# Rotate a perturbation vector by 90° in the harmonic frame (cos ↔ sin per harmonic),
# giving the orthogonal signal quadrature. Used to probe the squeezed quadrature of a
# phase-sensitive (degenerate) parametric amplifier; pair with phase_preserving_s11.
function rotate_quadrature(h_sys::HarmonicSystem, U::AbstractVector)
    rows = linearised_row_map(h_sys)
    inv_rows = Dict(v => k for (k, v) in rows)
    U2 = zeros(eltype(U), length(U))
    for (row, val) in enumerate(U)
        val == 0 && continue
        name, order, comp = inv_rows[row]
        if comp == :Cos
            U2[rows[(name, order, :Sin)]] = val
        elseif comp == :Sin
            U2[rows[(name, order, :Cos)]] = val
        else                                   # DC term has no quadrature partner
            U2[row] = val
        end
    end
    return U2
end

# Phase-preserving (signal-to-signal) reflection magnitude |S_ss| from the port-voltage
# responses to two orthogonal quadrature current drives of amplitude `δI0`: `V_cos` to a
# cos (in-phase, X) drive and `V_sin` to a sin (quadrature, P) drive. A degenerate
# parametric amplifier is phase-SENSITIVE — a single real-quadrature injection measures the
# amplified (or squeezed) quadrature, not the phase-preserving gain that nodal HB codes
# (e.g. JosephsonCircuits.jl) report. The reflected field per unit drive is r_X = 2V_cos/(Z0
# δI0) − 1 and r_P = 2V_sin/(Z0 δI0) − i (the P-drive's incident wave is 90° rotated, hence
# −i). Their real/imag parts form a 2×2 field-transfer matrix whose singular values σ± are
# the amplified/squeezed quadrature gains; the phase-preserving magnitude is the mean,
# |S_ss| = (σ₊ + σ₋)/2 (idler |S_si| = (σ₊ − σ₋)/2). Passively σ₊ = σ₋ → ordinary |S11|.
function phase_preserving_s11(V_cos::AbstractVector, V_sin::AbstractVector, Z0::Real, δI0::Real)
    map(zip(V_cos, V_sin)) do (vc, vs)
        r_X = 2vc/(Z0*δI0) - 1
        r_P = 2vs/(Z0*δI0) - im
        σ = LinearAlgebra.svdvals([real(r_X) real(r_P); imag(r_X) imag(r_P)])
        (σ[1] + σ[2]) / 2
    end
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

    f_re, f_im = compile_phasor(expression, Symbolics.unwrap.(linearised_vars(h_sys)),
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
