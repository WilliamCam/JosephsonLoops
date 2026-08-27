using Symbolics
using SymbolicUtils

function sample_collocation_grid!(residuals, Nt, res_expr, ω1, tvar)
    for n in 0:(Nt-1)
        phi_n = 2*pi*n/Nt
        res_at_coll = substitute(res_expr, Dict(tvar => phi_n / ω1))
        push!(residuals, res_at_coll ~ 0)
    end
end

function harmonic_solution(fourier_basis::FourierBasis, ω1::Num, ω2::Num, t::Num)
    indices = fourier_basis.fourier_indicies
    X = fourier_basis.dc_coeff
    cos_coeff, sin_coeff = fourier_basis.sin_coeffs, fourier_basis.cos_coeffs
    for n in 2:length(indices)
        pump_index, signal_index = indices[n]
        X += cos_coeff[n-1] * cos((pump_index * ω1 + signal_index * ω2)*t)
        X += sin_coeff[n-1] * sin((pump_index * ω1 + signal_index * ω2)*t)
    end
    return X
end

function harmonic_solution_symbolic_derrivative(fourier_basis::FourierBasis, ω1::Num, ω2::Num, tvar::Num)
    indices = fourier_basis.fourier_indicies
    N = length(indices)-1
    # must mirror harmonic_solution's quadrature convention (sin_coeffs on cos,
    # cos_coeffs on sin): the jacobian substitutes X and its derivatives into the same
    # equations, so both builders have to agree on which array carries which quadrature.
    cos_fourier, sin_fourier = fourier_basis.sin_coeffs, fourier_basis.cos_coeffs
    dcos_fourier, dsin_fourier = fourier_basis.d_sin_coeffs, fourier_basis.d_cos_coeffs
    d2cos_fourier, d2sin_fourier = fourier_basis.d2_sin_coeffs, fourier_basis.d2_cos_coeffs
    #TODO: Remove derrivative variable of DC term
    dc_of_t = (@variables dc(tvar))[1]
    d_dc_of_t = (@variables d_dc(t))[1]
    d2_dc_of_t = (@variables d2_dc(tvar))[1]
    cos_of_t = (@variables a(tvar)[1:N])[1]
    sin_of_t = (@variables b(tvar)[1:N])[1]
    d_cos_of_t = (@variables da(tvar)[1:N])[1]
    d_sin_of_t = (@variables db(tvar)[1:N])[1]
    d2_cos_of_t = (@variables d2a(tvar)[1:N])[1]
    d2_sin_of_t = (@variables d2b(tvar)[1:N])[1]
    lose_time_dependence = merge(
        Dict([dc_of_t => fourier_basis.dc_coeff, d_dc_of_t => fourier_basis.d_dc_coeff,
              d2_dc_of_t => fourier_basis.d2_dc_coeff]),
        Dict([var_t => sym for (var_t,sym) in zip(cos_of_t, cos_fourier)]),
        Dict([var_t => sym for (var_t,sym) in zip(sin_of_t, sin_fourier)]),
        Dict([var_t => sym for (var_t,sym) in zip(d_cos_of_t, dcos_fourier)]),
        Dict([var_t => sym for (var_t,sym) in zip(d_sin_of_t, dsin_fourier)]),
        Dict([var_t => sym for (var_t,sym) in zip(d2_cos_of_t, d2cos_fourier)]),
        Dict([var_t => sym for (var_t,sym) in zip(d2_sin_of_t, d2sin_fourier)])
    )
    #Retrieves symbolic expressions and collects all Num's for determiniation of jacobian
    X = dc_of_t # Start with the DC term A₁
    dX = d_dc_of_t
    d2X = d2_dc_of_t
    for n in 1:N
        pump_index, signal_index = indices[n+1] # ignore DC index
        @assert (pump_index, signal_index) !== (0,0) "Attempting to create harmonic from DC index"
        x_harmonic = cos_of_t[n] * cos((pump_index * ω1 + signal_index * ω2)*tvar)
        x_harmonic += sin_of_t[n] * sin((pump_index * ω1 + signal_index * ω2)*tvar)
        X += x_harmonic
        #dXdt
        dX_harmonic = expand_derivatives(D(x_harmonic))
        dX_harmonic_of_t = substitute(dX_harmonic, 
            Dict([
                D(cos_of_t[n])=>d_cos_of_t[n], 
                D(sin_of_t[n])=>d_sin_of_t[n]
                ])
            )
        dX += dX_harmonic_of_t
        d2X_harmonic_of_t = substitute(expand_derivatives(D(dX_harmonic_of_t)),
            Dict([
                D(cos_of_t[n])=>d_cos_of_t[n],
                D(sin_of_t[n])=>d_sin_of_t[n],
                D(d_cos_of_t[n])=>d2_cos_of_t[n],
                D(d_sin_of_t[n])=>d2_sin_of_t[n]
                ])
            )
        d2X += d2X_harmonic_of_t
    end
    @assert has_derivative(d2X)==false "Differential operator D(t) cannot be in $d2X -> replace with Symbolics variable e.g. 'dA[1]'"
    dX = substitute(dX, lose_time_dependence)
    d2X = substitute(d2X, lose_time_dependence)
    return dX, d2X
end

function harmonic_equation(eqs::Vector{Equation}, states::Vector{Num}, tvar::Num, ω::Union{Tuple{Num,Num}}, N::Int;
        jac=false, intermod_order = 0, commensurate = nothing)

    M = length(states)
    if M==1
        eqs = [eqs]
    end
    @assert (M == length(eqs)) "System does not have the same number of equations as state variables"
    @assert length(ω) <= 2 "maximum of two tones supported"

    ω1, ω2 = ω
    # isequal, not ==: Num == Num is symbolic when ω2 is a real frequency symbol
    single_tone = isequal(ω2, Num(0))

    basis_map = construct_fourier_basis(N, states, intermod_order = intermod_order, single_tone = single_tone)
    N_terms = length(basis_map[states[1]].fourier_indicies)

    if !single_tone
        @assert commensurate !== nothing "two-tone HB needs commensurate = (p, q) with ω1 = p·ω0, ω2 = q·ω0 (integers)"
        p, q = commensurate
        n_ints = [m * p + n * q for (m, n) in basis_map[states[1]].fourier_indicies[2:end]]
       
        Nt_two_tone = 0
        for cand in (2 * N_terms - 1):(2 * maximum(abs.(n_ints)) + 1)
            folded = [mod.(n_ints, cand); mod.(.-n_ints, cand)]  # where each slot's ± harmonics land
            if !any(iszero, folded) && allunique(folded)
                Nt_two_tone = cand
                break
            end
        end
        @assert Nt_two_tone > 0 "mixing products collide at every grid size — change (p, q) or the truncation"
    end
    X = Num[]
    residuals = Equation[]
    harmonic_eqs = eqs

    if jac
        d_residuals = Equation[]
        dvars = Num[]
        d2vars = Num[]
        vars = Num[]
        jac_subs = Dict{Num, Num}()
    end
    
    for k in 1:M
        cur_state = states[k]
        cur_fourier_basis = basis_map[cur_state]

        #ansatz
        harmonic_state = harmonic_solution(cur_fourier_basis, ω1, ω2, tvar)
        push!(X, harmonic_state)

        #derivatives
        if jac
            dX, d2X = harmonic_solution_symbolic_derrivative(cur_fourier_basis, ω1, ω2, tvar)
            jac_subs[Differential(tvar)(Differential(tvar)(states[k]))] = d2X
            jac_subs[Differential(tvar)(states[k])]                     = dX
            jac_subs[states[k]]                                         = harmonic_state
            push!(vars, cur_fourier_basis.dc_coeff)
            push!(dvars, cur_fourier_basis.d_dc_coeff)
            push!(d2vars, cur_fourier_basis.d2_dc_coeff)
            for n in 1:N_terms-1
                push!(vars, cur_fourier_basis.cos_coeffs[n], cur_fourier_basis.sin_coeffs[n])
                push!(dvars, cur_fourier_basis.d_cos_coeffs[n], cur_fourier_basis.d_sin_coeffs[n])
                push!(d2vars, cur_fourier_basis.d2_cos_coeffs[n], cur_fourier_basis.d2_sin_coeffs[n])
            end
        end
        dXdt, d2Xdt2 = get_derivatives(harmonic_state, tvar)
        harmonic_eqs = substitute(harmonic_eqs, Dict(Differential(tvar)(Differential(tvar)(cur_state))=>d2Xdt2))
        harmonic_eqs = substitute(harmonic_eqs, Dict(Differential(tvar)(cur_state)=>dXdt))
        harmonic_eqs = substitute(harmonic_eqs, Dict(cur_state=>harmonic_state))
    end
    if jac
        d_harmonic_eqs = substitute(eqs, jac_subs)
    end
    Nt = single_tone ? 2*N_terms-1 : Nt_two_tone
    ω_grid = single_tone ? ω1 : (1 // p) * ω1
    tone_pin = single_tone ? Dict{Num,Num}() : Dict(ω2 => (q // p) * ω1)
    for k in 1:M
        res_expr = harmonic_eqs[k].lhs - harmonic_eqs[k].rhs
        res_expr = single_tone ? Symbolics.simplify(res_expr) : substitute(res_expr, tone_pin)
        sample_collocation_grid!(residuals, Nt, res_expr, ω_grid, tvar)
        if jac
            d_res_expr = d_harmonic_eqs[k].lhs - d_harmonic_eqs[k].rhs
            d_res_expr = single_tone ? Symbolics.simplify(d_res_expr) : substitute(d_res_expr, tone_pin)
            sample_collocation_grid!(d_residuals, Nt, d_res_expr, ω_grid, tvar)
        end
    end
    if single_tone
        @named sys = NonlinearSystem(residuals)
    else
        projected = rotate_to_harmonic_frame(M, n_ints, Nt, residuals)
        @named sys = NonlinearSystem([expr ~ 0 for expr in projected])
    end
    if jac
        rotated_system = single_tone ? rotate_to_harmonic_frame(M, N, Nt, d_residuals) :
                                       rotate_to_harmonic_frame(M, n_ints, Nt, d_residuals)
        J0, J1, J2 = build_jacobians(rotated_system, vars, dvars, d2vars)
        return sys, X, basis_map, (J0, J1, J2)
    else
        return sys, X, basis_map, nothing
    end
end
