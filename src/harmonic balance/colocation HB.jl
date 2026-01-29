using Symbolics
using SymbolicUtils
using ModelingToolkit
using QuestBase


function harmonic_solution(N, tvar, wvar, Afourier, Bfourier; get_derivatives=false, dAfourier=nothing, dBfourier=nothing)
    if get_derivatives
        @assert all([dAfourier, dBfourier] .!= nothing) "Need derivatives of harmonic variables e.g. 'dA[1]' "
        dX = dAfourier[1]
        d2X = 0
        vars = [Afourier[1]]
        dvars = [dAfourier[1]]
    end
    #TODO rename DC term to avoid annoying indexing i.e. A₀ + A₁Sin(ωt) + B₁Cos(ωt) would be
    X = Afourier[1]  # Start with the DC term A₁
    for n in 1:N
        X += Afourier[n + 1] * cos(n * wvar * tvar) + Bfourier[n] * sin(n * wvar * tvar)
        if get_derivatives
            #dXdt
            dX += (dAfourier[n+1] + n*wvar*Bfourier[n]) * cos(n * wvar * tvar) + (dBfourier[n] - n*wvar*Afourier[n+1]) * sin(n * wvar * tvar)
            #d2Xdt2
            _cos_comp = (-(n*wvar)^2 * Afourier[n+1] - 2*n*wvar*dBfourier[n]) * cos(n * wvar * tvar)
            _sin_comp = (-(n*wvar)^2 * Bfourier[n] + 2*n*wvar*dAfourier[n+1]) * sin(n * wvar * tvar)
            d2X += _cos_comp + _sin_comp
            push!(dvars, dAfourier[n+1], dBfourier[1])
            push!(vars, Afourier[n+1], Bfourier[n])
        end
    end
    if get_derivatives
        return X, dX, d2X, vars, dvars
    else
        return X
    end
end

function harmonic_equation(eqs, states, tvar, wvar, N; jac=false)
    M = length(states)
    if M==1
        eqs = [eqs]
    end
    @assert (M == length(eqs)) "System does not have the same number of equations as state variables"
    @assert M < 26 "System of equations is too large"
    coeff_labels = 'A':'Z'
    X = Num[]
    harmonic_system, d_harmonic_system = Equation[], Equation[]
    harmonic_eqs, d_harmonic_eqs = eqs, eqs
    if jac
        dvars = []
        vars = []
    end
    # loop over each state varibale for multiple harmonic equations
    for k in 1:M
        cos_coeff_labels, sin_coeff_labels = Symbol(coeff_labels[2*k-1]), Symbol(coeff_labels[2*k])
        cos_coeffs = @variables $cos_coeff_labels[1:N+1]
        sin_coeffs = @variables $sin_coeff_labels[1:N]
        harmonic_state = harmonic_solution(N, tvar, wvar, cos_coeffs[1], sin_coeffs[1])
        #Linearisation via jacobian requires introduction of derivative harmonic variables i.e. A[1]'
        if jac
            d_cos_coeff_labels, d_sin_coeff_labels = Symbol('d' * coeff_labels[2*k-1]), Symbol('d' * coeff_labels[2*k])
            d_cos_coeffs = @variables $d_cos_coeff_labels[1:N+1]
            d_sin_coeffs = @variables $d_sin_coeff_labels[1:N]

            _, _dXdt, _d2Xdt2, _vars, _dvars = harmonic_solution(N, tvar, wvar, cos_coeffs[1], sin_coeffs[1], 
                get_derivatives=true, dAfourier = d_cos_coeffs[1], dBfourier = d_sin_coeffs[1]
            )
      
            d_harmonic_eqs = substitute(harmonic_eqs, Dict(
                Differential(tvar)(Differential(tvar)(states[k])) => _d2Xdt2,
                Differential(tvar)(states[k]) => _dXdt,
                states[k] => harmonic_state)
            )
            append!(vars, _vars)
            append!(dvars,_dvars)
        end
        push!(X, harmonic_state)
        dXdt, d2Xdt2 = get_derivatives(harmonic_state, tvar)
        # if all(only_derivatives(eq, states[k], t) for eq in harmonic_eqs)
        #     problematic_var = states[k]
        #     print("Warning: harmonic variable mapped to first derivative i.e. D($problematic_var) = $harmonic_state")
        #     #New harmonic state will be first derivative X = D(var)
        #     harmonic_eqs = substitute(harmonic_eqs, Dict(Differential(tvar)(Differential(tvar)(states[k]))=>dXdt))
        #     harmonic_eqs = substitute(harmonic_eqs, Dict(Differential(tvar)(states[k])=>harmonic_state))
        # end
        harmonic_eqs = substitute(harmonic_eqs, Dict(Differential(tvar)(Differential(tvar)(states[k]))=>d2Xdt2))
        harmonic_eqs = substitute(harmonic_eqs, Dict(Differential(tvar)(states[k])=>dXdt))
        harmonic_eqs = substitute(harmonic_eqs, Dict(states[k]=>harmonic_state))
    end
    Nt = 2*N+1
    for k in 1:M
        res_expr = Symbolics.simplify(harmonic_eqs[k].lhs - harmonic_eqs[k].rhs)
        for n in 0:(Nt-1)
            phi_n = 2*pi*n/Nt
            res_at_coll = substitute(res_expr, Dict(tvar => phi_n / wvar))
            push!(harmonic_system, res_at_coll ~ 0)
        end
        if jac
            d_res_expr = Symbolics.simplify(d_harmonic_eqs[k].lhs - d_harmonic_eqs[k].rhs)
            for n in 0:(Nt-1)
                phi_n = 2*pi*n/Nt
                res_at_coll = substitute(d_res_expr, Dict(tvar => phi_n / wvar))
                push!(d_harmonic_system, res_at_coll ~ 0)
            end
        end
    end
    @named sys = NonlinearSystem(harmonic_system)
    if jac
        rotated_system = rotate_to_harmonic_frame(N, Nt, d_harmonic_system)
        #TODO: Check orderiing for M>1 larger systems
        J0, J1 = build_jacobians(rotated_system, vars, dvars) 
        return sys, X, (J0, J1);
    else
        return sys, X
    end
end

