using Symbolics
using SymbolicUtils

function harmonic_solution(N, tvar, wvar, fourier_vars)
    @assert length(fourier_vars)==2 "Sin and Cos Fourier coefficients"
    Afourier, Bfourier = fourier_vars
    #TODO rename DC term to avoid annoying indexing i.e. A₀ + A₁Sin(ωt) + B₁Cos(ωt) would be
    X = Afourier[1]  # Start with the DC term A₁
    for n in 1:N
        X += Afourier[n + 1] * cos(n * wvar * tvar) + Bfourier[n] * sin(n * wvar * tvar)
    end
    return X
end

function fourier_coefficient_derrivatives(N, tvar, wvar, fourier_vars, d_fourier_vars)
    @assert length(fourier_vars)==2 "Sin and Cos Fourier coefficients"
    Afourier, Bfourier = fourier_vars
    dAfourier, dBfourier = d_fourier_vars
    #TODO: Remove derrivative variable of DC term
    cos_variables = (@variables a(t)[1:N+1])[1]
    sin_variables = (@variables b(t)[1:N])[1]
    d_cos_variables = (@variables da(t)[1:N+1])[1]
    d_sin_variables = (@variables db(t)[1:N])[1]
    sub_dict = merge(
        Dict([var_t => sym for (var_t,sym) in zip(cos_variables, Afourier)]),
        Dict([var_t => sym for (var_t,sym) in zip(sin_variables, Bfourier)]),
        Dict([var_t => sym for (var_t,sym) in zip(d_cos_variables, dAfourier)]),
        Dict([var_t => sym for (var_t,sym) in zip(d_sin_variables, dBfourier)])
    )
    #Retrieves symbolic expressions and collects all Num's for determiniation of jacobian
    X = cos_variables[1]  # Start with the DC term A₁
    dX = d_cos_variables[1]
    d2X = 0
    #TODO rename DC term to avoid annoying indexing i.e. A₀ + A₁Sin(ωt) + B₁Cos(ωt) would be
    for n in 1:N
        x_harmonic = cos_variables[n + 1] * cos(n * wvar * tvar) + sin_variables[n] * sin(n * wvar * tvar)
        X += x_harmonic
        #dXdt
        dX_harmonic = substitute(expand_derivatives(
            D(x_harmonic)), Dict([D(cos_variables[n+1])=>d_cos_variables[n+1],D(sin_variables[n])=>d_sin_variables[n]])
            )
        dX += dX_harmonic
        #d2Xdt2, simplify the expression and remove slow second order terms A''/B''
        d2X_harmonic = simplify(expand_derivatives(D(dX_harmonic)), expand=true, rewriter=remove_slow_d2_rewriter)
        d2X += d2X_harmonic
    end
    @assert has_derivative(d2X)==false "Differential operator D(t) cannot be in $d2X -> replace with Symbolics variable e.g. 'dA[1]'"
    dX = substitute(dX, sub_dict)
    d2X = substitute(d2X, sub_dict)
    return dX, d2X
end

function harmonic_equation(eqs, states, tvar, wvar, N; jac=false)
    M = length(states)
    if M==1
        eqs = [eqs]
    end
    @assert (M == length(eqs)) "System does not have the same number of equations as state variables"
    @assert M < 13 "System of equations is too large"
    coeff_labels = 'A':'Z'
    X = Num[]
    X_dt = Num[]   
    harmonic_system, d_harmonic_system = Equation[], Equation[]
    harmonic_eqs, d_harmonic_eqs = eqs, eqs
    if jac
        dvars = Num[]
        vars = Num[]
        dyn_subs = Dict{Num, Num}()
    end
    # loop over each state varibale for multiple harmonic equations
    variable_map = Dict{Tuple{SymbolicUtils.BasicSymbolic{Real}, Int, Symbol}, Num}()

    for k in 1:M
        cos_labels, sin_labels = Symbol(coeff_labels[2*k-1]), Symbol(coeff_labels[2*k])
        cos_sym_arr = @variables $cos_labels[1:N+1]
        sin_sym_arr = @variables $sin_labels[1:N]
        cos_vars, sin_vars = append!(cos_sym_arr, sin_sym_arr)

        #populate varibable map
        state_symbol = Symbolics.unwrap(states[k])
        variable_map[(state_symbol, 0, :Cos)] = cos_vars[1]
        for n in 1:N
            variable_map[(state_symbol, n, :Cos)] = cos_vars[n+1]
            variable_map[(state_symbol, n, :Sin)] = sin_vars[n]
        end

        harmonic_state = harmonic_solution(N, tvar, wvar, [cos_vars, sin_vars])

        push!(X, harmonic_state)
        dXdt, d2Xdt2 = get_derivatives(harmonic_state, tvar)
        push!(X_dt, dXdt)

        #Linearisation via jacobian requires introduction of derivative harmonic variables i.e. A[1]'
        if jac
            d_cos_labels, d_sin_labels = Symbol('d' * coeff_labels[2*k-1]), Symbol('d' * coeff_labels[2*k])
            d_cos_sym_arr = @variables $d_cos_labels[1:N+1]
            d_sin_sym_arr = @variables $d_sin_labels[1:N]
            d_cos_vars, d_sin_vars = append!(d_cos_sym_arr, d_sin_sym_arr)
        
            dXdt, d2Xdt2 = fourier_coefficient_derrivatives(N, tvar, wvar, [cos_vars, sin_vars] ,[d_cos_vars,d_sin_vars])
            dyn_subs[Differential(tvar)(Differential(tvar)(states[k]))] = d2Xdt2
            dyn_subs[Differential(tvar)(states[k])]                     = dXdt
            dyn_subs[states[k]]                                         = harmonic_state
            #Ordering of these variables is important for Jacobian formulation
            push!(vars, cos_vars[1]) #DC term
            push!(dvars, d_cos_vars[1]) #DC term ?
            for n in 1:N
                push!(vars, cos_vars[n+1], sin_vars[n])
                push!(dvars, d_cos_vars[n+1], d_sin_vars[n])
            end
        end

        harmonic_eqs = substitute(harmonic_eqs, Dict(Differential(tvar)(Differential(tvar)(states[k]))=>d2Xdt2))
        harmonic_eqs = substitute(harmonic_eqs, Dict(Differential(tvar)(states[k])=>dXdt))
        harmonic_eqs = substitute(harmonic_eqs, Dict(states[k]=>harmonic_state))
    end
    if jac
        d_harmonic_eqs = substitute(eqs, dyn_subs)
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
        rotated_system = rotate_to_harmonic_frame(M, N, Nt, d_harmonic_system)
        J0, J1 = build_jacobians(rotated_system, vars, dvars)
        return sys, X, variable_map, (J0, J1), X_dt
    else
        return sys, X, variable_map, X_dt
    end
end

