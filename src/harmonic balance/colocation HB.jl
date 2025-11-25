using Symbolics
using ModelingToolkit
using SymbolicUtils

function get_full_equations(model::ModelingToolkit.System, tvar::Num)
    eqs    = full_equations(model)
    states = unknowns(model)

    function var_is_in(vars::Vector, target_var::SymbolicUtils.BasicSymbolic{Real})
        ret = false
        for var in vars
            if isequal(var, target_var)
                ret = true
                break
            end
        end
        return ret
    end

    diff2vars   = Vector{Num}()   # “2nd derivative” helper vars
    diffvars    = Vector{Num}()   # corresponding 1st-derivative vars
    remove_idxs = Int[]

    # Find MTK-generated helper equations like z(t) ~ xˍt(t) etc.
    for (i, eq) in enumerate(eqs)
        vars = get_variables(eq.rhs)
        if length(vars) == 1 && var_is_in(states, vars[1])
            push!(diff2vars, vars[1])
            push!(diffvars, get_variables(eq.lhs)[1])
            push!(remove_idxs, i)
        end
    end

    # Remove those helper equations from the equation list
    for i in reverse(remove_idxs)
        deleteat!(eqs, i)
    end

    # Replace helper vars with actual derivatives
    for (i, var) in enumerate(diffvars)
        eqs = substitute(eqs, Dict(diff2vars[i] => Differential(tvar)(diffvars[i])))
    end

    # Remove helper vars from the state list
    empty!(remove_idxs)
    for (i, var) in enumerate(states)
        if var_is_in(diff2vars, var)
            push!(remove_idxs, i)
        end
    end
    for i in reverse(remove_idxs)
        deleteat!(states, i)
    end

    return eqs, states
end


function harmonic_solution(N, tvar, wvar, Afourier, Bfourier)
    X = Afourier[1]  # constant term A₀
    for n in 1:N
        X += Afourier[n + 1] * cos(n * wvar * tvar) +
             Bfourier[n]     * sin(n * wvar * tvar)
    end
    return X
end

function get_derivatives(X, t)
    D      = Differential(t)
    dXdt   = Symbolics.expand_derivatives(D(X))
    d2Xdt2 = Symbolics.expand_derivatives(D(dXdt))
    return dXdt, d2Xdt2
end

function harmonic_equation(eqs, states, tvar, wvar, N; dc_values = nothing)
    # normalise eqs to a vector
    if eqs isa Equation
        eqs = [eqs]
    end

    M = length(states)

    if M != length(eqs)
        error("System does not have the same number of equations as state variables.")
    end
    if M > 26
        error("System of equations is too large for single-letter coefficient labels (M > 26).")
    end

    coeff_labels    = 'A':'Z'
    X               = Num[]
    harmonic_system = Equation[]
    # Working copy so substitutions accumulate
    harmonic_eqs    = copy(eqs)

    #  1. Build Fourier series and substitute into ODEs 
    for k in 1:M
        cos_label = coeff_labels[2k - 1]  # e.g. 'A'
        sin_label = coeff_labels[2k]      # e.g. 'B'

        # Decide if DC is fixed for this state
        fix_dc = false
        dc_val = nothing
        if dc_values !== nothing
            dc_val = dc_values[k]
            if dc_val !== nothing
                fix_dc = true
            end
        end

        # --- Cosine coefficients (A_n) as scalars ---
        cos_coeffs = Num[]
        if fix_dc
            # A1..AN unknown, A0 fixed to dc_val
            for n in 1:N
                var_sym = Symbol(string(cos_label), "_", n)  # e.g. :A_1, :A_2, ...
                v, = @variables $(var_sym)
                push!(cos_coeffs, v)
            end
            a0 = dc_val
        else
            # A0..AN unknown
            for n in 0:N
                var_sym = Symbol(string(cos_label), "_", n)  # e.g. :A_0, :A_1, ...
                v, = @variables $(var_sym)
                push!(cos_coeffs, v)
            end
            a0 = cos_coeffs[1]  # A0
        end

        # --- Sine coefficients (B_n) as scalars ---
        sin_coeffs = Num[]
        for n in 1:N
            var_sym = Symbol(string(sin_label), "_", n)  # e.g. :B_1, :B_2, ...
            v, = @variables $(var_sym)
            push!(sin_coeffs, v)
        end

        # --- Build harmonic series for state k ---
        harmonic_state = a0
        for n in 1:N
            # For fix_dc: cos_coeffs[n] = A_n (A1..AN)
            # For !fix_dc: cos_coeffs[n+1] = A_n (A1..AN, since index1 is A0)
            ccoeff = fix_dc ? cos_coeffs[n] : cos_coeffs[n + 1]
            scoeff = sin_coeffs[n]
            harmonic_state += ccoeff * cos(n * wvar * tvar) +
                              scoeff * sin(n * wvar * tvar)
        end

        push!(X, harmonic_state)

        # Compute derivatives and substitute into the ODEs
        dXdt, d2Xdt2 = get_derivatives(harmonic_state, tvar)
        harmonic_eqs = substitute(
            harmonic_eqs,
            Dict(
                Differential(tvar)(Differential(tvar)(states[k])) => d2Xdt2,
                Differential(tvar)(states[k])                     => dXdt,
                states[k]                                         => harmonic_state,
            ),
        )
    end

    #  2. Collocation residuals
       Nt_total = 2N + 1  # full set of collocation points

    for k in 1:M
        res_expr = Symbolics.simplify(harmonic_eqs[k].lhs - harmonic_eqs[k].rhs)

        # Check if DC fixed for this state
        fix_dc_k = false
        if dc_values !== nothing
            if dc_values[k] !== nothing
                fix_dc_k = true
            end
        end

        # If DC fixed: skip n=0 collocation (that would duplicate DC info)
        if fix_dc_k
            collocation_indices = 1:(Nt_total - 1)   # 2N points
        else
            collocation_indices = 0:(Nt_total - 1)   # 2N+1 points
        end

        for n in collocation_indices
            phi_n = 2pi * n / Nt_total
            res_at_coll = substitute(res_expr, Dict(tvar => phi_n / wvar))
            push!(harmonic_system, res_at_coll ~ 0)
        end
    end

    return harmonic_system, X
end

function is_term(set, target_term)
    vars =
        if set isa Equation
            get_variables(set)
        elseif set isa SymbolicUtils.BasicSymbolic{Real}
            get_variables(set)
        elseif set isa Num
            get_variables(set)
        else
            set
        end
    ret = false
    for term in vars
        if isequal(term, target_term)
            ret = true
            break
        end
    end
    return ret
end
