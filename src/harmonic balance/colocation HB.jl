using Symbolics
using ModelingToolkit
using SymbolicUtils
using QuestBase

function get_full_equations(model::ModelingToolkit.System, tvar::Num)
    eqs = full_equations(model)
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

    diff2vars = Vector{Num}()
    diffvars = Vector{Num}()
    remove_idxs = Int[]
    for (i, eq) in enumerate(eqs)
        vars = get_variables(eq.rhs)
        if length(vars) == 1 && var_is_in(states, vars[1])
            push!(diff2vars, vars[1])
            push!(diffvars, get_variables(eq.lhs)[1])
            push!(remove_idxs, i)
        end
    end

    for i in reverse(remove_idxs)
        deleteat!(eqs, i)
    end

    for (i,var) in enumerate(diffvars)
        eqs = substitute(eqs, Dict(diff2vars[i]=>Differential(tvar)(diffvars[i])))
    end
    remove_idxs = Int[]
    for (i,var) in enumerate(states)
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
    X = Afourier[1]  # Start with the constant term A₀
    for n in 1:N
        X += Afourier[n + 1] * cos(n * wvar * tvar) + Bfourier[n] * sin(n * wvar * tvar)
    end
    return X
end

function get_derivatives(X, t)
    D = Differential(t)
    dXdt = Symbolics.expand_derivatives(D(X))
    d2Xdt2 = Symbolics.expand_derivatives(D(dXdt))
    return dXdt, d2Xdt2
end

function harmonic_equation(eqs, states, tvar, wvar, N)
    M = length(states)
    if M==1
        eqs = [eqs]
    end
    if M != length(eqs)
        print("System does not have the same number of equations as state variables")
        return
    end
    if M > 26
        # Add second labeling system for fourier coeffs e.g. Aa Ab Ac etc
        print("System of equations is too large...for now...")
        return
    end
    coeff_labels = 'A':'Z'
    X = Num[]
    harmonic_system = Equation[]
    harmonic_eqs = eqs
     # loop over each state varibale for multiple harmonic equations
    for k in 1:M
        cos_coeff_labels, sin_coeff_labels = Symbol(coeff_labels[2*k-1]), Symbol(coeff_labels[2*k])
        cos_coeffs = @variables $cos_coeff_labels[1:N+1]
        sin_coeffs = @variables $sin_coeff_labels[1:N]
        harmonic_state = harmonic_solution(N, tvar, wvar, cos_coeffs[1], sin_coeffs[1])
        push!(X, harmonic_state)
        dXdt, d2Xdt2 = get_derivatives(harmonic_state, tvar)
        harmonic_eqs = substitute(harmonic_eqs, Dict(Differential(tvar)(Differential(tvar)(states[k]))=>d2Xdt2))
        harmonic_eqs = substitute(harmonic_eqs, Dict(Differential(tvar)(states[k])=>dXdt))
        harmonic_eqs = substitute(harmonic_eqs, Dict(states[k]=>harmonic_state))
    end
    Nt = 2*N+1
    for k in 1:M
        res_expr = Symbolics.simplify(harmonic_eqs[k].lhs - harmonic_eqs[k].rhs)
        for n in 0:(Nt-1)
            # phase angle (numeric) -> phi_n = 2π*n/Nt
            phi_n = 2*pi*n/Nt
            # substitute tvar -> phi_n / wvar to evaluate residual at that phase
            # Note: this causes terms like cos(m*w*t) to become cos(m*phi_n), purely numeric
            res_at_coll = substitute(res_expr, Dict(tvar => phi_n / wvar))
            # Optionally simplify trig of numeric arguments
            #res_at_coll = QuestBase.trig_reduce(Symbolics.expand(res_at_coll))
            push!(harmonic_system, res_at_coll ~ 0)
        end
    end
    @mtkcompile ns = NonlinearSystem(harmonic_system)
    return ns, X
end

function is_term(set, target_term)
    if typeof(set) == Equation
        vars = get_variables(set)
    elseif typeof(set) == SymbolicUtils.BasicSymbolic{Real}
        vars = get_variables(set)
    elseif typeof(set) == Num
        vars = get_variables(set)
    else
        vars = set
    end
    ret = false
    for term in vars
        if isequal(term, target_term)
            ret = true
            break
        else
            ret = false
        end
    end
    return ret
end


#### Example Usage

# @variables t x(t) # declare constant variables and a function x(t)
# @parameters  α ω ω0 F η 
# diff_eq = Differential(t)(Differential(t)(x)) + ω0^2*x + α*x^3 + η*Differential(t)(x)*x^2 - F*cos(ω*t) ~ 0

# Nharmonics = 3

# harmonic_sys, harmonic_states = harmonic_equation(diff_eq, x, t, ω, 3)
# harmonic_sys
# @named ns = NonlinearSystem(harmonic_sys)
# sys = structural_simplify(ns)

# N = 300
# ω_vec = range(0.9,1.2,N)
# solution=[]

# for i in 1:1:N
#     ps = [α => 0.05, ω0 => 1.0, F => 10, η => 0.1, ω=> ω_vec[i]]
#     prob = NonlinearProblem(sys,zeros(2*Nharmonics+1), ps)
#     sol = solve(prob)
#     push!(solution,  sol[ns.A[1]] + sqrt(sol[ns.A[4]]^2+sol[ns.B[3]]^2))
# end

# plot(solution)

