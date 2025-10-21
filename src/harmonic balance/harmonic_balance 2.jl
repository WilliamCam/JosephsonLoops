using Symbolics, DifferentialEquations, HarmonicBalance, ModelingToolkit

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
     # loop over each state varibale for multiple harmonic equations
    for k in 1:M
        cos_coeff_labels, sin_coeff_labels = Symbol(coeff_labels[2*k-1]), Symbol(coeff_labels[2*k])
        cos_coeffs = @variables $cos_coeff_labels[1:N+1]
        sin_coeffs = @variables $sin_coeff_labels[1:N]
        harmonic_state = harmonic_solution(N, tvar, wvar, cos_coeffs[1], sin_coeffs[1])
        push!(X, harmonic_state)
        dXdt, d2Xdt2 = get_derivatives(harmonic_state, tvar)
        harmonic_eq = substitute(eqs[k], Dict(Differential(tvar)(Differential(tvar)(states[k]))=>d2Xdt2))
        harmonic_eq = substitute(harmonic_eq, Dict(Differential(tvar)(states[k])=>dXdt))
        harmonic_eq = substitute(harmonic_eq, Dict(states[k]=>harmonic_state))
            for n in 1:(N+1)
                if n==1
                    push!(harmonic_system, 0~sum(collect_static_terms(harmonic_eq, tvar)))
                else
                    push!(harmonic_system, 0~sum(collect_terms(harmonic_eq, cos((n-1)*wvar*tvar))))
                    push!(harmonic_system, 0~sum(collect_terms(harmonic_eq, sin((n-1)*wvar*tvar))))
                end
            end
    end
    return harmonic_system, X
end

function is_term(set, target_term)
    if typeof(set) == Equation
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

function collect_terms(eq::Equation, target_term::Num)
    terms = Num[]
    #set rhs = 0
    expr = HarmonicBalance.trig_reduce(Symbolics.expand(eq.lhs))
    num_to_sym = expr~0
    expr = Symbolics.arguments(num_to_sym.lhs)
    for term in expr
        if isequal(Symbolics.coeff(term, target_term), Num(0))
            continue
        else
            push!(terms, Num(Symbolics.coeff(term, target_term)))
        end
    end
    return terms
end

function collect_static_terms(eq::Equation, indvar::Num)
    terms=Num[]
    expr = HarmonicBalance.trig_reduce(Symbolics.expand(eq.lhs))
    num_to_sym = expr~0
    expr = Symbolics.arguments(num_to_sym.lhs)
    for term in expr
        if is_term(Num(term), indvar)
            continue
        else
            push!(terms,term)
        end
    end
    return terms
end

##### Example Usage

# @variables t x(t) # declare constant variables and a function x(t)
# @parameters  α ω ω0 F η 
# diff_eq = Differential(t)(Differential(t)(x)) + ω0^2*x + α*x^3 + η*Differential(t)(x)*x^2 - F*cos(ω*t) ~ 0

# Nharmonics = 3

# harmonic_sys, harmonic_states = harmonic_equation(diff_eq, x, t, ω, 3)

# @named ns = NonlinearSystem(harmonic_sys)

# sys = structural_simplify(ns)

# N = 300
# ω_vec = range(0.9,1.2,N)
# solution=[]

# for i in 1:1:N
#     ps = [α => 0.05, ω0 => 1.0, F => 0.01, η => 0.1, ω=> ω_vec[i]]
#     prob = NonlinearProblem(sys,zeros(2*Nharmonics+1), ps)
#     sol = solve(prob)
#     push!(solution,  sol[ns.A[1]] + sqrt(sol[ns.A[4]]^2+sol[ns.B[3]]^2))
# end

# plot(ω_vec, solution)

