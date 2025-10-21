using Symbolics, DifferentialEquations, HarmonicBalance, ModelingToolkit

@variables t x(t) # declare constant variables and a function x(t)
@parameters α ω0 F η ω 
diff_eq = d(x,t,2) + ω0^2*x + α*x^3 + η*d(x,t)*x^2 -F*cos(ω*t) ~ 0
D = Differential(t)

Nharmonics = 3

@variables A[1:Nharmonics+1] B[1:Nharmonics]

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

test = get_derivatives(X, t)

function harmonic_equation(eq, xvar, tvar, X, wvar, N)
    harmonic_system = Equation[]
    dXdt, d2Xdt2 = get_derivatives(X, tvar)
    harmonic_eq = substitute(eq, Dict(Differential(tvar)(Differential(tvar)(xvar))=>d2Xdt2))
    harmonic_eq = substitute(harmonic_eq, Dict(Differential(tvar)(xvar)=>dXdt))
    harmonic_eq = substitute(harmonic_eq, Dict(xvar=>X))

    for n in 1:(N+1)
        if n==1
            push!(harmonic_system, 0~sum(collect_static_terms(harmonic_eq, tvar)))
        else
            push!(harmonic_system, 0~sum(collect_terms(harmonic_eq, cos((n-1)*wvar*tvar))))
            push!(harmonic_system, 0~sum(collect_terms(harmonic_eq, sin((n-1)*wvar*tvar))))
        end
    end
    return harmonic_system
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

X = harmonic_solution(Nharmonics, t, ω, A, B)

harmonic_sys = harmonic_equation(diff_eq, x, t, X, ω, Nharmonics)


@named ns = NonlinearSystem(harmonic_sys)

guesses = [A[1]=>0.0
A[2]=>1.0
B[1]=>1.0]
for n in 2:Nharmonics
    push!(guesses, A[n+1]=>0.0)
    push!(guesses, B[n]=>0.0)
end


sys = structural_simplify(ns)

N = 100
ω_vec = range(0.9,1.2,N)
solution=[]

for i in 1:1:N
    ps = [α => 0.05, ω0 => 1.0, F => 0.01, η => 0.1, ω=> ω_vec[i]]
    prob = NonlinearProblem(sys, guesses, ps)
    sol = solve(prob)
    push!(solution,  sol[A[1]] + sqrt(sol[A[4]]^2+sol[B[3]]^2))
end

plot(ω_vec, solution)
