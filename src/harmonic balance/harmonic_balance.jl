using Symbolics, DifferentialEquations, HarmonicBalance, ModelingToolkit

@variables t x(t) # declare constant variables and a function x(t)
@parameters α ω0 F η ω 
diff_eq = d(x,t,2) + ω0^2*x + α*x^3 + η*d(x,t)*x^2 ~ F*cos(ω*t)

Nharmonics = 3

@variables A0 A1 A2 B0 B1 B2

X = A1*cos(ω*t) + A2*cos(2*ω*t) + B1*sin(ω*t) + B2*sin(2*ω*t)

dx = -A1*ω*sin(ω*t) -2*A2*ω*sin(2*ω*t) + B1*ω*cos(ω*t) + 2*B2*ω*cos(2*ω*t)
d2dx = -A1*ω^2*cos(ω*t) - 4*A2*ω^2*cos(2*ω*t) - B1*ω^2*sin(ω*t) - 4*B2*ω^2*sin(2*ω*t)

harmonic_eq = d2dx + ω0^2*X + α*X^3 + η*dx*X^2 - F*cos(ω*t) ~ 0

simplify(harmonic_eq)

using Symbolics

test= HarmonicBalance.trig_reduce(Symbolics.expand(harmonic_eq.lhs))
function is_term(set, target_term)
    if typeof(set) == Equation
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


NL_eq1 = sum(collect_terms(harmonic_eq, cos(ω*t)))
NL_eq2 = sum(collect_terms(harmonic_eq, sin(ω*t)))
NL_eq3 = sum(collect_terms(harmonic_eq, cos(2*ω*t)))
NL_eq4 = sum(collect_terms(harmonic_eq, sin(2*ω*t)))

eqs = [0~NL_eq1
       0~NL_eq2
       0~NL_eq3
       0~NL_eq4]

@named ns = NonlinearSystem(eqs)

guesses = [A1=>1.0, A2=>0.0, B1=>1.0, B2=>0.0]

ps = [α => 1., ω0 => 1.0, F => 0.01, η => 0.1, ω=>1.1]

sys = structural_simplify(ns)

ω_vec = range(0.9,1.2,100)
solution=[]

for i in 1:1:100
    ps = [α => 0.02, ω0 => 1.0, F => 0.01, η => 0.1, ω=> ω_vec[i]]
    prob = NonlinearProblem(sys, guesses, ps)
    sol = solve(prob)
    push!(solution, sqrt(sol.u[1]^2+sol.u[3]^2))
end

plot(solution)