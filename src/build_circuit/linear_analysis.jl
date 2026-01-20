using Symbolics
using SymbolicUtils
using ModelingToolkit
using QuestBase
using DifferentialEquations
using Plots
include("../harmonic balance/colocation HB.jl")

#Linear analysis adapted from Kosata 2022 thesis

#set up example differential equation in time domain. We will use duffing oscillator example from p.59
@variables t x(t)
@parameters  α ω ω0 F γ η
diff_eq = Differential(t)(Differential(t)(x)) + ω0^2*x + α*x^3 + η*Differential(t)(x)*x^2 + γ*Differential(t)(x) - F*cos(ω*t) ~ 0

#number of harmonics to 1 just to test
Nharmonics = 1

#Collocation HB (this is different from Kosata 2022 / HarmonicBalance.jl)
harmonic_sys, harmonic_states = harmonic_equation(diff_eq, x, t, ω, Nharmonics)
harmonic_sys
@named ns = NonlinearSystem(harmonic_sys)
sys_nns = toggle_namespacing(ns, false)
#determine the system jacobian using MTK (Must be done before system tearing)
@time jac = calculate_jacobian(sys_nns)
#system tearing
sys =  mtkcompile(ns)
sys
#Solve HB system for frequencies (Large signal pump)
N = 100
ω_vec = range(0.9,1.2,N) #pump freqs
solution=[]
state_storage = []
u0 = zeros(2*Nharmonics+1)
for i in 1:1:N
    ps = [α => 1.0, ω0 => 1.0, F => 0.01, η => 0.1, ω=> ω_vec[i], γ=>1.0e-3]
    prob = NonlinearProblem(sys,u0, ps)
    sol = solve(prob)
    # update u0 for continuation, for this example this gives us the high amplitude branch
    #u0 = sol.u
    push!(solution,  sol[ns.A[1]] + sqrt(sol[ns.A[2]]^2+sol[ns.B[1]]^2))
    push!(state_storage,  sol.u)
end
plot(ω_vec ,solution)

#we will linearize around a single pump tone
j = 26
ωp = ω_vec[j]
ps = [α => 1.0, ω0 => 1.0, F => 0.01, η => 0.1, ω=> ωp, γ=>1.0e-3]

# can either solve the problem again for U₀ at pump tone, or re-use solution from previous loop
# prob = NonlinearProblem(sys,zeros(3), ps)
# sol = solve(prob)
#working_point = sol.u
#vars=unknowns(sys_nns)
#_sub_rules = Dict(vars.=>sol[vars])

#set the working point U₀
working_point = state_storage[j]
vars=unknowns(sys)
# we do this because MTK does not preserve the ordering of varaibles i.e. [A[1], A[2], B[1], etc]
_sub_rules = Dict(vars.=>working_point)

#J₁ is approximated as δJ/δω
jac_1 = Symbolics.derivative.(jac, ω)

N=500
#small signal frequency vector
Ω = range(0.9,1.2,N)

#update jacobian with numeric values
sub_rules=merge(Dict(ps), _sub_rules, Dict(ω=>ωp))
J₀ = substitute(jac, (sub_rules))
J₁ = substitute(jac_1, sub_rules)

#create a amplitude perturbation for A[2] and B[1] i.e. a perturbation at Ω
perturb = [1e-3; 0.0; 1e-3]

#Linear analysis solve loop
out = []
for i in 1:1:N
    # this matrix is precicesly eq. (5.12) of Kosata 2022 simplified for a single harmonic (N=1)
    mat = J₀ + (Ω[i] - ωp)*J₁
    #Julia linear algebra matrix solver M = A \ B
    resp = mat \ perturb
    #push!(out, 2*(resp[1]^2+resp[3]^2))
    u,v = resp[1], resp[3]
    push!(out, u^2 + v^2)
end
#plot our solution
plot(Ω, abs.(out))


