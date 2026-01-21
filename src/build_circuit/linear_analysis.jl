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
harmonic_sys, harmonic_states, hsys = harmonic_equation(diff_eq, x, t, ω, Nharmonics, jac=true)
@named ns = NonlinearSystem(harmonic_sys)
@named lin_sys = NonlinearSystem(hsys)
sys_nns = toggle_namespacing(lin_sys, false)
#system tearing
sys =  mtkcompile(ns)
sys
#Solve HB system for frequencies (Large signal pump)
N = 200
ω_vec = range(0.8,1.2,N) #pump freqs
solution=[]
state_storage = []
u0 = zeros(2*Nharmonics+1)
for i in 1:1:N
    ps = [α => 1.0, ω0 => 1.0, F => 0.01, η => 1.0e-1, ω=> ω_vec[i], γ=>1.0e-3]
    prob = NonlinearProblem(sys,u0, ps)
    sol = solve(prob)
    # update u0 for continuation, for this example this gives us the high amplitude branch
    u0 = sol.u
    push!(solution,  sol[ns.A[1]] + sqrt(sol[ns.A[2]]^2+sol[ns.B[1]]^2))
    push!(state_storage,  sol.u)
end
plot(ω_vec ,solution)

#we will linearize around a single pump tone
j = 150
ωp = ω_vec[j]
ps = [α => 1.0, ω0 => 1.0, F => 0.01, η => 1.0e-1, γ=>1.0e-3]

#set the working point U₀
working_point = state_storage[j]
vars=unknowns(sys)
# we do this because MTK does not preserve the ordering of varaibles i.e. [A[1], A[2], B[1], etc]
_sub_rules = Dict(vars.=>working_point)

#J₁ is approximated as δJ/δω
function build_gamma(H, N)
    # Use Rational{Int} for exact symbolic cancellation
    Γ = Matrix{Num}(undef, 2H + 1, N)
    for j in 1:N
        for n in 1:H
            phase = n * (j - 1) * (2π / N)
            # Use Rational for the lead fraction
            Γ[n, j] = Num((2//N) * cos(phase))
        end
        Γ[H+1, j] = Num(1//N) # Exact 1/3
        for n in 1:H
            phase = n * (j - 1) * (2π / N)
            Γ[H+1+n, j] = Num((2//N) * sin(phase))
        end
    end
    return Γ
end

Γ = build_gamma(1, 3) 

unknowns(lin_sys)
vars_dot = [unknowns(lin_sys)[4], unknowns(lin_sys)[6], unknowns(lin_sys)[1]]

Rfreq = Γ*[h.lhs for h in hsys]
Rfreq = simplify.(Rfreq)
_jac = Symbolics.jacobian(Rfreq, unknowns(ns))
jac = Num.(simplify(substitute(_jac, Dict(vars_dot .=> 0))))

_jac1 = Symbolics.jacobian(Rfreq, vars_dot)

N=800
#small signal frequency vector
Ω = range(0.8,1.2,N)

#update jacobian with numeric values

sub_rules=merge(Dict(ps), _sub_rules, Dict(ω=>ωp))
J₀ = substitute(jac, sub_rules)
J₁ = substitute(_jac1, sub_rules)

#create a amplitude perturbation for A[2] and B[1] i.e. a perturbation at Ω
perturb = zeros(2*Nharmonics+1)
perturb[1] = 1.0e-3
#Linear analysis solve loop
out = []
out2=[]
for i in 1:1:N
    Δ = Ω[i] - ωp
    # this matrix is precicesly eq. (5.12) of Kosata 2022 simplified for a single harmonic (N=1)
    mat = J₀ + 1im*(Ω[i] - ωp)*J₁
    #Julia linear algebra matrix solver M = A \ B
    resp = mat \ perturb
    u,v = resp[1], resp[3]

    signal_amp = abs((u + im*v) / 2)
    idler_amp = abs((conj(u) - im*conj(v)) / 2)
    push!(out, signal_amp)
    push!(out2, idler_amp)
end
#plot our solution
plot(Ω, abs.(out)) 
plot!(Ω, abs.(out2)) 

# Print the diagonal of your numeric J0 when F=0
println("J0 Diagonals: ", [J₀[1,1], J₀[2,2], J₀[3,3]])

# 1. Check the Coriolis coupling in J1
# For Variables [A2, A1, B1] (Cos, DC, Sin)
# J1[1, 3] should be -2*ωp and J1[3, 1] should be 2*ωp
println("Coriolis check: ", J₁[1, 3], " and ", J₁[3, 1])

# 2. Check the Mass diagonal in J1
# It should be 1.0 (or very close to it)
println("Mass diagonal: ", J₁[1, 1], " and ", J₁[3, 3])