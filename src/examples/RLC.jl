#In this example we model an parallel RLC circuit driven by a voltage source

#Load JosephsonLoops package
using JosephsonLoops 
const jls = JosephsonLoops

### Write Netlist
#Define componeents in each loop
loops = [
["I1", "R1"],
["R1", "C1", "J1"]
]

circuit = jls.process_netlist(loops)

#cirucit model ODAE system and initial condition vector are created.
model, u0 = jls.build_circuit(circuit)

# we set the values of circuit parameters, for any parameters not specified; default values will be assigned.
    ps = [
        jls.loop1.sys.ω => 100e6*2*pi
        jls.loop1.sys.I => 1e-12
        jls.C1.sys.C => 100.0e-15
        jls.J1.sys.L => 1000.0e-12
        jls.J1.sys.C=> 1000.0e-15
        jls.J1.sys.I0 => 1e-6
        jls.R1.sys.R => 50.0
        jls.J1.sys.R => 1.0
        jls.loop2.sys.Φₑ => 0.0
    ]
using DifferentialEquations

#specify transient window for solver
tspan = (0.0, 1e-6)
using DifferentialEquations
# Create dictionary of initial conditions
u0 = Dict(u => 0.0 for u in unknowns(model))
prob = DAEProblem(model, u0, tspan, ps)
#transient circuit analysis
sol = jls.tsolve(model, u0, tspan, ps, alg = Rodas5())
jls.tplot(sol, jls.R1, units = "amps")

# we can pass any arguments known to the problem interface of DifferentialEquations.jl, to achieve better results

#saveat force the integrator to step at certain time points
saveat = LinRange(10.0, 100.0, 1000)

#specify a different solving algorithim
alg = Rodas5()

@time sol = scsim.tsolve(model, u_initial, tspan, ps; saveat = saveat, alg = alg, abstol = 1e-6)
scsim.tplot(sol, scsim.R1, units = "amps")

x = scsim.ensemble_fsolve(model, u0, tspan, (0.1, 10.0), ps, scsim.loop1, scsim.R1, units = "amps", Ntraj = 500)   

using Plots
plot(x.u)


#Harmonic Balance

include("../harmonic balance/colocation HB.jl")

eqs, states = get_full_equations(model, jls.t)

harmonic_sys, harmonic_states = harmonic_equation(eqs, states, jls.t, jls.loop1.sys.ω, 3)

@named ns = NonlinearSystem(harmonic_sys)

sys = structural_simplify(ns)

I₀ = 1e-6
R₀ = 5.0
Id = 0.05e-6
ωc = sqrt(2*pi *I₀/(jls.Φ₀*1000.0e-15))/(2*pi)

1/(2*pi*sqrt(1000.0e-12*1000.0e-15))

N = 300
ω_vec = 2*pi*(3:0.01:6)*1e9
solution1=Float64[]
solution2=Float64[]
u0_prev = zeros(18)
for i in 1:1:length(ω_vec)
    ps = [
        jls.loop1.sys.ω => ω_vec[i]
        jls.loop1.sys.I => Id
        jls.C1.sys.C => 100.0e-15
        jls.J1.sys.L => 1000.0e-12
        jls.J1.sys.C=> 1000.0e-15
        jls.J1.sys.I0 => 0.3e-6
        jls.R1.sys.R => 50.0
        jls.J1.sys.R => 1000.0
        jls.loop2.sys.Φₑ => 0.0
    ]
    prob = NonlinearProblem(sys,u0_prev, ps)
    sol = solve(prob)
    u0_prev = sol
    push!(solution1, sqrt(sol[ns.A[2]]^2+sol[ns.B[1]]^2))
end
using Plots

Ii = @. Id  - jls.Φ₀/(2*pi)*ω_vec*solution1/50.0
Vi = @. jls.Φ₀/(2*pi)*ω_vec*solution1

ai = 0.5*(Vi+50*Ii)/50
bi = 0.5*(Vi-50*conj.(Ii))/50

plot(ω_vec/(2*pi), 10*log10.(abs2.(ai./bi)))

using BifurcationKit

bif_par = jls.loop1.sys.ω
p_start = [
        jls.loop1.sys.ω => ω_vec[1]
        jls.loop1.sys.I => Id
        jls.C1.sys.C => 100.0e-15
        jls.J1.sys.L => 1000.0e-12
        jls.J1.sys.C=> 1000.0e-15
        jls.J1.sys.I0 => 0.3e-6
        jls.R1.sys.R => 50.0
        jls.J1.sys.R => 1000.0
        jls.loop2.sys.Φₑ => 0.0
]
u0_guess = sol
plot_var = sys.A[2]

bprob = BifurcationProblem(sys,
    u0_guess,
    p_start,
    bif_par;
    plot_var = plot_var,
    jac = false)

p_span = (ω_vec[1], ω_vec[end])
opts_br = ContinuationPar(nev = 2,
    p_min = p_span[1],
    p_max = p_span[2])

bf = bifurcationdiagram(bprob, PALC(), 2)