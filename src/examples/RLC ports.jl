#In this example we model an parallel RLC circuit driven by a voltage source

#Load JosephsonLoops package
using JosephsonLoops 
const jls = JosephsonLoops

loops = [
["I1", "R1"],
["R1", "C1", "J1"]
]

ext_flux = [false, true]


circuit = jls.process_netlist(loops, ext_flux=ext_flux)

#cirucit model ODAE system and initial condition vector are created.
model, u0, guesses = jls.build_circuit(circuit)


# we set the values of circuit parameters, for any parameters not specified; default values will be assigned.
    ps = [
        jls.I1.ω => 100e6*2*pi
        jls.I1.I => 1e-12
        jls.R1.R => 50.0
        jls.C1.C => 100.0e-15
        jls.J1.C=> 1000.0e-15
        jls.J1.I0 => 1e-6
        jls.J1.R => 1.0
        jls.Φₑ2.Φₑ => 0.0
    ]

#specify transient window for solver
tspan = (0.0, 1e-6)
#transient circuit analysis
sol = jls.tsolve(model,u0, ps, tspan, guesses=guesses)
jls.tplot(sol, jls.C1, units = "amps")

# we can pass any arguments known to the problem interface of DifferentialEquations.jl, to achieve better results
#saveat force the integrator to step at certain time points
saveat = LinRange(0.0, 1e-6, 1000)

#specify a different solving algorithim
alg = Rodas5()

sol = jls.tsolve(model, u0, ps, tspan;
    guesses=guesses, saveat = saveat, solver_opts = alg, abstol = 1e-6)
jls.tplot(sol, jls.C1, units = "amps")

x = scsim.ensemble_fsolve(model, u0, tspan, (0.1, 10.0), ps, scsim.loop1, scsim.R1, units = "amps", Ntraj = 500)   

using Plots
plot(x.u)


#Harmonic Balance
include("../harmonic balance/colocation HB.jl")
eqs, states = jls.get_full_equations(model, jls.t)

harmonic_system, harmonic_states = jls.harmonic_equation(eqs, states, jls.t, jls.I1.ω, 3)

@named ns = NonlinearSystem(harmonic_system[1:end-1])
sys =  mtkcompile(ns)

I₀ = 1e-6
R₀ = 5000.0
Id = 0.05e-6
ωc = sqrt(2*pi *I₀/(jls.Φ₀*1000.0e-15))/(2*pi)
1/sqrt(1000.0e-15*jls.Φ₀/(2*pi*I₀))/(2*pi)


1/(2*pi*sqrt(1000.0e-12*1000.0e-15))

N = 100
ω_vec = 2*pi*(8.5:0.005:9.5)*1e9
solution1=Float64[]
solution2=Float64[]
u0_prev = zeros(12)
for i in 1:1:length(ω_vec)
    ps = [
        jls.I1.ω => ω_vec[i]
        jls.I1.I => 2.5e-9
        jls.C1.C => 100.0e-15
        jls.J1.C=> 1000.0e-15
        jls.J1.I0 => I₀
        jls.R1.R => 5000.0
        jls.J1.R => 5000.0
    ]
    prob = NonlinearProblem(sys,zeros(12), ps)
    sol = solve(prob)
    u0_prev = sol
    push!(solution1, sol[ns.A[1]]+ sqrt(sol[ns.A[2]]^2+sol[ns.B[1]]^2))
    #push!(solution2, sol[ns.C[1]]+ sqrt(sol[ns.C[2]]^2+sol[ns.D[1]]^2))
end
using Plots

Ii = @. Id  - jls.Φ₀/(2*pi)*ω_vec*solution1/50.0
Vi = @. jls.Φ₀/(2*pi)*ω_vec*solution1

ai = 0.5*(Vi+50*Ii)/50
bi = 0.5*(Vi-50*conj.(Ii))/50

plot(ω_vec/(2*pi), ((bi./ai)))
#using BifurcationKit

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
        jls.loop2.sys.Φₑ => 0.5
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