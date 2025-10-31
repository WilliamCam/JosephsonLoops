# In this example we model a DC SQUID driven by a current source

#Load JLoop package
using JosephsonLoops 
const jls = JosephsonLoops

### Write Netlist
#Define componeents in each loop
loops = [
["J1"],
["L2","C1"],
["C1","R1"],
["R1","I1"]]

coupling = [(1,2)]

circuit = jls.process_netlist(loops, mutalInd = coupling)

#cirucit model ODAE system and initial condition vector are created.
model, u0 = jls.build_circuit(circuit)

# we set the values of circuit parameters, for any parameters not specified; default values will be assigned.

I₀ = 1.0e-6
R₀ = 5.0
Φ₀ = jls.Φ₀

βc  = 2*pi/Φ₀ * I₀ * R₀^2
βL = 2*pi/Φ₀ * I₀
fd = 100e6

ps = [
    jls.loop4.sys.ω => 2*pi*fd
    jls.loop4.sys.I => 1.0*I₀
    jls.J1.sys.I0 => I₀
    jls.J1.sys.R => R₀
    jls.J1.sys.C => 0.01/βc
    jls.J1.sys.L => 2.0/βL
    jls.R1.sys.R => 50.0
    jls.C1.sys.C => 2.0/βc
    jls.L2.sys.L => 100.0/βL
    jls.M12.sys.L => 8.0/βL
    jls.loop1.sys.Φₑ => 0.5*Φ₀

]

tspan = (0.0, 1e-6)
saveat = LinRange(tspan[2]/10.0, tspan[2], 10000)
sol = jls.tsolve(model, u0, tspan, ps; saveat = saveat)
jls.tplot(sol, jls.R1, units="volts")

## Parameter Sweeps
using Plots
Φspan = (0.0, 2.0*Φ₀)
ensemble_sol = jls.ensemble_parameter_sweep(model, u0, tspan, Φspan, ps, jls.loop1.sys.Φₑ, jls.R1, saveat = saveat)
plot(ensemble_sol.u)


Ispan = (0.0, 2*I₀)
ensemble_sol = jls.ensemble_parameter_sweep(model, u0, tspan, Ispan, ps, jls.loop4.sys.I, jls.R1, saveat = saveat)
plot(ensemble_sol.u)


### Harmonic Balance ######
using ModelingToolkit

eqs, states = jls.get_full_equations(model, jls.t)

Nharmonics = 3
harmonic_sys, harmonic_states = jls.harmonic_equation(eqs, states, jls.t, jls.loop4.sys.ω, 3)


@named ns = NonlinearSystem(harmonic_sys)

sys = structural_simplify(ns)

N = 300
ω_vec = range(0.8,1.2,N)*2*pi*fd
solution=[]

for i in 1:1:N
    ps = [
    jls.loop4.sys.ω => ω_vec[i]
    jls.loop4.sys.I => 1.0*I₀
    jls.J1.sys.I0 => I₀
    jls.J1.sys.R => R₀
    jls.J1.sys.C => 0.01/βc
    jls.J1.sys.L => 2.0/βL
    jls.R1.sys.R => 50.0
    jls.C1.sys.C => 2.0/βc
    jls.L2.sys.L => 100.0/βL
    jls.M12.sys.L => 8.0/βL
    jls.loop1.sys.Φₑ => 0.5*Φ₀]
    prob = NonlinearProblem(sys,zeros(2*Nharmonics+1), ps)
    sol = solve(prob)
    push!(solution,  sol[ns.E[1]] + sqrt(sol[ns.E[3]]^2+sol[ns.F[1]]^2))
end

using Plots
plot(ω_vec, solution)

bif_par = jls.loop4.sys.ω
p_start = [
    jls.loop4.sys.ω => ω_vec[1]
    jls.loop4.sys.I => 1.0*I₀
    jls.J1.sys.I0 => I₀
    jls.J1.sys.R => R₀
    jls.J1.sys.C => 0.01/βc
    jls.J1.sys.L => 2.0/βL
    jls.R1.sys.R => 50.0
    jls.C1.sys.C => 2.0/βc
    jls.L2.sys.L => 100.0/βL
    jls.M12.sys.L => 8.0/βL
    jls.loop1.sys.Φₑ => 0.5*Φ₀]
u0_guess = Dict(unknowns(sys) .=> 0.0)
plot_var = ns.E[1] + sqrt(ns.E[3]^2+ns.F[1]^2)


