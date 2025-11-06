using JosephsonLoops 
const jls = JosephsonLoops

loops = [
["J1","L1"],
["L2","C1"],
["C1","R1"],
["R1","I1"]]

coupling = [(1,2)]
ext_flux = [true, false, false, false]
circuit = jls.process_netlist(loops, mutual_coupling = coupling, ext_flux = ext_flux)

model = jls.build_circuit(circuit)  

I₀ = 1.0e-6
R₀ = 5.0
Φ₀ = jls.Φ₀

βc  = 2*pi/Φ₀ * I₀ * R₀^2
βL = 2*pi/Φ₀ * I₀
fdrive = 100e6

ps = [
    jls.I1.ω => 2*pi*fdrive
    jls.I1.I => 1.0*I₀
    jls.J1.I0 => I₀
    jls.J1.R => R₀
    jls.J1.C => 0.01/βc
    jls.R1.R => 50.0
    jls.C1.C => 2.0/βc
    jls.L1.L => 2.0/βc
    jls.L2.L => 100.0/βL
    jls.M12.L => 8.0/βL

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
harmonic_sys, harmonic_states = jls.harmonic_equation(eqs, states, jls.t, jls.I1.ω, 3)


@named ns = NonlinearSystem(harmonic_sys)

sys = structural_simplify(ns)

N = 300
ω_vec = range(0.8,1.2,N)*2*pi*fdrive
solution=[]

for i in 1:1:N
    ps = [
    jls.I1.ω => ω_vec[i]
    jls.I1.I => 0.6*I₀
    jls.J1.I0 => I₀
    jls.J1.R => R₀
    jls.J1.C => 0.01/βc
    jls.R1.R => 50.0
    jls.C1.C => 2.0/βc
    jls.L1.L => 2.0/βc
    jls.L2.L => 100.0/βL
    jls.M12.L => 8.0/βL
    jls.Φₑ1.Φₑ => 0.5*Φ₀
    ]
    prob = NonlinearProblem(sys,zeros(19), ps)
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


