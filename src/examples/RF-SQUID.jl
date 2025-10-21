# In this example we model a DC SQUID driven by a current source

#Load JLoop package
using JosephsonLoops 
const jls = JosephsonLoops

#open the circuit netlist file named RF-SQUID.jld2, to create a new circuit netlist we can use jls.new_netlist()
jls.open_file("src//examples//RF-SQUID.jld2")


#cirucit model ODAE system and initial condition vector are created.
model, u0 = jls.build_circuit()

# we set the values of circuit parameters, for any parameters not specified; default values will be assigned.

I₀ = 1.0e-6
R₀ = 5.0
Φ₀ = jls.Φ₀

βc  = 2*pi/Φ₀ * I₀ * R₀^2
βL = 2*pi/Φ₀ * I₀

ps = [
    jls.loop4.sys.ω => 2*pi*100.0e+6
    jls.loop4.sys.I => 0.8*I₀
    jls.J1.sys.I0 => I₀
    jls.J1.sys.R => R₀
    jls.J1.sys.C => 0.01/βc
    jls.J1.sys.L => 1.0/βc
    jls.R1.sys.R => 50.0
    jls.C1.sys.C => 2.0/βc
    jls.L1.sys.L => 2.0/βL
    jls.L2.sys.L => 100.0/βL
    jls.M12.sys.L => 8.0/βL
    jls.loop1.sys.Φₑ => 0.5*Φ₀

]

tspan = (0.0, 1e-6)
saveat = LinRange(tspan[2]/10.0, tspan[2], 10000)
sol = jls.tsolve(model, u0, tspan, ps; saveat = saveat)
jls.tplot(sol, jls.R1, units="volts")

## Parameter Sweeps
Φspan = (0.0, 2.0*Φ₀)
ensemble_sol = jls.ensemble_parameter_sweep(model, u0, tspan, Φspan, ps, jls.loop1.sys.Φₑ, jls.R1, saveat = saveat)
plot(ensemble_sol.u)


Ispan = (0.0, 2*I₀)
ensemble_sol = jls.ensemble_parameter_sweep(model, u0, tspan, Ispan, ps, jls.loop4.sys.I, jls.R1, saveat = saveat)
plot(ensemble_sol.u)
