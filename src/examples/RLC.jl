#In this example we model an parallel RLC circuit driven by a voltage source

#Load scsim package

using JLoop

#open the circuit netlist file named RLC.jld2, to create a new circuit netlist we can use scsim.new_netlist()
JLoop.open_file("src//examples//RLC.jld2")

#cirucit model ODE system and initial condition vector are created.
model, u0 = JLoop.build_circuit()

# we set the values of circuit parameters, for any parameters not specified; default values will be assigned.
ps = [
    JLoop.loop1.sys.ω => 1.0
    JLoop.loop1.sys.I => 1.0
    JLoop.C1.sys.C => 1.0e-1
    JLoop.L1.sys.L => 1.0e-1
    JLoop.R1.sys.R => 10.0
]
using DifferentialEquations

#specify transient window for solver
tspan = (0.0, 100.0)

#transient circuit analysis
sol = JLoop.tsolve(model, u0, tspan, ps, alg = Rodas5())
JLoop.tplot(sol, JLoop.R1, units = "amps")

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


## Can we achieve this using harmonic balance ? ## Yes! TODO: Transformation function to put ODESystem into diff_eq form
eqs = equations(model)
using Symbolics

obs = observed(model)
sub1 = Symbolics.solve_for(observed(model)[4], scsim.loop3.sys.iₘ)
sub2 = Symbolics.solve_for(obs[3], scsim.loop2.sys.iₘ)
sub3 = Symbolics.solve_for(obs[1], scsim.loop1.sys.iₘ)

eqs1 = substitute(eqs, Dict(scsim.loop3.sys.iₘ => sub1))
eqs2 = substitute(eqs1, Dict(scsim.loop2.sys.iₘ => sub2))
eqs3 = substitute(eqs2, Dict(scsim.loop1.sys.iₘ => sub3))

sub4 =  Symbolics.solve_for(obs[2], scsim.R1.sys.i)
eqs4 = substitute(eqs3, Dict(scsim.R1.sys.i => sub4))

sub5 = Symbolics.solve_for(eqs4[4], scsim.C1.sys.i)
eqs5 = substitute(eqs4, Dict(scsim.C1.sys.i => sub5))

hb_diff_eq = simplify(eqs5[1])

hb_diff_eq

using HarmonicBalance

diff_eq = DifferentialEquation(
    d(scsim.C1.sys.θ, scsim.t, 2)*(2.067833848e-15*scsim.C1.sys.C*scsim.L1.sys.L*scsim.R1.sys.R) + 2.067833848e-15*scsim.R1.sys.R*scsim.C1.sys.θ + 2.067833848e-15*scsim.L1.sys.L*d(scsim.C1.sys.θ,scsim.t) ~ 6.283185307179586*scsim.L1.sys.L*scsim.R1.sys.R*scsim.loop1.sys.I*cos(scsim.loop1.sys.ω*scsim.t),
    scsim.C1.sys.θ  
)

# diff_eq = DifferentialEquation(
#     d(scsim.C1.sys.θ, scsim.t, 2) + scsim.C1.sys.θ + scsim.L1.sys.L*d(scsim.C1.sys.θ,scsim.t) ~ scsim.loop1.sys.I*cos(scsim.loop1.sys.ω*scsim.t),
#     scsim.C1.sys.θ  
#)
add_harmonic!(diff_eq, scsim.C1.sys.θ, scsim.loop1.sys.ω) 
harmonic_eq = get_harmonic_equations(diff_eq)

varied = scsim.loop1.sys.ω => LinRange(900.0, 1100.0, 100) # range of parameter values

fixed = (scsim.L1.sys.L => 1.0e-3, scsim.loop1.sys.I => 1.0e-6, scsim.R1.sys.R => 50.0, scsim.C1.sys.C=>1.0e-3)

solutions = get_steady_states(harmonic_eq, varied, fixed)

plt = plot_1D_solutions(solutions, x="loop1₊ω", y="sqrt(u1^2 + v1^2)");

x = [solutions[i][1][scsim.loop1.sys.ω] for i in range(1,100)]

y = [sqrt(solutions[i][1][HarmonicBalance.u1]^2 + solutions[i][1][HarmonicBalance.v1]^2) for i in range(1,100)]

plot(real(x),abs.(y))


##############################Testing HarmonicBalance.jl to see if it can handle algebraic expressions

diff_eq = DifferentialEquation(
    [d(scsim.C1.sys.θ, scsim.t, 2)*(2.067833848e-15*scsim.C1.sys.C)  ~ 6.283185307179586*scsim.C1.sys.i,
    (scsim.L1.sys.L*(scsim.loop1.sys.I*cos(scsim.loop1.sys.ω*scsim.t) + (-3.2910597840193497e-16*d(scsim.C1.sys.θ, scsim.t)) / scsim.R1.sys.R) - 3.2910597840193497e-16*scsim.C1.sys.θ) / scsim.L1.sys.L],
    [scsim.C1.sys.θ, scsim.C1.sys.i]  
)

# diff_eq = DifferentialEquation(
#     d(scsim.C1.sys.θ, scsim.t, 2) + scsim.C1.sys.θ + scsim.L1.sys.L*d(scsim.C1.sys.θ,scsim.t) ~ scsim.loop1.sys.I*cos(scsim.loop1.sys.ω*scsim.t),
#     scsim.C1.sys.θ  
#)
add_harmonic!(diff_eq, scsim.C1.sys.θ, scsim.loop1.sys.ω) 
harmonic_eq = get_harmonic_equations(diff_eq)

varied = scsim.loop1.sys.ω => LinRange(900.0, 1100.0, 100) # range of parameter values

fixed = (scsim.L1.sys.L => 1.0e-3, scsim.loop1.sys.I => 1.0e-6, scsim.R1.sys.R => 50.0, scsim.C1.sys.C=>1.0e-3)

solutions = get_steady_states(harmonic_eq, varied, fixed)

plt = plot_1D_solutions(solutions, x="loop1₊ω", y="sqrt(u1^2 + v1^2)");

x = [solutions[i][1][scsim.loop1.sys.ω] for i in range(1,100)]

y = [sqrt(solutions[i][1][HarmonicBalance.u1]^2 + solutions[i][1][HarmonicBalance.v1]^2) for i in range(1,100)]

plot(real(x),abs.(y))