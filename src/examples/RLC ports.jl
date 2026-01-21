#In this example we model an parallel RLC circuit driven by a voltage source

#Load JosephsonLoops package
using JosephsonLoops 
const jls = JosephsonLoops

loops = [
["P1", "C1", "J1"]
]

circuit = jls.process_netlist(loops)
#cirucit model ODAE system and initial condition vector are created.
model, u0, guesses = jls.build_circuit(circuit)

# we set the values `of circuit parameters, for any parameters not specified; default values will be assigned.
    ps = [
        jls.P1.ω => 100e6*2*pi
        jls.P1.I => 1e-12
        jls.P1.Rport.R => 50.0
        jls.C1.C => 100.0e-15
        jls.J1.C=> 1000.0e-15
        jls.J1.I0 => 1e-6
        jls.J1.R => 1.0
    ]

#specify transient window for solver
tspan = (0.0, 1e-6)
#transient circuit analysis
sol = jls.tsolve(model,u0, ps, tspan, guesses=guesses)
jls.tplot(sol, jls.C1, model, units = "amps")

# we can pass any arguments known to the problem interface of DifferentialEquations.jl, to achieve better results
#saveat force the integrator to step at certain time points
saveat = LinRange(0.0, 1e-6, 1000)

#specify a different solving algorithim
alg = Rodas5()

sol = jls.tsolve(model, u0, ps, tspan;
    guesses=guesses, saveat = saveat, solver_opts = alg, abstol = 1e-6)
jls.tplot(sol, jls.C1, units = "amps")

#Harmonic Balance
using ModelingToolkit
eqs, states = jls.get_full_equations(model, jls.t)
 
harmonic_system, harmonic_states = jls.harmonic_equation(eqs, states, jls.t, jls.P1.Isrc.ω, 1)

heqs = harmonic_system[1:end-1]

#[0~Icos-ns.F[1]]
@named ns = NonlinearSystem(heqs)
sys_nns = toggle_namespacing(ns, false)
@time jac = calculate_jacobian(sys_nns)
sys =  mtkcompile(ns)
ns


# I₀ = 1e-6
# R₀ = 50.0
# Id = 0.05e-6
# ωc = sqrt(2*pi *I₀/(jls.Φ₀*1000.0e-15))/(2*pi)
# Ic = jls.Φ₀/(2*pi*1000.0e-12)

# 1/(2*pi*sqrt(1000.0e-12*1000.0e-15))

# N = 100
# ω_vec = 2*pi*(4.5:0.001:5.0)*1e9
# solution1=Float64[]
# solution2=Float64[]
# solution3=Float64[]
# u0_prev = zeros(6)
# for i in 1:1:length(ω_vec)
#     ps = [
#         jls.P1.Isrc.ω => ω_vec[i]
#         jls.P1.Isrc.I => 0.00565e-6
#         jls.C1.C => 100.0e-15
#         jls.J1.C=> 1000.0e-15
#         jls.J1.I0 => Ic
#         jls.P1.Rsrc.R => 50.0
#         jls.J1.R => 1e9
#     ]
#     prob = NonlinearProblem(sys,zeros(6), ps)
#     sol = solve(prob)
#     #u0_prev = sol
#     push!(solution1,  sqrt(sol[ns.E[2]]^2+sol[ns.F[1]]^2))
#     push!(solution2,  sqrt(sol[ns.G[2]]^2+sol[ns.H[1]]^2))
# end

ω_vec, solution1, solution2 = @time hbsweep(sys, jls, ns)
using Plots
Ii = @. abs(0.00565e-6 - (solution1))
Vi = @. abs(jls.Φ₀/(2*pi)*solution2)

ai = 0.5*(Vi+50.0*Ii)/sqrt(50.0)
bi = 0.5*(Vi-50.0*(Ii))/sqrt(50.0)

p = plot(ω_vec/(2*pi), (bi./ai))
display(p)

#Linear Analysis

ps = [
    jls.P1.Isrc.ω => 4.75e9*2*pi
    jls.P1.Isrc.I => 0.00565e-6
    jls.C1.C => 100.0e-15
    jls.J1.C=> 1000.0e-15
    jls.J1.I0 => Ic
    jls.P1.Rsrc.R => 50.0
    jls.J1.R => 1e9
]
prob = NonlinearProblem(sys, zeros(6), ps)
@time sol = solve(prob)
op_point = sol
unknowns(sys_nns)