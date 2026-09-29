using JosephsonLoops
using Plots

# ---- circuit ----------------------------------------------------------------------
loops = [
    ["I1", "R1", "Lin"],
    ["Lt", "J1", "J2"],
    ["J2", "Idc"],
    ["Idc", "R2"], 
]

circuit = process_netlist(loops, mutual_coupling=[(1,2)], ext_flux=[false,true,false,false])
squid, u0, guesses = build_circuit(circuit)

# ---- normalisation and parameters --------------------------------------------------
I₀ = 10e-6  # critical current of a 1000 pH junction
R₀ = 10.0e3
ωc = R₀*I₀/(Φ₀/2π)
Z0 = 50.0
f_pump = 500e6
1/(2*pi/Φ₀ * I₀)


ps = Dict(
    squid.I1.ω => 2π*f_pump/ωc,
    squid.I1.I => 10e-9/I₀,
    squid.R1.r => 50.0/R₀,
    squid.R2.r => 50.0/R₀,
    squid.Idc.I       => 3.0,
    squid.Idc.ω       => 0.0,
    squid.Lin.βL      => 400,
    squid.J1.βc       => 1000.0e-15*R₀*ωc,
    squid.J1.r        => 6.0/R₀,
    squid.J2.βc       => 1000.0e-15*R₀*ωc,
    squid.J2.r        => 6.0/R₀,
    squid.Lt.βL       => 1.0,
    squid.M12.βM       => 2.0,
    squid.Φₑ2.Φₑ => 0.5*2π
)


#timescale τ is normalised                       
tspan = (0.0, 1e-8) .* ωc


using DifferentialEquations
tsol  = tsolve(squid, guesses, ps, tspan; guesses = guesses,solver_opts = Rodas5())
p_td  = plot(tsol[squid.R2.i][1:end] .* I₀, legend = false,
                 xlabel = "Sample", ylabel = "I(C1) (A)", title = "Transient at 100 MHz")