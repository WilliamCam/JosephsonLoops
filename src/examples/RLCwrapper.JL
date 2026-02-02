using ModelingToolkit
using JosephsonLoops
using Plots
using Symbolics
using DifferentialEquations

#need to include harmonic balance files on repl/terminal

const jls = JosephsonLoops

loops = [
    ["P1", "C1", "J1"]
]

circuit = jls.process_netlist(loops)
model, u0, guesses = jls.build_circuit(circuit)

ps = Dict(
    jls.P1.Isrc.ω => 100e6*2*pi,
    jls.P1.Isrc.I => 1e-12,
    jls.P1.Rsrc.R => 50.0,
    jls.C1.C => 100.0e-15,
    jls.J1.C => 1000.0e-15,
    jls.J1.I0 => 1e-6,
    jls.J1.R => 1.0
)
# Constants
Ic = 1e-6 
Φ₀ = 2.067833848e-15

#time domain simulation wrapper updates- todo
#linear analysis wrapper - todo



#hb Setup
#ask will how he set up includes/ throwing compilation error here when i try to include HB files separately or from Josephson Loops

h_prob = jls.HarmonicProblem(model, jls.P1.Isrc.ω; N=1, tearing = true)

# --- 4. Frequency Sweep ---
# Define the sweep range (4.5 to 5.0 GHz)
ω_vec = 2*pi*(4.5:0.001:5.0)*1e9

# Create a copy of parameters for the sweep
sweep_params = copy(ps)
sweep_params[jls.P1.Isrc.I] = 0.00565e-6 # Specific current for this test

sweep_params[jls.J1.R] = 1e9             # High resistance (open circuit approximation)
# Solve the harmonic balance problem over the frequency sweep
sweep_res = jls.solve_sweep(h_prob, sweep_params, jls.P1.Isrc.ω => ω_vec)




#Testing and debugging as of now. Able to get P1, C1 , J1 but not P.Rsrc.i 

current_p = jls.get_harmonic_coeffs(h_prob, model, sweep_res, ps, "P1₊i")
thetap = jls.get_harmonic_coeffs(h_prob, model, sweep_res, ps, "P1₊dθ")


phase = jls.get_harmonic_phase(h_prob, model, sweep_res, ps, "P1")

# solution1 = get_harmonic_coeffs(h_prob, model, sweep_res, ps, "P1")
# solution2 = get_harmonic_coeffs(h_prob, model, sweep_res, ps, "C1")

Ii = @. abs(0.00565e-6 - current_p)

Vi = @. abs(jls.Φ₀ / (2*pi) * thetap) 

# Calculate Power Waves (a = incident, b = reflected)
Z0 = 50.0
ai = @. 0.5 * (Vi + Z0 * Ii) / sqrt(Z0)
bi = @. 0.5 * (Vi - Z0 * Ii) / sqrt(Z0)

#4. Plot S11-
p = plot(ω_vec / (2*pi), abs.(bi ./ ai), 
    title = "S11 Scattering Parameter",
    xlabel = "Frequency (Hz)", 
    ylabel = "|S11|", 
    label = "S11", 
    lw = 2)

display(p)