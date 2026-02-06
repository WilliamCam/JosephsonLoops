using JosephsonLoops
using ModelingToolkit
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

# --- 4. Frequency Sweep ---
# Define the sweep range (4.5 to 10.0 GHz)
ω_vec = 2*pi*(8:0.001:10.0)*1e9

# Create a copy of parameters for the sweep
sweep_params = copy(ps)
sweep_params[jls.P1.Isrc.I] = 0.00565e-6 # Specific current for this test
sweep_params[jls.J1.R] = 1e9   # High resistance (open circuit approximation)

#create the HarmonicProblem and HarmonicSweepProblem
h_prob = jls.HarmonicProblem(model, jls.P1.Isrc.ω, sweep_params ;tearing = true, N=1)
sweep_prob = jls.HarmonicSweepProblem(h_prob, jls.P1.Isrc.ω, ω_vec)
sweep_params[jls.J1.R] = 1e9             # High resistance (open circuit approximation)
# Solve the harmonic balance problem over the frequency sweep

# Solve the harmonic balance problem over the frequency sweep and the fundamental harmonic balance problem

sweep_res = jls.solve(sweep_prob)
h_prob_res = jls.solve(h_prob)


#Testing =======

current_p_mag = jls.get_harmonic(h_prob, sweep_res, "P1₊i", 1; output_type=:magnitude)
theta_p_mag = jls.get_harmonic(h_prob, sweep_res, "P1₊dθ", 1; output_type=:magnitude)

Ii = @. (0.00565e-6 - current_p_mag)

Vi = @. (jls.Φ₀ / (2*pi) * theta_p_mag) 

# Calculate Power Waves (a = incident, b = reflected)
Z0 = 50.0
#ask will if we need abs. here or not and ask about
ai = @. 0.5 * (Vi + Z0 * Ii) / sqrt(Z0)
bi = @. 0.5 * (Vi - Z0 * Ii) / sqrt(Z0)

plot(ω_vec/(2*pi), (bi./ai))