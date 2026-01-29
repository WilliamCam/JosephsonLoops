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

#hb setup
h_prob = HarmonicProblem(model, jls.P1.Isrc.ω; N=1, tearing = true)

# --- 4. Frequency Sweep ---
# Define the sweep range (4.5 to 5.0 GHz)
ω_vec = 2*pi*(4.5:0.001:5.0)*1e9

# Create a copy of parameters for the sweep
sweep_params = copy(ps)
sweep_params[jls.P1.Isrc.I] = 0.00565e-6 # Specific current for this test
sweep_params[jls.J1.R] = 1e9             # High resistance (open circuit approximation)

# UPDATED CALL: Matches `solve_sweep(prob, base_params, sweep_pair)`
sweep_res = solve_sweep(h_prob, sweep_params, jls.P1.Isrc.ω => ω_vec)

thetap = get_harmonic_coeffs(model, sweep_res, "P1")
    
    