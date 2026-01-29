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


println("Running Frequency Sweep...")

# UPDATED CALL: Matches `solve_sweep(prob, base_params, sweep_pair)`
sweep_res = solve_sweep(h_prob, sweep_params, jls.P1.Isrc.ω => ω_vec)



#-----Confused about this entire section- need help from will

# Get variables from the system
sys_vars = unknowns(h_prob.sys)
# Helper to find variable by string name
find_var(name) = first(filter(v -> occursin(name, string(v)), sys_vars))

# Extract Phase Coefficients (A=Cos, B=Sin) for N=1
# A[2] corresponds to Cos(wt), B[1] corresponds to Sin(wt)
A_fund = sweep_res.results[find_var("A[2]")] 
B_fund = sweep_res.results[find_var("B[1]")] 


# solution2 = Phase Amplitude
solution2 = @. sqrt(A_fund^2 + B_fund^2)

# Physics: Calculate Voltage Amplitude |V| = (Φ₀/2π) * ω * |φ|
# We multiply by sweep_vals (ω) because V is the derivative of flux/phase
Vi = @. abs(Φ₀/(2*pi) * sweep_res.sweep_vals * solution2)

# solution1 = Current drawn by the resistor (Vi / R_src)
# This allows us to calculate Ii = |I_src - I_R|
solution1 = Vi ./ 50.0 


Ii = @. abs(0.00565e-6 - (solution1))

# Scattering Parameters
ai = @. 0.5*(Vi + 50.0*Ii)/sqrt(50.0)
bi = @. 0.5*(Vi - 50.0*Ii)/sqrt(50.0)

# Plotting
p = plot(ω_vec/(2*pi), abs.(bi ./ ai), 
         title="Reflection Coefficient (S11)", 
         xlabel="Frequency (Hz)", ylabel="|b/a|", 
         lw=2, legend=false)

display(p)
savefig("rlc_harmonic_sweep.png")