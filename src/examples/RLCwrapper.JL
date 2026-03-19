using JosephsonLoops

const jls = JosephsonLoops
loops = [
    ["P1", "C1", "J1"]
]

@time circuit = jls.process_netlist(loops)
@time model, u0, guesses = jls.build_circuit(circuit)
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

#time domain simulation 
tspan = (0.0, 1e-6)
tsol = jls.tsolve(model, u0, ps, tspan; guesses=guesses)
p1 = jls.plot(tsol[jls.C1.i][end-400:end], title = "Transient Time Plot", xlabel = "t", ylabel = "I_C1")


#hb Setup
# Define the sweep range (8 to 10.0 GHz)
ω_vec = 2*pi*(8:0.001:10.0)*1e9

# Create a copy of parameters for the sweep
sweep_params = copy(ps) 
sweep_params[jls.P1.Isrc.I] = 0.00565e-6 # Specific current for this test
sweep_params[jls.J1.R] = 1e9  
sweep_params[jls.J1.R] = 1e9    

#create the HarmonicProblem and HarmonicSweepProblem
@time h_prob = jls.HarmonicProblem(model, jls.P1.Isrc.ω, sweep_params ;tearing = true, N=2)
@time sweep_prob = jls.HarmonicSweepProblem(h_prob, jls.P1.Isrc.ω, ω_vec)


# Solve the harmonic balance problem over the frequency sweep and the fundamental harmonic balance problem
@time sweep_res = jls.solve(sweep_prob)
@time h_prob_res = jls.solve(h_prob)


@time current_p_mag = jls.get_harmonic(h_prob, sweep_res, "J1₊i", 2; output_type=:magnitude)
@time theta_p_mag = jls.get_harmonic(h_prob, sweep_res, "J1₊θ", 2; output_type=:magnitude)

# obs variables C1.i and J1.i
@time c1_i_mag = jls.get_harmonic(h_prob, sweep_res, "C1₊i", 1; output_type=:magnitude)
@time j1_i_mag = jls.get_harmonic(h_prob, sweep_res, "J1₊i", 1; output_type=:magnitude)

Ii = @. (0.00565e-6 - current_p_mag)
Vi = @. (jls.Φ₀ / (2*pi) * theta_p_mag) 

# Calculate Power Waves (a = incident, b = reflected)
Z0 = 50.0
ai = @. 0.5 * (Vi + Z0 * Ii) / sqrt(Z0)
bi = @. 0.5 * (Vi - Z0 * Ii) / sqrt(Z0)

p = jls.plot(ω_vec/(2*pi), 20*log10.(abs.(bi./ai)), xlabel="Frequency (Hz)", ylabel="S11 (dB)", title="RLC S-Parameter", lw=2)
