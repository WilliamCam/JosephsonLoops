using Symbolics
using SymbolicUtils
using ModelingToolkit
using QuestBase
using DifferentialEquations
using Plots
using JosephsonLoops
using LinearAlgebra
#Linear analysis adapted from Kosata 2022 thesis

#set up example differential equation in time domain. We will use duffing oscillator example from p.59
@variables t x(t)
@parameters  α ω ω0 F γ η
diff_eq = Differential(t)(Differential(t)(x)) + ω0^2*x + α*x^3 + η*Differential(t)(x)*x^2 + γ*Differential(t)(x) - F*cos(ω*t) ~ 0

D = Differential(t)
# Duffing equation: x'' + 2γx' + ω0^2x + αx^3 = F*cos(ωt)
# Note: Using small damping and nonlinearity for testing
diff_eq = D(D(x)) + 2*γ*D(x) + ω0^2*x + α*x^3 ~ F*cos(ω*t)

@named sys = ODESystem([diff_eq], t)

# Parameters
# ω0 = 1.0 (natural freq)
# α = 0.1 (nonlinearity)
# F = 0.5 (drive amplitude)
# γ = 0.05 (damping)
flux_vals = Dict(ω0 => 1.0, α => 0.1, F => 0.5, γ => 0.05)

# 2. Setup Harmonic Problem with Linearization
# We'll sweep frequency ω around the resonance (ω0 = 1.0)
ω_sweep_range = 0.5:0.01:1.5
param_sweep = ω => ω_sweep_range

# Create HarmonicProblem dealing with 1 harmonic (N=1)
# IMPORTANT: Set linear=true to compute Jacobians
h_prob = HarmonicProblem(sys, ω, flux_vals; N=1, linear=true)

# 3. Solve the nonlinear sweep
sweep_prob = HarmonicSweepProblem(h_prob, ω, ω_sweep_range)
res = solve(sweep_prob)

# Plot the nonlinear response amplitude
# Reconstruct x(t) amplitude from harmonic coefficients
# x(t) ≈ A0 + A1*cos(ωt) + B1*sin(ωt)
# Amplitude ≈ sqrt(A1^2 + B1^2) (ignoring DC and higher harmonics for now)

# Get variable names for harmonics
A1_sym = h_prob.variable_map[("x", 1, :Cos)]
B1_sym = h_prob.variable_map[("x", 1, :Sin)]

amp = sqrt.(res.results[A1_sym].^2 .+ res.results[B1_sym].^2)
plot(res.sweep_vals, amp, label="Nonlinear Response", xlabel="Frequency ω", ylabel="Amplitude", title="Duffing Oscillator Sweep")
savefig("duffing_nonlinear_sweep.png")

# 4. Calculate Linear Response
# Let's pick an operating point near resonance, e.g., ω = 1.0
op_freq = 1.0
println("Calculating linear response at ω = $op_freq")

# Define probe frequencies for linear response (small signal analysis)
# We look at response to perturbations at frequencies Ω around the drive frequency
probe_freqs = 0.8:0.01:1.2

# Define perturbation vector
# We perturb the cosine component of the first harmonic (A1)
# The system variables order determines the index
n_vars = length(h_prob.sys_vars)
perturb_vec = zeros(n_vars)

# Find index of A1
A1_idx = findfirst(isequal(A1_sym), h_prob.sys_vars)
perturb_vec[A1_idx] = 1.0 

# Calculate linear response
lin_resp = linear_response(h_prob, res, op_freq, probe_freqs, perturb_vec)

# Plot Linear Response Magnitude
# We can just look at the magnitude of the response vector
resp_mag = [norm(lin_resp[:, i]) for i in 1:length(probe_freqs)]

plot(probe_freqs, resp_mag, label="Linear Response Magnitude", xlabel="Probe Frequency Ω", ylabel="|Response|", title="Linear Response at ω=$op_freq")
savefig("duffing_linear_response.png")

println("Test complete. Plots saved.")

