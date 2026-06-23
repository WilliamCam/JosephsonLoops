
using Revise
using JosephsonLoops
using Symbolics
using ModelingToolkit
using BenchmarkTools

const jls = JosephsonLoops
loops = [
    ["P1", "C1", "J1"]
]

@time circuit = jls.process_netlist(loops)
@time model, u0, guesses = jls.build_circuit(circuit)
ps = Dict(
    jls.P1.Isrc.ω => 100e6*2*pi,
    jls.P1.Isrc.I => 0.00565e-6,
    jls.P1.Rsrc.R => 50.0,
    jls.C1.C => 100.0e-15,
    jls.J1.C => 1000.0e-15,
    jls.J1.I0 => jls.Φ₀/(2π*1000.0e-12),
    jls.J1.R => 1e12
)

#time domain simulation 
#tspan = (0.0, 1e-6)
#tsol = jls.tsolve(model, guesses, ps, tspan; guesses=guesses)
#p1 = jls.plot(tsol[jls.C1.i][end-400:end], title = "Transient Time Plot", xlabel = "t", ylabel = "I_C1")

#hb Setup
# Define the sweep range (8 to 10.0 GHz)
ω_vec = collect(2*pi*(4.5:0.001:5.0)*1e9)


sweep_params = delete!(ps, jls.P1.Isrc.ω)

sys = jls.HarmonicSystem(model, jls.P1.Isrc.ω, 16)
prob = jls.HarmonicProblem(sys, ω_vec, sweep_params)

result = jls.solve!(prob)
out = prob.result.solution[jls.P1.Isrc.ω]

current = jls.get_output(prob, result, "C1₊i", 1)
theta_p_mag = jls.get_output(prob, result, "P1₊dθ",  1)

Ii = @. (0.00565e-6 - current)
Vi = @. (jls.Φ₀ / (2*pi) * real.(theta_p_mag)) 
# Calculate Power Waves (a = incident, b = reflected)
Z0 = 50.0
ai = @. 0.5 * (Vi + Z0 * Ii) / sqrt(Z0)
bi = @. 0.5 * (Vi - Z0 * Ii) / sqrt(Z0)
p = jls.plot(ω_vec/(2*pi), 20*log10.(abs.(bi./ai)), xlabel="Frequency (Hz)", ylabel="S11 (dB)", title="RLC S-Parameter", lw=2)

#  Linear (small-signal) analysis around a pump tone — mirrors the JPA example in MIT's
#  JosephsonCircuits.jl (port ∥ 50Ω, series Cc=100fF, junction Lj=1nH ∥ Cj=1000fF).
sys = jls.HarmonicSystem(model, jls.P1.Isrc.ω, 2, determine_jacobian=true)

# Linearised responses are ordered by the jacobian's `vars` ordering — [DC, cos₁, sin₁,
# cos₂, sin₂] per state — NOT by unknowns(system). The row map gives the bookkeeping:
jls.linearised_row_map(sys)

# Pump tone at the MIT example's frequency. Two conventions to mind when comparing:
# (1) JosephsonCircuits' source `current=Ip` is a one-sided spectral amplitude — physical
#     peak current is 2·Ip, so their 5.65 nA ≡ 11.3 nA here;
# (2) the MIT junction is lossless (Lj ∥ Cj) while our RCSJ junction carries the 10 kΩ
#     shunt, raising the critical pump further. Drive at 15 nA, where the linearised
#     gain peaks (+13.3 dB at the degenerate point).
ωp = 2*pi*4.75001e9
jpa_params = copy(sweep_params)
jpa_params[jls.P1.Isrc.I] = 11.3e-9

# 1 nA test signal on the port current source: U = -∂F/∂I locates the source equation
# rows, quadratures, signs and equation scalings automatically.
δI = 1.0e-10
pert = jls.source_perturbation_vector(sys, jls.P1.Isrc.I, jpa_params; amplitude=δI)

# Pump working point via downward continuation (5.0 GHz → ωp): the softening junction
# pulls its resonance toward the pump, so approaching from above lands the driven branch.
ω_down = collect(range(2*pi*5.0e9, ωp, 120))
pump_prob = jls.HarmonicProblem(sys, ω_down, jpa_params)
jls.solve!(pump_prob)
U₀ = real.(pump_prob.result.solution[jls.P1.Isrc.ω][:, end])

Ω_vec = collect(2*pi*(4.5:0.001:5.0)*1e9)
lin_prob = jls.HarmonicProblem(sys, Ω_vec, jpa_params; U₀=U₀, linear_response = (ωp, pert))
lin_res = jls.solve!(lin_prob)

# Port voltage response to the sin (P-quadrature) source drive. V = Φ₀/2π · dθ.
Z0 = 50.0
V_sin = (jls.Φ₀ / (2*pi)) .* jls.get_output(sys, lin_prob, lin_res, "P1₊dθ", 1)

# Single-quadrature S11 at the Norton port (Isrc ∥ Rsrc ∥ DUT share one node pair, so
# V_DUT = V_port and I_DUT = δI − V/Z0):  S11 = 2V/(Z0·δI) − 1, with the sin drive's δI
# phasor i·δI. At the degenerate point this measures the AMPLIFIED quadrature (~+19 dB).
S11_amp = @. 2V_sin / (Z0 * (im * δI)) - 1

# A degenerate JPA is phase-SENSITIVE: the single quadrature above is the amplified (or
# squeezed) gain, not the phase-PRESERVING signal gain |S_ss| that nodal HB codes
# (JosephsonCircuits.jl) report. Recover |S_ss| by also driving the orthogonal (cos)
# quadrature and combining: |S_ss| = (σ₊ + σ₋)/2 of the 2×2 quadrature field-transfer map.
# This is the curve to compare against JosephsonCircuits' S(0,0) (≈ +13.3 dB here).
pert_cos = jls.rotate_quadrature(sys, pert)
lin_cos  = jls.HarmonicProblem(sys, Ω_vec, jpa_params; U₀=U₀, linear_response=(ωp, pert_cos))
jls.solve!(lin_cos)
V_cos = (jls.Φ₀ / (2*pi)) .* jls.get_output(sys, lin_cos, lin_cos.result, "P1₊dθ", 1)
S_ss  = jls.phase_preserving_s11(V_cos, V_sin, Z0, δI)

p_amp = jls.plot(Ω_vec/(2*pi*1e9), 20*log10.(abs.(S11_amp)),
    xlabel="Frequency (GHz)", ylabel="Gain (dB)",
    title="Amplified quadrature (phase-sensitive)",
    lw=2, legend=false)

p_ss = jls.plot(Ω_vec/(2*pi*1e9), 20*log10.(S_ss),
    xlabel="Frequency (GHz)", ylabel="Gain (dB)",
    title="Phase-preserving |S_ss| (matches JosephsonCircuits.jl)",
    lw=2, ls=:dash, legend=false)