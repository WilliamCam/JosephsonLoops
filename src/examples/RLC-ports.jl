
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

sys = jls.HarmonicSystem(model, jls.P1.Isrc.ω, 2)
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

# The whole fix: build the correct input perturbation. source_perturbation_vector locates the
# source drive (U = −∂F/∂I — rows, quadratures, signs, scalings all automatic) and combines its
# two quadratures into the complex injection U_cos − i·U_sin, a pure signal sideband (e^{+iΩt}).
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

# Upper-sideband port-voltage response (V = Φ₀/2π · dθ) to the unit signal drive.
Z0 = 50.0
V_sig = (jls.Φ₀ / (2*pi)) .* jls.get_output(sys, lin_prob, lin_res, "P1₊dθ", 1)

# Phase-PRESERVING reflection at the Norton port (Isrc ∥ Rsrc ∥ DUT share one node pair, so
# V_DUT = V_port): the unit signal sideband carries incident current 2δI, so S(0,0) =
# V_sig/(Z0·δI) − 1. This is COMPLEX — it keeps the reflection phase — and is the signal-to-
# signal scattering that JosephsonCircuits.jl reports. |S(0,0)| ≈ +13.3 dB at the degenerate
# point; with the pump off it collapses to the ordinary linear S11. (A single real-quadrature
# drive would instead read the amplified/squeezed quadrature of the phase-sensitive amplifier.)
S11 = @. V_sig / (Z0 * δI) - 1

p_mag = jls.plot(Ω_vec/(2*pi*1e9), 20*log10.(abs.(S11)),
    xlabel="Frequency (GHz)", ylabel="|S₁₁| (dB)",
    title=" S₁₁ — magnitude (matches JosephsonCircuits.jl)",
    lw=2, legend=false)

