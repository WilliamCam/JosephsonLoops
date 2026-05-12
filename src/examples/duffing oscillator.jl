using Symbolics
using SymbolicUtils
using ModelingToolkit
using QuestBase
using DifferentialEquations
using Plots
using JosephsonLoops
using Revise 
const jls = JosephsonLoops
# Linear analysis adapted from Kosata 2022 thesis. Mirrors build_circuit/linear_analysis.jl
# but goes through the HarmonicSystem / LinearizedProblem wrappers.

# set up example differential equation in time domain. Duffing oscillator from p.59
@variables t x(t)
@parameters  α ω ω0 F γ η
diff_eq = Differential(t)(Differential(t)(x)) + ω0^2*x + α*x^3 + η*Differential(t)(x)*x^2 + γ*Differential(t)(x) - F*cos(ω*t) ~ 0

Nharmonics = 1

# Wrap the Duffing equation in an ODESystem so HarmonicSystem can consume it
@named duffing_sys = System([diff_eq], t)
duffing_model = mtkcompile(duffing_sys)

# Build the harmonic system with linearization data populated (J0, J1)
h_sys = jls.HarmonicSystem(duffing_model, ω; N=Nharmonics, linear=true)

# Harmonic-coefficient symbols from the variable_map
A_dc = h_sys.variable_map[("x", 0, :Cos)]
A1   = h_sys.variable_map[("x", 1, :Cos)]
B1   = h_sys.variable_map[("x", 1, :Sin)]

# Solve HB system for frequencies (large-signal pump)
N = 200
ω_vec = range(0.8, 1.2, N)
ps = Dict(α => 1.0, ω0 => 1.0, F => 0.01, η => 1.0e-1, γ => 1.0e-3)
h_prob = jls.HarmonicProblem(h_sys, ps; sweep_var=ω, sweep_vals=ω_vec)
h_result = jls.solve(h_prob)

# Pump amplitude vs ω (DC + fundamental magnitude, matches prototype)
solution = h_result.results[A_dc] .+ sqrt.(h_result.results[A1].^2 .+ h_result.results[B1].^2)
plot(ω_vec, solution)

# Linearize around a single pump tone at j-th sweep point (warm-start lands on high branch)
j = 70
ωp = ω_vec[j]

# Working point: every harmonic coefficient evaluated at ωp
U₀ = Dict{Num, Float64}(v => vals[j] for (v, vals) in h_result.results)

# Small-signal frequency sweep
N = 800
Ω = range(0.8, 1.2, N)

# Perturb B[1] (sin@1) at amplitude 1e-3 — same as prototype's perturb[3] = 1e-3
lin_prob = jls.LinearizedProblem(h_sys, ps, ωp, U₀, Ω;
    perturb_var="x", perturb_order=1, perturb_component=:Sin, perturb_amplitude=1.0e-3)
lin_result = jls.solve(lin_prob)

# Small-signal phasor: u = B[1] response, v = A[2] response, signal = |(u + iv)/2|
u_resp = lin_result.results[B1]
v_resp = lin_result.results[A1]
out = abs.((u_resp .+ 1im .* v_resp) ./ 2)

plot(Ω, abs.(out))
