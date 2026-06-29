
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
I₀ = jls.Φ₀/(2π*1000.0e-12)
R₀ = 100e3
ωc = R₀*I₀/(jls.Φ₀/2π)

ps = Dict(
    jls.P1.Isrc.ω => 100e6*2*pi/ωc,
    jls.P1.Isrc.I => 0.00565e-6/I₀,
    jls.P1.Rsrc.r => 50.0/R₀,
    jls.C1.βc => 100.0e-15*R₀*ωc,
    jls.J1.βc => 1000.0e-15*R₀*ωc,
    jls.J1.r => 1.0
)

#time domain simulation 
tspan = (0.0, 1e-6).*ωc
tsol = jls.tsolve(model, guesses, ps, tspan; guesses=guesses)
p1 = jls.plot(tsol[jls.C1.i][end-400:end].*I₀, title = "Transient Time Plot", xlabel = "t", ylabel = "I_C1")

#hb Setup
# Define the sweep range (8 to 10.0 GHz)
ω_vec = collect(2*pi*(4.5:0.001:5.0)*1e9/ωc)

sweep_params = delete!(ps, jls.P1.Isrc.ω)

sys = jls.HarmonicSystem(model, jls.P1.Isrc.ω, 2)
prob = jls.HarmonicProblem(sys, ω_vec, sweep_params)

result = jls.solve!(prob)
out = prob.result.solution[jls.P1.Isrc.ω]

current = jls.get_output(prob, result, "C1₊i", 1)
theta_p_mag = jls.get_output(prob, result, "P1₊dφ",  1)

eq = jls.get_harmonic_expression(prob, "C1₊i", 1)

Ii = @. (0.00565e-6/I₀ - current)*I₀
Vi = @. real.(theta_p_mag)*R₀*I₀
# Calculate Power Waves (a = incident, b = reflected)
Z0 = 50.0
ai = @. 0.5 * (Vi + Z0 * Ii) / sqrt(Z0)
bi = @. 0.5 * (Vi - Z0 * Ii) / sqrt(Z0)
p = jls.plot(ω_vec/(2*pi), 20*log10.(abs.(bi./ai)), xlabel="Frequency (Hz)", ylabel="S11 (dB)", title="RLC S-Parameter", lw=2)

#  Linear (small-signal) analysis around a pump tone — mirrors the JPA example in MIT's
#  JosephsonCircuits.jl (port ∥ 50Ω, series Cc=100fF, junction Lj=1nH ∥ Cj=1000fF).
sys = jls.HarmonicSystem(model, jls.P1.Isrc.ω, 2, determine_jacobian=true)

ωp = 2*pi*4.75001e9./ωc
δI = 11.3e-9/I₀

jpa_params = copy(sweep_params)
jpa_params[jls.P1.Isrc.I] = δI

pert = jls.perturbation_response(sys, jls.P1.Isrc.I, jpa_params)

ω_down = collect(range(2*pi*5.0e9/ωc, ωp, 120))
pump_prob = jls.HarmonicProblem(sys, ω_down, jpa_params)
jls.solve!(pump_prob)
U₀ = real.(pump_prob.result.solution[jls.P1.Isrc.ω][:, end])

Ω_vec = collect(2*pi*(4.5:0.001:5.0)*1e9)/ωc
lin_prob = jls.HarmonicProblem(sys, Ω_vec, jpa_params; U₀=U₀, linear_response = (ωp, pert))
lin_res = jls.solve!(lin_prob)

Z0 = 50.0
V_sig = I₀*R₀.*jls.get_output(sys, lin_prob, lin_res, "P1₊dφ", 1)
using Plots
plot(abs.(V_sig))
S11 = @. V_sig / (Z0/R₀ * I₀) - 1

p_mag = jls.plot(Ω_vec/(2*pi*1e9), 20*log10.(abs.(S11)),
    xlabel="Frequency (GHz)", ylabel="|S₁₁| (dB)",
    title="JPA gain — JosephsonLoops vs JosephsonCircuits.jl",
    lw=2, label="JosephsonLoops", legend=:topright)



# Overlay JosephsonCircuits.jl. First generate the CSV once, from the MIT project:
#   cd ../JosephsonCircuits-MIT/JosephsonCircuits.jl
#   julia --project=. ../../JosephsonLoops/mit_jpa_export.jl
using DelimitedFiles
mit_csv = joinpath(@__DIR__, "..", "..", "mit_jpa.csv")
if isfile(mit_csv)
    mit = readdlm(mit_csv, ',')
    jls.plot!(p_mag, mit[:, 1], mit[:, 2], lw=2, ls=:dash, label="JosephsonCircuits.jl")
else
    @warn "mit_jpa.csv not found — run mit_jpa_export.jl in the JosephsonCircuits project first."
end


