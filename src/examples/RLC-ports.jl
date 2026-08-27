using JosephsonLoops
using Revise
using Symbolics
using ModelingToolkit
using BenchmarkTools

const jls = JosephsonLoops
loops = [
    ["P1", "C1", "J1"]
]

circuit = jls.process_netlist(loops)
model, u0, guesses = jls.build_circuit(circuit)
I₀ = jls.Φ₀/(2π*1000.0e-12)
R₀ = 10.0e3
ωc = R₀*I₀/(jls.Φ₀/2π)

ps = Dict(
    jls.P1.source.ω => 100e6*2*pi/ωc,
    jls.P1.source.I => 0.00565e-6/I₀,
    jls.P1.Rₙ.r => 50.0/R₀,
    jls.C1.βc => 100.0e-15*R₀*ωc,
    jls.J1.βc => 1000.0e-15*R₀*ωc,
    jls.J1.r => 1e8/R₀,
    jls.J1.α => 1.0
)

#time domain simulation 
tspan = (0.0, 1e-6).*ωc
tsol = jls.tsolve(model, guesses, ps, tspan; guesses=guesses)
p1 = jls.plot(tsol[jls.C1.i][end-400:end].*I₀, title = "Transient Time Plot", xlabel = "t", ylabel = "I_C1")

#hb Setup
# Define the sweep range (8 to 10.0 GHz)
ω_vec = collect(2*pi*(4.5:0.001:5.0)*1e9/ωc)

sweep_params = delete!(ps, jls.P1.source.ω)

sys = jls.HarmonicSystem(model, jls.P1.source.ω, 2)
prob = jls.HarmonicProblem(sys, sweep_params, parameter_sweep = [jls.P1.source.ω=>ω_vec])

result = jls.solve!(prob)

current = jls.get_solution(prob, jls.P1.i, 1)
theta_p_mag =  jls.get_solution(prob, jls.P1.dφ, 1)


Ii = @. current*I₀
Vi = @. (theta_p_mag)*R₀*I₀
# Calculate Power Waves (a = incident, b = reflected)
Z0 = 50.0
ai = @. 0.5 * (Vi + Z0 * Ii) / sqrt(Z0)
bi = @. 0.5 * (Vi - Z0 * Ii) / sqrt(Z0)
p = jls.plot(ω_vec/(2*pi), 20*log10.(abs.(bi./ai)), xlabel="Frequency (Hz)", ylabel="S11 (dB)", title="RLC S-Parameter", lw=2)

sys = jls.HarmonicSystem(model, jls.P1.source.ω, 2, determine_jacobian=true)

s11_expr = jls.get_HB_scattering_matrix(model, '1', '1')[1]
test = jls.get_solution(prob, s11_expr, (1,0))


# Linearised (JPA) analysis: hold the pump fixed and sweep a small signal across it.
# The pump tone has to be back in the parameter dict (the steady-state sweep above
# deleted it), and Ω_vec is the signal grid the response is evaluated on.
Ω_vec = ω_vec
lin_params = copy(sweep_params)
lin_params[jls.P1.source.ω] = 2π*4.75001e9/ωc     # MIT README JPA pump
lin_params[jls.P1.source.I] = 11.3e-9/I₀          # their one-sided Ip = 5.65 nA ≡ our 11.3 nA
δU_jpa = jls.perturbation_response(sys, jls.P1.source.I, lin_params, amplitude = lin_params[jls.P1.source.I])
jpa = jls.LinearisedProblem(sys, lin_params, δU_jpa, Ω_vec)
jls.solve!(jpa)
I_jpa = jls.get_solution(jpa, jls.P1.i, 1) .* I₀
V_jpa = jls.get_solution(jpa, jls.P1.dφ, 1) .* R₀ .* I₀
a_jpa = @. 0.5*(V_jpa + Z0*I_jpa)/sqrt(Z0)
b_jpa = @. 0.5*(V_jpa - Z0*I_jpa)/sqrt(Z0)
gain_dB = 20*log10.(abs.(b_jpa./a_jpa))
p_gain = jls.plot(Ω_vec.*ωc./(2π*1e9), gain_dB,
    xlabel="Frequency (GHz)", ylabel="|S₁₁| (dB)",
    title="JPA gain — JosephsonLoops vs JosephsonCircuits.jl",
    lw=2, label="JosephsonLoops", legend=:topright)


using DelimitedFiles
mit_csv = joinpath(pkgdir(jls), "mit_jpa.csv") 
if isfile(mit_csv)
    mit = readdlm(mit_csv, ',')
    jls.plot!(p_gain, mit[:, 1], mit[:, 2], lw=2, ls=:dash, label="JosephsonCircuits.jl")
    Δ = gain_dB .- Float64.(mit[:, 2])
    fs = Ω_vec .* ωc ./ (2π*1e9)
    bw = fs[gain_dB .>= maximum(gain_dB) - 3]
    println("ours: peak ", round(maximum(gain_dB), digits=3), " dB @ ",
            round(fs[argmax(gain_dB)], digits=4), " GHz, 3 dB width ",
            round(1000*(maximum(bw)-minimum(bw)), digits=1), " MHz")
    println("MIT:  peak ", round(maximum(Float64.(mit[:, 2])), digits=3), " dB")
    println("band max|Δ| = ", round(maximum(abs.(Δ)), digits=3),
            " dB, median |Δ| = ", round(sort(abs.(Δ))[end÷2+1], sigdigits=3), " dB")
else
    @warn "mit_jpa.csv not found — run mit_jpa_export.jl with the JosephsonCircuits-MIT project first."
end
# regenerates the figure used in the docs
jls.savefig(p_gain, joinpath(pkgdir(jls), "docs", "images", "jpa-single-tone.png"))
p_gain 