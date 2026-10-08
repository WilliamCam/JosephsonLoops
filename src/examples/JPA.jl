

using JosephsonLoops
using ModelingToolkit
using Plots

# ---- circuit ----------------------------------------------------------------------
loops = [["P1", "C1", "J1"]]
circuit = process_netlist(loops)
jpa, u0, guesses = build_circuit(circuit)

# ---- normalisation ----------------------------------------------------------------
I₀ = Φ₀/(2π*1000.0e-12)          # critical current of a 1000 pH junction, so J1.α = 1
R₀ = 10.0e3
ωc = R₀*I₀/(Φ₀/2π)
Z0 = 50.0
GHz(f) = 2π*f*1e9/ωc             # frequency in GHz to the normalised angular frequency


circuit_ps = Dict(
    jpa.P1.Rₙ.r => 50.0/R₀,          # 50 ohm port
    jpa.C1.βc   => 100.0e-15*R₀*ωc,   # 100 fF coupling capacitor
    jpa.J1.βc   => 1000.0e-15*R₀*ωc,  # 1000 fF junction capacitance
    jpa.J1.r    => 1.0e8/R₀,          # junction shunt resistance
    jpa.J1.α    => 1.0,               # critical current, in units of I₀
)

# probe frequencies, shared by both parts
f_vec = collect(4.5:0.001:5.0)
Ω_vec = GHz.(f_vec)

s11 = get_HB_scattering_matrix(jpa, '1', '1')[1]

# ===================================================================================
# part 1: one pump at 4.75001 GHz
# ===================================================================================
ps1 = merge(circuit_ps, Dict(
    jpa.P1.source.ω => GHz(4.75001),
    jpa.P1.source.I => 11.3e-9/I₀,
))

sys1 = HarmonicSystem(jpa, jpa.P1.source.ω, 2, determine_jacobian = true)

prob1 = HarmonicProblem(sys1, ps1, parameter_sweep = [jpa.P1.source.ω => Ω_vec])
solve!(prob1)
S11_pump = get_solution(prob1, s11, 1)

#---- small signal: analysis -----------
δU1  = perturbation_response(sys1, jpa.P1.source.I, ps1, amplitude = ps1[jpa.P1.source.I])
lin1 = LinearisedProblem(sys1, ps1, δU1, Ω_vec)
solve!(lin1)
gain1 = 20 .* log10.(abs.(get_solution(lin1, s11, 1)))


# ===================================================================================
# part 2: two pumps at 4.65001 GHz and 4.85001 GHz through the same port
# ===================================================================================

f1, f2 = 4.65001, 4.85001
sys2 = HarmonicSystem(jpa, (jpa.P1.source.ω, jpa.P1.source.ω₂), 2,
                      determine_jacobian = true, intermod_order = 3, tones = (f1*1e9, f2*1e9))

Ip = 2*1.7*0.00565e-6/I₀           # 1.7 x 5.65 nA per pump in the reference, one sided
ps2 = merge(circuit_ps, Dict(
    jpa.P1.source.ω  => GHz(f1),
    jpa.P1.source.I  => Ip,
    jpa.P1.source.ω₂ => (97/93)*GHz(f1),
    jpa.P1.source.I₂ => Ip,
))


U = zeros(length(unknowns(sys2.system)))
for frac in (0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 1.0)
    p = merge(ps2, Dict(jpa.P1.source.I => frac*Ip, jpa.P1.source.I₂ => frac*Ip))
    step = HarmonicProblem(sys2, p, U₀ = U)
    solve!(step)
    global U = real.(step.result.solution)
end

# ---- small signal gain around the two pump working point -----------------------------
δU2  = perturbation_response(sys2, jpa.P1.source.I, ps2, amplitude = ps2[jpa.P1.source.I])
lin2 = LinearisedProblem(sys2, ps2, δU2, Ω_vec, U₀ = U)
solve!(lin2)
gain2 = 20 .* log10.(abs.(get_solution(lin2, s11, 1)))
report("two pumps", gain2)

# ---- figure -------------------------------------------------------------------------
p = plot(f_vec, gain1, lw = 2, label = "one pump, 4.75 GHz",
         xlabel = "Frequency (GHz)", ylabel = "|S₁₁| (dB)", title = "JPA gain", legend = :topright)
plot!(p, f_vec, gain2, lw = 2, label = "two pumps, 4.65 and 4.85 GHz")
savefig(p, joinpath(pkgdir(JosephsonLoops), "docs", "images", "jpa.png"))
display(p)
