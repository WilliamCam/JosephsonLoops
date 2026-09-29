# Josephson parametric amplifier, pumped once and pumped twice.
#
# The circuit is a 50 ohm port, a 100 fF coupling capacitor and a 1000 pH junction shunted by
# 1000 fF, all in one loop. It is the amplifier from the JosephsonCircuits.jl documentation,
# so the numbers printed here can be compared with that package directly.
#
# Part 1 drives it with one pump just above its resonance and reads the small signal gain
# seen by a weak probe swept across the pump. Part 2 drives the same circuit with two pumps,
# both entering through the single port as the second tone of its current source, and reads
# the gain again. Both parts follow the same two steps:
#
#   HarmonicSystem      the Fourier ansatz, the collocation grid and the jacobians, built once
#   LinearisedProblem   the probe response around the pumped working point, which gives S11
#
# The two pump case adds a working point ramp, because a strongly driven circuit has more
# than one steady state and a cold solve can land on the trivial one.


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

# The reference quotes one sided spectral amplitudes, which are half the peak amplitude of
# an I*sin(ωt) source, so its 5.65 nA single pump is 11.3 nA here.
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

# reflection at port 1 as a symbolic expression in the port voltage and current: the same
# expression serves the steady state sweep and the linearised response
s11 = get_HB_scattering_matrix(jpa, '1', '1')[1]

function report(label, gain_dB)
    k = argmax(gain_dB)
    inband = f_vec[gain_dB .>= gain_dB[k] - 3]
    println(label, ": peak ", round(gain_dB[k], digits = 3), " dB at ", f_vec[k], " GHz, 3 dB bandwidth ",
            round((inband[end] - inband[1])*1e3, digits = 1), " MHz")
end

# ===================================================================================
# part 1: one pump at 4.75001 GHz
# ===================================================================================
ps1 = merge(circuit_ps, Dict(
    jpa.P1.source.ω => GHz(4.75001),
    jpa.P1.source.I => 11.3e-9/I₀,
))

# two harmonics of the pump: DC, the fundamental and the second harmonic. N = 3 moves the
# peak from 13.12 dB to 13.30 dB, within 0.006 dB of the reference.
sys1 = HarmonicSystem(jpa, jpa.P1.source.ω, 2, determine_jacobian = true)

# ---- steady state: the reflection of the pump itself, swept over the pump frequency ----
# The swept parameter is removed from the fixed parameters; delete! mutates, so work on a copy.
sweep_ps = delete!(copy(ps1), jpa.P1.source.ω)
prob1 = HarmonicProblem(sys1, sweep_ps, parameter_sweep = [jpa.P1.source.ω => Ω_vec])
solve!(prob1)
S11_pump = get_solution(prob1, s11, 1)      # order 1 is the pump frequency itself

# ---- small signal: the gain seen by a weak probe around the fixed pump ------------
δU1  = perturbation_response(sys1, jpa.P1.source.I, ps1, amplitude = ps1[jpa.P1.source.I])
lin1 = LinearisedProblem(sys1, ps1, δU1, Ω_vec)
solve!(lin1)
gain1 = 20 .* log10.(abs.(get_solution(lin1, s11, 1)))
report("one pump", gain1)

# ---- time domain check of the same working point -----------------------------------
# The pumped circuit is integrated from rest at the pump frequency and the amplitude of the
# port current in the last periods is compared with the harmonic balance fundamental.
periods = 300
T = 2π/ps1[jpa.P1.source.ω]
tsol = tsolve(jpa, guesses, ps1, (0.0, periods*T); guesses = guesses, saveat = T/64)
I_td = maximum(abs.(tsol[jpa.P1.i][end-64*10:end])) * I₀
j_pump = argmin(abs.(f_vec .- 4.75001))
I_hb = abs(get_solution(prob1, jpa.P1.i, 1)[j_pump]) * I₀
println("port current at the pump: harmonic balance ", round(I_hb*1e9, digits = 3), " nA, time domain ",
        round(I_td*1e9, digits = 3), " nA, difference ", round(100*abs(I_td - I_hb)/I_hb, digits = 2), " percent")

p_td = plot(tsol.t ./ ωc .* 1e9, tsol[jpa.P1.i] .* I₀ .* 1e9, label = false,
            xlabel = "t (ns)", ylabel = "port current (nA)", title = "Time domain, pump at 4.75 GHz")
savefig(p_td, joinpath(pkgdir(JosephsonLoops), "docs", "images", "jpa-time-domain.png"))
display(p_td)

# ===================================================================================
# part 2: two pumps at 4.65001 GHz and 4.85001 GHz through the same port
# ===================================================================================
# The pump ratio 4.65:4.85 is 93:97, so both tones are harmonics of one base frequency near
# 50 MHz and a one dimensional collocation grid is exact. The backend reports the ratio it
# chose and moves the second tone by 430 Hz to make it exact; ω₂ is set to that value so the
# parameter and the grid agree. intermod_order = 3 admits the mixing products m*ω1 + n*ω2 with
# m + |n| <= 3, which is where the conversion between the two pumps takes place.
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

# ---- working point: ramp both pumps together, carrying the solution forward -----------
# Each step starts from the previous converged state, which keeps the solve on the driven
# branch. A HarmonicProblem without a sweep solves one point; its result is the state vector.
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
