# Nonlinear RLC resonator: the same circuit solved twice, in the time domain and by
# harmonic balance, so the two solvers can be checked against each other.
#
# A current source drives a 50 ohm shunt resistor in parallel with a series branch of a
# coupling capacitor and a Josephson junction. The junction carries its own capacitance, so
# the branch has a series resonance where the capacitor cancels the junction inductance and
# a parallel resonance where the junction and its capacitance cancel each other. Both fall
# inside the swept band.
#
# The drive is small, so the junction stays in its linear regime and the two solvers must
# agree. This is the simplest check that a harmonic balance result is trustworthy: solve the
# same circuit as an initial value problem, let the transient die away, and measure the
# amplitude that remains.

using JosephsonLoops
const jls = JosephsonLoops

# ---- circuit ----------------------------------------------------------------------
# Loop 1 is the source and the resistor. Loop 2 is the resistor, the capacitor and the
# junction, so the resistor is the branch the two loops share: the source sees the resistor
# in parallel with the series capacitor and junction.
loops = [["I1", "R1"],
         ["R1", "C1", "J1"]]
circuit = jls.process_netlist(loops)
model, u0, guesses = jls.build_circuit(circuit)

# ---- normalisation ----------------------------------------------------------------
Lj = 1000.0e-12                      # junction inductance, so α = 1 below
I₀ = jls.Φ₀/(2π*Lj)                  # 0.33 uA
R₀ = 50.0
ωc = R₀*I₀/(jls.Φ₀/2π)

f_drive = 4.8e9                      # the series resonance of the branch
ps = Dict(
    jls.I1.ω  => 2π*f_drive/ωc,
    jls.I1.I  => 1.0e-9/I₀,          # 1 nA, small enough to stay linear
    jls.R1.r  => 50.0/R₀,
    jls.C1.βc => 100.0e-15*R₀*ωc,
    jls.J1.βc => 1000.0e-15*R₀*ωc,
    jls.J1.r  => 1.0e8/R₀,
    jls.J1.α  => 1.0,
)

# ---- harmonic balance: sweep the drive frequency -----------------------------------
# The swept parameter has to leave the fixed parameter dict, and solve! continues from the
# previous point, so the sweep follows one branch.
sys = jls.HarmonicSystem(model, jls.I1.ω, 2)

f_vec = collect(4.5:0.002:5.2)
ω_vec = collect(2π .* f_vec .* 1e9 ./ ωc)
sweep_ps = delete!(copy(ps), jls.I1.ω)
prob = jls.HarmonicProblem(sys, sweep_ps, parameter_sweep = [jls.I1.ω => ω_vec])
jls.solve!(prob)

# first harmonic of the current through the capacitor, in amps
I_C1 = jls.get_solution(prob, jls.C1.i, 1) .* I₀
amp_hb = abs.(I_C1)

# ---- time domain: the same circuit at one frequency --------------------------------
# tspan is in normalised time, so a physical duration is multiplied by ωc. The last tenth
# of the window is taken as the steady state.
tspan = (0.0, 60.0e-9) .* ωc
tsol = jls.tsolve(model, guesses, ps, tspan; guesses = guesses, saveat = LinRange(tspan..., 20_000))
i_td = tsol[jls.C1.i] .* I₀
tail = i_td[(9*length(i_td)) ÷ 10 : end]   # the last tenth, after the transient has died
amp_td = (maximum(tail) - minimum(tail))/2

k = argmin(abs.(f_vec .- f_drive/1e9))
println("at ", f_vec[k], " GHz: harmonic balance ", round(amp_hb[k]*1e9, digits = 4),
        " nA, time domain ", round(amp_td*1e9, digits = 4), " nA, difference ",
        round(100*abs(amp_hb[k] - amp_td)/amp_td, digits = 2), " percent")
println("series resonance at ", f_vec[argmax(amp_hb)], " GHz, where the capacitor cancels the junction ",
        "inductance and the branch carries the most current; parallel resonance at ", f_vec[argmin(amp_hb)],
        " GHz, where the junction and its own capacitance cancel and the branch carries the least")

# ---- figure -------------------------------------------------------------------------
p = jls.plot(f_vec, amp_hb .* 1e9, lw = 2, label = "harmonic balance",
             xlabel = "Drive frequency (GHz)", ylabel = "|I(C1)| at the drive (nA)",
             title = "Nonlinear RLC resonator", legend = :topleft, size = (760, 420))
jls.scatter!(p, [f_vec[k]], [amp_td*1e9], color = :black, ms = 5,
             label = "time domain at $(f_drive/1e9) GHz")
display(p)
jls.savefig(p, joinpath(pkgdir(jls), "docs", "images", "rlc-resonator.png"))
