# One loop, one capacitor, one inductor: is the inductor's sign convention right?
#
# This is the smallest circuit that can answer the question, and it needs no junction, no
# flux and no harmonic balance subtleties. A port drives a capacitor and an inductor in
# series. Such a loop has one series resonance, at f0 = 1/(2π√(LC)), where the capacitor's
# reactance cancels the inductor's and the branch carries the whole source current.
#
# The component library defines an inductor as
#
#     in.Φ ~ βL*i
#
# while every Branch element (junction, capacitor, resistor) takes its phase as φ ~ out.Φ
# with i ~ in.iₘ - out.iₘ. If those two parities disagree then the inductor delivers its flux
# to the loop with the opposite sign to everything else, and the loop behaves as though the
# inductance were negative: the two reactances add instead of cancelling, and the resonance
# never happens.
#
# The test needs no change to the library, because negating the parameter is exactly the same
# thing as flipping the sign in the equation. So the same circuit is solved twice, once with
# βL = +L/Lj and once with βL = -L/Lj, and both are compared with the textbook series
# resonance. Whichever curve matches identifies the correct convention.
#
# This matters beyond tidiness. Under the convention that does not match, any loop holding
# both an inductor and a junction settles at the wrong operating point, which is what the
# rf-SQUID coupler's half flux quantum offset compensates for.

using JosephsonLoops
using Symbolics
const jls = JosephsonLoops

# ---- circuit: port, capacitor and inductor in one loop ---------------------------------
loops = [["P1", "C1", "L1"]]
circuit = jls.process_netlist(loops)
model, u0, guesses = jls.build_circuit(circuit)

# ---- parameters ------------------------------------------------------------------------
L  = 10.0e-9
C  = 10.0e-15
Z0 = 50.0
f0 = 1/(2π*sqrt(L*C))                # 15.9 GHz, the series resonance

I₀ = jls.Φ₀/(2π*L)                   # so that βL = 1 for this inductor
R₀ = Z0
ωc = R₀*I₀/(jls.Φ₀/2π)
I_drive = 1.0e-9/I₀

f_vec = collect(10.0:0.05:22.0)
ω_vec = collect(2π .* f_vec .* 1e9 ./ ωc)

make_ps(βL_sign) = Dict{Num,Float64}(
    jls.P1.source.I => I_drive,
    jls.P1.Rₙ.r     => Z0/R₀,
    jls.C1.βc       => C*R₀*ωc,
    jls.L1.βL       => βL_sign * L/(jls.Φ₀/(2π*I₀)),
)

# ---- the textbook answer ----------------------------------------------------------------
# The port is a current source in parallel with Z0, driving the series branch, so the current
# in the loop is I_source * Z0/(Z0 + Z). With the inductor entering the loop backwards the
# branch impedance is -jωL + 1/(jωC) instead of jωL + 1/(jωC): the two reactances then add,
# so the magnitude has no zero and the resonance disappears.
Z_series(f)  = im*2π*f*L + 1/(im*2π*f*C)
Z_flipped(f) = -im*2π*f*L + 1/(im*2π*f*C)
loop_current(Zf, f) = abs(Z0 / (Z0 + Zf(f)))          # normalised to the source current

exact_series  = [loop_current(Z_series,  f*1e9) for f in f_vec]
exact_flipped = [loop_current(Z_flipped, f*1e9) for f in f_vec]

# ---- solve the same circuit under each convention ---------------------------------------
sys = jls.HarmonicSystem(model, jls.P1.source.ω, 1)

sweep_pos = delete!(copy(make_ps(+1.0)), jls.P1.source.ω)
prob_pos  = jls.HarmonicProblem(sys, sweep_pos, parameter_sweep = [jls.P1.source.ω => ω_vec])
jls.solve!(prob_pos)
i_pos = abs.(jls.get_solution(prob_pos, jls.C1.i, 1)) ./ I_drive

sweep_neg = delete!(copy(make_ps(-1.0)), jls.P1.source.ω)
prob_neg  = jls.HarmonicProblem(sys, sweep_neg, parameter_sweep = [jls.P1.source.ω => ω_vec])
jls.solve!(prob_neg)
i_neg = abs.(jls.get_solution(prob_neg, jls.C1.i, 1)) ./ I_drive

# ---- which one is the physical series resonance? ----------------------------------------
err(sim, exact) = maximum(abs.(sim .- exact)) / maximum(exact)
println("series resonance of L = ", L*1e9, " nH with C = ", C*1e15, " fF is at ",
        round(f0/1e9, digits = 2), " GHz")
println("βL = +L/Lj  (the library as written): peak ", round(maximum(i_pos), sigdigits = 4),
        " of the source current at ", f_vec[argmax(i_pos)], " GHz")
println("βL = -L/Lj  (the sign flipped)      : peak ", round(maximum(i_neg), sigdigits = 4),
        " of the source current at ", f_vec[argmax(i_neg)], " GHz")
println("against the textbook series resonance, worst deviation over the band:")
println("    βL = +L/Lj : ", round(100*err(i_pos, exact_series), sigdigits = 3), " percent")
println("    βL = -L/Lj : ", round(100*err(i_neg, exact_series), sigdigits = 3), " percent")
println("against the textbook curve for an inductor entering backwards:")
println("    βL = +L/Lj : ", round(100*err(i_pos, exact_flipped), sigdigits = 3), " percent")
println("    βL = -L/Lj : ", round(100*err(i_neg, exact_flipped), sigdigits = 3), " percent")

# ---- figure -------------------------------------------------------------------------------
p = jls.plot(f_vec, exact_series, lw = 3, color = :black, label = "textbook series resonance",
             xlabel = "Frequency (GHz)", ylabel = "|I| in the loop / source current",
             title = "One loop, one capacitor, one inductor", legend = :topright)
jls.plot!(p, f_vec, exact_flipped, lw = 3, ls = :dot, color = :gray,
          label = "textbook, inductor entering backwards")
jls.plot!(p, f_vec, i_pos, lw = 2, ls = :dash, label = "JosephsonLoops, βL = +L/Lj")
jls.plot!(p, f_vec, i_neg, lw = 2, ls = :dash, label = "JosephsonLoops, βL = -L/Lj")
display(p)
jls.savefig(p, joinpath(pkgdir(jls), "docs", "images", "inductor-sign.png"))
