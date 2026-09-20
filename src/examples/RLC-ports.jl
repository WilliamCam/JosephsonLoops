# Josephson parametric amplifier: a junction shunted by a coupling capacitor, driven and
# read out through a single 50 ohm port. This is the single tone benchmark circuit.
#
# The example runs the two harmonic balance problems the package provides, on one circuit:
#
#   HarmonicProblem   the circuit's own steady state response to its drive, swept over the
#                     drive frequency, which gives the resonance of the loaded junction;
#   LinearisedProblem the response to a weak probe on top of a fixed pump, which is what an
#                     amplifier's gain is.
#
# At a pump of 11.3 nA just above the resonance the reflected probe comes back amplified.
# The reference curve is in mit_jpa.csv, written by mit_jpa_export.jl, and is overlaid below
# if that file is present.
#
# The script also integrates the same circuit in the time domain, which is a sanity check that
# the model solves as an initial value problem. RLC.jl makes that comparison quantitative.

using JosephsonLoops
using Symbolics
using DelimitedFiles
const jls = JosephsonLoops

# ---- circuit ----------------------------------------------------------------------
loops = [["P1", "C1", "J1"]]
circuit = jls.process_netlist(loops)
model, u0, guesses = jls.build_circuit(circuit)

# ---- normalisation and parameters --------------------------------------------------
I₀ = jls.Φ₀/(2π*1000.0e-12)          # critical current of a 1000 pH junction, so α = 1
R₀ = 10.0e3
ωc = R₀*I₀/(jls.Φ₀/2π)
Z0 = 50.0
# The 10 kHz offset from 4.75 GHz is the reference's own pump frequency. The signal sweep
# below steps in 1 MHz, so 4.750 GHz is exactly a grid point, and a pump sitting on one
# makes signal, pump and idler degenerate there.
f_pump = 4.75001e9
# The reference drives with 5.65 nA as a one sided spectral amplitude, which is half the
# peak amplitude of this source, so the same drive is 11.3 nA here.
I_pump = 11.3e-9/I₀

ps = Dict{Num,Float64}(
    jls.P1.source.ω => 2π*f_pump/ωc,
    jls.P1.source.I => I_pump,
    jls.P1.Rₙ.r     => 50.0/R₀,
    jls.C1.βc       => 100.0e-15*R₀*ωc,
    jls.J1.βc       => 1000.0e-15*R₀*ωc,
    jls.J1.r        => 1.0e8/R₀,
    jls.J1.α        => 1.0,
)

# ---- time domain check ----------------------------------------------------------------

ps_td = merge(ps, Dict(jls.P1.source.ω => 2π*100e6/ωc,
                       jls.P1.source.I => 0.00565e-6/I₀))
tspan = (0.0, 1e-6) .* ωc
tsol  = jls.tsolve(model, guesses, ps_td, tspan; guesses = guesses)
p_td  = jls.plot(tsol[jls.C1.i][end-400:end] .* I₀, legend = false,
                 xlabel = "Sample", ylabel = "I(C1) (A)", title = "Transient at 100 MHz")
display(p_td)

f_vec = collect(4.5:0.001:5.0)
ω_vec = collect(2π .* f_vec .* 1e9 ./ ωc)

# ---- the harmonic system, built once ------------------------------------------------
# determine_jacobian is needed by the linearised problem below, and the same system serves
# the steady state sweep, so it is built once. This is the expensive symbolic step.
@time sys = jls.HarmonicSystem(model, jls.P1.source.ω, 2, determine_jacobian = true)

# ---- steady state: sweep the drive frequency ----------------------------------------

# The swept parameter must be absent from the fixed dict, and solve! continues from the
# previous point. S11 comes from the package's own port expression rather than by hand.
sweep_ps = delete!(copy(ps), jls.P1.source.ω)
prob = jls.HarmonicProblem(sys, sweep_ps, parameter_sweep = [jls.P1.source.ω => ω_vec])
jls.solve!(prob)

# The package builds the port's reflection expression, so the wave algebra is not repeated
# here.
s11_expr  = jls.get_HB_scattering_matrix(model, '1', '1')[1]
S11_drive = jls.get_solution(prob, s11_expr, (1, 0))

# The only loss in this circuit is the junction's 100 MΩ shunt, so a driven one port must
# reflect everything. |S11| = 0 dB across the band is therefore a conservation check on the
# steady state solve, not a feature of the circuit.
println("steady state sweep: |S11| departs from 0 dB by at most ",
        round(maximum(abs.(20*log10.(abs.(S11_drive)))), sigdigits = 2), " dB across the band")

# the resonance itself shows in how hard the drive swings the junction
φ1 = abs.(jls.get_solution(prob, jls.J1.φ, 1))
println("junction phase swing peaks at ", round(maximum(φ1), sigdigits = 4), " rad at ",
        f_vec[argmax(φ1)], " GHz, the resonance of the loaded junction")

# ---- small signal: gain around the pump ----------------------------------------------
# The drive derivative gives the injection vector. The amplitude must match the pump, or the
# incident wave is undercounted and |S11| comes out above 0 dB even with the pump off.
δU  = jls.perturbation_response(sys, jls.P1.source.I, ps, amplitude = ps[jls.P1.source.I])
lin = jls.LinearisedProblem(sys, ps, δU, ω_vec)
jls.solve!(lin)

I_sig = jls.get_solution(lin, jls.P1.i, 1)  .* I₀
V_sig = jls.get_solution(lin, jls.P1.dφ, 1) .* (R₀*I₀)
a = @. 0.5*(V_sig + Z0*I_sig)/sqrt(Z0)
b = @. 0.5*(V_sig - Z0*I_sig)/sqrt(Z0)
gain_dB = 20 .* log10.(abs.(b ./ a))

bw = f_vec[gain_dB .>= maximum(gain_dB) - 3]
println("JosephsonLoops: peak ", round(maximum(gain_dB), digits = 3), " dB @ ",
        round(f_vec[argmax(gain_dB)], digits = 4), " GHz, 3 dB width ",
        round(1000*(maximum(bw) - minimum(bw)), digits = 1), " MHz")

# ---- figure ------------------------------------------------------------------------
p_gain = jls.plot(f_vec, gain_dB, lw = 2, label = false,
                  xlabel = "Frequency (GHz)", ylabel = "|S₁₁| (dB)",
                  title = "JPA gain")
# pkgdir resolves to the repository root whatever the working directory is
mit_csv = joinpath(pkgdir(jls), "mit_jpa.csv")
if isfile(mit_csv)
    mit = readdlm(mit_csv, ',')
    Δ = gain_dB .- Float64.(mit[:, 2])
    println("JosephsonCircuits.jl: peak ", round(maximum(Float64.(mit[:, 2])), digits = 3), " dB")
    println("band max|Δ| = ", round(maximum(abs.(Δ)), digits = 3), " dB, median |Δ| = ",
            round(sort(abs.(Δ))[end÷2+1], sigdigits = 3), " dB")
else
    @warn "reference not found at $mit_csv: generate it with mit_jpa_export.jl"
end
display(p_gain)
jls.savefig(p_gain, joinpath(pkgdir(jls), "docs", "images", "jpa-single-tone.png"))
