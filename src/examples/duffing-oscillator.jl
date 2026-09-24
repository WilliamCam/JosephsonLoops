# Driven Duffing oscillator, solved by collocation harmonic balance.
#
#     ẍ + γẋ + ω₀²x + αx³ + ηẋx² = F cos(ωt)
#
# Nothing here is a circuit. The harmonic balance backend takes any ModelingToolkit system,
# so the same code path that solves a Josephson circuit solves this. It is the smallest
# demonstration of the machinery, with no netlist in the way.
#
# The example produces two plots:
#
#   1. the amplitude response, the DC component plus the size of the fundamental, against
#      the drive frequency. Because `solve!` continues from the previous point, the sweep
#      follows one branch of the response rather than jumping between them.
#   2. the small signal response around one point of that sweep. A weak probe is applied on
#      top of the drive and the oscillator's response is read at each probe frequency, which
#      is the same linearisation the amplifier examples use for gain.
#
# Reference: Kosata, "Harmonic balance methods for nonlinear oscillators" (2022), chapter 5.

using JosephsonLoops
using ModelingToolkit
using Symbolics

# ---- the oscillator as a ModelingToolkit system ------------------------------------
@independent_variables t
@variables x(t)
@parameters α ω ω₀ F γ η
D = Differential(t)
duffing = D(D(x)) + ω₀^2*x + α*x^3 + η*D(x)*x^2 + γ*D(x) - F*cos(ω*t) ~ 0
@named duffing_sys = System([duffing], t)
model = mtkcompile(duffing_sys)

# one harmonic of the drive; determine_jacobian builds what the linearised problem needs
sys = HarmonicSystem(model, ω, 1, determine_jacobian = true)

ω_vec = collect(range(0.8, 1.2, 200))
# the swept parameter still needs a value, because a parameter declared with @parameters
# carries no default and the problem is built before the sweep starts
ps = Dict{Num,Float64}(α => 1.0, ω₀ => 1.0, F => 0.01, η => 1.0e-1, γ => 1.0e-3,
                       ω => first(ω_vec))

# ---- 1. amplitude response against drive frequency ---------------------------------
prob = HarmonicProblem(sys, ps, parameter_sweep = [ω => ω_vec])
solve!(prob)

# Both nonlinearities are odd, x³ and ẋx², and the drive is a single tone, so the response
# contains only odd harmonics and carries no DC component. The amplitude is the size of
# the fundamental.
amplitude = abs.(get_solution(prob, x, 1))

println("amplitude response: peak ", round(maximum(amplitude), sigdigits = 4), " at ω = ",
        round(ω_vec[argmax(amplitude)], digits = 4), ", from ", round(minimum(amplitude), sigdigits = 3),
        " at the edges of the sweep")

# ---- 2. small signal response around one point of that sweep -----------------------
# The working point is taken from the sweep, so the probe sees the oscillator as the drive
# has left it. The probe enters through the drive amplitude, the same way a pump does in the
# circuit examples.
j  = 70
ωp = ω_vec[j]
ps_wp = merge(ps, Dict(ω => ωp))
U₀ = real.(prob.result.solution[:, j])

δU = perturbation_response(sys, F, ps_wp, amplitude = 1.0e-3)
Ω  = collect(range(0.8, 1.2, 800))
lin = LinearisedProblem(sys, ps_wp, δU, Ω, U₀ = U₀)
solve!(lin)

# the sideband amplitude at Ω is |A + iB|/2
response = abs.(get_solution(lin, x, 1)) ./ 2

println("small signal around ω = ", round(ωp, digits = 4), ": response peaks at Ω = ",
        round(Ω[argmax(response)], digits = 4), ", ",
        round(maximum(response)/minimum(response), digits = 1), " times its smallest value")

# ---- figure -------------------------------------------------------------------------
p1 = plot(ω_vec, amplitude, lw = 2, label = false,
              xlabel = "Drive frequency ω", ylabel = "|fundamental|",
              title = "Amplitude response, F = $(ps[F])")
vline!(p1, [ωp], ls = :dash, color = :black, label = "working point for the panel below")
p2 = plot(Ω, response, lw = 2, label = false,
              xlabel = "Probe frequency Ω", ylabel = "|x(Ω)| / 2",
              title = "Small signal response around ω = $(round(ωp, digits = 3))")
p = plot(p1, p2, layout = (2, 1), size = (760, 680), left_margin = 8Plots.mm)
display(p)
