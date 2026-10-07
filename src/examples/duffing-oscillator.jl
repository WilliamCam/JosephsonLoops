# Driven Duffing oscillator, analyzed with harmonic balance and BifurcationKit.
#
#     ẍ + γẋ + ω₀²x + αx³ + ηẋx² = F cos(ωt)
#
# Nothing here is a circuit. The harmonic balance backend takes any ModelingToolkit system,
# so the same code path that solves a Josephson circuit solves this. It is the smallest
# demonstration of the machinery, with no netlist in the way.
#
# BifurcationKit continues the compiled harmonic-balance residual of the driven system, so
# the branches are periodic responses rather than equilibria of a source-free equation.

using BifurcationKit
using JosephsonLoops
using ModelingToolkit
using NonlinearSolve
using LinearAlgebra: norm
import Plots
using Symbolics

# ---- the oscillator as a ModelingToolkit system ------------------------------------
@independent_variables t
@variables x(t)
@parameters α ω ω₀ F γ η
D = Differential(t)
duffing = D(D(x)) + ω₀^2*x + α*x^3 + η*D(x)*x^2 + γ*D(x) - F*cos(ω*t) ~ 0
@named duffing_sys = System([duffing], t)
model = mtkcompile(duffing_sys)

sys = HarmonicSystem(model, ω, 1, determine_jacobian = true)

ω_vec = collect(range(0.8, 1.2, 200))

ps = Dict{Num,Float64}(α => 1.0, ω₀ => 1.0, F => 0.01, η => 1.0e-1, γ => 1.0e-3,
                       ω => first(ω_vec))

# ---- 1. amplitude response against drive frequency ---------------------------------
prob = HarmonicProblem(sys, ps, parameter_sweep = [ω => ω_vec])
JosephsonLoops.solve!(prob)

# ---- bifurcation analysis of the driven harmonic-balance equations ------------------
# Continue the compiled residual from three converged roots at one frequency. The distinct
# initial guesses are needed because a frequency sweep follows only one of the coexisting
# solutions. The frequency is updated in the MTK parameter vector used by the residual.
ω_branch = 1.1
ω_index = findfirst(isequal(ω), parameters(sys.system))

seed_parameters = deepcopy(prob.problem.p)
seed_parameters.tunable[ω_index] = ω_branch
continuation_parameters = deepcopy(seed_parameters)

hb_residual! = (out, u, drive_frequency) -> begin
    continuation_parameters.tunable[ω_index] = drive_frequency
    prob.problem.f(out, u, continuation_parameters)
    out
end

continuation_options = ContinuationPar(
    p_min = first(ω_vec), p_max = last(ω_vec), ds = 0.0004, dsmin = 1e-7,
    dsmax = 0.002, max_steps = 3000, detect_bifurcation = 3, nev = 3,
)
hb_branches = []
seed_amplitudes = Float64[]

for guess_amplitude in (0.03, 0.3, 0.5)
    guess = zeros(size(prob.result.solution, 1))
    guess[1] = guess_amplitude
    root = NonlinearSolve.solve(
        NonlinearProblem(prob.problem.f, guess, seed_parameters),
        NewtonRaphson(),
    )
    residual = similar(root.u)
    prob.problem.f(residual, root.u, seed_parameters)
    push!(seed_amplitudes, hypot(root.u[1], root.u[2]))
    push!(hb_branches, continuation(
        BifurcationProblem(
            hb_residual!,
            root.u,
            ω_branch;
            inplace = true,
            record_from_solution = (u, p; kwargs...) -> hypot(u[1], u[2]),
        ),
        PALC(),
        continuation_options;
        bothside = true,
        plot = false,
        verbosity = 0,
    ))
end
println("Distinct driven HB seed amplitudes at ω = ", ω_branch, ": ",
        round.(seed_amplitudes, sigdigits = 5))
println("Driven harmonic-balance branch points: ",
        [[(point.type, round(point.param, digits = 6)) for point in branch.specialpoint]
         for branch in hb_branches])

# The ordinary Duffing nonlinearities are odd and the drive is a single tone, so the
# harmonic-balance response has no DC component. Its plotted amplitude is the fundamental.
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
JosephsonLoops.solve!(lin)

# the sideband amplitude at Ω is |A + iB|/2
response = abs.(get_solution(lin, x, 1)) ./ 2

println("small signal around ω = ", round(ωp, digits = 4), ": response peaks at Ω = ",
        round(Ω[argmax(response)], digits = 4), ", ",
        round(maximum(response)/minimum(response), digits = 1), " times its smallest value")

# ---- figure -------------------------------------------------------------------------
p1 = Plots.plot(ω_vec, amplitude, lw = 2, label = "frequency sweep",
          xlabel = "Drive frequency ω", ylabel = "|fundamental|",
          title = "Driven response and HB continuation, F = $(ps[F])")
for (i, branch) in enumerate(hb_branches)
    Plots.plot!(p1, branch.param, branch.branch.x, lw = 2,
                label = "continuation from A = $(round(seed_amplitudes[i], sigdigits = 3))")
end
Plots.vline!(p1, [ωp], ls = :dash, color = :black, label = "working point for the panel below")
p2 = Plots.plot(Ω, response, lw = 2, label = false,
              xlabel = "Probe frequency Ω", ylabel = "|x(Ω)| / 2",
              title = "Small signal response around ω = $(round(ωp, digits = 3))")
p = Plots.plot(p1, p2, layout = (2, 1), size = (760, 720), left_margin = 8Plots.mm)
display(p)
