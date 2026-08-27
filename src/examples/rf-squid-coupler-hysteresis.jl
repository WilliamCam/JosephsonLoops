# Hysteretic rf-SQUID coupler (betaL > 1) — a capability demo AGAINST
# JosephsonCircuits.jl, on the circuit of arXiv:2408.07861 Fig. 8.
#
# For betaL = L_loop/LJ0 > 1 the rf-SQUID is HYSTERETIC: in a window around each
# half-integer flux two stable junction states coexist, and the state you observe
# depends on sweep direction. Analytically the folds sit where cos(phi*) = -1/betaL:
# for betaL = 2, phi* = 2pi/3, f_fold = (phi* + betaL*sin(phi*))/2pi = 0.391 / 0.609
# around each half-integer of effective flux.
#
# JosephsonLoops traces the physical loop because the working point is CONTINUED —
# each flux step starts from the previous solution, so the solver stays on-branch
# until the fold and then switches, differently on the up and down sweeps.
# JosephsonCircuits.jl's hbsolve cold-starts every call (x0 = nothing is hardcoded),
# so it has no state memory and CANNOT show hysteresis — see the companion probe
# mit_rf_squid_hysteresis.jl, which runs the identical up/down sweep there.
#
# The sweep spans -0.25..1.25 Phi0 (effective 0.25..1.75 with the parity offset),
# also testing flux excursions beyond one quantum.

using JosephsonLoops
using Symbolics
using ModelingToolkit
using Revise
const jls = JosephsonLoops

loops = [["P1", "L1"], ["L1", "J1", "L2"], ["L2", "P2"]]
circuit = jls.process_netlist(loops, ext_flux = [false, true, false])
model, u0, guesses = jls.build_circuit(circuit)

Ic  = 5.0e-6
I₀  = Ic
R₀  = 50.0
ωc  = R₀*I₀/(jls.Φ₀/2π)
f_probe = 5.0e9
Ω_probe = 2π*f_probe/ωc

βL = 2.0                           # > 1: hysteretic
make_ps(f_ext) = Dict(
    jls.P1.source.ω => Ω_probe,
    jls.P1.source.I => 0.0,
    jls.P2.source.ω => Ω_probe,
    jls.P2.source.I => 0.0,
    jls.P1.Rₙ.r     => 50.0/R₀,
    jls.P2.Rₙ.r     => 50.0/R₀,
    jls.J1.r        => 1.0e6/R₀,
    jls.J1.βc       => 1.0e-15*R₀*ωc,
    jls.J1.α        => 1.0,
    jls.L1.βL       => 0.5*βL,
    jls.L2.βL       => 0.5*βL,
    jls.Φₑ2.Φₑ      => 2π*(f_ext + 0.5),   # same parity compensation as the coupler example
)

@time sys_hb = jls.HarmonicSystem(model, jls.P1.source.ω, 1, determine_jacobian = true)
sysu = jls.unknowns(sys_hb.system)
δU = jls.perturbation_response(sys_hb, jls.P1.source.I, make_ps(0.0), amplitude = 1.0)
Z0 = 50.0

# At betaL > 1 no single solver is reliable (measured: LM falls into a least-squares
# trap JᵀR = 0, R ≠ 0 that even flux-ramping cannot escape, while Newton — which fails
# at betaL ≤ 1.2 — reaches machine precision at 1.5/2.0). So: solver cascade, and only
# a point with a tiny true residual is ever accepted into the continuation — dragging
# one unconverged state forward poisons every later warm start.
function solve_wp(ps, U0)
    for alg in (jls.NonlinearSolve.NewtonRaphson(), jls.NonlinearSolve.LevenbergMarquardt(),
                jls.NonlinearSolve.TrustRegion())
        prob = jls.NonlinearProblem(sys_hb.system, merge(Dict(sysu .=> U0), ps))
        sol = jls.ModelingToolkit.solve(prob, alg; maxiters = 2000)
        maximum(abs.(sol.resid)) < 1e-9 && return real.(sol.u), true
    end
    return U0, false
end

function sweep_s21!(f_path, Uwp)
    out = fill(NaN, length(f_path))
    for (jf, f_ext) in enumerate(f_path)
        ps = make_ps(f_ext)
        U, ok = solve_wp(ps, Uwp)
        ok || continue                     # leave a gap, keep the last GOOD state as guess
        Uwp = U                            # continuation: next point starts on this branch

        lin = jls.LinearisedProblem(sys_hb, ps, δU, [Ω_probe], U₀ = Uwp)
        jls.solve!(lin)
        V₁ = jls.get_solution(lin, jls.P1.dφ, 1)[1] * R₀ * I₀
        I₁ = jls.get_solution(lin, jls.P1.i, 1)[1] * I₀
        V₂ = jls.get_solution(lin, jls.P2.dφ, 1)[1] * R₀ * I₀
        a₁ = 0.5*(V₁ + Z0*I₁)/sqrt(Z0)
        out[jf] = 20*log10(abs(V₂/sqrt(Z0)) / abs(a₁))
    end
    return out, Uwp
end

f_up   = collect(-0.25:0.005:1.25)
f_down = reverse(f_up)

S21_up, Uend = sweep_s21!(f_up, fill(0.0, length(sysu)))
S21_down, _  = sweep_s21!(f_down, Uend)   # down-sweep continues from the up-sweep's end state

# where the two branches disagree = the measured hysteresis windows
# (expected for betaL = 2: f in [-0.109, 0.109] and [0.891, 1.109], from cosφ* = -1/βL)
println("unconverged points: ", count(isnan, S21_up), " up, ", count(isnan, S21_down), " down")
Δ = abs.(S21_up .- reverse(S21_down))
open_pts = f_up[.!isnan.(Δ) .& (Δ .> 1.0)]
println("hysteresis opens over ", length(open_pts), " flux points",
        isempty(open_pts) ? "" : " (f = $(minimum(open_pts)) .. $(maximum(open_pts)))")

p = jls.plot(xlabel = "Norm_Ext_Flux", ylabel = "dB(S(2,1))",
             title = "rf-SQUID coupler, βL = $(βL) — hysteresis",
             ylim = (-80, 0), legend = :bottomright)
jls.plot!(p, f_up, S21_up, lw = 2, label = "sweep up")
jls.plot!(p, reverse(f_down), reverse(S21_down), lw = 2, ls = :dash, label = "sweep down")
display(p)
# regenerates the figure used in the docs
jls.savefig(p, joinpath(pkgdir(jls), "docs", "images", "rf-squid-hysteresis.png"))
