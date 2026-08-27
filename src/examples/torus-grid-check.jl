# Torus-grid check: build the doubly-pumped docs JPA (intermod 0, N = 2) on BOTH
# two-tone grids and compare them end to end. Run top to bottom in a fresh REPL.
#
# What to look for, in order:
#  1. Grid dispatch: the first build prints "commensurate two-tone grid: ω1:ω2 = 93:97",
#     the second prints "no tones given — using the 2D torus grid".
#  2. Working points: every ramp step should report a small max residual (< ~1e-9).
#     A "Stalled"/"MaxIters" retcode with a LARGE residual means that grid's working
#     point is garbage and its gain curve is meaningless (this was the torus failure
#     mode before the projection-weight zero-snapping fix in utils.jl).
#  3. Cross-residuals: each grid's root plugged into the OTHER grid's equations.
#     Small both ways = the two formulations agree and any gain difference is
#     truncation/snap-sized. Large = the equation sets genuinely differ (bug).
#  4. Gain curves: peaks should sit within ~0.1 dB of each other (the tone-2 snap is
#     430 Hz, so physical differences between the grids are far below plot resolution).

using JosephsonLoops
using Symbolics
using ModelingToolkit
using Revise
const jls = JosephsonLoops

loops = [["P1", "C1", "J1"]]
circuit = jls.process_netlist(loops)
model, u0, guesses = jls.build_circuit(circuit)

I₀ = jls.Φ₀/(2π*1000.0e-12)
R₀ = 10.0e3
ωc = R₀*I₀/(jls.Φ₀/2π)
Ip = 2 * 1.7 * 0.00565e-6 / I₀
base_ps(ω2val) = Dict(
    jls.P1.source.ω  => 2π*4.65001e9/ωc,
    jls.P1.source.I  => Ip,
    jls.P1.source.ω₂ => ω2val,
    jls.P1.source.I₂ => Ip,
    jls.P1.Rₙ.r      => 50.0/R₀,
    jls.C1.βc        => 100.0e-15*R₀*ωc,
    jls.J1.βc        => 1000.0e-15*R₀*ωc,
    jls.J1.r         => 1e8/R₀,
    jls.J1.α         => 1.0,
)
ps_comm  = base_ps((97/93)*2π*4.65001e9/ωc)   # the snapped value the 1D grid pins anyway
ps_torus = base_ps(2π*4.85001e9/ωc)           # torus needs no snap: exact docs value
Ω_vec = collect(2π*(4.5:0.001:5.0)*1e9/ωc)
fs = collect(4.5:0.001:5.0)

@time sys_comm = jls.HarmonicSystem(model, (jls.P1.source.ω, jls.P1.source.ω₂), 2,
    determine_jacobian = true, intermod_order = 0, tones = (4.65001e9, 4.85001e9))
@time sys_torus = jls.HarmonicSystem(model, (jls.P1.source.ω, jls.P1.source.ω₂), 2,
    determine_jacobian = true, intermod_order = 0)

# -- working points, with convergence made visible ------------------------------------
function ramp(sys, ps, name)
    sysu = jls.unknowns(sys.system)
    U = fill(0.0, length(sysu))
    println(name, " ramp:")
    for frac in (0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 1.0)
        p = copy(ps)
        p[jls.P1.source.I]  = frac * Ip
        p[jls.P1.source.I₂] = frac * Ip
        prob = jls.NonlinearProblem(sys.system, merge(Dict(sysu .=> U), p))
        sol = jls.ModelingToolkit.solve(prob, jls.NonlinearSolve.LevenbergMarquardt(); maxiters = 2000)
        U = real.(sol.u)
        println("  ", frac, ": ", sol.retcode, "  maxresid = ", round(maximum(abs.(sol.resid)), sigdigits=3))
    end
    return Dict(zip(sysu, U)), U
end
wp_comm, U_comm  = ramp(sys_comm,  ps_comm,  "commensurate")
wp_torus, U_torus = ramp(sys_torus, ps_torus, "torus")

# -- cross-residuals: each root in the other system's equations -----------------------
function resid_at(sys, ps, wp)
    sysu = jls.unknowns(sys.system)
    prob = jls.NonlinearProblem(sys.system, merge(Dict(u => get(wp, u, 0.0) for u in sysu), ps))
    return maximum(abs.(prob.f(prob.u0, prob.p)))
end
println("cross-residuals (own root should be tiny; other root small = same physics):")
println("  comm @ own root:   ", round(resid_at(sys_comm,  ps_comm,  wp_comm),  sigdigits=3))
println("  comm @ torus root: ", round(resid_at(sys_comm,  ps_comm,  wp_torus), sigdigits=3))
println("  torus @ own root:  ", round(resid_at(sys_torus, ps_torus, wp_torus), sigdigits=3))
println("  torus @ comm root: ", round(resid_at(sys_torus, ps_torus, wp_comm),  sigdigits=3))

# -- linearised gain on both grids ----------------------------------------------------
function gain_curve(sys, ps, Upump)
    δU  = jls.perturbation_response(sys, jls.P1.source.I, ps, amplitude = ps[jls.P1.source.I])
    lin = jls.LinearisedProblem(sys, ps, δU, Ω_vec, U₀ = Upump)
    jls.solve!(lin)
    I_sig = jls.get_solution(lin, jls.P1.i, 1) .* I₀
    V_sig = jls.get_solution(lin, jls.P1.dφ, 1) .* R₀ .* I₀
    Z0 = 50.0
    a = @. 0.5*(V_sig + Z0*I_sig)/sqrt(Z0)
    b = @. 0.5*(V_sig - Z0*I_sig)/sqrt(Z0)
    return 20 .* log10.(abs.(b ./ a))
end
g_comm  = gain_curve(sys_comm,  ps_comm,  U_comm)
g_torus = gain_curve(sys_torus, ps_torus, U_torus)

println("commensurate: peak ", round(maximum(g_comm),  digits=3), " dB @ ", fs[argmax(g_comm)],  " GHz")
println("torus:        peak ", round(maximum(g_torus), digits=3), " dB @ ", fs[argmax(g_torus)], " GHz")
D = g_torus .- g_comm
println("grid vs grid: max|Δ| = ", round(maximum(abs.(D)), sigdigits=3),
        " dB, median |Δ| = ", round(sort(abs.(D))[end÷2+1], sigdigits=3), " dB")

p = jls.plot(fs, g_comm, lw = 2, label = "commensurate (93:97)",
             xlabel = "Frequency (GHz)", ylabel = "|S₁₁| (dB)",
             title = "Two-tone grid check — intermod 0")
jls.plot!(p, fs, g_torus, lw = 2, ls = :dash, label = "2D torus")
display(p)
