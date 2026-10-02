 
using JosephsonLoops, ModelingToolkit, Plots
@parameters ωHB

# One loop inductance, two junctions, DC bias, external flux, and voltage readout.
loops = [["Lt", "J1", "J2"], ["J2", "Idc"], ["Idc", "R2"]]
squid, _, guesses = build_circuit(process_netlist(loops, ext_flux=[true, false, false]))

Φ0 = JosephsonLoops.Φ₀
Ic, Rj, Lsq = 6.3e-6, 6.0, 200e-12
ωc = Rj * Ic / (Φ0 / 2π)
βC_paper = 1.0 # Assumed: the paper gives βC ≈ 1 as optimum guidance, not device Cj.
Cj_assumed = βC_paper * Φ0 / (π * Ic * Rj^2)

ps = Dict{Num,Float64}(
    squid.Idc.I => 1.8, squid.Idc.ω => 0.0,
    squid.Lt.βL => 2π * Lsq * Ic / Φ0, # Model uses flux units Φ0/(2π); Lt is total Lsq.
    squid.J1.βc => 2βC_paper, squid.J2.βc => 2βC_paper,
    squid.J1.r => 1.0, squid.J2.r => 1.0, # Rj/Rj
    squid.R2.r => 1.0, # Assumed readout load Rload = Rj; not specified in the paper.
    squid.Φₑ1.Φₑ => π / 2, # Φ0/4 in model units.
)

# Average voltage after transients for small static flux offsets estimates the
# quasistatic responsivity without resolving the GHz carrier as a probe.
δφ = 1e-3
function mean_voltage(φₑ)
    p = merge(ps, Dict{Num,Float64}(squid.Φₑ1.Φₑ => φₑ))
    sol = tsolve(squid, guesses, p, (0.0, 10e-9) .* ωc; guesses)
    first_fit = cld(length(sol.t), 2)
    v = p[squid.R2.r] .* sol[squid.R2.i] .* (Rj * Ic)
    return sum(v[first_fit:end]) / length(v[first_fit:end])
end

V₋, V₊ = mean_voltage(π / 2 - δφ), mean_voltage(π / 2 + δφ)
model_μV_per_Φ0 = abs(V₊ - V₋) / (2δφ) * 2π * 1e6
theory_μV_per_Φ0 = Rj / Lsq * Φ0 * 1e6
βL_paper = Lsq * Ic / Φ0


# Linearise about the running state and sweep upper sidebands from the carrier
# through 1 GHz offset. This is flux-to-voltage responsivity in V/Φ₀, not power gain.
tspan = (0.0, 10e-9) .* ωc
tsol = tsolve(squid, guesses, ps, tspan; guesses)
function late_phase_rate(t, φ)
    i = cld(length(t), 2)
    x, y = t[i:end], φ[i:end]
    xc, yc = x .- sum(x) / length(x), y .- sum(y) / length(y)
    sum(xc .* yc) / sum(abs2, xc)
end
ω0_guess = (abs(late_phase_rate(tsol.t, tsol[squid.J1.φ])) +
            abs(late_phase_rate(tsol.t, tsol[squid.J2.φ]))) / 2

sys = HarmonicSystem(squid, ωHB, 2;
    autonomous=true,
    drifting_states=Dict(squid.J1.φ => 1, squid.J2.φ => -1),
    determine_jacobian=true)
prob = HarmonicProblem(sys, ps)
iω0 = JosephsonLoops.var_index(unknowns(sys.system), sys.autonomous_frequency)
prob.U₀[iω0] = ω0_guess
U = real.(JosephsonLoops.solve!(prob))
ω0 = U[iω0]
residual = maximum(abs.(prob.problem.f(U, prob.problem.p)))
@assert residual < 1e-8 "Autonomous HB did not converge"

δU = perturbation_response(sys, squid.Φₑ1.Φₑ, ps, amplitude=1.0)
f_gain = collect(range(0.0, 1.0e9, length=101))
Ω_gain = ω0 .+ 2π .* f_gain ./ ωc
lin = LinearisedProblem(sys, ps, δU, Ω_gain; U₀=U)
JosephsonLoops.solve!(lin)
V_per_Φ0 = 2π * Rj * Ic .* get_solution(lin, squid.R2.r * squid.R2.i, 1)
@assert all(isfinite, real.(V_per_Φ0)) && all(isfinite, imag.(V_per_Φ0)) "Non-finite flux response"

p_gain = plot(f_gain ./ 1e9, 20 .* log10.(abs.(V_per_Φ0 ./ 1e-6)), lw=2,
    xlabel="Upper-sideband offset from carrier (GHz)",
    ylabel="Flux-to-voltage response (dB re 1 μV/Φ₀)",
    title="DC-SQUID small-signal response, 0–1 GHz offset", label=false)
display(plot(p_static, p_gain, layout=(2, 1), size=(800, 700)))
