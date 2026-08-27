# RF SQUID coupler (two-port), from "Modeling flux-quantizing Josephson junction
# circuits in Keysight ADS" (arXiv:2408.07861), Sec. IV.2.1 / Fig. 8:

# Plot: dB(S21) at a fixed 5 GHz probe vs normalised external flux (0..1 Φ₀),
# one curve per βL = L_loop/LJ0 ∈ 0.6..1.0. Expected (their Fig. 8b): ~-40 dB
# baseline, transmission -> 0 dB at Φext = 0.5·Φ₀, βL-dependent nulls near 0.35/0.65.
#
# No pump: the working point is the flux-biased DC state (P1.source.I = 0) and the
# probe enters only through the linearised response. βL ≤ 1 keeps the rf-SQUID
# non-hysteretic; the flux sweep still carries the working point forward as a
# continuation so the steep region around 0.5·Φ₀ stays on-branch.

using JosephsonLoops
using Symbolics
using ModelingToolkit
using Revise
const jls = JosephsonLoops

# mesh loops: port loop 1 | SQUID loop (flux-threaded) | port loop 2
loops = [["P1", "L1"], ["L1", "J1", "L2"], ["L2", "P2"]]
circuit = jls.process_netlist(loops, ext_flux = [false, true, false])
model, u0, guesses = jls.build_circuit(circuit)

Ic  = 5.0e-6                       # paper junction
I₀  = Ic                           # -> J1.α = 1 and LJ0·(2π I₀/Φ₀) = 1
R₀  = 50.0                         # -> port r = 1
ωc  = R₀*I₀/(jls.Φ₀/2π)
LJ0 = jls.Φ₀/(2π*Ic)               # 65.8 pH (paper's 32.9 pH is the 0.5·LJ0 per-side value at βL=1)
f_probe = 5.0e9
Ω_probe = 2π*f_probe/ωc

# with I₀ = Ic the nondimensional inductance of each side is exactly 0.5·βL
make_ps(βL, f_ext) = Dict(
    jls.P1.source.ω => Ω_probe,
    jls.P1.source.I => 0.0,        # no pump — flux-biased passive working point
    jls.P2.source.ω => Ω_probe,
    jls.P2.source.I => 0.0,        # P2's source default is I = 1.0: MUST be zeroed
    jls.P1.Rₙ.r     => 50.0/R₀,
    jls.P2.Rₙ.r     => 50.0/R₀,
    jls.J1.r        => 1.0e6/R₀,
    jls.J1.βc       => 1.0e-15*R₀*ωc,   # negligible junction capacitance (model needs βc > 0)
    jls.J1.α        => 1.0,
    jls.L1.βL       => 0.5*βL,
    jls.L2.βL       => 0.5*βL,
    # loop phases are in radians: 1 Φ₀ = 2π. The +0.5 compensates a half-quantum
    # offset in the mesh branch parity (at Φₑ = 0 the junction sits at φ = π, not 0:
    # verified against the paper's Fig. 8b values AND an analytic ABCD calculation of
    # the Pi network — every curve maps under f -> f + 0.5). Library-level sign
    # convention question for Will (junction i ~ in.iₘ - out.iₘ vs φ ~ out.Φ parity).
    jls.Φₑ2.Φₑ      => 2π*(f_ext + 0.5),
)

@time sys_hb = jls.HarmonicSystem(model, jls.P1.source.ω, 1, determine_jacobian = true)
sysu = jls.unknowns(sys_hb.system)

βL_vals = 0.6:0.1:1.0
f_vals  = collect(0.0:0.005:1.0)
Z0 = 50.0
S21_dB = zeros(length(f_vals), length(βL_vals))

for (jβ, βL) in enumerate(βL_vals)
    # probe drive vector: source enters linearly, so one δU per βL is enough and the
    # amplitude cancels in the ratio b₂/a₁
    δU = jls.perturbation_response(sys_hb, jls.P1.source.I, make_ps(βL, 0.0), amplitude = 1.0)
    Uwp = fill(0.0, length(sysu))
    for (jf, f_ext) in enumerate(f_vals)
        ps = make_ps(βL, f_ext)
        prob = jls.NonlinearProblem(sys_hb.system, merge(Dict(sysu .=> Uwp), ps))
        sol = jls.ModelingToolkit.solve(prob, jls.NonlinearSolve.LevenbergMarquardt(); maxiters = 2000)
        Uwp = real.(sol.u)

        lin = jls.LinearisedProblem(sys_hb, ps, δU, [Ω_probe], U₀ = Uwp)
        jls.solve!(lin)

        V₁ = jls.get_solution(lin, jls.P1.dφ, 1)[1] * R₀ * I₀
        I₁ = jls.get_solution(lin, jls.P1.i, 1)[1] * I₀
        V₂ = jls.get_solution(lin, jls.P2.dφ, 1)[1] * R₀ * I₀
        a₁ = 0.5*(V₁ + Z0*I₁)/sqrt(Z0)
        
        S21_dB[jf, jβ] = 20*log10(abs(V₂/sqrt(Z0)) / abs(a₁))
    end
    println("βL = ", βL, ": S21 range [", round(minimum(S21_dB[:, jβ]), digits=1), ", ",
            round(maximum(S21_dB[:, jβ]), digits=1), "] dB, peak @ f_ext = ",
            f_vals[argmax(S21_dB[:, jβ])])
end

# axes framed exactly like the paper's Fig. 8b
p = jls.plot(xlabel = "Norm_Ext_Flux", ylabel = "dB(S(2,1))",
             title = "RF SQUID coupler — 5 GHz probe",
             xlim = (0.0, 1.0), xticks = 0.0:0.1:1.0,
             ylim = (-80, 0), yticks = -80:20:0, legend = :topright)
for (jβ, βL) in enumerate(βL_vals)
    jls.plot!(p, f_vals, S21_dB[:, jβ], lw = 2, label = "βL = $(βL)")
end
display(p)
