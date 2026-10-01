using JosephsonLoops
using ModelingToolkit
using Plots

# ---- circuit ----------------------------------------------------------------------
loops = [
    ["I1", "R1", "Lin"],
    ["Lt", "J1", "J2"],
    ["J2", "Idc"],
    ["Idc", "R2"],
]

circuit = process_netlist(loops, mutual_coupling=[(1, 2)], ext_flux=[false, true, false, false])
squid, u0, guesses = build_circuit(circuit)

# ---- normalisation and bias ---------------------------------------------------------
I₀ = 1e-6
R₀ = 10.0e3
ωc = R₀ * I₀ / (Φ₀ / 2π)
Z0 = 50.0
f_probe = collect(0.1:0.1:10.0) # MHz offset from the running-state carrier
Ω_probe_offset = 2π .* f_probe .* 1e6 ./ ωc

ps = Dict(
    squid.I1.ω       => first(Ω_probe_offset),
    squid.I1.I       => 0.0,
    squid.R1.r       => Z0 / R₀,
    squid.R2.r       => Z0 / R₀,
    squid.Idc.I      => 2.2, # 2.2 I₀: voltage-running state for this flux and circuit
    squid.Idc.ω      => 0.0,
    squid.Lin.βL     => 400.0,
    squid.J1.βc      => 1000.0e-15 * R₀ * ωc,
    squid.J1.r       => 6.0 / R₀,
    squid.J2.βc      => 1000.0e-15 * R₀ * ωc,
    squid.J2.r       => 6.0 / R₀,
    squid.Lt.βL      => 1.0,
    squid.M12.βM     => 2.0,
    squid.Φₑ2.Φₑ    => 0.5 * 2π,
)

# Estimate the Josephson carrier from the late-time phase slopes. The two junction
# phases wind in opposite directions in this netlist, but share the same fundamental.
tspan = (0.0, 4e-9) .* ωc
tsol = tsolve(squid, guesses, ps, tspan; guesses=guesses)

function phase_slope(t, phase)
    first_idx = cld(length(t), 2)
    t_fit = t[first_idx:end]
    phase_fit = phase[first_idx:end]
    t_mean = sum(t_fit) / length(t_fit)
    phase_mean = sum(phase_fit) / length(phase_fit)
    sum((t_fit .- t_mean) .* (phase_fit .- phase_mean)) /
        sum((t_fit .- t_mean).^2)
end

ωJ1_td = phase_slope(tsol.t, tsol[squid.J1.φ])
ωJ2_td = phase_slope(tsol.t, tsol[squid.J2.φ])
ω0_guess = (abs(ωJ1_td) + abs(ωJ2_td)) / 2
println("Time-domain phase-rate estimates: ",
        round(ωJ1_td * ωc / (2π) / 1e9, digits=4), " and ",
        round(ωJ2_td * ωc / (2π) / 1e9, digits=4), " GHz")

# ---- autonomous HB of the voltage-running state ------------------------------------
# Each selected phase is represented as winding*ω₀*t plus a periodic Fourier correction.
# The sign records the direction of phase advance in the circuit's branch convention.
sys = HarmonicSystem(
    squid, squid.I1.ω, 2;
    autonomous=true,
    drifting_states=Dict(squid.J1.φ => 1, squid.J2.φ => -1),
    determine_jacobian=true,
)
prob = HarmonicProblem(sys, ps)
ω0_index = JosephsonLoops.var_index(unknowns(sys.system), sys.autonomous_frequency)
prob.U₀[ω0_index] = ω0_guess
sol = JosephsonLoops.solve!(prob)

U_wp = real.(sol)
ω0 = U_wp[ω0_index]
hb_residual = prob.problem.f(U_wp, prob.problem.p)
max_hb_residual = maximum(abs.(hb_residual))
f0 = ω0 * ωc / (2π)
println("Autonomous HB carrier: ", round(f0 / 1e9, digits=4), " GHz")
println("Maximum autonomous-HB residual: ", max_hb_residual)
@assert ω0 > 0 "Autonomous HB returned a non-positive carrier frequency"
@assert max_hb_residual < 1e-8 "Autonomous HB did not converge to a consistent solution"

# ---- MHz small-signal response about the running state ------------------------------
# LinearisedProblem takes absolute probe frequencies. Shift the probe offsets by the
# solved carrier, and set the reference tone to that same carrier.
ps_run = merge(ps, Dict(squid.I1.ω => ω0))
δU = perturbation_response(sys, squid.I1.I, ps_run, amplitude=1.0)
Ω_probe = ω0 .+ Ω_probe_offset
lin = LinearisedProblem(sys, ps_run, δU, Ω_probe, U₀=U_wp)
JosephsonLoops.solve!(lin)

# The resistor relation D(φ) = r*i gives its voltage in the model's normalised units.
vout = get_solution(lin, squid.R2.r * squid.R2.i, 1)
Ztrans = R₀ .* vout
Ztrans_dBΩ = 20 .* log10.(abs.(Ztrans))

println("MHz probe sweep: transimpedance magnitude ",
        round(minimum(abs.(Ztrans)), sigdigits=4), " to ",
        round(maximum(abs.(Ztrans)), sigdigits=4), " Ω")
println("Response finite across sweep: ", all(isfinite, real.(Ztrans)) &&
        all(isfinite, imag.(Ztrans)))
@assert all(isfinite, real.(Ztrans)) && all(isfinite, imag.(Ztrans)) "Non-finite small-signal response"

p1 = plot(f_probe, Ztrans_dBΩ, lw=2, label=false,
          xlabel="Probe offset from carrier (MHz)", ylabel="|Vout / Iin| (dBΩ)",
          title="Running DC-SQUID small-signal transimpedance")
p2 = plot(f_probe, angle.(Ztrans) .* 180 / π, lw=2, label=false,
          xlabel="Probe offset from carrier (MHz)", ylabel="Phase (degrees)")
display(plot(p1, p2, layout=(2, 1), size=(760, 650)))
