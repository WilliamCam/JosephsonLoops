# Lumped nonlinear Purcell-filter readout model: pump-dependent S11 phase.
#
# This is a useful circuit-level design problem because a Purcell filter couples a readout
# resonator to a strongly nonlinear, dissipative mode. Harmonic balance resolves the driven
# operating point; LinearisedProblem then predicts the weak-probe reflection around it.
# The comparison target is the measured readout-frequency/amplitude-dependent S11 phase in
# Fig. 4 of "Photon-noise-tolerant dispersive readout of a superconducting qubit using a
# nonlinear Purcell filter", arXiv:2309.04315. The paper reports filter-mode targets
# f_f = 9.791 GHz, κ_f/2π = 310 MHz and α_f/2π = -120 MHz.
#
# Netlist layout follows the JPA example's process_netlist/build_circuit pattern. This
# reduced two-loop model couples a readout RLC loop to a JJ-based nonlinear filter loop.
# Qubit states are represented only by a small effective readout-mode detuning; this is not
# a quantum qubit model. Source current is in normalized circuit units and is not calibrated
# to the paper's Ωro. Compare trends and spectral shape, not absolute pulse powers or
# readout fidelities. The full experiment also has finite excitation and photon-noise effects
# outside this classical lumped model.

using JosephsonLoops
using ModelingToolkit
import Plots

# ---- coupled lumped circuit ---------------------------------------------------------
# Readout loop: P1-C1-L1. Nonlinear filter loop: L1-C2-J1, coupled through shared L1.
loops = [["P1", "C1", "L1"], ["L1", "C2", "J1"]]
circuit = process_netlist(loops)
filter, _, _ = build_circuit(circuit)

# ---- paper-based mode estimates and normalization -----------------------------------
const f_filter = 9.791e9
const κ_filter = 0.31e9
const α_filter = -120e6
const e_charge = 1.602176634e-19
const h_planck = 6.62607015e-34
const Φbar = Φ₀ / (2π)

# EC/h ≈ |α| gives the filter-mode capacitance; the plasma-frequency estimate sets Ic.
const C_filter = e_charge^2 / (2 * h_planck * abs(α_filter))
const C_junction = 10e-15
const C_mode = C_filter + C_junction
const I_c = C_mode * Φbar * (2π * (f_filter + abs(α_filter)))^2
const I₀ = I_c
const R₀ = 50.0
const ωc = R₀ * I₀ / Φbar
ωnorm(f) = 2π * f / ωc
pF(C) = C * R₀ * ωc
βL(L) = L * I₀ / Φbar

const C_readout = 100e-15
const f_readout = 9.8e9
const L_readout = 2.64e-9
const R_filter = 1 / (2π * κ_filter * C_filter)

base_ps = Dict(
    filter.P1.Rₙ.r     => 1.0,                     # 50 Ω port with R₀ = 50 Ω
    filter.P1.source.ω => ωnorm(f_readout),
    filter.P1.source.I => 0.01,
    filter.C1.βc       => pF(C_readout),
    filter.L1.βL       => βL(L_readout),
    filter.C2.βc       => pF(C_filter),
    filter.J1.βc       => pF(C_junction),
    filter.J1.α        => 1.0,
    filter.J1.r        => R_filter / R₀,
)

# Treat the qubit-state dispersive shift as a small change to the effective readout mode.
# The ±1 MHz shift is illustrative; the pump remains fixed at the readout frequency.
const state_shift = 1e6
readout_states = [
    ("|g⟩", f_readout - state_shift),
    ("|e⟩", f_readout + state_shift),
]
pump_levels = [0.03, 0.3, 1.0]
f_probe = collect(9.7995:0.000001:9.8005) # GHz; 1 kHz spacing resolves the reflection notch
Ω_probe = ωnorm.(f_probe .* 1e9)

sys = HarmonicSystem(filter, filter.P1.source.ω, 2, determine_jacobian = true)
s11 = get_HB_scattering_matrix(filter, '1', '1')[1]

function unwrap_phase_degrees(phase)
    unwrapped = similar(phase)
    unwrapped[1] = phase[1]
    offset = 0.0
    for i in 2:length(phase)
        delta = phase[i] - phase[i - 1]
        delta > 180 && (offset -= 360)
        delta < -180 && (offset += 360)
        unwrapped[i] = phase[i] + offset
    end
    return unwrapped
end

# ---- driven solutions and probe reflection ------------------------------------------
function main()
    phase_spectra = Dict{Tuple{String,Float64},Vector{Float64}}()
    s11_spectra = Dict{Tuple{String,Float64},Vector{ComplexF64}}()
    for (state, state_frequency) in readout_states
        # L ∝ 1/f² for the effective readout frequency at fixed capacitance.
        ps_state = merge(base_ps, Dict(
            filter.L1.βL => βL(L_readout * (f_readout / state_frequency)^2),
        ))
        U₀ = zeros(length(unknowns(sys.system)))

        for level in pump_levels
            ps = merge(ps_state, Dict(filter.P1.source.I => level))
            working_point = HarmonicProblem(sys, ps, U₀ = U₀)
            JosephsonLoops.solve!(working_point)
            U₀ = real.(working_point.result.solution)

            δU = perturbation_response(
                sys,
                filter.P1.source.I,
                ps,
                amplitude = level,
            )
            linear = LinearisedProblem(sys, ps, δU, Ω_probe, U₀ = U₀)
            JosephsonLoops.solve!(linear)
            S11 = get_solution(linear, s11, 1)
            @assert all(isfinite, S11) "Non-finite reflection spectrum for state $state at drive $level"

            phase = unwrap_phase_degrees(rad2deg.(angle.(S11)))
            phase_spectra[(state, level)] = phase
            s11_spectra[(state, level)] = S11
            peak_index = argmax(abs.(S11))
            dip_index = argmin(abs.(S11))
            println("state $state, normalized drive $level: unwrapped phase span ",
                    "$(round(last(phase) - first(phase), digits = 2))°; peak |S₁₁| ",
                    "$(round(abs(S11[peak_index]), digits = 3)) at ",
                    "$(round(f_probe[peak_index], digits = 5)) GHz; minimum ",
                    "$(round(abs(S11[dip_index]), digits = 3)) at ",
                    "$(round(f_probe[dip_index], digits = 5)) GHz")
        end
    end

    # ---- plot against the experiment's phase-spectrum observable ---------------------
    p = Plots.plot(
        xlabel = "Probe frequency (GHz)",
        ylabel = "arg(S₁₁) (degrees)",
        title = "Two-mode lumped nonlinear Purcell filter",
        legend = :outerright,
    )
    pmag = Plots.plot(
        xlabel = "Probe frequency (GHz)",
        ylabel = "20 log₁₀ |S₁₁| (dB)",
        title = "Weak-probe reflection magnitude",
        legend = :outerright,
    )
    for (state, _) in readout_states
        for level in pump_levels
            Plots.plot!(
                p,
                f_probe,
                phase_spectra[(state, level)],
                lw = 2,
                label = "$state, drive = $level",
            )
            Plots.plot!(
                pmag,
                f_probe,
                20 .* log10.(max.(abs.(s11_spectra[(state, level)]), eps(Float64))),
                lw = 2,
                label = "$state, drive = $level",
            )
        end
    end
    Plots.vline!(p, [f_filter / 1e9], ls = :dot, color = :black, label = "filter-mode target")
    display(Plots.plot(p, pmag, layout = (2, 1), size = (900, 800)))

    println("Reference filter targets: f = $(f_filter / 1e9) GHz, κ/2π = ",
            "$(κ_filter / 1e6) MHz, α/2π = $(α_filter / 1e6) MHz")
    println("Lumped estimates: Cfilter = $(round(C_filter * 1e15, digits = 1)) fF, ",
            "Ic = $(round(I_c * 1e9, digits = 1)) nA, Lreadout = ",
            "$(round(L_readout * 1e9, digits = 2)) nH, Rfilter = ",
            "$(round(R_filter, digits = 1)) Ω")
end

main()
