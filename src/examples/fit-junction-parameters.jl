# Fit junction and coupling-capacitor parameters from a pumped JPA S11 spectrum.

using JosephsonLoops
using ModelingToolkit
using Optimization
using OptimizationOptimJL
using Plots
using LinearAlgebra: norm

jpa, _, _ = build_circuit(process_netlist([["P1", "C1", "J1"]]))
I₀ = Φ₀ / (2π * 1000e-12)
R₀ = 10e3
ωc = R₀ * I₀ / (Φ₀ / (2π))
GHz(f) = 2π * f * 1e9 / ωc

βc_coupling_nominal = 100e-15 * R₀ * ωc
βc_junction_nominal = 1000e-15 * R₀ * ωc
pump_current = 11.3e-9 / I₀
pump_frequency = 4.75001
probe_frequency = collect(range(4.65, 4.85, length = 201))
probe_Ω = GHz.(probe_frequency)
working_frequency = collect(4.50:0.005:5.00)
working_Ω = GHz.(working_frequency)

hb_system = HarmonicSystem(jpa, jpa.P1.source.ω, 2; determine_jacobian = true)
s11 = get_HB_scattering_matrix(jpa, '1', '1')[1]
fixed_parameters = Dict(
    jpa.P1.Rₙ.r => 50 / R₀,
    jpa.P1.source.I => pump_current,
    jpa.P1.source.ω => GHz(pump_frequency),
    jpa.J1.r => 1e8 / R₀,
)

# Use a frequency continuation sweep to find the pumped branch near the pump tone.
seed_parameters = copy(fixed_parameters)
delete!(seed_parameters, jpa.P1.source.ω)
seed_problem = HarmonicProblem(
    hb_system,
    merge(
        seed_parameters,
        Dict(
            jpa.C1.βc => βc_coupling_nominal,
            jpa.J1.α => 1.0,
            jpa.J1.βc => βc_junction_nominal,
        ),
    );
    parameter_sweep = [jpa.P1.source.ω => working_Ω],
)
JosephsonLoops.solve!(seed_problem)
seed_index = argmin(abs.(working_frequency .- pump_frequency))
seed = real.(seed_problem.result.solution[:, seed_index])

function predict_s11(α, junction_capacitance_scale, coupling_capacitance_scale)
    parameters = merge(
        fixed_parameters,
        Dict(
            jpa.C1.βc => coupling_capacitance_scale * βc_coupling_nominal,
            jpa.J1.α => Float64(α),
            jpa.J1.βc => junction_capacitance_scale * βc_junction_nominal,
        ),
    )
    working_problem = HarmonicProblem(hb_system, parameters; U₀ = seed)
    JosephsonLoops.solve!(working_problem)
    working_point = real.(working_problem.result.solution)
    modulation = perturbation_response(
        hb_system,
        jpa.P1.source.I,
        parameters;
        amplitude = pump_current,
    )
    linear_problem = LinearisedProblem(
        hb_system,
        parameters,
        modulation,
        probe_Ω;
        U₀ = working_point,
    )
    JosephsonLoops.solve!(linear_problem)
    spectrum = get_solution(linear_problem, s11, 1)
    all(isfinite, spectrum) || error("Non-finite S11 for α=$α, Cj=$junction_capacitance_scale, Cc=$coupling_capacitance_scale")
    return spectrum
end

# Synthetic target and the optimizer's initial prediction; all three scales are fitted.
α_true, junction_capacitance_true, coupling_capacitance_true = 1.01, 1.01, 1.01
α_initial, junction_capacitance_initial, coupling_capacitance_initial = 0.98, 0.98, 0.98
target_s11 = predict_s11(α_true, junction_capacitance_true, coupling_capacitance_true)
initial_s11 = predict_s11(α_initial, junction_capacitance_initial, coupling_capacitance_initial)

function s11_loss(log_α, log_junction_capacitance, log_coupling_capacitance)
    prediction = predict_s11(
        exp(log_α),
        exp(log_junction_capacitance),
        exp(log_coupling_capacitance),
    )
    return 1e6 * sum(abs2, prediction .- target_s11) / length(target_s11)
end

@register_symbolic s11_loss(log_α, log_junction_capacitance, log_coupling_capacitance)
@variables log_α log_junction_capacitance log_coupling_capacitance
@named fit_system = OptimizationSystem(
    s11_loss(log_α, log_junction_capacitance, log_coupling_capacitance),
    [log_α, log_junction_capacitance, log_coupling_capacitance],
    [],
)
fit_system = complete(fit_system)
fit_problem = OptimizationProblem(
    fit_system,
    [
        log_α => log(α_initial),
        log_junction_capacitance => log(junction_capacitance_initial),
        log_coupling_capacitance => log(coupling_capacitance_initial),
    ],
)
fit_solution = Optimization.solve(
    fit_problem,
    NelderMead(initial_simplex = OptimizationOptimJL.Optim.AffineSimplexer(a = 0.1, b = 0.5));
    maxiters = 500,
    g_abstol = 1e-20,
)
α_fit, junction_capacitance_fit, coupling_capacitance_fit = exp.(fit_solution.u)
fit_s11 = predict_s11(α_fit, junction_capacitance_fit, coupling_capacitance_fit)
relative_error = norm(fit_s11 - target_s11) / norm(target_s11)

println("Optimization result: $(fit_solution.retcode)")
println("Truth:   α=$α_true, junction C scale=$junction_capacitance_true, coupling C scale=$coupling_capacitance_true")
println("Initial: α=$α_initial, junction C scale=$junction_capacitance_initial, coupling C scale=$coupling_capacitance_initial")
println("Fitted:  α=$(round(α_fit, sigdigits = 6)), junction C scale=$(round(junction_capacitance_fit, sigdigits = 6)), coupling C scale=$(round(coupling_capacitance_fit, sigdigits = 6))")
println("Target peak gain: $(round(maximum(20 .* log10.(abs.(target_s11))), digits = 2)) dB")
println("Initial relative S11 error: $(round(norm(initial_s11 - target_s11) / norm(target_s11), sigdigits = 5))")
println("Fitted relative S11 error: $(round(relative_error, sigdigits = 5))")
@assert(
    isfinite(relative_error) && relative_error < 1e-3 &&
        abs(α_fit - α_true) / α_true < 0.02 &&
        abs(junction_capacitance_fit - junction_capacitance_true) / junction_capacitance_true < 0.02 &&
        abs(coupling_capacitance_fit - coupling_capacitance_true) / coupling_capacitance_true < 0.02,
    "JPA fit failed its synthetic parameter-recovery tolerances",
)

p_gain = plot(
    probe_frequency,
    20 .* log10.(abs.(target_s11));
    linewidth = 2.5,
    label = "Synthetic target",
    xlabel = "Probe frequency (GHz)",
    ylabel = "20 log₁₀ |S₁₁| (dB)",
    title = "Pumped JPA reflection gain",
    grid = true,
)
plot!(p_gain, probe_frequency, 20 .* log10.(abs.(initial_s11)); linewidth = 2, linestyle = :dot, label = "Initial guess")
plot!(p_gain, probe_frequency, 20 .* log10.(abs.(fit_s11)); linewidth = 2, linestyle = :dash, label = "Optimized fit")

p_phase = plot(
    probe_frequency,
    rad2deg.(angle.(target_s11));
    linewidth = 2.5,
    label = "Synthetic target",
    xlabel = "Probe frequency (GHz)",
    ylabel = "Phase(S₁₁) (degrees)",
    title = "Reflection phase",
    grid = true,
)
plot!(p_phase, probe_frequency, rad2deg.(angle.(initial_s11)); linewidth = 2, linestyle = :dot, label = "Initial guess")
plot!(p_phase, probe_frequency, rad2deg.(angle.(fit_s11)); linewidth = 2, linestyle = :dash, label = "Optimized fit")

p_residual = plot(
    probe_frequency,
    abs.(fit_s11 .- target_s11);
    linewidth = 2,
    color = :firebrick,
    label = false,
    xlabel = "Probe frequency (GHz)",
    ylabel = "|S₁₁ fit - target|",
    title = "Optimized complex residual",
    grid = true,
)

diagnostic_plot = plot(
    p_gain,
    p_phase,
    p_residual;
    layout = (3, 1),
    size = (900, 1050),
    plot_title = "JPA S₁₁ fit: initial and optimized spectra",
    left_margin = 8Plots.mm,
)
image_path = joinpath(@__DIR__, "..", "..", "docs", "images", "fit-junction-parameters.png")
savefig(diagnostic_plot, image_path)
display(diagnostic_plot)
println("S₁₁ comparison plot saved to $image_path")
