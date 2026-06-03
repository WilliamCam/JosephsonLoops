
using Revise
using JosephsonLoops
using Symbolics
using ModelingToolkit
using BenchmarkTools

const jls = JosephsonLoops
loops = [
    ["P1", "C1", "J1"]
]

@time circuit = jls.process_netlist(loops)
@time model, u0, guesses = jls.build_circuit(circuit)
ps = Dict(
    jls.P1.Isrc.ω => 100e6*2*pi,
    jls.P1.Isrc.I => 0.00565e-6,
    jls.P1.Rsrc.R => 50.0,
    jls.C1.C => 100.0e-15,
    jls.J1.C => 1000.0e-15,
    jls.J1.I0 => jls.Φ₀/(2*pi*1.0e-9),
    jls.J1.R => 50.0
)

#time domain simulation 
tspan = (0.0, 1e-6)
tsol = jls.tsolve(model, guesses, ps, tspan; guesses=guesses)
p1 = jls.plot(tsol[jls.C1.i][end-400:end], title = "Transient Time Plot", xlabel = "t", ylabel = "I_C1")

#hb Setup
# Define the sweep range (8 to 10.0 GHz)
ω_vec = collect(2*pi*(1:0.1:30.0)*1e9)


sweep_params = delete!(ps, jls.P1.Isrc.ω)

sys = jls.HarmonicSystem(model, jls.P1.Isrc.ω, 1)
prob = jls.HarmonicProblem(sys, ω_vec, sweep_params)

result = jls.solve!(prob)
out = prob.result.solution[jls.P1.Isrc.ω]

current = jls.get_output(prob, result, "C1₊i", 1)
theta_p_mag = jls.get_output(prob, result, "P1₊dθ",  1)


Ii = @. (0.00565e-6 - current)
Vi = @. (jls.Φ₀ / (2*pi) * real.(theta_p_mag)) 
# Calculate Power Waves (a = incident, b = reflected)
Z0 = 50.0
ai = @. 0.5 * (Vi + Z0 * Ii) / sqrt(Z0)
bi = @. 0.5 * (Vi - Z0 * Ii) / sqrt(Z0)
p = jls.plot(ω_vec/(2*pi), 20*log10.(abs.(bi./ai)), xlabel="Frequency (Hz)", ylabel="S11 (dB)", title="RLC S-Parameter", lw=2)

# Build HarmonicSystem once (expensive symbolic expansion, shared across problems)
resp = zeros(size(sys.jacobian[1],1))
resp[11] = 1.0
resp[12] = 1.0
lin_prob = jls.HarmonicProblem(sys, ω_vec, sweep_params, linear_response = (2*pi*4.75001*1e9, resp))
# Solve
sweep_res = jls.solve!(lin_prob)
out = lin_prob.result.solution[jls.P1.Isrc.ω]
plot(out[12,:])


#  JPA linearization (matches MIT JosephsonCircuits.jl JPA example)
# Pump tone (MIT example)
ωp = 2*pi*4.75001e9
Ip = 0.00565e-6
jpa_params[jls.P1.Isrc.I] = Ip

# HarmonicSystem with linearization data; tearing=false keeps all 12 vars as unknowns
@time h_sys_lin = jls.HarmonicSystem(model, jls.P1.Isrc.ω; N=1, linear=true, tearing=false)

# Single-point nonlinear HB at the pump frequency to get the working point

pump_prob = jls.HarmonicProblem(h_sys_lin, jpa_params;
    sweep_var=jls.P1.Isrc.ω, sweep_vals=[ωp])
pump_res = jls.solve(pump_prob)

# U₀: every harmonic coefficient at the pump
U₀ = Dict{Num, Float64}(v => vals[1] for (v, vals) in pump_res.results)

# Small-signal frequency sweep (4.5–10 GHz)
Ω_vals = 2*pi*(4.5:0.001:10.0)*1e9

# Linearize: unit kick on the port-current cos coefficient at the fundamental
lin_prob = jls.LinearizedProblem(h_sys_lin, jpa_params, ωp, U₀, Ω_vals;
    perturb_var="P1₊i", perturb_order=1, perturb_component=:Cos, perturb_amplitude=1.0)
lin_result = jls.solve(lin_prob)

# S11 from the small-signal current/phase response (mirrors the nonlinear S11 above)
P1i_A  = h_sys_lin.variable_map[("P1₊i",  1, :Cos)]
P1i_B  = h_sys_lin.variable_map[("P1₊i",  1, :Sin)]
P1dθ_A = h_sys_lin.variable_map[("P1₊dθ", 1, :Cos)]
P1dθ_B = h_sys_lin.variable_map[("P1₊dθ", 1, :Sin)]
current_p = lin_result.results[P1i_A]  .+ 1im .* lin_result.results[P1i_B]
theta_p   = lin_result.results[P1dθ_A] .+ 1im .* lin_result.results[P1dθ_B]

Ii = @. (1.0 - current_p)
Vi = @. (jls.Φ₀ / (2*pi) * theta_p)
Z0 = 50.0
ai = @. 0.5 * (Vi + Z0 * Ii) / sqrt(Z0)
bi = @. 0.5 * (Vi - Z0 * Ii) / sqrt(Z0)
jls.plot(Ω_vals/(2*pi*1e9), 20*log10.(abs.(bi./ai)),
    xlabel="Frequency (GHz)", ylabel="S11 (dB)",
    title="JPA linearized S11 (LinearizedProblem)", lw=4, label="JosephsonLoops")


arr1 = [1, 2, 3]
arr2 = [:a, :b]

# Creates an iterator of tuples
combinations = Iterators.product(arr1, arr2)

for (x, y) in combinations
    println("Combination: $x, $y")
end