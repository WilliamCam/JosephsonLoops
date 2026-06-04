
using Symbolics
using SymbolicUtils
using ModelingToolkit
using QuestBase
using DifferentialEquations
using Plots
using JosephsonLoops
using Revise
const jls = JosephsonLoops
# Driven Duffing oscillator (Kosata 2022, p.59) solved by collocation harmonic balance.

# Time-domain equation
@variables t x(t)
@parameters α ω ω0 F γ η
diff_eq = Differential(t)(Differential(t)(x)) + ω0^2*x + α*x^3 + η*Differential(t)(x)*x^2 + γ*Differential(t)(x) - F*cos(ω*t) ~ 0

Nharmonics = 1

# Wrap the Duffing equation in an ODESystem so HarmonicSystem can consume it
@named duffing_sys = System([diff_eq], t)
duffing_model = mtkcompile(duffing_sys)

# Build the harmonic system with linearization data populated (J0, J1)
h_sys = jls.HarmonicSystem(duffing_model, ω, Nharmonics; determine_jacobian=true)

# Large-signal pump: sweep ω and solve the harmonic balance system. The same h_sys feeds
# the linear response below; determine_jacobian=true (above) is what builds its Jacobians.
ω_vec = collect(range(0.8, 1.2, 200))
ps = Dict(α => 1.0, ω0 => 1.0, F => 0.01, η => 1.0e-1, γ => 1.0e-3, ω => first(ω_vec))
h_prob = jls.HarmonicProblem(h_sys, ω_vec, ps)
sweep_res = jls.solve!(h_prob)

# Pump amplitude vs ω: DC + fundamental magnitude. get_output takes the variable name
# directly, so there's no need to pull coefficient symbols out of the variable_map.
x_dc   = jls.get_output(h_prob, sweep_res, "x", 0)
x_fund = jls.get_output(h_prob, sweep_res, "x", 1)
solution = real.(x_dc) .+ abs.(x_fund)
plot(ω_vec, solution)

#  Small-signal linearisation around a pump tone (Kosata 2022 eq. 5.12)
j  = 70
ωp = ω_vec[j]

# Probe: unit kick (1e-3) on x's sin@1 coefficient. For a single state the HB coefficients
# are ordered [DC, cos₁, sin₁], so the sin@1 entry is index 3.
perturbation = zeros(Float64, length(unknowns(h_sys.system)))
perturbation[3] = 1.0e-3

# Warm-start the working-point solve from the j-th sweep point (lands on the high branch)
U₀ = sweep_res.solution[ω][:, j]
Ω = collect(range(0.8, 1.2, 800))
lin_prob = jls.HarmonicProblem(h_sys, Ω, ps; linear_response=(ωp, perturbation))
jls.solve!(lin_prob)

# Response is stored as real/imag slices: solution[ω][coeff, Ω, (1=real, 2=imag)].
sol   = lin_prob.result.solution[ω]
u_resp = sol[3, :]  # sin@1 response
v_resp = sol[2, :]  # cos@1 response
out = abs.((u_resp .+ 1im .* v_resp) ./ 2)
plot(Ω, out)
