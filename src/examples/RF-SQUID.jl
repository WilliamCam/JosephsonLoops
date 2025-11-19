# Example: Coupled Josephson circuit with mutual inductance and HB

using JosephsonLoops
using ModelingToolkit
using Symbolics
using DifferentialEquations
using Plots

const jls = JosephsonLoops

#########################################################
# 1. Build circuit from netlist
#########################################################

loops = [
    ["J1", "L1"],
    ["L2", "C1"],
    ["C1", "R1"],
    ["R1", "I1"],
]

coupling = [(1, 2)]
ext_flux = [true, false, false, false]

circuit = jls.process_netlist(loops; mutual_coupling = coupling, ext_flux = ext_flux)

# circuit model ODAE system and initial condition vector are created.
model, u0, guesses = jls.build_circuit(circuit)
sys = model

#########################################################
# 2. Fix cyclic/self-referential guesses
#########################################################

function sanitize_guesses(guesses, u0)
    g = Dict{Any, Any}()
    for (k, v) in guesses
        try
            if v === k || v == k
                if haskey(u0, k)
                    g[k] = float(u0[k])
                else
                    g[k] = 0.0
                end
            else
                g[k] = v
            end
        catch
            g[k] = v
        end
    end
    return g
end

guesses = sanitize_guesses(guesses, u0)

#########################################################
# 3. Backwards-compatible loop aliases (if needed)
#########################################################

try
    if !isdefined(jls, :loop1) && isdefined(jls, :I1)
        Core.eval(jls, :(loop1 = I1))
    end
catch e
    @warn "Failed to create jls.loop1 alias: $e"
end

try
    if !isdefined(jls, :loop2)
        if isdefined(jls, Symbol("Φₑ2"))
            Core.eval(jls, :(loop2 = Φₑ2))
        elseif isdefined(jls, Symbol("Φₑ1"))
            Core.eval(jls, :(loop2 = Φₑ1))
        end
    end
catch e
    @warn "Failed to create jls.loop2 alias: $e"
end

#########################################################
# 4. Parameters
#########################################################

I₀ = 1.0e-6
R₀ = 5.0
Φ₀ = jls.Φ₀

βc  = 2π / Φ₀ * I₀ * R₀^2
βL  = 2π / Φ₀ * I₀
fdrive = 100e6

ps = [
    jls.I1.ω  => 2π * fdrive,
    jls.I1.I  => 1.0 * I₀,
    jls.J1.I0 => I₀,
    jls.J1.R  => R₀,
    jls.J1.C  => 0.01 / βc,
    jls.R1.R  => 50.0,
    jls.C1.C  => 2.0 / βc,
    jls.L1.L  => 2.0 / βc,
    jls.L2.L  => 100.0 / βL,
    jls.M12.L => 8.0 / βL,
]

p_dict = Dict(ps)  # convenient mutable form

#########################################################
# 5. Time-domain transient using ODEProblem
#########################################################

tspan = (0.0, 1e-6)
saveat = range(tspan[2] / 10.0, tspan[2], length = 10_000)

prob = ODEProblem(
    sys,
    merge(Dict(u0), p_dict),
    tspan;
    guesses = guesses,
)

sol = solve(prob, Rodas5(); saveat = saveat)

# Example transient plots (currents/voltages on R1)
p1 = plot(sol.t, sol[jls.R1.i],
          xlabel = "t (s)", ylabel = "I_R1 (A)",
          title = "R1 current (transient)")

# If your R1 has a voltage variable, you can also do:
# p2 = plot(sol.t, sol[jls.R1.v],
#           xlabel = "t (s)", ylabel = "V_R1 (V)",
#           title = "R1 voltage (transient)")

display(p1)
# display(p2)  # uncomment if you plot voltage too

#########################################################
# 6. Parameter sweep over external flux Φₑ1
#########################################################

# For a single external flux on loop 1, JosephsonLoops will usually
# create a component Φₑ1 with a parameter Φₑ:
Φ_param = jls.Φₑ1.Φₑ

Φspan = (0.0, 2.0 * Φ₀)
Φ_vals = range(Φspan[1], Φspan[2], length = 80)

R1_final_vs_Φ = Float64[]

for Φ in Φ_vals
    pΦ = copy(p_dict)
    pΦ[Φ_param] = Φ
    probΦ = ODEProblem(sys, merge(Dict(u0), pΦ), tspan; guesses = guesses)
    solΦ = solve(probΦ, Rodas5(); saveat = saveat)
    push!(R1_final_vs_Φ, solΦ[jls.R1.i][end])  # R1 current at final time
end

pΦ = plot(Φ_vals ./ Φ₀, R1_final_vs_Φ,
          xlabel = "Φₑ1 / Φ₀",
          ylabel = "I_R1(t_end) (A)",
          title = "R1 current vs external flux")
display(pΦ)

#########################################################
# 7. Parameter sweep over drive current I1.I
#########################################################

Ispan = (0.0, 2.0 * I₀)
I_vals = range(Ispan[1], Ispan[2], length = 80)

R1_final_vs_I = Float64[]

for Id in I_vals
    pI = copy(p_dict)
    pI[jls.I1.I] = Id
    probI = ODEProblem(sys, merge(Dict(u0), pI), tspan; guesses = guesses)
    solI = solve(probI, Rodas5(); saveat = saveat)
    push!(R1_final_vs_I, solI[jls.R1.i][end])
end

pI = plot(I_vals ./ I₀, R1_final_vs_I,
          xlabel = "I_drive / I₀",
          ylabel = "I_R1(t_end) (A)",
          title = "R1 current vs drive amplitude")
display(pI)

#########################################################
# 8. Harmonic Balance (HB) – RLC-style
#########################################################

include("../harmonic balance/colocation HB.jl")

eqs, states = jls.get_full_equations(model, jls.t)

Nharmonics = 3
harmonic_sys, harmonic_states =
    jls.harmonic_equation(eqs, states, jls.t, jls.I1.ω, Nharmonics)

@named ns = NonlinearSystem(harmonic_sys)

hb_sys = structural_simplify(ns;
                              fully_determined = false,
                              check_consistency = false)

state_syms = collect(unknowns(hb_sys))

N = 300
ω_vec = range(0.8, 1.2, length = N) .* (2π * fdrive)
solution = Float64[]

u0_prev = zeros(length(state_syms))

for ω in ω_vec
    psHB = Dict(
        jls.I1.ω  => ω,
        jls.I1.I  => 0.6 * I₀,
        jls.J1.I0 => I₀,
        jls.J1.R  => R₀,
        jls.J1.C  => 0.01 / βc,
        jls.R1.R  => 50.0,
        jls.C1.C  => 2.0 / βc,
        jls.L1.L  => 2.0 / βc,
        jls.L2.L  => 100.0 / βL,
        jls.M12.L => 8.0 / βL,
        Φ_param   => 0.5 * Φ₀,   # bias flux
    )

    state_guess = isempty(state_syms) ? Dict() : Dict(state_syms .=> u0_prev)
    prob_map    = merge(state_guess, psHB)

    hb_prob = NonlinearProblem(hb_sys, prob_map;
                               allow_incomplete = true,
                               check_length = false)

    hb_sol = solve(hb_prob)

    u0_prev .= hb_sol.u

    # This combination is what you had originally;
    # assumes ns.E and ns.F exist in the HB system.
    push!(solution,
          hb_sol[ns.E[1]] +
          sqrt(hb_sol[ns.E[3]]^2 + hb_sol[ns.F[1]]^2))
end

pHB = plot(ω_vec ./ (2π),
           solution,
           xlabel = "Frequency (Hz)",
           ylabel = "Amplitude (arb.)",
           title = "HB response of coupled Josephson circuit")
display(pHB)


# bif_par = jls.loop4.sys.ω
# p_start = [
#     jls.loop4.sys.ω => ω_vec[1]
#     jls.loop4.sys.I => 1.0*I₀
#     jls.J1.sys.I0 => I₀
#     jls.J1.sys.R => R₀
#     jls.J1.sys.C => 0.01/βc
#     jls.J1.sys.L => 2.0/βL
#     jls.R1.sys.R => 50.0
#     jls.C1.sys.C => 2.0/βc
#     jls.L2.sys.L => 100.0/βL
#     jls.M12.sys.L => 8.0/βL
#     jls.loop1.sys.Φₑ => 0.5*Φ₀]
# u0_guess = Dict(unknowns(sys) .=> 0.0)
# plot_var = ns.E[1] + sqrt(ns.E[3]^2+ns.F[1]^2)


