# In this example e model a parallel RLC circuit driven by a voltage source

using JosephsonLoops
using ModelingToolkit
using Symbolics
using DifferentialEquations
using Plots
using NonlinearSolve




const jls = JosephsonLoops
#build circuit
loops = [
    ["I1", "R1"],
    ["R1", "C1", "J1"],
]


ext_flux = [true, true]

circuit = jls.process_netlist(loops; ext_flux = ext_flux)

# circuit model ODAE system and initial condition vector are created.
model, u0, guesses = jls.build_circuit(circuit)
sys = model

# 2. Fix cyclic/self-referential guesses

guesses = jls.sanitize_guesses(guesses, u0)



# 3. Create backward-compatible aliases loop1 / loop2


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

# 4. Set parameter values
p = [
    jls.I1.ω  => 2pi * 8.3e9,
    jls.I1.I  => 0.5e-6,      # Reduced to 10 nA for small signal linear response
    jls.C1.C  => 100.0e-15,
    jls.J1.C  => 1000.0e-15,
    jls.J1.I0 => 1e-6,
    jls.R1.R  => 50.0,      # Increased to 50 Ohm for higher Q
    jls.J1.R  => 50.0,    # Increased to 10 kOhm to remove internal damping
    jls.Φₑ2.Φₑ => 0.5,
]


# 5. Time-domain simulation via ODEProblem


tspan = (0.0, 1e-6)

prob = ODEProblem(
    sys,
    merge(Dict(u0), Dict(p)),
    tspan;
    guesses = guesses,
)

sol = solve(prob, Rodas5())

p1 = plot(sol[jls.C1.i][end-400:end], title = "Transient Time Plot", xlabel = "t", ylabel = "I_C1")
savefig(p1, "transient_plot.png")
p_dc = jls.set_param(p, jls.I1.I, 0.0)

tspan_dc = (0.0, 5e-6)


prob_dc = ODEProblem(
    sys,
    merge(Dict(u0), Dict(p_dc)),
    tspan_dc;
    guesses = guesses,
)

sol_dc= solve(prob_dc, Rodas5())


# 7. Ensemble solving example (parameter/noise exploration)
# Removed commented out code.



# 8. Harmonic Balance setup


include("../harmonic balance/colocation HB.jl")

# Extract 2nd-order equations and their state variables
eqs, states = jls.get_full_equations(model, jls.t)

# Build DC values for each state from the DC time-domain solution.
# `sol_dc[st]` is the time series for that state; take the final value as the DC operating point.
dc_values = [sol_dc[st][end] for st in states]

# Build harmonic balance system around the DC point.
# Assumes harmonic_equation(eqs, states, t, ω, N; dc_values=...) is implemented.

harmonic_sys, harmonic_states =
    jls.harmonic_equation(eqs, states, jls.t, jls.I1.ω, 3; dc_values = dc_values)

@named ns = NonlinearSystem(harmonic_sys)
hb_sys = structural_simplify(ns; fully_determined = true, check_consistency = true)

state_syms = collect(unknowns(hb_sys))


# 9. Sweep over drive frequency using HB
#work in progress/ confused about this section 
I₀ = 1e-6
R₀ = 5.0
Id = 0.05e-6
ωc = sqrt(2pi * I₀ / (jls.Φ₀ * 1000.0e-15)) / (2pi)

_ = 1 / (2pi * sqrt(1000.0e-12 * 1000.0e-15))  # keep original line

# frequency grid (rad/s)
# Resonant frequency approx 8.3 GHz. Scanning 1-15 GHz.
ω_vec = 2pi * (1:0.2:15) * 1e9
Nω    = length(ω_vec)

# preallocate responses
solution1 = zeros(Nω)
solution2 = zeros(Nω)
solution_gain = zeros(Nω)

# continuation seed: one entry per unknown in HB system
u0_prev = zeros(length(state_syms))

for (i, drive_freq) in enumerate(ω_vec)
    # parameter set at this drive frequency
    ps = merge(Dict(p), Dict(
        jls.I1.ω  => drive_freq,
        jls.I1.I  => 10e-9,      # Match the small signal drive
        jls.J1.R  => 10000.0,    # High internal resistance
    ))

    # use previous HB solution as initial guess (pseudo-continuation)
    state_guess = isempty(state_syms) ? Dict() : Dict(state_syms .=> u0_prev)
    prob_map    = merge(state_guess, ps)

    hb_prob = NonlinearProblem(hb_sys, prob_map;
                               allow_incomplete = true,
                               check_length     = false)

    hb_sol = solve(hb_prob)

    # update continuation seed
    u0_prev .= hb_sol.u

    #jj curent amp and fundamental phase amp
    A1 = hb_sol[ns.A_1]
    B1 = hb_sol[ns.B_1]
    ampϕ = sqrt(A1^2 + B1^2) 

    I0_val = ps[jls.J1.I0]
    ampI = I0_val * ampϕ #need to take derivative of this, find current whats flowing through components
    solution1[i]   = ampϕ
    solution2[i] = ampI

end

freqs = collect(ω_vec) ./ (2pi)

# sanity checks
@show length(freqs)
@show length(solution1)
@show length(solution2)

# phase amplitude vs frequency
p2 = plot(
    freqs, solution1;
    xlabel = "Frequency (Hz)",
    ylabel = "Phase amplitude (rad)",
    title  = "JJ phase response (fundamental)",
    label  = "phase amplitude",
    lw     = 2,
)
savefig(p2, "frequency_response_phase.png")

# overlay current amplitude on a new figure (or the same, up to you)
p3 = plot(
    freqs, solution2;
    xlabel = "Frequency (Hz)",
    ylabel = "Current amplitude (A)",
    title  = "JJ current response (approx.)",
    label  = "current amplitude",
    lw     = 2,
)
savefig(p3, "frequency_response_current.png")

