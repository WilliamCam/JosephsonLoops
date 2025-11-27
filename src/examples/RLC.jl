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


function sanitize_guesses(guesses, u0)
    g = Dict{Any, Any}()
    for (k, v) in guesses
        try
            # if the guess maps a variable to itself (or is identical), replace it with a numeric guess
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
            # if comparison fails for some symbolic type, fallback to the original value
            g[k] = v
        end
    end
    return g
end

guesses = sanitize_guesses(guesses, u0)

function set_param(pairs, key, val)
    out = Pair{Num,Float64}[]
    for (k, v) in pairs
        if isequal(k, key)
            push!(out, key => val)
        else
            push!(out, k => v)
        end
    end
    return out
end



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
    jls.I1.ω  => 100e6 * 2pi,
    jls.I1.I  => 1e-4,
    jls.C1.C  => 100.0e-15,
    jls.J1.C  => 1000.0e-15,
    jls.J1.I0 => 1e-6,
    jls.R1.R  => 50.0,
    jls.J1.R  => 1.0,
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

plot(sol[jls.C1.i], title = "Transient Time Plot", xlabel = "t", ylabel = "I_C1")


# 7. Ensemble solving example (parameter/noise exploration)

# x = jls.ensemble_fsolve(
#     model,
#     u0,
#     tspan,
#     (0.1, 10.0),
#     p,
#     jls.loop1,
#     jls.R1;
#     units   "amps",
#     Ntraj  = 500,
# )

# plot(x.u, title = "Ensemble trajectories")


# 8. Harmonic Balance setup (DC via algebraic steady state, no time integration)


include("../harmonic balance/colocation HB.jl")

# second-order equations and their state variables
eqs, states = jls.get_full_equations(model, jls.t)

# Build a DC system by forcing all time derivatives to zero:
D  = Differential(jls.t)
dc_eqs = eqs
for st in states
    dc_eqs = substitute(dc_eqs, Dict(
        D(D(st)) => 0,   # d²x/dt² = 0
        D(st)    => 0,   # dx/dt   = 0
    ))
end

# Nonlinear algebraic system: dc_eqs(states, params) = 0
@named dc_ns = NonlinearSystem(dc_eqs, states, parameters(model))

dc_sys = structural_simplify(dc_ns;
    fully_determined = true,
    check_consistency = false,
)
dc_state_syms = collect(unknowns(dc_sys))

# 9. Sweep over drive frequency using DC → HB at each step

I₀ = 1e-4
R₀ = 5.0
Id = 0.05e-6
ωc = sqrt(2pi * I₀ / (jls.Φ₀ * 1000.0e-15)) / (2pi)

_ = 1 / (2pi * sqrt(1000.0e-12 * 1000.0e-15))  # keep original line

# frequency grid: 1–20 GHz (linear here; change to log/whatever if you like)
freqs = range(1e9, 20e9; length = 50)   # Hz
ω_vec = 2pi .* freqs                        # rad/s
Nω    = length(ω_vec)

solution1 = zeros(Nω)   # phase amplitude
solution2 = zeros(Nω)   # current amplitude

# continuation seeds for DC and HB
dc_guess = zeros(length(dc_state_syms))   # initial guess for DC solve
hb_guess = nothing                 # we’ll set this after first HB solve



for (i, ωd) in enumerate(ω_vec)
    # DC solve at this frequency: AC drive OFF (I1.I = 0)
    ps_dc = Dict(
        jls.I1.ω  => ωd,          # might drop out of DC eqs, but fine
        jls.I1.I  => 0.0,         # NO AC drive for DC operating point
        jls.C1.C  => 100.0e-15,
        jls.J1.C  => 1000.0e-15,
        jls.J1.I0 => 1e-6,
        jls.R1.R  => 50.0,
        jls.J1.R  => 1.0,
        jls.Φₑ2.Φₑ => 0.5,
    )

    dc_state_guess = isempty(dc_state_syms) ? Dict() : Dict(dc_state_syms .=> dc_guess)
    dc_prob = NonlinearProblem(dc_sys, merge(dc_state_guess, ps_dc, Dict(jls.t => 0.0)))
    dc_sol  = solve(dc_prob)

    # 
    dc_guess .= dc_sol.u

    
    dc_values = Float64[]
    for st in states
        val = try
            dc_sol[st]
        catch
            try
                dc_sol[ModelingToolkit.observed(dc_sys, st)]
            catch
                0.0
            end
        end
        push!(dc_values, float(val))
    end

    #  build HB system around this DC point, with AC drive ON ---
    ps_hb = Dict(
        jls.I1.ω  => ωd,
        jls.I1.I  =>-1e-6,       # AC drive amplitude
        jls.C1.C  => 100.0e-15,
        jls.J1.C  => 1000.0e-15,
        jls.J1.I0 => 1e-6,
        jls.R1.R  => 50.0,
        jls.J1.R  => 1.0,
        jls.Φₑ2.Φₑ => 0.5,
    )

    Nh = 3  # number of harmonics
    harmonic_sys, harmonic_states =
        jls.harmonic_equation(eqs, states, jls.t, jls.I1.ω, Nh;
                              dc_values = dc_values)

    @named hb_ns = NonlinearSystem(harmonic_sys) 
    hb_ns = structural_simplify(hb_ns;
                                      fully_determined = true,
                                      check_consistency = true)

    state_syms = collect(unknowns(hb_ns))

    # continuation guess for HB:
    if hb_guess === nothing
        hb_guess = zeros(length(state_syms))
    end

    hb_guess_dict = Dict(state_syms .=> hb_guess)
    prob_map      = merge(hb_guess_dict, ps_hb)

    hb_prob = NonlinearProblem(hb_ns, prob_map;
                               allow_incomplete = false,
                               check_length     = false)

    hb_sol = solve(hb_prob)

    hb_guess .= hb_sol.u  # update HB continuation seed

    # extract JJ phase & current amplitudes from fundamental 
    A1 = hb_sol[hb_ns.A_1]
    B1 = hb_sol[hb_ns.B_1]
    ampϕ = sqrt(A1^2 + B1^2)      # fundamental phase amplitude

    I0_val = ps_hb[jls.J1.I0]
    ampI   = I0_val * ampϕ        # crude current amplitude estimate

    solution1[i] = ampϕ
    solution2[i] = ampI
end
# frequency axis in GHz for plotting
freqs_GHz = freqs ./ 1e9

@show length(freqs_GHz)
@show length(solution1)
@show length(solution2)

# phase amplitude vs frequency
plot(
    freqs_GHz, solution1;
    xlabel = "Frequency (GHz)",
    ylabel = "Phase amplitude (rad)",
    title  = "JJ phase response (fundamental)",
    label  = "phase amplitude",
    lw     = 2,
)

# current amplitude vs frequency
plot(
    freqs_GHz, solution2;
    xlabel = "Frequency (GHz)",
    ylabel = "Current amplitude (A)",
    title  = "JJ current response (approx.)",
    label  = "current amplitude",
    lw     = 2,
)




# # 11. Bifurcation analysis with BifurcationKit

# # Choose continuation parameter (drive frequency)


# # Parameter vector for starting point: use original jls symbols,
# # *not* sys.I1.I etc., because hb_sys does not expose those as fields.
# bif_par = jls.loop1.ω
# p_start = [
#     jls.loop1.ω => ω_vec[1]
#     jls.loop1.I => Id
#     jls.C1.C    => 100.0e-15
#     jls.J1.C    => 1000.0e-15
#     jls.J1.I0   => 0.3e-6
#     jls.R1.R    => 50.0
#     jls.J1.R    => 1000.0
#     jls.loop2.Φₑ => 0.0
# ]


# # Use the last HB solution as initial guess for continuation
# u0_guess = copy(u0_prev)

# # Pick an observable for plotting along the branch.
# # We just use the second unknown here; adjust if you want a specific state.
# plot_var = state_syms[2]

# bprob = BifurcationProblem(
#     hb_sys,
#     u0_guess,
#     p_start,
#     bif_par;
#     plot_var = plot_var,
#     jac = false,
# )

# p_span = (ω_vec[1], ω_vec[end])

# opts_br = ContinuationPar(
#     nev   = 2,
#     p_min = p_span[1],
#     p_max = p_span[2],
# )

# bf = bifurcationdiagram(bprob, PALC(), 2, opts_br)
# --- final display section ---
