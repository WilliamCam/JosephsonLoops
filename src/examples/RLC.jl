# In this example we model a parallel RLC circuit driven by a voltage source
# Load JosephsonLoops package
using JosephsonLoops
using DifferentialEquations
using Plots

const jls = JosephsonLoops

loops = [
    ["I1", "R1"], #linear
    ["R1", "C1", "J1"], #nonlinear
]
ext_flux = [false, true]

circuit = jls.process_netlist(loops; ext_flux = ext_flux)

# circuit model ODAE system and initial condition vector are created.
model, u0, guesses = jls.build_circuit(circuit)

# 4. Set parameter values
p = [
    jls.I1.ω  => 2pi * 8.3e9,
    jls.I1.I  => 0.5e-6,     
    jls.C1.C  => 100.0e-15,
    jls.J1.C  => 1000.0e-15,
    jls.J1.I0 => 1e-6,
    jls.R1.R  => 50.0,     
    jls.J1.R  => 50.0,    
    jls.Φₑ2.Φₑ => 0.5,
]

# 5. Time-domain simulation via ODEProblem
tspan = (0.0, 1e-6)
prob = ODEProblem(
    model,
    merge(Dict(u0), Dict(p)),
    tspan;
    guesses = guesses,
)

sol = solve(prob, Rodas5())

p1 = plot(sol[jls.C1.i][end-400:end], title = "Transient Time Plot", xlabel = "t", ylabel = "I_C1")
savefig(p1, "transient_plot.png")


# Turn off AC drive to find DC point
p_dc = jls.set_param(p, jls.I1.I, 0.0)
tspan_dc = (0.0, 5e-6)

prob_dc = ODEProblem(
    model,
    merge(Dict(u0), Dict(p_dc)),
    tspan_dc;
    guesses = guesses,
)

sol_dc = solve(prob_dc, Rodas5())



# 7. Ensemble solving example (parameter/noise exploration)

x = jls.ensemble_fsolve(model, u0, tspan, (0.1, 10.0), p, jls.loop1, jls.R1, units = "amps", Ntraj = 500)   
p2 = plot(x.u)
savefig(p2, "ensemble_fsolve_plot.png")





# 8. Harmonic Balance setup
include("../harmonic balance/colocation HB.jl")

# Extract 2nd-order equations and their state variables
eqs, states = jls.get_full_equations(model, jls.t)

# 1. Get all fixed DC values from the transient simulation
dc_values = Vector{Any}([sol_dc[st][end] for st in states])


dc_values[1] = nothing

harmonic_sys, harmonic_states =
    jls.harmonic_equation(eqs, states, jls.t, jls.I1.ω, 3; dc_values = dc_values)

@named ns = NonlinearSystem(harmonic_sys)
hb_sys = structural_simplify(ns; fully_determined = true, check_consistency = true)

state_syms = collect(unknowns(hb_sys))

I₀ = 1e-6
R₀ = 5.0

Id_vec = range(5e-8, 1e-6; length = 10) 
ωc = sqrt(2pi * I₀ / (jls.Φ₀ * 1000.0e-15)) / (2pi)
ω_vec = 2pi * (1:0.2:15) * 1e9
Nω    = length(ω_vec)

# preallocate responses
solution1 = zeros(length(Id_vec), Nω)
solution2 = zeros(length(Id_vec), Nω)

for (k, Id_val) in enumerate(Id_vec)

    println("Sweeping Id = $Id_val")

    # reset continuation seed for each new drive current loop
    u0_prev = zeros(length(state_syms))

    for (i, drive_freq) in enumerate(ω_vec)

        ps = merge(Dict(p), Dict(
            jls.I1.ω => drive_freq,
            jls.I1.I => Id_val,
            jls.J1.R => 10000.0,
        ))

        # Base guess from previous step (continuation)
        current_guess = isempty(state_syms) ? Dict() : Dict(state_syms .=> u0_prev)
        prob_map = merge(current_guess, ps)

        hb_prob = NonlinearProblem(
            hb_sys,
            prob_map;
            allow_incomplete = true,
            check_length     = false,
        )
        hb_sol = solve(hb_prob)
        u0_prev .= hb_sol.u # Update continuation vector
        # Extract amplitudes
        A1 = hb_sol[ns.A_1]
        B1 = hb_sol[ns.B_1]
        ampϕ = sqrt(A1^2 + B1^2)

        I0_val = ps[jls.J1.I0]
        ampI   = I0_val * ampϕ
       
        solution1[k, i] = ampϕ
        solution2[k, i] = ampI  
    end
end
freqs = collect(ω_vec) ./ (2pi)

p3 = plot(xlabel="Frequency (Hz)", ylabel="ampϕ", title="HB sweep vs frequency")
for (k, Id_val) in enumerate(Id_vec)
    plot!(p3, freqs, solution1[k, :], label="Id=$(Id_val)")
end
display(p3)

# Create the heatmap
p_heatmap = heatmap(
    freqs,              # x-axis: Frequency (Hz)
    Id_vec,             # y-axis: Drive Current (A)
    log10.(solution1),          # z-axis: Amplitude Matrix
    xlabel = "Frequency (Hz)",
    ylabel = "Drive Current (A)",
    title = "HB Amplitude Heatmap",
    color = :viridis    
)

display(p_heatmap)
savefig(p_heatmap, "hb_heatmap.png")

freqs_GHz = collect(ω_vec) ./ (2pi * 1e9)
# 1 Reshape vectors for matrix broadcasting
Id_col = reshape(Id_vec, :, 1)   # Column Vector (N x 1)
w_row  = reshape(ω_vec, 1, :)    # Row Vector    (1 x M)

# 2. Calculate Voltage (Vi) and Net Current (Ii) at the port

Vi = @. (jls.Φ₀ / (2*pi)) * w_row * solution1
Ii = @. Id_col - (jls.Φ₀ / (2*pi)) * w_row * solution1 / 50.0
ai = @. 0.5 * (Vi + 50 * Ii) / 50.0
bi = @. 0.5 * (Vi - 50 * conj(Ii)) / 50.0 
S11_mag = @. abs(bi / ai)

# PLOTTING
p_s11 = plot(
    xlabel = "Frequency (GHz)",      # Changed to GHz
    ylabel = "S11 Log Scale", 
    title  = "S11 Resonance Dips",
    legend = :bottomright,
    yscale = :log10,                 
    ylims  = (1e-3, 10),            # Fix limits: 0.001 to 1.1
    yticks = [1, 0.1, 0.01, 0.001],  # Explicit powers of 10
)

for (k, Id_val) in enumerate(Id_vec)
    # Convert Current to µA for cleaner legend labels
    Id_uA = round(Id_val * 1e6, sigdigits=2)
    
    plot!(p_s11, freqs_GHz, S11_mag[k, :], 
          label = "Id = $(Id_uA) µA", 
          lw = 1.5) # Thicker lines
end

display(p_s11)
savefig(p_s11, "S11.png")



# using BifurcationKit 

bif_par = jls.loop1.sys.ω
p_start = [
        jls.loop1.sys.ω => ω_vec[1]
        jls.loop1.sys.I => Id
        jls.C1.sys.C => 100.0e-15
        jls.J1.sys.L => 1000.0e-12
        jls.J1.sys.C=> 1000.0e-15
        jls.J1.sys.I0 => 0.3e-6
        jls.R1.sys.R => 50.0
        jls.J1.sys.R => 1000.0
        jls.loop2.sys.Φₑ => 0.5
]
u0_guess = sol
plot_var = sys.A[2]

bprob = BifurcationProblem(sys,
    u0_guess,
    p_start,
    bif_par;
    plot_var = plot_var,
    jac = false)

p_span = (ω_vec[1], ω_vec[end])
opts_br = ContinuationPar(nev = 2,
    p_min = p_span[1],
    p_max = p_span[2])

bf = bifurcationdiagram(bprob, PALC(), 2)
