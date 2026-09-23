# Josephson parametric amplifier: a junction connected by a coupling capacitor, driven and
# read out through a single 50 ohm port.

using JosephsonLoops
using Plots

# ---- circuit ----------------------------------------------------------------------
loops = [["P1", "C1", "J1"]]
circuit = process_netlist(loops)
jpa, u0, guesses = build_circuit(circuit)

# ---- normalisation and parameters --------------------------------------------------
I₀ = Φ₀/(2π*1000.0e-12)          # critical current of a 1000 pH junction
R₀ = 10.0e3
ωc = R₀*I₀/(Φ₀/2π)
Z0 = 50.0
f_pump = 4.75001e9
I_pump = 11.3e-9/I₀

ps = Dict(
    jpa.P1.source.ω => 2π*f_pump/ωc,
    jpa.P1.source.I => I_pump,
    jpa.P1.Rₙ.r     => 50.0/R₀,
    jpa.C1.βc       => 100.0e-15*R₀*ωc,
    jpa.J1.βc       => 1000.0e-15*R₀*ωc,
    jpa.J1.r        => 1.0e8/R₀
)

# ---- time domain solution ----------------------------------------------------------------

ps_td = merge(ps, Dict(jpa.P1.source.ω => 2π*100e6/ωc,
                       jpa.P1.source.I => 0.00565e-6/I₀))

#timescale τ is normalised                       
tspan = (0.0, 1e-6) .* ωc

tsol  = tsolve(jpa, guesses, ps_td, tspan; guesses = guesses)
p_td  = plot(tsol[jpa.C1.i][end-400:end] .* I₀, legend = false,
                 xlabel = "Sample", ylabel = "I(C1) (A)", title = "Transient at 100 MHz")

f_vec = collect(4.5:0.001:5.0)
ω_vec = collect(2π .* f_vec .* 1e9 ./ ωc)

# ---- Construct harmonic system ------------------------------------------------
#harmonic system object is built once and contains all symbolic information about the system
sys = HarmonicSystem(jpa, jpa.P1.source.ω, 2, determine_jacobian = true)
#from here can perform both nonlinear harmonic balance, and linearisation

# ---- steady state harmonic balance ----------------------------------------
prob = HarmonicProblem(sys, ps, parameter_sweep = [jpa.P1.source.ω => ω_vec])
solve!(prob)

# As the package is fully symbolic we can extract any desired expression, 
# as a utility we include some simple scattering parameter formulas
s11_expr  = get_HB_scattering_matrix(jpa, '1', '1')[1]
# get solution can take any symbolic input expression as long as all varaibles are defined within HarmonicSystem,
# it will output the numeric solution of the input formula applied to the problem assuming solve! has been called
S11_drive = get_solution(prob, s11_expr, (1, 0))


# ---- small signal: gain around the pump ----------------------------------------------
# We can instead linearise about the pump tone to observe the small signal gain

δU  = perturbation_response(sys, jpa.P1.source.I, ps, amplitude = ps[jpa.P1.source.I])
lin = LinearisedProblem(sys, ps, δU, ω_vec)
solve!(lin)


S11 = get_solution(lin, s11_expr, 1)

gain_dB = 20*log10.(abs.(S11))
plot(f_vec, gain_dB, lw = 2, label = false,
                  xlabel = "Frequency (GHz)", ylabel = "|S₁₁| (dB)",
                  title = "JPA gain")


