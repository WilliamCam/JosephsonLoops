# Doubly pumped Josephson parametric amplifier: the same circuit as RLC-ports.jl, driven by
# two pumps through the one port using the current source's second tone. This is the two tone
# benchmark circuit.
#
# The tone ratio decides the collocation grid, and that choice is made in the backend: pass
# the two frequencies as plain numbers through `tones` and a commensurate one dimensional
# grid is used when the ratio is close to a small rational. This example's ratio is 93:97,
# so it always takes that grid. Mixing products need `intermod_order` of at least 3, since
# the conversion between the two pumps lives in those slots.
#
# The pump amplitude is ramped by hand because both tones have to rise together and a
# parameter sweep moves one parameter per axis.

using JosephsonLoops
using Symbolics
using ModelingToolkit
using NonlinearSolve
using DelimitedFiles
const jls = JosephsonLoops

loops = [["P1", "C1", "J1"]]
circuit = jls.process_netlist(loops)
model, u0, guesses = jls.build_circuit(circuit)

I₀ = jls.Φ₀/(2π*1000.0e-12)          # α=1 → Lj = 1000 pH
R₀ = 10.0e3
ωc = R₀*I₀/(jls.Φ₀/2π)

Ip = 2 * 1.7 * 0.00565e-6 / I₀        # docs drive: 0.00565 uA * 1.7 per pump. Their current is
                                      # a one sided amplitude, so the value here is twice it.
ps = Dict{Num,Float64}(
    jls.P1.source.ω  => 2π*4.65001e9/ωc,             # docs pump 1 = 93·ω0 (declared)
    jls.P1.source.I  => Ip,
    jls.P1.source.ω₂ => (97/93)*2π*4.65001e9/ωc,     # docs pump 2 = 97·ω0 (= 4.85001 GHz
    jls.P1.source.I₂ => Ip,                          #  to within 430 Hz, ratio exact)
    jls.P1.Rₙ.r      => 50.0/R₀,
    jls.C1.βc        => 100.0e-15*R₀*ωc,
    jls.J1.βc        => 1000.0e-15*R₀*ωc,
    jls.J1.r         => 1e8/R₀,
    jls.J1.α         => 1.0,
)
Ω_vec = collect(2π*(4.5:0.001:5.0)*1e9/ωc)     # the docs grid

# the commensurate integers are derived in the backend from the tone ratio:
# 4.65:4.85 rationalizes to 93:97, snapping tone 2 by 430 Hz (reported via @info)
@time sys_mt = jls.HarmonicSystem(model, (jls.P1.source.ω, jls.P1.source.ω₂), 2,
    determine_jacobian = true, intermod_order = 3, tones = (4.65001e9, 4.85001e9))

# pump working point: an amplitude ramp with Levenberg-Marquardt, carrying the solution forward
sysu = jls.unknowns(sys_mt.system)
Upump = fill(0.0, length(sysu))
for frac in (0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 1.0)
    p = copy(ps)
    p[jls.P1.source.I]  = frac * Ip
    p[jls.P1.source.I₂] = frac * Ip
    prob = NonlinearProblem(sys_mt.system, merge(Dict(sysu .=> Upump), p))
    sol = solve(prob, LevenbergMarquardt(); maxiters = 2000)
    # only a converged point is carried forward: one unconverged state poisons every later
    # warm start, and the residual is the only thing that says so
    resid = maximum(abs.(sol.resid))
    resid < 1e-9 || error("pump ramp did not converge at frac = $frac, residual $resid")
    global Upump = real.(sol.u)
    println("ramp ", frac, ": residual ", round(resid, sigdigits = 2))
end

δU  = jls.perturbation_response(sys_mt, jls.P1.source.I, ps, amplitude = ps[jls.P1.source.I])
lin = jls.LinearisedProblem(sys_mt, ps, δU, Ω_vec, U₀ = Upump)
jls.solve!(lin)

I_sig = jls.get_solution(lin, jls.P1.i, 1) .* I₀
V_sig = jls.get_solution(lin, jls.P1.dφ, 1) .* R₀ .* I₀
Z0 = 50.0
a = @. 0.5*(V_sig + Z0*I_sig)/sqrt(Z0)
b = @. 0.5*(V_sig - Z0*I_sig)/sqrt(Z0)
gain = 20 .* log10.(abs.(b ./ a))

fs = Ω_vec .* ωc ./ (2π*1e9)          # the same grid the solve used, in GHz
p_mt = jls.plot(fs, gain, lw = 2, label = false,
    xlabel = "Frequency (GHz)", ylabel = "|S₁₁| (dB)",
    title = "Doubly pumped JPA gain")

println("JosephsonLoops: peak ", round(maximum(gain), digits=3), " dB @ ", fs[argmax(gain)], " GHz")

# pkgdir resolves to the repository root whatever the working directory is. A path relative
# to @__DIR__ silently misses when this code is pasted into a REPL.
mit_csv = joinpath(pkgdir(jls), "mit_jpa_multitone.csv")
if isfile(mit_csv)
    mit = readdlm(mit_csv, ',')
    Δ = gain .- Float64.(mit[:,2])
    println("JosephsonCircuits.jl: peak ", round(maximum(Float64.(mit[:,2])), digits=3), " dB @ ",
            mit[argmax(Float64.(mit[:,2])), 1], " GHz")
    println("band max|Δ| = ", round(maximum(abs.(Δ)), digits=3), " dB, median |Δ| = ",
            round(sort(abs.(Δ))[end÷2+1], sigdigits=3), " dB")
else
    @warn "reference not found at $mit_csv: generate it with mit_jpa_multitone_export.jl"
end
display(p_mt)
# regenerates the figure used in the docs
jls.savefig(p_mt, joinpath(pkgdir(jls), "docs", "images", "jpa-doubly-pumped.png"))

