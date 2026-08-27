
using JosephsonLoops
using Symbolics
using ModelingToolkit
using DelimitedFiles
using Revise 
const jls = JosephsonLoops

loops = [["P1", "C1", "J1"]]
circuit = jls.process_netlist(loops)
model, u0, guesses = jls.build_circuit(circuit)

I₀ = jls.Φ₀/(2π*1000.0e-12)          # α=1 → Lj = 1000 pH
R₀ = 10.0e3
ωc = R₀*I₀/(jls.Φ₀/2π)

Ip = 2 * 1.7 * 0.00565e-6 / I₀        # docs drive: Ip = 0.00565 uA * 1.7 per pump
                                      # (their one-sided amplitude = 2x our I*sin)
ps = Dict(
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

@time sys_mt = jls.HarmonicSystem(model, (jls.P1.source.ω, jls.P1.source.ω₂), 2,
    determine_jacobian = true, intermod_order = 3, commensurate = (93, 97))

# pump working point: amplitude ramp with LM, carrying the solution forward
sysu = jls.unknowns(sys_mt.system)
Upump = fill(0.0, length(sysu))
for frac in (0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 1.0)
    p = copy(ps)
    p[jls.P1.source.I]  = frac * Ip
    p[jls.P1.source.I₂] = frac * Ip
    prob = jls.NonlinearProblem(sys_mt.system, merge(Dict(sysu .=> Upump), p))
    sol = jls.ModelingToolkit.solve(prob, jls.NonlinearSolve.LevenbergMarquardt(); maxiters = 2000)
    global Upump = real.(sol.u)
    println("ramp ", frac, ": retcode = ", sol.retcode)
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

fs = collect(4.5:0.001:5.0)
p_mt = jls.plot(fs, gain, lw=2, label="JosephsonLoops",
    xlabel="Frequency (GHz)", ylabel="|S₁₁| (dB)", title="Doubly-pumped JPA — benchmark")

# pkgdir resolves to the repo root regardless of the REPL's working directory —
# @__DIR__-relative paths silently miss when this code is pasted into a REPL.
mit_csv = joinpath(pkgdir(jls), "mit_jpa_multitone.csv")
if isfile(mit_csv)
    mit = readdlm(mit_csv, ',')
    jls.plot!(p_mt, mit[:,1], mit[:,2], lw=2, ls=:dash, label="JosephsonCircuits.jl")
    Δ = gain .- Float64.(mit[:,2])
    println("ours: peak ", round(maximum(gain), digits=3), " dB @ ", fs[argmax(gain)], " GHz")
    println("MIT:  peak ", round(maximum(Float64.(mit[:,2])), digits=3), " dB @ ",
            mit[argmax(Float64.(mit[:,2])), 1], " GHz")
    println("band max|Δ| = ", round(maximum(abs.(Δ)), digits=3), " dB, median |Δ| = ",
            round(sort(abs.(Δ))[end÷2+1], sigdigits=3), " dB")
else
    @warn "reference not found at $mit_csv — run mit_jpa_multitone_export.jl with the JosephsonCircuits-MIT project"
end
display(p_mt)

