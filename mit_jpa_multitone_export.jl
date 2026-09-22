# Run this WITH THE JosephsonCircuits PROJECT to export the doubly pumped JPA gain curve to
# the CSV that multitone-jpa-benchmark.jl reads for the overlay:
#   julia --project=../JosephsonCircuits-MIT/JosephsonCircuits.jl mit_jpa_multitone_export.jl
#
# Same circuit as mit_jpa_export.jl, driven by two pumps through the one port instead of one.
# Conventions: their `current` is a one-sided amplitude, so our I*cos amplitude is twice it,
# and gain is the phase-preserving S(0,0) reflection.
using JosephsonCircuits
using DelimitedFiles

@variables R Cc Lj Cj
circuit = [
    ("P1","1","0",1),
    ("R1","1","0",R),
    ("C1","1","2",Cc),
    ("Lj1","2","0",Lj),
    ("C2","2","0",Cj)]
circuitdefs = Dict(Lj => 1000.0e-12, Cc => 100.0e-15, Cj => 1000.0e-15, R => 50.0)

ws = 2*pi*(4.5:0.001:5.0)*1e9            # same sweep as the example's Ω_vec
wp = (2*pi*4.65001e9, 2*pi*4.85001e9)    # the two pumps, ratio 93:97
Ip = 1.7 * 0.00565e-6                    # per pump
sources = [(mode=(1,0), port=1, current=Ip), (mode=(0,1), port=1, current=Ip)]

jpa = hbsolve(ws, wp, sources, (8,8), (16,16), circuit, circuitdefs)
S = jpa.linearized.S(outputmode=(0,0), outputport=1, inputmode=(0,0), inputport=1, freqindex=:)
gain_dB = 10 .* log10.(abs2.(S))

out = joinpath(@__DIR__, "mit_jpa_multitone.csv")
writedlm(out, [collect(ws) ./ (2π*1e9)  gain_dB], ',')
println("wrote $out  ($(length(gain_dB)) points)  peak $(round(maximum(gain_dB), digits=3)) dB @ ",
        round(ws[argmax(gain_dB)]/(2π*1e9), digits=4), " GHz")
