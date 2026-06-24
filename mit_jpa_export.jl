# Run this WITH THE JosephsonCircuits PROJECT to export the MIT JPA gain curve to a CSV
# that RLC-ports.jl reads for the overlay:
#   cd ../JosephsonCircuits-MIT/JosephsonCircuits.jl
#   julia --project=. ../../JosephsonLoops/mit_jpa_export.jl
#
# This is the lossless README JPA (matches our J1.R=1e12 example). Conventions: their
# current=Ip is a one-sided amplitude (≡ our 2·Ip), and gain is the phase-preserving S(0,0).
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

ws = 2*pi*(4.5:0.001:5.0)*1e9            # same sweep as RLC-ports.jl's Ω_vec
wp = (2*pi*4.75001*1e9,)
Ip = 0.00565e-6
sources = [(mode=(1,), port=1, current=Ip)]

jpa = hbsolve(ws, wp, sources, (8,), (16,), circuit, circuitdefs)
S = jpa.linearized.S(outputmode=(0,), outputport=1, inputmode=(0,), inputport=1, freqindex=:)
gain_dB = 10 .* log10.(abs2.(S))

out = raw"C:\Users\Divyank Sharma\OneDrive\Desktop\thesis\JosephsonLoops\mit_jpa.csv"
writedlm(out, [collect(ws) ./ (2π*1e9)  gain_dB], ',')
println("wrote $out  ($(length(gain_dB)) points)  peak $(round(maximum(gain_dB), digits=2)) dB @ ",
        round(ws[argmax(gain_dB)]/(2π*1e9), digits=4), " GHz")
