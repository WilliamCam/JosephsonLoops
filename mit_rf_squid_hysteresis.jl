# Run WITH THE JosephsonCircuits PROJECT — companion probe to
# src/examples/rf-squid-coupler-hysteresis.jl:
#   julia --project=..\JosephsonCircuits-MIT\JosephsonCircuits.jl mit_rf_squid_hysteresis.jl
#
# Same betaL = 2 coupler, same up/down flux sweep. hbsolve cold-starts its Newton
# solve on every call (x0 = nothing is hardcoded in the wrapper), so it has no branch
# memory: up- and down-sweeps must return IDENTICAL curves (or fail to converge in
# the bistable window) — i.e. the physical hysteresis of a betaL > 1 rf-SQUID cannot
# be reproduced through the supported API.
using JosephsonCircuits
using Printf

const phi0 = 2.067833848e-15
Ic  = 5.0e-6
LJ0 = phi0/(2*pi*Ic)
betaL = 2.0

@variables R Ls Lj Lb Kb Rsh
circuit = [
    ("P1","1","0",1),
    ("R1","1","0",R),
    ("L1","1","0",Ls),
    ("Lj1","1","2",Lj),
    ("Rsh1","1","2",Rsh),
    ("L2","2","0",Ls),
    ("P2","2","0",2),
    ("R2","2","0",R),
    ("L3","3","0",Lb),
    ("K1","L1","L3",Kb),
    ("P3","3","0",3),
    ("R3","3","0",1000.0),
]
Lsval = 0.5*betaL*LJ0
Kval = 0.1; Lbval = 100.0e-12
M = Kval*sqrt(Lsval*Lbval)
circuitdefs = Dict(R=>50.0, Ls=>Lsval, Lj=>LJ0, Lb=>Lbval, Kb=>Kval, Rsh=>1.0e6)

ws = [2*pi*5.0e9]
wp = (2*pi*1.0e9,)

f_up = collect(-0.25:0.005:1.25)
f_down = reverse(f_up)

function s21_at(f)
    sources = [(mode=(0,), port=3, current=f*phi0/M)]
    try
        sol = hbsolve(ws, wp, sources, (4,), (4,), circuit, circuitdefs,
                      dc = true, threewavemixing = true, fourwavemixing = true)
        S21 = sol.linearized.S(outputmode=(0,), outputport=2,
                               inputmode=(0,), inputport=1, freqindex=:)
        return 10*log10(abs2(S21[1]))
    catch
        return NaN
    end
end

up   = [s21_at(f) for f in f_up]
down = [s21_at(f) for f in f_down]

d = abs.(up .- reverse(down))
nvalid = count(!isnan, d)
@printf("max |up - down| = %.3g dB over %d valid points (NaN failures: %d up, %d down)\n",
        maximum(filter(!isnan, d)), nvalid, count(isnan, up), count(isnan, down))
println("identical curves (no hysteresis possible: cold start each point) or failures — either")
println("way the betaL > 1 branch structure is not accessible through the supported API.")
for f in (-0.109, 0.0, 0.109, 0.391, 0.5, 0.609, 1.0)
    j = argmin(abs.(f_up .- f))
    @printf("  f=%6.3f  up %8.2f dB   down %8.2f dB\n", f_up[j], up[j], reverse(down)[j])
end
