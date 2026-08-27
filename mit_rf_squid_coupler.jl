# Run WITH THE JosephsonCircuits PROJECT to generate the rf-SQUID coupler reference
# (arXiv:2408.07861 Fig. 8b) from JosephsonCircuits.jl:
#   julia --project=..\JosephsonCircuits-MIT\JosephsonCircuits.jl mit_rf_squid_coupler.jl
#
# VERDICT (probed 2026-08-27): JosephsonCircuits.jl CAN solve this circuit. Their
# component set has no flux-bias element (their README: "3WM and flux-biasing are
# relatively untested"), so external flux must be emulated with a DC current source
# (mode=(0,)) at a bias port, K-coupled into one loop inductor, following their
# flux-pumped JPA README pattern (dc=true, threewavemixing=true). With that idiom all
# points converge, including deep frustration (betaL=1, f=0.5), and the values match
# the ADS paper and the analytic Pi network exactly:
#   f=0 baselines -46.7/-42.7/-39.7 dB and f=0.5 peaks -34.6/-23.6/-0.0 dB
#   for betaL = 0.6/0.8/1.0.
#
# Writes mit_rf_squid_coupler.csv: column 1 = Norm_Ext_Flux, then one S21(dB) column
# per betaL in BETAS. Overlay against src/examples/rf-squid-coupler.jl.
using JosephsonCircuits
using DelimitedFiles
using Printf

const phi0 = 2.067833848e-15
Ic  = 5.0e-6
LJ0 = phi0/(2*pi*Ic)              # 65.8 pH; paper's per-side L = 0.5*betaL*LJ0

@variables R Ls Lj Lb Kb Rsh
circuit = [
    ("P1","1","0",1),
    ("R1","1","0",R),
    ("L1","1","0",Ls),
    ("Lj1","1","2",Lj),
    ("Rsh1","1","2",Rsh),         # paper's 1 MOhm junction shunt
    ("L2","2","0",Ls),
    ("P2","2","0",2),
    ("R2","2","0",R),
    ("L3","3","0",Lb),            # flux-bias coil: loop flux = M*Idc, M = Kb*sqrt(Ls*Lb)
    ("K1","L1","L3",Kb),
    ("P3","3","0",3),
    ("R3","3","0",1000.0),        # their README idiom: bias applied across a resistive port
]

ws = [2*pi*5.0e9]                 # the paper's fixed 5 GHz probe
wp = (2*pi*1.0e9,)                # dummy strong-tone frequency; only the DC source is nonzero

BETAS = collect(0.6:0.1:1.0)
fvals = collect(0.0:0.005:1.0)
Kval  = 0.1
Lbval = 100.0e-12

S21_dB = fill(NaN, length(fvals), length(BETAS))
for (jb, beta) in enumerate(BETAS)
    Lsval = 0.5*beta*LJ0
    M = Kval*sqrt(Lsval*Lbval)
    circuitdefs = Dict(R=>50.0, Ls=>Lsval, Lj=>LJ0, Lb=>Lbval, Kb=>Kval, Rsh=>1.0e6)
    for (jf, f) in enumerate(fvals)
        sources = [(mode=(0,), port=3, current=f*phi0/M)]
        sol = hbsolve(ws, wp, sources, (4,), (4,), circuit, circuitdefs,
                      dc = true, threewavemixing = true, fourwavemixing = true)
        S21 = sol.linearized.S(outputmode=(0,), outputport=2,
                               inputmode=(0,), inputport=1, freqindex=:)
        S21_dB[jf, jb] = 10*log10(abs2(S21[1]))
    end
    @printf("betaL = %.1f: f=0 -> %7.2f dB, f=0.5 -> %7.2f dB, min %7.2f dB\n",
            beta, S21_dB[1, jb], S21_dB[fvals .== 0.5, jb][1], minimum(S21_dB[:, jb]))
end

writedlm(joinpath(@__DIR__, "mit_rf_squid_coupler.csv"), hcat(fvals, S21_dB), ',')
println("wrote mit_rf_squid_coupler.csv")
