using JosephsonLoops
using Revise
using Symbolics
using ModelingToolkit
using BenchmarkTools

const jls = JosephsonLoops
loops = [
    ["P1", "J1"],
    ["J1", "C1", "L2"],
    ["L2", "C2"],
    ["C2", "C3", "Rl"]
]

circuit = jls.process_netlist(loops)
model, u0, guesses = jls.build_circuit(circuit)
I₀ = jls.Φ₀/(2π*1000.0e-12)
R₀ = 10.0e3
ωc = R₀*I₀/(jls.Φ₀/2π)