"""
    JosephsonLoops

Time domain and harmonic balance simulation of lumped element superconducting circuits
containing Josephson junctions, written in loop currents.

A circuit is entered as a list of loops and assembled into an acausal ModelingToolkit model
from a small component library. That one model is either integrated as an initial value
problem with [`tsolve`](@ref) or turned into a harmonic balance problem: the Fourier ansatz is
substituted symbolically, the residuals are sampled on a collocation grid, and the Fourier
coefficients are solved for with the full nonlinearity in place. Around such a working point
the response to a weak probe is a linear problem, which is how amplifier gain and scattering
parameters are computed.

The package is aimed at strongly pumped nonlinear circuits with a modest number of junctions,
such as parametric amplifiers, flux tunable couplers and SQUIDs, and at any driven nonlinear
oscillator that can be written as a differential equation. It is not yet optimised for large
junction arrays.

The typical flow is

    process_netlist -> build_circuit -> HarmonicSystem -> HarmonicProblem -> solve! -> get_solution
                                                       \\-> LinearisedProblem -> solve! -> get_solution

with `tsolve` taking the built model directly for a time domain solve. `get_solution` accepts
any symbolic expression of the model's variables and parameters, not only a state.

See the readme for the loop current formulation, units and normalisation, and worked
examples.
"""
module JosephsonLoops

# ModelingToolkit provides the acausal component models and the circuit assembly; Symbolics
# and SymbolicUtils carry the Fourier ansatz, the collocation residuals and the jacobians of
# the linearised problem, all of which stay symbolic until they are compiled; NonlinearSolve
# solves the harmonic balance system; DifferentialEquations integrates the same model in the
# time domain; Plots is loaded here so that the examples can plot without a second import.
# These are named explicitly rather than pulled in through other packages, which keeps the
# precompilation time down.
using ModelingToolkit, Plots, DifferentialEquations, Symbolics, DataStructures, LinearAlgebra, NonlinearSolve
using SymbolicUtils
using NonlinearSolve
using BenchmarkTools
using StaticArrays

# The files below follow the analysis flow. A netlist of loops is parsed and assembled into a
# ModelingToolkit model from the component library; that model is integrated in the time
# domain, or its equations are turned into a harmonic balance system by substituting the
# Fourier ansatz and sampling a collocation grid; a solved problem is read back through
# symbolic expressions; and the linearisation around a working point gives the small signal
# response used for gain and scattering parameters.

# Circuit model: the acausal components, the netlist parser and the time domain solver.
include("build_circuit/component_library.jl")   # Loop connector, Branch, the components, the port and Φ₀
include("build_circuit/circuit_model.jl")       # process_netlist and build_circuit
include("build_circuit/utils.jl")               # tsolve and the ensemble sweep helpers

# Harmonic balance: the Fourier basis, the collocation grid and the steady state solve.
include("harmonic balance/harmonic_system.jl")  # HarmonicSystem, HarmonicProblem, solve!, pretty printing
include("harmonic balance/utils.jl")            # Fourier basis, projection onto the harmonic frame, jacobians, S parameters
include("harmonic balance/symbolic_rules.jl")   # rewrite rules that drop slow derivatives when linearising
include("harmonic balance/get_phasor.jl")       # get_solution: any expression of the model, read on a solved problem
include("harmonic balance/colocation.jl")       # collocation sampling and the harmonic equations, one and two tone

# Small signal: the linearised response around a working point.
include("linearisation/linear_system.jl")       # LinearisedProblem, perturbation_response, solve!, get_solution

@doc """
    Φ₀

The magnetic flux quantum in weber, `h/(2e) = 2.067833848e-15`. The package works in
normalised units where fluxes are measured in `Φ₀/2π`, so this constant is what converts a
junction inductance to a critical current, `I₀ = Φ₀/(2π*Lj)`, and sets the characteristic
frequency `ωc = R₀*I₀/(Φ₀/2π)`.
""" Φ₀

# The public API. Everything else in the package is internal and may change.
export process_netlist, build_circuit,                                # netlist to model
    tsolve,                                                           # time domain
    HarmonicSystem, HarmonicProblem, LinearisedProblem, solve!,       # harmonic balance and small signal
    get_solution, perturbation_response, get_HB_scattering_matrix,    # reading results, probes, S parameters
    Φ₀

end # module JosephsonLoops
