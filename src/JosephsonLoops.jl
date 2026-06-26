module JosephsonLoops

#Specific packages loading for reduce precompilation time
using ModelingToolkit, Plots, DifferentialEquations, Symbolics, DataStructures, LinearAlgebra, NonlinearSolve

#Internal API
include("build_circuit/component_library.jl")
include("build_circuit/circuit_model.jl")
include("build_circuit/utils.jl")
include("harmonic balance/harmonic_system.jl")   # defines HarmonicSystem/HarmonicProblem structs
include("harmonic balance/utils.jl")             # perturbation_response annotates ::HarmonicSystem, so load after the structs
include("harmonic balance/get_phasor.jl")
include("harmonic balance/colocation.jl")
end # module JLoop
