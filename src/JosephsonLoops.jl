
module JosephsonLoops

using ModelingToolkit, Plots, DifferentialEquations, Symbolics, DataStructures, LinearAlgebra, NonlinearSolve
using SymbolicUtils
using NonlinearSolve
using BenchmarkTools
using StaticArrays


include("build_circuit/component_library.jl")  
include("build_circuit/circuit_model.jl")       
include("build_circuit/utils.jl")             


include("harmonic balance/harmonic_system.jl")  
include("harmonic balance/utils.jl")            
include("harmonic balance/symbolic_rules.jl")   
include("harmonic balance/get_phasor.jl")       
include("harmonic balance/colocation.jl")       


include("linearisation/linear_system.jl")       

@doc """
    Φ₀

The magnetic flux quantum in weber, `h/(2e) = 2.067833848e-15`. The package works in
normalised units where fluxes are measured in `Φ₀/2π`, so this constant is what converts a
junction inductance to a critical current, `I₀ = Φ₀/(2π*Lj)`, and sets the characteristic
frequency `ωc = R₀*I₀/(Φ₀/2π)`.
""" Φ₀


export process_netlist, build_circuit,                                
    tsolve,                                                           
    HarmonicSystem, HarmonicProblem, LinearisedProblem, solve!,       
    get_solution, perturbation_response, get_HB_scattering_matrix,    
    Φ₀
end 
