using Symbolics
using SymbolicUtils
using ModelingToolkit
using QuestBase
using DifferentialEquations
using Plots
using JosephsonLoops
#Linear analysis adapted from Kosata 2022 thesis

#set up example differential equation in time domain. We will use duffing oscillator example from p.59
@variables t x(t)
@parameters  α ω ω0 F γ η
diff_eq = Differential(t)(Differential(t)(x)) + ω0^2*x + α*x^3 + η*Differential(t)(x)*x^2 + γ*Differential(t)(x) - F*cos(ω*t) ~ 0


