using Symbolics
using SymbolicUtils


struct HarmonicProblem
    sys  # The algebraic NonlinearSystem
    N::Int
    sys_vars::Vector{Num}
    ωvar::Num
end

struct HarmonicSweepResult
    sweep_var::Num
    sweep_vals::Vector{Float64}
    results::Dict{Num, Vector{Float64}} # Maps symbolic variables to result vectors
    
end

function var_is_in(vars::Vector, target_var::SymbolicUtils.BasicSymbolic{Real})
    ret = false
    for var in vars
        if isequal(var, target_var)
            ret = true
            break
        end
    end
return ret
end

function var_is_in(vars::Vector, target_var::Num)
    ret = false
    for var in vars
        if isequal(var, target_var)
            ret = true
            break
        end
    end
return ret
end

function var_index(vars::Vector, target_var::SymbolicUtils.BasicSymbolic{Real})
    return findfirst(x->isequal(x, target_var),vars)
end

function get_HB_scattering_matrix(model::System,i::Char,j::Char)
    port_i_sym = Symbol('P'*i)
    port_j_sym = Symbol('P'*j)
    @assert (hasproperty(model,port_i_sym)) "Error: Port $i not defined in variable map"
    @assert (hasproperty(model,port_j_sym)) "Error: Port $j not defined in variable map"
    N_ports = maximum([parse(Int,port_index) for port_index in [i,j]])
    @assert (N_ports < 3) "Maximum of 2 port networks supported"
    a = zeros(Num,1,N_ports)
    b = zeros(Num,1,N_ports)
    #TODO: assert port reference impedances are equal
    for k in N_ports
        port_k_sym = Symbol('P'*string(k))
        port_k = getproperty(model, port_k_sym)
        a[k] = 0.5/sqrt(port_k.Rsrc.R)*(port_k.dθ*Φ₀/2π+port_k.Rsrc.R*port_k.i)
        b[k] = 0.5/sqrt(port_k.Rsrc.R)*(port_k.dθ*Φ₀/2π-port_k.Rsrc.R*port_k.i)
    end
    return b / a
end

function only_derivatives(expr, var, tvar)
    # The maximum order of derivative to check against.
    # We must substitute all known derivative terms of `var`.
    sub_rules = Dict(
        Differential(tvar)(var) => 0,
        Differential(tvar)(Differential(tvar)(var)) => 0,
    )
    expr_sub = substitute(expr, sub_rules)
    remaining_vars = get_variables(expr_sub)
    return !var_is_in(remaining_vars, var)
end

function get_derivatives(X, t)
    D = Differential(t)
    dXdt = Symbolics.expand_derivatives(D(X))
    d2Xdt2 = Symbolics.expand_derivatives(D(dXdt))
    return dXdt, d2Xdt2
end

function get_full_equations(model::ModelingToolkit.System, tvar::Num)
    eqs = full_equations(model)
    states = unknowns(model)

    diff2vars = Vector{Num}()
    diffvars = Vector{Num}()
    remove_idxs = Int[]
    for (i, eq) in enumerate(eqs)
        vars = get_variables(eq.rhs)
        if length(vars) == 1 && var_is_in(states, vars[1])
            push!(diff2vars, vars[1])
            push!(diffvars, get_variables(eq.lhs)[1])
            push!(remove_idxs, i)
        end
    end

    for i in reverse(remove_idxs)
        deleteat!(eqs, i)
    end

    for (i,var) in enumerate(diffvars)
        eqs = substitute(eqs, Dict(diff2vars[i]=>Differential(tvar)(diffvars[i])))
    end
    remove_idxs = Int[]
    for (i,var) in enumerate(states)
        if var_is_in(diff2vars, var)
            push!(remove_idxs, i)   
        end
    end
    for i in reverse(remove_idxs)
        deleteat!(states, i)
    end
    return eqs, states
end

function is_term(set, target_term)
    if typeof(set) == Equation
        vars = get_variables(set)
    elseif typeof(set) == SymbolicUtils.BasicSymbolic{Real}
        vars = get_variables(set)
    elseif typeof(set) == Num
        vars = get_variables(set)
    else
        vars = set
    end
    ret = false
    for term in vars
        if isequal(term, target_term)
            ret = true
            break
        else
            ret = false
        end
    end
    return ret
end

# function hbsweep(sys, jls, ns)
#     I₀ = 1e-6
#     R₀ = 50.0
#     Id = 0.05e-6
#     ωc = sqrt(2*pi *I₀/(jls.Φ₀*1000.0e-15))/(2*pi)
#     Ic = jls.Φ₀/(2*pi*1000.0e-12)
#     ω_vec = 2*pi*(4.5:0.001:5.0)*1e9
#     N = length(ω_vec)
#     solution1 = Vector{Float64}(undef, N)
#     solution2 = Vector{Float64}(undef, N)

#     ps = [
#         jls.P1.Isrc.ω => ω_vec[1]
#         jls.P1.Isrc.I => 0.00565e-6
#         jls.C1.C      => 100.0e-15
#         jls.J1.C      => 1000.0e-15
#         jls.J1.I0     => Ic
#         jls.P1.Rsrc.R => 50.0
#         jls.J1.R      => 1e9
#     ]
#     u0_vals = zeros(6)
#     u0_map = unknowns(sys) .=> u0_vals
#     prob = NonlinearProblem(sys, u0_map, ps)
#     for i in 1:N
#         # Update only frequency using remake
#         new_prob = remake(prob, p = [jls.P1.Isrc.ω => ω_vec[i]])
        
#         sol = solve(new_prob)
        
#         # Calculate magnitudes (Harmonic Ansatz: sqrt(A^2 + B^2))
#         solution1[i] = sqrt(sol[ns.C[1]]^2 + sol[ns.D[1]]^2) #check if push! is better
#         solution2[i] = sqrt(sol[ns.A[1]]^2 + sol[ns.B[1]]^2)
#     end

#     return ω_vec, solution1, solution2
# end 
function build_jacobians(rotated_system, vars, dvars)
    #TODO check ordering
    _jac = Symbolics.jacobian(rotated_system, vars)
    jac_0 = Num.(simplify(substitute(_jac, Dict(dvars .=> 0))))
    jac_1 = Symbolics.jacobian(rotated_system, dvars)
    return jac_0, jac_1
end
function rotate_to_harmonic_frame(N, Nt, harmonic_system)
    Γ = Matrix{Num}(undef, 2N + 1, Nt)
    # Currently only works for one harmonic ansatz. e.g. M=1
    for j in 1:Nt
        # 1. Place DC at the first index [1, j]
        Γ[1, j] = Num(1//Nt)
        # 2. Place Cosines at indices [2 to N+1]
        for n in 1:N
            phase = n * (j - 1) * (2π / Nt)
            Γ[n + 1, j] = Num((2//Nt) * cos(phase))
        end
        # 3. Place Sines at indices [N+2 to 2N+1]
        for n in 1:N
            phase = n * (j - 1) * (2π / Nt)
            Γ[N + 1 + n, j] = Num((2//Nt) * sin(phase))
        end
    end
    
    # Perform the rotation/projection
    rotated_system = Γ * [equation.lhs for equation in harmonic_system]
    
    return simplify.(rotated_system)
end





"""
    HarmonicProblem(sys, omega_pair::Pair; N::Int=1)

Constructs a harmonic balance problem from a time-domain dynamical system.

This function transforms a differential equation system (likely an ODESystem) into a system of nonlinear algebraic equations representing the harmonic coefficients. It automatically identifies or defines the independent variable (time) and expands the system variables into their harmonic series representations up to order N

# Arguments
- `sys`: The system model (typically a `ModelingToolkit.ODESystem`) containing the differential equations.
- `omega_pair::Pair`: A pair defining the fundamental frequency variable and its fixed value (e.g., `ω => 2.0`).

# Keywords
- `N::Int=1`: The number of harmonics to include in the expansion (truncation order). Higher values increase accuracy but increase computational cost.
compile::Bool=false`: Whether to compile and tear the resulting nonlinear system for performance.
# Returns
- `HarmonicProblem`: A struct containing the expanded nonlinear system (`complete_sys`), the harmonic order, and the frequency definitions.

# Details
If the generated system is over-determined (more equations than variables), the function automatically truncates the equation set to match the number of unknowns.
"""
function HarmonicProblem(sys, ωvar::Num; tearing::Bool=true, N::Int=1)
    # 1. Handle Time Variable
    tvar = ModelingToolkit.get_iv(sys) #put _ in tvar and wvar
    tvar = Num(tvar)

    eqs, states = get_full_equations(sys, tvar)

    nonlinear_sys, _ = harmonic_equation(eqs, states, tvar, ωvar, N) # wouldnt we need an omega value here??
    sys_eqs = equations(nonlinear_sys)
    sys_vars = unknowns(nonlinear_sys)
    
    if length(sys_eqs) > length(sys_vars)
        n_drop = length(sys_eqs) - length(sys_vars)
         @warn "System is overdetermined: $(length(sys_eqs)) equations for $(length(sys_vars)) variables. " 
             "Dropping the last equation(s). Caution: This behavior depends on variable order."
        sys_eqs = sys_eqs[1:end-n_drop]
    end
        
    @named nonlinear_sys = NonlinearSystem(sys_eqs, sys_vars, parameters(sys))
    if tearing
        complete_sys = mtkcompile(nonlinear_sys)
    else 
        complete_sys = nonlinear_sys
    end
  

    return HarmonicProblem(complete_sys,N,  unknowns(complete_sys), ωvar)
end

"""
    solve_sweep(prob::HarmonicProblem, base_params, sweep_pair) -> HarmonicSweepResult

Performs a parameter sweep on the harmonic problem using zero-order continuation.

This function structurally simplifies the harmonic system and solves it repeatedly across a range of parameter values. It uses the solution from the previous step as the initial guess for the current step to ensure convergence along the solution branch.

# Arguments
- `prob::HarmonicProblem`: The harmonic problem struct created by `HarmonicProblem`.
- `base_params`: A collection (Dict or Vector) of fixed parameter values required to fully define the system.
- `sweep_pair`: A pair where the first element is the symbolic parameter to vary and the second is an iterable of values (e.g., `k => 0.0:0.1:5.0`).

# Returns
- `HarmonicSweepResult`: A struct containing:
    - The swept variable name.
    - The vector of swept values.
    - A Dictionary mapping system variables (Num) to vectors of their computed values across the sweep.
    - The original problem definition.

# Details
The function automatically initializes unknown variables to `0.001` for the first solve. For subsequent steps, it uses `remake` on the `NonlinearProblem` to update parameters and initial guesses efficiently.
"""

function solve_sweep(prob::HarmonicProblem, base_params, sweep_pair)
    sweep_var = first(sweep_pair)
    sweep_vals = last(sweep_pair)
    sys = prob.sys

    # Setup Parameters
    current_params = Dict(base_params)
   
    
    current_params[sweep_var] = first(sweep_vals)
    
    # Use prob.sys_vars if available, or fetch from system
    system_unknowns = hasproperty(prob, :sys_vars) ? prob.sys_vars : unknowns(sys)
    
    # Initial guess
    u0_guess = [v => 0.001 for v in system_unknowns]
    
    # Define Problem ONCE
    nl_prob = NonlinearProblem(sys, u0_guess, current_params)
    
    results = Dict{Num, Vector{Float64}}()
    for v in system_unknowns
        results[v] = Float64[]
        sizehint!(results[v], length(sweep_vals))
    end
    
    println("Sweeping $(sweep_var) over $(length(sweep_vals)) points...")

    # Initialize last_u with the numeric initial guess
    last_u = nl_prob.u0
  
    for val in sweep_vals
        # Continuation: Update parameter and use previous solution (last_u) as guess
        nl_prob = remake(nl_prob; u0 = last_u, p = [sweep_var => val])
        
        sol = solve(nl_prob)
        last_u = sol.u
        

        for (i, v) in enumerate(system_unknowns)
            push!(results[v], sol.u[i])
        end
    end

    return HarmonicSweepResult(sweep_var, collect(sweep_vals), results)
end