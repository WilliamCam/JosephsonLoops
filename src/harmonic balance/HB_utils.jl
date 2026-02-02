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

This function transforms a differential equation system  into a system of nonlinear algebraic equations representing the harmonic 
coefficients. It automatically identifies or defines the independent variable (time) and expands the system variables into their harmonic series representations
up to order N

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

    nonlinear_sys, _ = harmonic_equation(eqs, states, tvar, ωvar, N) 
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
    solve_sweep(prob::HarmonicProblem, params, sweep_params) -> HarmonicSweepResult

Performs a parameter sweep on the harmonic problem using zero-order continuation.

This function structurally simplifies the harmonic system and solves it repeatedly across a range of parameter values. It uses the solution from the previous 
step as the initial guess for the current step to ensure convergence along the solution branch.

# Arguments
- `prob::HarmonicProblem`: The harmonic problem struct created by `HarmonicProblem`.
- `params`: A collection (Dict or Vector) of fixed parameter values required to fully define the system.
- `sweep_params`: A pair where the first element is the symbolic parameter to vary and the second is an iterable of values (e.g., `k => 0.0:0.1:5.0`).

# Returns
- `HarmonicSweepResult`: A struct containing:
    - The swept variable name.
    - The vector of swept values.
    - A Dictionary mapping system variables (Num) to vectors of their computed values across the sweep.
  

# Details
The function automatically initializes unknown variables to `0.001` for the first solve. For subsequent steps, it uses `remake` on the `NonlinearProblem` to update parameters and initial guesses efficiently.
"""

function solve_sweep(prob::HarmonicProblem, params, sweep_params)
    sweep_var = first(sweep_params)
    sweep_vals = last(sweep_params)
    sys = prob.sys

    # Setup Parameters
    current_params = Dict(params)
   
    
    current_params[sweep_var] = first(sweep_vals)
    
    
    system_unknowns = hasproperty(prob, :sys_vars) ? prob.sys_vars : unknowns(sys)
    
    # Initial guess
    u0_guess = [v => 0 for v in system_unknowns]
    u0_guess = [v => 0 for v in system_unknowns]
    
    # Define Problem 
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


"""
    flatten_harmonic_equation(expr, obs_dict::Dict, known_vars::Set{String})

Recursively substitutes observed variables into a symbolic expression until only known system variables or parameters remain.

This function resolves dependency chains (e.g., `F` depends on `D`, `D` depends on `A`) by iteratively looking up unknowns in `obs_dict`. It handles namespace matching (e.g., matching `sys₊D` to `D`) and defaults unresolved variables to `0.0` after a fixed number of iterations.

# Arguments
- `expr`: The symbolic expression to flatten.
- `obs_dict::Dict`: A dictionary mapping observed variable names (String) to their symbolic equations (RHS).
- `known_vars::Set{String}`: A set of variable names (strings) that are considered "fundamental" (e.g., system states or parameters) and should not be substituted.

# Returns
- A symbolic expression containing only variables found in `known_vars` (or constants).
"""

#  Recursively resolves observed variables (e.g. F -> D -> A)

function flatten_harmonic_equation(expr, obs_dict::Dict, known_vars::Set{String})
    current_expr = expr
    obs_keys = collect(keys(obs_dict))

    for _ in 1:20 
        vars_in_expr = Symbolics.get_variables(current_expr)
        unknowns = filter(v -> string(v) ∉ known_vars, vars_in_expr)
        
        if isempty(unknowns); return current_expr; end
        
        sub_rules = Dict()
        for u in unknowns
            u_str = string(u)
            # Find definition handling namespaces (e.g. "sys₊D" matches "D")
            idx = findfirst(k -> k == u_str || endswith(k, "₊"*u_str), obs_keys)
            
            if idx !== nothing
                key = obs_keys[idx]
                sub_rules[u] = obs_dict[key]
            else
                sub_rules[u] = 0.0 
            end
        end
        current_expr = substitute(current_expr, sub_rules)
    end
    return current_expr
end

"""
    evaluate_harmonic_sweep(func, h_prob, sweep_res, params)

Evaluates a compiled harmonic function across all points in a simulation sweep.

This helper function iterates through the sweep results, updating the sweep parameter (if present) and the system state vector `u`, then invokes the compiled function `func` for each point.

# Arguments
- `func`: A compiled Julia function (usually generated by `build_function`) that takes `(u, p)` as arguments.
- `h_prob`: The `HarmonicProblem` struct containing system definitions and variable orderings.
- `sweep_res`: The result object from the sweep, containing `sweep_vals` and `results`.
- `params`: A dictionary or indexable object containing base parameter values.

# Returns
- `Vector{Float64}`: An array of evaluated results corresponding to each point in `sweep_res.sweep_vals`.
"""

function evaluate_harmonic_sweep(func, h_prob, sweep_res, params)
    n_points = length(sweep_res.sweep_vals)
    out_vals = zeros(Float64, n_points)
    
    param_syms = parameters(h_prob.sys)
    p_vals_base = [params[p] for p in param_syms]
    sweep_idx = findfirst(isequal(sweep_res.sweep_var), param_syms)

    for i in 1:n_points
        u_val = [sweep_res.results[v][i] for v in h_prob.sys_vars]
        if sweep_idx !== nothing
            p_vals_base[sweep_idx] = sweep_res.sweep_vals[i]
        end
        
      
        out_vals[i] = Base.invokelatest(func, u_val, p_vals_base)
    end
    return out_vals
end

"""
    fetch_harmonic_vector(target_key::String, h_prob::HarmonicProblem, sweep_res, params, obs_dict, known_vars)

Retrieves the data vector for a specific harmonic variable, checking both direct results and observed equations.

The function attempts to resolve the `target_key` in the following order:
1. Checks if `target_key` exists directly in `sweep_res.results`.
2. If not found, checks `obs_dict` (observed variables). If found, it flattens the observed equation, compiles it into a function, and evaluates it across the sweep.
3. If the key is not found in either, returns a vector of zeros.

# Arguments
- `target_key::String`: The name or substring of the variable to retrieve (e.g., "u", "sys₊u").
- `h_prob`: The `HarmonicProblem` context.
- `sweep_res`: The simulation sweep results.
- `params`: Base parameter values.
- `obs_dict`: Dictionary of observed equations.
- `known_vars`: Set of fundamental system variables.

# Returns
- `Vector{Float64}`: The numerical values of the target variable across the sweep.
"""
#currently not returning all variables from observed equations i think?

function fetch_harmonic_vector(target_key::String, h_prob::HarmonicProblem, sweep_res, params, obs_dict, known_vars)
    
    # A. Check Results (Fast Path)
    res_keys = collect(keys(sweep_res.results))
    res_idx = findfirst(k -> occursin(target_key, string(k)), res_keys)
    if res_idx !== nothing
        return sweep_res.results[res_keys[res_idx]]
    end

    # B. Check Observed Equations (Slow Path)
    # FIX: Collect keys into a Vector so 'findfirst' works
    obs_keys = collect(keys(obs_dict))
    match_idx = findfirst(k -> occursin(target_key, k), obs_keys)
    
    if match_idx !== nothing
        key = obs_keys[match_idx]
        
        # 1. Get RHS and Flatten
        raw_rhs = obs_dict[key]
        flat_rhs = flatten_harmonic_equation(raw_rhs, obs_dict, known_vars)
        
        # 2. Compile Function (Val{true} returns an expression we can eval)
        func_ex = build_function(flat_rhs, h_prob.sys_vars, parameters(h_prob.sys), expression=Val{true})
        func = eval(func_ex)
        
        # 3. Evaluate safely
        return evaluate_harmonic_sweep(func, h_prob, sweep_res, params)
    end

    # C. Not Found (Assume 0.0)
    return zeros(Float64, length(sweep_res.sweep_vals))
end


"""
    get_harmonic_coeffs(h_prob::HarmonicProblem, model, sweep_res, params, var_name::String)

Extracts and computes the complex harmonic coefficients for a specific physical variable.

This function maps a physical variable name (e.g., "position") to its corresponding harmonic components (Cos/Sin terms) based on the variable's index in the system. It retrieves these components from the results or observed equations and combines them into a complex phasor.

# Arguments
- `h_prob::HarmonicProblem`: The problem definition containing the system structure.
- `model`: The symbolic model (usually ODESystem) defining unknowns.
- `sweep_res`: The output data from the parameter sweep.
- `params`: Parameter values for evaluation.
- `var_name::String`: The string representation of the physical variable to query (e.g., "x").

# Returns
- `Vector{ComplexF64}`: A vector of complex numbers representing the harmonic coefficient for each point in the sweep.
"""
function get_harmonic_coeffs(h_prob::HarmonicProblem , model, sweep_res, params, var_name::String)

    # 1. Map Physical Name -> Harmonic Keys
    sys_states = unknowns(model)
    k = findfirst(s -> occursin(var_name, string(s)), sys_states)
    if k === nothing; error("Variable '$var_name' not found."); end

    alphabet = 'A':'Z'
    cos_key = "$(alphabet[2*k - 1])[2]"
    sin_key = "$(alphabet[2*k])[1]"
    
    println("Variable '$var_name' (Index $k) maps to Harmonic Vars: $cos_key, $sin_key")

    # 2. Prepare Lookup Data
    obs_dict = Dict(string(eq.lhs) => eq.rhs for eq in observed(h_prob.sys))
    known_vars = Set(string.([h_prob.sys_vars; parameters(h_prob.sys)]))
    
    # 3. Fetch Components
    A_vec = fetch_harmonic_vector(cos_key, h_prob, sweep_res, params, obs_dict, known_vars)
    B_vec = fetch_harmonic_vector(sin_key, h_prob, sweep_res, params, obs_dict, known_vars)

    return @. A_vec - im*B_vec
end


#double check formulas with will 

"""
    get_harmonic_magnitude(h_prob, model, sweep_res, params, var_name)

Wrapper that gets the complex phasor and returns its Peak Magnitude.
Magnitude = sqrt(A^2 + B^2)
"""
function get_harmonic_magnitude(h_prob, model, sweep_res, params, var_name::String)
    # 1. Get the complex phasor (A - iB)
    phasor = get_harmonic_coeffs(h_prob, model, sweep_res, params, var_name)
    
    # 2. Return Magnitude (Element-wise)
    return abs.(phasor)
end

"""
    get_harmonic_phase(h_prob, model, sweep_res, params, var_name)

Wrapper that gets the complex phasor and returns its Phase in Radians.
Phase = atan(-B / A)
"""
function get_harmonic_phase(h_prob, model, sweep_res, params, var_name::String)
    phasor = get_harmonic_coeffs(h_prob, model, sweep_res, params, var_name)
    return angle.(phasor)
end