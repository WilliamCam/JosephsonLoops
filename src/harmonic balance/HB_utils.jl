using Symbolics
using SymbolicUtils
using NonlinearSolve

struct HarmonicProblem
    sys  # The algebraic NonlinearSystem
    N::Int
    sys_vars::Vector{Num}
    ωvar::Num
    params::Dict
    variable_map::Dict{Tuple{String, Int, Symbol}, Num}
    observed::Vector{Equation}
end

struct HarmonicSweepProblem
    prob::HarmonicProblem
    sweep_var::Num
    sweep_vals::AbstractVector
end

struct HarmonicSweepResult
    sweep_var::Num
    sweep_vals::AbstractVector
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
function HarmonicProblem(sys, ωvar::Num, params; tearing::Bool=true, N::Int=1) #jac::Bool=false)
    # 1. Handle Time Variable
    tvar = ModelingToolkit.get_iv(sys) #put _ in tvar and wvars
    tvar = Num(tvar)

    eqs, states = get_full_equations(sys, tvar)

    nonlinear_sys, _, variable_map = harmonic_equation(eqs, states, tvar, ωvar, N) 
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
  

    return HarmonicProblem(complete_sys, N, unknowns(complete_sys), ωvar, params, variable_map, observed(sys))
end



function solve(prob::HarmonicProblem)
    u0_guess = fill(0.0, length(prob.sys_vars))
    current_params = prob.params
    #to remove deprecated error from NonlinearProblem, we have to pass merge(dict, u0)
    combined_args = merge(Dict(prob.sys_vars .=> u0_guess), current_params)
    nl_prob = NonlinearProblem(prob.sys, combined_args)
    return ModelingToolkit.solve(nl_prob)
end

function solve(sweepprob::HarmonicSweepProblem)
    prob = sweepprob.prob
    sweep_var = sweepprob.sweep_var
    sweep_vals = sweepprob.sweep_vals
    sys = prob.sys

    # Setup Parameters
    current_params = copy(prob.params)
   
    current_params[sweep_var] = first(sweep_vals)
    
    
    system_unknowns = prob.sys_vars
    
    # Initial guess
    u0_guess = fill(0.0, length(prob.sys_vars))

    # Define Problem - got rid of deprecated error
    combined_args = merge(Dict(prob.sys_vars .=> u0_guess), current_params)
    nl_prob = NonlinearProblem(sys, combined_args)
    
    results = Dict{Num, Vector{Float64}}()
    for v in system_unknowns
        results[v] = Float64[]
        sizehint!(results[v], length(sweep_vals))
    end
    println("Sweeping $(sweep_var) over $(length(sweep_vals)) points...")
    # Initialize last_u with the numeric initial guess
    last_u = nl_prob.u0
  
    for  val in sweep_vals     
        # Continuation: use previous solution (last_u)
        nl_prob = remake(nl_prob; u0 = last_u, p = [sweep_var => val])
        sol = ModelingToolkit.solve(nl_prob)
        last_u = sol.u
        for (i, v) in enumerate(prob.sys_vars)
            push!(results[v], sol.u[i])
        end
    end
    return HarmonicSweepResult(sweep_var, collect(sweep_vals), results)
end


"""
    get_harmonic(h_prob::HarmonicProblem, sweep_res, var_name::String, order::Int; output_type::Symbol=:complex)

Extracts harmonic components. Checks variable_map first, then falls back to observed equations.
output_type: :complex (default), :magnitude, :phase
"""
function get_harmonic(h_prob::HarmonicProblem, sweep_res, var_name::String, order::Int; output_type::Symbol=:complex)
    
    # Build observed dict from harmonic system (for vars like D, E, F that are in terms of A, B, C)
    obs_dict = Dict(string(eq.lhs) => eq.rhs for eq in observed(h_prob.sys))
    params = h_prob.params

    # Try to get from variable_map first
    if order == 0
        key = (var_name, 0, :Cos)
        if haskey(h_prob.variable_map, key)
            sym = h_prob.variable_map[key]
            val_vec = fetch_harmonic_vector(string(sym), h_prob, sweep_res, params, obs_dict)
            phasor = complex.(val_vec)
        else
            # Not in map - must be an observed variable from original system
            # Try to find it in the time-domain observed equations (h_prob.observed)
            phasor = complex.(reconstruct_from_observed(h_prob, sweep_res, var_name, 0, :Cos))
        end
    else
        key_cos = (var_name, order, :Cos)
        key_sin = (var_name, order, :Sin)
        
        if haskey(h_prob.variable_map, key_cos) && haskey(h_prob.variable_map, key_sin)
            sym_cos = h_prob.variable_map[key_cos]
            sym_sin = h_prob.variable_map[key_sin]
            
            A_vec = fetch_harmonic_vector(string(sym_cos), h_prob, sweep_res, params, obs_dict)
            B_vec = fetch_harmonic_vector(string(sym_sin), h_prob, sweep_res, params, obs_dict)
            
            phasor = @. (A_vec - im*B_vec)
        else
            # Not in map - reconstruct from time-domain observed
            A_vec = reconstruct_from_observed(h_prob, sweep_res, var_name, order, :Cos)
            B_vec = reconstruct_from_observed(h_prob, sweep_res, var_name, order, :Sin)
            phasor = @. (A_vec - im*B_vec)
        end
    end

    if output_type == :magnitude
        return abs.(phasor)
    elseif output_type == :phase
        return angle.(phasor)
    else
        return phasor
    end
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


function fetch_harmonic_vector(target_key::String, h_prob::HarmonicProblem, sweep_res, params, obs_dict)
    
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
        
        # Get RHS (already in harmonic form)
        raw_rhs = obs_dict[key]
        
        # Check if expression contains unknown harmonic coefficients (like G[1])
        vars_in_rhs = Symbolics.get_variables(raw_rhs)
        known_syms = Set(h_prob.sys_vars)
        param_syms = Set(parameters(h_prob.sys))
        
        unknown_vars = filter(v -> v ∉ known_syms && v ∉ param_syms, vars_in_rhs)
        
        if !isempty(unknown_vars)
            # Recursively substitute unknown harmonic coefficients
            for _ in 1:10
                vars_in_expr = Symbolics.get_variables(raw_rhs)
                unknowns = filter(v -> v ∉ known_syms && v ∉ param_syms, vars_in_expr)
                
                if isempty(unknowns)
                    break
                end
                
                subs_dict = Dict()
                for u in unknowns
                    u_str = string(u)
                    if haskey(obs_dict, u_str)
                        subs_dict[u] = obs_dict[u_str]
                    end
                end
                
                if isempty(subs_dict)
                    break
                end
                
                raw_rhs = substitute(raw_rhs, subs_dict)
                raw_rhs = Symbolics.simplify(raw_rhs)
            end
        end
        
        # Compile Function
        func_ex = build_function(raw_rhs, h_prob.sys_vars, parameters(h_prob.sys), expression=Val{true})
        func = eval(func_ex)
        
        # Evaluate
        return evaluate_harmonic_sweep(func, h_prob, sweep_res, params)
    end

    # C. Not Found (Assume 0.0)
    return zeros(Float64, length(sweep_res.sweep_vals))
end


"""
    reconstruct_from_observed(h_prob, sweep_res, var_name, order, component)

Reconstructs a harmonic component from time-domain observed equations.
This is used for variables like C1₊i that are observed in the original system
but not included in the harmonic variable_map after tearing.
"""

function reconstruct_from_observed(h_prob::HarmonicProblem, sweep_res, var_name::String, order::Int, component::Symbol)
    # Find the observed equation for this variable in the original (time-domain) system
    target_eq = nothing
    for eq in h_prob.observed
        lhs_str = string(eq.lhs)
        # Match "C1₊i(t)" to "C1₊i"
        if replace(lhs_str, "(t)" => "") == var_name
            target_eq = eq
            break
        end
    end
    
    if target_eq === nothing
        error("Variable $var_name not found in variable_map or observed equations")
    end
    
    # Extract time variable
    rhs_expr = target_eq.rhs
    lhs_term = first(Symbolics.get_variables(target_eq.lhs))
    tvar = Symbolics.arguments(Symbolics.unwrap(lhs_term))[1]
    
    # First, recursively substitute any observed variable dependencies (e.g., J1₊i in C1₊i equation)
    obs_dict = Dict(replace(string(eq.lhs), "(t)" => "") => eq.rhs for eq in h_prob.observed)
    
    # Recursively expand observed dependencies including their derivatives
    for _ in 1:10
        vars_in_rhs = Symbolics.get_variables(rhs_expr)
        has_obs = false
        obs_subs = Dict()
        
        for v in vars_in_rhs
            v_name = replace(string(v), "(t)" => "")
            if haskey(obs_dict, v_name)
                # Only substitute the base variable - let state substitution handle derivatives
                obs_subs[v] = obs_dict[v_name]
                has_obs = true
            end
        end
        
        if !has_obs
            break
        end
        
        # Just substitute, don't expand derivatives yet
        rhs_expr = substitute(rhs_expr, obs_subs)
    end
    
    N = h_prob.N
    ωvar = h_prob.ωvar
    
    # Build substitution map: state_var(t) => harmonic expansion
    state_subs = Dict()
    
    # Group harmonic coefficients by state variable name
    grouped_vars = Dict{String, Dict}()
    for (key, sym) in h_prob.variable_map
        (name, n, type) = key
        if !haskey(grouped_vars, name)
            grouped_vars[name] = Dict()
        end
        grouped_vars[name][(n, type)] = sym
    end
    
    # For each state variable,  build its harmonic expansion
    for (name, coeffs) in grouped_vars
        A0 = get(coeffs, (0, :Cos), 0.0)
        
        expand_expr = A0
        dX_expr = 0
        d2X_expr = 0
        
        for n in 1:N
            An = get(coeffs, (n, :Cos), 0.0)
            Bn = get(coeffs, (n, :Sin), 0.0)
            
            phase = n * ωvar * tvar
            expand_expr += An * cos(phase) + Bn * sin(phase)
            dX_expr += -An * n * ωvar * sin(phase) + Bn * n * ωvar * cos(phase)
            d2X_expr += -An * (n*ωvar)^2 * cos(phase) - Bn * (n*ωvar)^2 * sin(phase)
        end
        
        # Find matching variable in rhs_expr and substitute
        vars_in_rhs = Symbolics.get_variables(rhs_expr)
        for v in vars_in_rhs
            v_str = replace(string(v), "(t)" => "")
            
            # Match base variable (e.g., "J1₊θ" or "sys₊J1₊θ")
            if v_str == name || endswith(v_str, "₊"*name)
                state_subs[v] = expand_expr
                state_subs[Differential(tvar)(v)] = dX_expr
                state_subs[Differential(tvar)(Differential(tvar)(v))] = d2X_expr
            # Match derivative notation (e.g., "J1₊θˍt" represents D(J1₊θ))
            elseif v_str == name * "ˍt" || endswith(v_str, "₊" * name * "ˍt")
                state_subs[v] = dX_expr
            end
        end
    end
    
    # Substitute state expansions into observed equation
    expanded_rhs = substitute(rhs_expr, state_subs)
    
    # Expand any remaining derivatives (from the state substitutions)
    # Use multiple passes to ensure all nested derivatives are expanded
    for _ in 1:3
        expanded_rhs = Symbolics.expand_derivatives(expanded_rhs)
        expanded_rhs = Symbolics.simplify(expanded_rhs)
    end
    
    # Extract harmonic component via DFT
    Nt = 2*N+1
    dft_val = 0
    for k in 0:(Nt-1)
        phi_k = 2*pi*k/Nt
        val_k = substitute(expanded_rhs, Dict(tvar => phi_k / ωvar))
        
        if component == :Cos && order == 0
            dft_val += val_k
        elseif component == :Cos
            dft_val += val_k * cos(order * phi_k)
        elseif component == :Sin
            dft_val += val_k * sin(order * phi_k)
        end
    end
    
    factor = (order == 0) ? (1/Nt) : (2/Nt)
    final_expr = dft_val * factor
    
    # Check if expression contains harmonic coefficients from observed equations (like D[1])
    vars_in_final = Symbolics.get_variables(final_expr)
    known_syms = Set(h_prob.sys_vars)
    param_syms = Set(parameters(h_prob.sys))
    
    # Find unknown variables (harmonic coefficients from observed equations)
    unknown_vars = filter(v -> v ∉ known_syms && v ∉ param_syms, vars_in_final)
    
    if !isempty(unknown_vars)
        
        # Build observed dict from harmonic system
        obs_dict = Dict(string(eq.lhs) => eq.rhs for eq in observed(h_prob.sys))
        
        # Recursively substitute unknown harmonic coefficients
        for _ in 1:10
            vars_in_expr = Symbolics.get_variables(final_expr)
            unknowns = filter(v -> v ∉ known_syms && v ∉ param_syms, vars_in_expr)
            
            if isempty(unknowns)
                break
            end
            
            subs_dict = Dict()
            for u in unknowns
                u_str = string(u)
                if haskey(obs_dict, u_str)
                    subs_dict[u] = obs_dict[u_str]
                end
            end
            
            if isempty(subs_dict)
                break
            end
            
            final_expr = substitute(final_expr, subs_dict)
            final_expr = Symbolics.simplify(final_expr)
        end
    end
    

    
    # Compile and evaluate
    func_ex = build_function(final_expr, h_prob.sys_vars, parameters(h_prob.sys), expression=Val{true})
    func = eval(func_ex)
    
    return evaluate_harmonic_sweep(func, h_prob, sweep_res, h_prob.params)
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


function fetch_harmonic_vector(target_key::String, h_prob::HarmonicProblem, sweep_res, params, obs_dict)
    
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
        
        # Get RHS (already in harmonic form)
        raw_rhs = obs_dict[key]
        
        # Check if expression contains unknown harmonic coefficients (like G[1])
        vars_in_rhs = Symbolics.get_variables(raw_rhs)
        known_syms = Set(h_prob.sys_vars)
        param_syms = Set(parameters(h_prob.sys))
        
        unknown_vars = filter(v -> v ∉ known_syms && v ∉ param_syms, vars_in_rhs)
        
        if !isempty(unknown_vars)
            # Recursively substitute unknown harmonic coefficients
            for _ in 1:10
                vars_in_expr = Symbolics.get_variables(raw_rhs)
                unknowns = filter(v -> v ∉ known_syms && v ∉ param_syms, vars_in_expr)
                
                if isempty(unknowns)
                    break
                end
                
                subs_dict = Dict()
                for u in unknowns
                    u_str = string(u)
                    if haskey(obs_dict, u_str)
                        subs_dict[u] = obs_dict[u_str]
                    end
                end
                
                if isempty(subs_dict)
                    break
                end
                
                raw_rhs = substitute(raw_rhs, subs_dict)
                raw_rhs = Symbolics.simplify(raw_rhs)
            end
        end
        
        # Compile Function
        func_ex = build_function(raw_rhs, h_prob.sys_vars, parameters(h_prob.sys), expression=Val{true})
        func = eval(func_ex)
        
        # Evaluate
        return evaluate_harmonic_sweep(func, h_prob, sweep_res, params)
    end

    # C. Not Found (Assume 0.0)
    return zeros(Float64, length(sweep_res.sweep_vals))
end

