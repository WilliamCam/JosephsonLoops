using Symbolics
using SymbolicUtils
using NonlinearSolve

struct HarmonicSystem
    complete_sys::ModelingToolkit.AbstractSystem
    ω_var::Num
    t_var::Num          # independent variable from the original ODE system
    N::Int
    variable_map::Dict{Tuple{String, Int, Symbol}, Num}
    observed_equations::Dict{Any, Any}
end

struct HarmonicProblem
    harmonic_system::HarmonicSystem
    params::Dict
    sweep_var::Union{Num, Nothing}
    sweep_vals::Union{AbstractVector, Nothing}
    u0::Vector{Float64}
end

struct HarmonicResult
    sweep_var::Num
    sweep_vals::AbstractVector
    results::Dict{Num, Vector{Float64}}
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

function rotate_to_harmonic_frame(M, N, Nt, harmonic_system)
    # M: Number of state variables (dimension of state vector)
    block_rows = 2N + 1
    # Total Matrix Dimensions
    total_rows = M * block_rows
    total_cols = M * Nt
    # Initialize the large block-diagonal matrix
    Γ_total = zeros(Num, total_rows, total_cols)
    # Generate the single-variable operator (the logic you already have)
    Γ_single = Matrix{Num}(undef, block_rows, Nt)
    for j in 1:Nt
        # DC Term
        Γ_single[1, j] = Num(1//Nt)
        # AC Terms
        for n in 1:N
            phase = n * (j - 1) * (2π / Nt)
            Γ_single[n + 1, j] = Num((2//Nt) * cos(phase))      # Cosine
            Γ_single[N + 1 + n, j] = Num((2//Nt) * sin(phase))  # Sine
        end
    end
    # Place the single-variable operator into the block diagonal
    for d in 1:M
        row_range = (d-1)*block_rows + 1 : d*block_rows
        col_range = (d-1)*Nt + 1 : d*Nt
        Γ_total[row_range, col_range] .= Γ_single
    end
    rotated_system = Γ_total * [equation.lhs for equation in harmonic_system]    
    return simplify.(rotated_system)
end

"""
    HarmonicSystem(sys, ωvar; tearing=true, N=1) -> HarmonicSystem

Construct a `HarmonicSystem` by applying the harmonic balance method to a
time-domain ODE system.

Each state variable `x(t)` is expanded into a truncated Fourier series
    x(t) = A₀ + Σₙ₌₁ᴺ [Aₙ cos(nωt) + Bₙ sin(nωt)]
and the ODE residual is projected onto each harmonic basis function, converting
the differential equations into a set of nonlinear algebraic equations in the
Fourier coefficients {A₀, Aₙ, Bₙ}.

# Arguments
- `sys`: A `ModelingToolkit.ODESystem` (or compatible system) containing the
  time-domain differential equations.
- `ωvar::Num`: Symbolic variable representing the fundamental angular frequency ω.
  Must already be declared as a parameter in `sys`.

# Keywords
- `N::Int=1`: Harmonic truncation order. `N=1` retains only the fundamental; larger
  values improve accuracy at the cost of a larger nonlinear system.
- `tearing::Bool=true`: Whether to apply structural simplification (`mtkcompile`) to
  the resulting `NonlinearSystem`. Tearing reduces the number of active unknowns and
  typically speeds up the solve.

# Returns
A `HarmonicSystem` containing the compiled nonlinear system, frequency metadata,
the variable map, and the original observed equations.

# Details
- The independent variable (time) is inferred automatically via
  `ModelingToolkit.get_iv(sys)`.
- If the expanded system is over-determined (more equations than unknowns — which
  can occur when symmetry constraints force certain coefficients to zero), trailing
  equations are dropped with a warning. The ordering is determined by
  `harmonic_equation`, so inspect results carefully in this case.
- Observed equations from `sys` are stored verbatim in `HarmonicSystem.observed_equations`
  so that non-state quantities (e.g. branch currents, voltages across elements) can be
  reconstructed from the HB solution without re-solving the system.
"""
function HarmonicSystem(sys, ωvar::Num; tearing::Bool=true, N::Int=1)
    tvar = Num(ModelingToolkit.get_iv(sys))
    eqs, states = get_full_equations(sys, tvar)

    nonlinear_sys, _, variable_map = harmonic_equation(eqs, states, tvar, ωvar, N)
    
    sys_eqs = equations(nonlinear_sys)
    sys_vars = unknowns(nonlinear_sys)
    
    if length(sys_eqs) > length(sys_vars)
        n_drop = length(sys_eqs) - length(sys_vars)
        @warn "System is overdetermined: $(length(sys_eqs)) equations for $(length(sys_vars)) variables. " *
              "Dropping the last equation(s). Caution: This behavior depends on variable order."
        sys_eqs = sys_eqs[1:end-n_drop]
    end
        
    @named nonlinear_sys = NonlinearSystem(sys_eqs, sys_vars, parameters(sys))
    complete_sys = tearing ? mtkcompile(nonlinear_sys) : nonlinear_sys
    
    # Store original observed equations to enable reconstruction of non-state variables (eg. C1.i)
    observed_eqs = Dict{Any, Any}(eq.lhs => eq.rhs for eq in observed(sys))
    
    return HarmonicSystem(complete_sys, ωvar, tvar, N, variable_map, observed_eqs)
end

function HarmonicProblem(harmonic_sys::HarmonicSystem, params; sweep_var::Union{Num, Nothing}=nothing, sweep_vals::Union{AbstractVector, Nothing}=nothing)
    u0 = fill(0.0, length(unknowns(harmonic_sys.complete_sys)))
    return HarmonicProblem(harmonic_sys, params, sweep_var, sweep_vals, u0)
end
"""
    solve(prob::HarmonicProblem; kwargs...) -> NonlinearSolution | HarmonicResult

Solve the harmonic balance system. If `sweep_var` and `sweep_vals` are set on the
problem, performs a parameter sweep using continuation (warm-starting each point from
the previous solution). Otherwise, solves at a single parameter point.

All keyword arguments are forwarded to `ModelingToolkit.solve` (e.g. `alg`, `abstol`).

# Returns
- **No sweep**: A `NonlinearSolution` whose harmonic coefficients can be read via `sol[var]`.
- **With sweep**: A `HarmonicResult` containing the swept parameter, sweep values, and
  result arrays for each HB unknown.
"""
function solve(prob::HarmonicProblem; kwargs...)
    hsys = prob.harmonic_system
    sys = hsys.complete_sys
    system_unknowns = unknowns(sys)

    # no sweep
    if prob.sweep_var === nothing || prob.sweep_vals === nothing
        combined_args = merge(Dict(system_unknowns .=> prob.u0), prob.params)
        nl_prob = NonlinearProblem(sys, combined_args)
        return ModelingToolkit.solve(nl_prob; kwargs...)
    end
    #Sweep solve 
    sweep_var = prob.sweep_var
    sweep_vals = prob.sweep_vals

    # Remove sweep variable from fixed params to avoid conflicts
    current_params = copy(prob.params)
    if haskey(current_params, sweep_var)
        delete!(current_params, sweep_var)
    end
    sweep_var_uw = Symbolics.unwrap(sweep_var)
    if haskey(current_params, sweep_var_uw)
        delete!(current_params, sweep_var_uw)
    end

    current_params[sweep_var] = first(sweep_vals)
    combined_args = merge(Dict(system_unknowns .=> prob.u0), current_params)
    nl_prob = NonlinearProblem(sys, combined_args)

    # Results containers
    results = Dict{Num, Vector{Float64}}()
    for v in system_unknowns
        results[v] = Float64[]
        sizehint!(results[v], length(sweep_vals))
    end

    println("Sweeping $(sweep_var) over $(length(sweep_vals)) points...")
    last_u = nl_prob.u0

    for val in sweep_vals
        nl_prob = remake(nl_prob; u0 = last_u, p = [sweep_var => val])
        sol = ModelingToolkit.solve(nl_prob, kwargs...)
        last_u = sol.u
        for (i, v) in enumerate(system_unknowns)
            push!(results[v], sol.u[i])
        end
    end
    return HarmonicResult(sweep_var, collect(sweep_vals), results)
end