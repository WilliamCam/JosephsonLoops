# Extracts numeric phasor A - im B from Harmonic balance result for some order in ω
using ModelingToolkit
#include("./utils.jl")
#  Shared Helpers

function get_phasor(h_prob::Union{HarmonicProblem, LinearizedProblem}, result::HarmonicResult, var_name::String, order::Int = 1)
    #TODO: parameterised or updated call to harmonic_expression to avoid recomputing output expression
    expr = get_harmonic_expression(h_prob, var_name, order)
    return apply_harmonic_expression(h_prob, expr[1], expr[2])(result)
end

"""Strip `(t)` suffix and whitespace from a symbolic variable name."""
clean_name(x) = replace(string(x), "(t)" => "", " " => "")

"""Find a matching key in a variable_map by exact or fuzzy name match."""
function find_varmap_key(variable_map, vname, order, component)
    haskey(variable_map, (vname, order, component)) && return (vname, order, component)
    for k in keys(variable_map)
        k[2] == order && k[3] == component && occursin(vname, k[1]) && return k
    end
    return nothing
end

"""Unified symbolic lookup across one or more Dicts. Tries hash, isequal, then string match."""
function _find_sym(sym::Num, dicts...)
    sym_uw = Symbolics.unwrap(sym)
    sym_str = clean_name(sym)
    for d in dicts
        haskey(d, sym) && return d[sym]
        for (k, v) in d
            (isequal(k, sym) || isequal(Symbolics.unwrap(k), sym_uw)) && return v
        end
        for (k, v) in d
            clean_name(k) == sym_str && return v
        end
    end
    return nothing
end

function _get_coeff_expr(sym::Num, h_prob::HarmonicProblem)
    sys = h_prob.harmonic_system.complete_sys
    # If sym is a system unknown, return it directly — it will be substituted numerically
    for u in unknowns(sys)
        isequal(u, sym) && return sym
    end
    # isequal can fail for array-indexed vars (different underlying array objects after
    # mtkcompile).  Retry with string comparison, returning the COMPILED system's symbol
    # so that build_function can map it to the correct input position.
    sym_str = clean_name(sym)
    for u in unknowns(sys)
        clean_name(u) == sym_str && return Num(Symbolics.unwrap(u))
    end
    # sym was eliminated by tearing — find it in observed equations
    for eq in observed(sys)
        isequal(eq.lhs, sym) && return Num(Symbolics.unwrap(eq.rhs))
    end
    for eq in observed(sys)
        clean_name(eq.lhs) == sym_str && return Num(Symbolics.unwrap(eq.rhs))
    end
    # Not found anywhere — return the symbol itself (may be zero-forced by structure)
    @warn "_get_coeff_expr: could not resolve $(clean_name(sym)) — returning as-is. unknowns=$(clean_name.(unknowns(sys))), observed_lhs=$(clean_name.(eq.lhs for eq in observed(sys)))"
    return sym
end

function get_harmonic_expression(h_prob::HarmonicProblem, var_name::String, order::Int)
    var_clean = clean_name(var_name)
    vmap = h_prob.harmonic_system.variable_map
    #Build symbolic phasor expression using raw vmap symbols 
    if order == 0
        key = find_varmap_key(vmap, var_clean, 0, :Cos)
        phasor_re = key !== nothing ? vmap[key] : reconstruct_from_observed(h_prob, var_clean, 0, :Cos)
        phasor_im = Num(0.0)
    else
        key_cos = find_varmap_key(vmap, var_clean, order, :Cos)
        key_sin = find_varmap_key(vmap, var_clean, order, :Sin)
        if key_cos !== nothing && key_sin !== nothing
            A_expr = vmap[key_cos]
            B_expr = vmap[key_sin]
        else
            A_expr = reconstruct_from_observed(h_prob, var_clean, order, :Cos)
            B_expr = reconstruct_from_observed(h_prob, var_clean, order, :Sin)
        end
        phasor_re = A_expr
        phasor_im = B_expr
    end
    return phasor_re, phasor_im
end

function apply_harmonic_expression(h_prob::HarmonicProblem, expression_re::Num, expression_im::Num)
    vmap = h_prob.harmonic_system.variable_map
    sys = h_prob.harmonic_system.complete_sys
    all_vmap_syms = collect(values(vmap))
    _params     = collect(parameters(sys))
    all_vars_uw   = vcat(Symbolics.unwrap.(all_vmap_syms), Symbolics.unwrap.(_params))

    #Compile build_function 
    f_re = Symbolics.build_function(Symbolics.unwrap(expression_re), all_vars_uw; expression=Val{false})
    f_im = Symbolics.build_function(Symbolics.unwrap(expression_im), all_vars_uw; expression=Val{false})

    # Return closure for fast numeric evaluation ─
    n_vmap   = length(all_vmap_syms)
    n_params = length(_params)

    # Pre-compute fixed parameter values
    param_vals = zeros(Float64, n_params)
    for (j, p) in enumerate(_params)
        val = _find_sym(Num(p), h_prob.params)
        val !== nothing && (param_vals[j] = Float64(val))
    end

    return function(sweep_res::HarmonicResult)
        n_points = length(sweep_res.sweep_vals)
        input_vec = zeros(Float64, n_vmap + n_params)
        input_vec[n_vmap+1:end] .= param_vals
        vmap_vecs = Vector{Union{Nothing, Vector{Float64}}}(undef, n_vmap)
        for (j, sym) in enumerate(all_vmap_syms)
            vmap_vecs[j] = _find_sym(sym, sweep_res.results)
        end
        # Locate sweep parameter index in param section
        sweep_j = nothing
        sweep_str = clean_name(sweep_res.sweep_var)
        for (j, p) in enumerate(_params)
            if clean_name(p) == sweep_str
                sweep_j = j; break
            end
        end
        phasor = Vector{ComplexF64}(undef, n_points)
        for i in 1:n_points
            for j in 1:n_vmap
                input_vec[j] = vmap_vecs[j] !== nothing ? vmap_vecs[j][i] : 0.0
            end
            sweep_j !== nothing && (input_vec[n_vmap + sweep_j] = sweep_res.sweep_vals[i])
            phasor[i] = complex(f_re(input_vec), f_im(input_vec))
        end
        return phasor
    end
end
function reconstruct_from_observed(h_prob::HarmonicProblem, var_name::String, order::Int, component::Symbol)
    sys = h_prob.harmonic_system.complete_sys
    t_sym = h_prob.harmonic_system.t_var
    tvar_uw = Symbolics.unwrap(t_sym)

    obs_dict = Dict{Any, Any}()
    for (k, v) in h_prob.harmonic_system.observed_equations
        obs_dict[Symbolics.unwrap(k)] = Symbolics.unwrap(v)
    end
    for eq in observed(sys)
        obs_dict[Symbolics.unwrap(eq.lhs)] = Symbolics.unwrap(eq.rhs)
    end
    # Find target variable
    var_name_clean = clean_name(var_name)
    target_sym = nothing
    for k in keys(obs_dict)
        k_str = clean_name(k)
        if k_str == var_name_clean || occursin(var_name_clean, k_str)
            target_sym = k
            break
        end
    end
    target_sym === nothing && error("Variable $var_name not found in observed equations. Available: $([clean_name(k) for k in keys(obs_dict)])")
    println("[reconstruct_from_observed] target = \"$var_name\"  order=$order  component=$component")
    println("  matched symbol: $target_sym")
    println("  defining RHS:   $(obs_dict[target_sym])")

    # Start with the RHS of the observed equation
    expr = obs_dict[target_sym]
    states_set = Set(Symbolics.unwrap(v) for v in unknowns(sys))
    params_set = Set(Symbolics.unwrap(p) for p in parameters(sys))

    # Flatten observed-in-observed dependencies
    for _ in 1:20
        vars_in_expr = Symbolics.get_variables(expr)
        unknown_vars = filter(v -> !in(v, states_set) && !in(v, params_set) && !isequal(v, tvar_uw), vars_in_expr)
        isempty(unknown_vars) && break
        subs = Dict()
        for u in unknown_vars
            if SymbolicUtils.iscall(u) && SymbolicUtils.operation(u) isa Differential
                base_var = SymbolicUtils.arguments(u)[1]
                haskey(obs_dict, base_var) && (subs[u] = Symbolics.expand_derivatives(Differential(tvar_uw)(obs_dict[base_var])))
            elseif haskey(obs_dict, u)
                subs[u] = obs_dict[u]
            end
        end
        isempty(subs) && break
        expr = Symbolics.substitute(expr, subs)
        expr = Symbolics.expand_derivatives(expr)
    end
    println("  flattened expr (only states + params + t remain):\n    $expr")

    # Identify original ODE state variables still present in expr
    orig_state_names = Set(k[1] for k in keys(h_prob.harmonic_system.variable_map))
    orig_state_map = Dict{Any, Tuple{String, Int}}()
    for v in Symbolics.get_variables(expr)
        in(v, states_set) && continue
        in(v, params_set) && continue
        isequal(v, tvar_uw) && continue
        v_str = clean_name(v)
        found = false
        for name in orig_state_names
            if v_str == name
                orig_state_map[v] = (name, 0); found = true; break
            elseif v_str == name * "ˍt" || endswith(v_str, "₊" * name * "ˍt")
                orig_state_map[v] = (name, 1); found = true; break
            elseif v_str == name * "ˍtt" || endswith(v_str, "₊" * name * "ˍtt")
                orig_state_map[v] = (name, 2); found = true; break
            end
            if !found && SymbolicUtils.iscall(v) && SymbolicUtils.operation(v) isa Differential
                base = SymbolicUtils.arguments(v)[1]
                if clean_name(base) == name
                    orig_state_map[v] = (name, 1); break
                end
            end
        end
    end

    vmap = h_prob.harmonic_system.variable_map
    ωvar_sym = h_prob.harmonic_system.ω_var
    N = h_prob.harmonic_system.N

    println("  states to expand into the harmonic ansatz: $(collect(values(orig_state_map)))")

    # Build symbolic ansatz substitution
    # Use _get_coeff_expr to pre-resolve any vmap symbols that were eliminatedby mtkcompile 
    ansatz_subs = Dict{Any, Any}()
    for (orig_var, (orig_name, deriv_order)) in orig_state_map
        x_t = dx_t = ddx_t = Num(0.0)

        key_dc = find_varmap_key(vmap, orig_name, 0, :Cos)
        key_dc !== nothing && (x_t += _get_coeff_expr(vmap[key_dc], h_prob))

        for n in 1:N
            key_c = find_varmap_key(vmap, orig_name, n, :Cos)
            key_s = find_varmap_key(vmap, orig_name, n, :Sin)
            An = key_c !== nothing ? _get_coeff_expr(vmap[key_c], h_prob) : Num(0.0)
            Bn = key_s !== nothing ? _get_coeff_expr(vmap[key_s], h_prob) : Num(0.0)
            nω = n * ωvar_sym
            x_t   +=  An * cos(nω * t_sym) + Bn * sin(nω * t_sym)
            dx_t  += -An * nω * sin(nω * t_sym) + Bn * nω * cos(nω * t_sym)
            ddx_t += -An * nω^2 * cos(nω * t_sym) - Bn * nω^2 * sin(nω * t_sym)
        end

        ansatz_subs[orig_var] = Symbolics.unwrap(
            deriv_order == 0 ? x_t : deriv_order == 1 ? dx_t : ddx_t)
    end

    # Substitute symbolic ansatz into expression
    val_t = Symbolics.substitute(expr, ansatz_subs)
    val_t = Symbolics.expand_derivatives(val_t)
    println("  after ansatz substitution + derivative expansion:\n    $val_t")

    #  Isolate the target harmonic component
    remaining = Symbolics.get_variables(val_t)
    if any(isequal(v, tvar_uw) for v in remaining)
        val_expanded = Symbolics.expand(val_t)
        isolate = Dict{Any, Any}()
        for k in 1:(2*N+1)
            cos_k = Symbolics.unwrap(cos(Num(k) * ωvar_sym * t_sym))
            sin_k = Symbolics.unwrap(sin(Num(k) * ωvar_sym * t_sym))
            if k == order
                isolate[cos_k] = (component == :Cos ? 1.0 : 0.0)
                isolate[sin_k] = (component == :Sin ? 1.0 : 0.0)
            else
                isolate[cos_k] = 0.0
                isolate[sin_k] = 0.0
            end
        end
        final_expr = Symbolics.substitute(Symbolics.unwrap(val_expanded), isolate)
    else
        final_expr = val_t
    end
    println("  final isolated expression: $final_expr\n")

    return Num(final_expr)
end
  