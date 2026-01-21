using Symbolics
using SymbolicUtils

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

function hbsweep(sys, jls, ns)
    I₀ = 1e-6
    R₀ = 50.0
    Id = 0.05e-6
    ωc = sqrt(2*pi *I₀/(jls.Φ₀*1000.0e-15))/(2*pi)
    Ic = jls.Φ₀/(2*pi*1000.0e-12)
    ω_vec = 2*pi*(4.5:0.001:5.0)*1e9
    N = length(ω_vec)
    solution1 = Vector{Float64}(undef, N)
    solution2 = Vector{Float64}(undef, N)

    ps = [
        jls.P1.Isrc.ω => ω_vec[1]
        jls.P1.Isrc.I => 0.00565e-6
        jls.C1.C      => 100.0e-15
        jls.J1.C      => 1000.0e-15
        jls.J1.I0     => Ic
        jls.P1.Rsrc.R => 50.0
        jls.J1.R      => 1e9
    ]
    u0_vals = zeros(6)
    u0_map = unknowns(sys) .=> u0_vals
    prob = NonlinearProblem(sys, u0_map, ps)
    for i in 1:N
        # Update only frequency using remake
        new_prob = remake(prob, p = [jls.P1.Isrc.ω => ω_vec[i]])
        
        sol = solve(new_prob)
        
        # Calculate magnitudes (Harmonic Ansatz: sqrt(A^2 + B^2))
        solution1[i] = sqrt(sol[ns.C[1]]^2 + sol[ns.D[1]]^2)
        solution2[i] = sqrt(sol[ns.A[1]]^2 + sol[ns.B[1]]^2)
    end

    return ω_vec, solution1, solution2
end
    