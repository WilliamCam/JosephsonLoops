
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
    # copy: unknowns(model)/full_equations(model) return references into the model, and the
    # deleteat! calls below would otherwise mutate it (corrupting it for the next call).
    eqs = copy(full_equations(model))
    states = copy(unknowns(model))

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
    # Columns of J0 follow `vars` and columns of J1 follow `dvars`; harmonic_equation
    # assembles both in the same order, so column i of J1 is d/d(vars[i]').
    _jac = Symbolics.jacobian(rotated_system, vars)
    jac_0 = Num.((substitute(_jac, Dict(dvars .=> 0))))
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
    # Single-variable Fourier projection of the collocation residuals. Rows are
    # interleaved [DC, cos₁, sin₁, cos₂, sin₂, ...] so that equation rows carry the same
    # ordering as the `vars` jacobian columns ([DC, cos₁, sin₁, ...] per state):
    # one index map then addresses both the perturbation vector and the response vector.
    Γ_single = Matrix{Num}(undef, block_rows, Nt)
    for j in 1:Nt
        # DC Term
        Γ_single[1, j] = Num(1//Nt)
        # AC Terms
        for n in 1:N
            phase = n * (j - 1) * (2π / Nt)
            Γ_single[2n, j] = Num((2//Nt) * cos(phase))          # Cosine
            Γ_single[2n + 1, j] = Num((2//Nt) * sin(phase))      # Sine
        end
    end
    # Place the single-variable operator into the block diagonal
    for d in 1:M
        row_range = (d-1)*block_rows + 1 : d*block_rows
        col_range = (d-1)*Nt + 1 : d*Nt
        Γ_total[row_range, col_range] .= Γ_single
    end
    #ordering should be preserved as equations in colocation.jl are created in this order
    rotated_system = Γ_total * [equation.lhs for equation in harmonic_system]
    return (rotated_system)
end

#Creates system pertubation response in harmonic frame
function perturbation_response(h_sys::HarmonicSystem, source_param::Num, parameters::Dict; amplitude::Float64 = 1.0)
    t = Num(ModelingToolkit.get_iv(h_sys.time_domain_system))
    eqs, _states = get_full_equations(h_sys.time_domain_system, t)
    N, ω = h_sys.N, h_sys.ω
    Nt = 2N + 1
    U = zeros(ComplexF64, length(eqs) * Nt)        # physical (real, single-quadrature) drive

    fixed_params = copy(parameters)
    delete!(fixed_params, ω)
    system_response_symbolic = Symbolics.jacobian([eq.lhs - eq.rhs for eq in eqs], [source_param])
    system_response = Symbolics.substitute(system_response_symbolic, fixed_params)
    @assert any(x -> !isequal(x, Num(0)), system_response) "Parameter $source_param not found in time domain system"  

    zero_harmonics = Dict{Any, Any}()
    for n in 1:(2N + 1)
        zero_harmonics[Symbolics.unwrap(cos(Num(n) * ω * t))] = 0.0
        zero_harmonics[Symbolics.unwrap(sin(Num(n) * ω * t))] = 0.0
    end

    for (k, eq) in enumerate(system_response)
        isequal(Symbolics.simplify(eq), Num(0)) && continue
        base = (k - 1) * Nt
        U[base + 1] -= amplitude * Symbolics.substitute(eq, zero_harmonics)
        for n in 1:N
            c_cos = Symbolics.coeff(eq, cos(Num(n) * ω * t))
            c_sin = Symbolics.coeff(eq, sin(Num(n) * ω * t))
            if c_cos == 0.0 && c_sin == 0.0
                continue
            end
            U[base + 2n + 1] -= amplitude * (c_cos - im * c_sin)

            U[base + 2n] -= amplitude * (c_sin + im * c_cos)
        end
    end
    return U
end

using Symbolics, SymbolicUtils
using SymbolicUtils.Rewriters


@variables t
@variables A(t), B(t)
D2 = Differential(t)^2

# Define individual rules
# Rule 1: A'' / B'' => 0
# Rule 2: B'' / A'' => 0
rules = Chain([
    @rule(~c * D2(~f) / D2(~g) => 0),
    @rule(~c * D2(~g) / D2(~f) => 0)
])

# Example expression containing both types of terms
expr = (D2(A)/D2(B)) + (D2(B)/D2(A)) + A(t)

# Apply the rules
simplified_expr = Postwalk(rules)(expr)
# Result: A(t)
