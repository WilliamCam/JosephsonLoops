
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

function has_derivative(equation)
    return any(eq -> !isempty(Symbolics.filterchildren(Symbolics.is_derivative, eq)), equation)
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
        a[k] = 0.5/sqrt(port_k.Rₙ.r)*(port_k.dφ+port_k.Rₙ.r*port_k.i)
        b[k] = 0.5/sqrt(port_k.Rₙ.r)*(port_k.dφ-port_k.Rₙ.r*port_k.i)
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

function get_full_equations(model::ModelingToolkit.System)
    # copy: unknowns(model)/full_equations(model) return references into the model, and the
    # deleteat! calls below would otherwise mutate it (corrupting it for the next call).
    eqs = copy(full_equations(model))
    states = copy(unknowns(model))
    tvar = ModelingToolkit.get_iv(model)

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
    return eqs, states, diffvars, diff2vars
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


function build_jacobians(rotated_system, vars, dvars, d2vars = nothing)
    
    _jac = Symbolics.jacobian(rotated_system, vars)
    zeroed = d2vars === nothing ? Dict(dvars .=> 0) : Dict(vcat(dvars, d2vars) .=> 0)
    jac_0 = Num.((substitute(_jac, zeroed)))
    jac_1 = Symbolics.jacobian(rotated_system, dvars)
    d2vars === nothing && return jac_0, jac_1
    jac_2 = Symbolics.jacobian(rotated_system, d2vars)
    return jac_0, jac_1, jac_2
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

function rotate_to_harmonic_frame(M, n_ints::Vector{Int}, Nt::Int, harmonic_system)
    block_rows = 2 * length(n_ints) + 1
    total_rows = M * block_rows
    total_cols = M * Nt
    Γ_total = zeros(Num, total_rows, total_cols)
    Γ_single = Matrix{Num}(undef, block_rows, Nt)
    for j in 1:Nt
        Γ_single[1, j] = Num(1//Nt)
        for (k, nk) in enumerate(n_ints)
            phase = nk * (j - 1) * (2π / Nt)
            Γ_single[2k, j]     = Num((2//Nt) * cos(phase))
            Γ_single[2k + 1, j] = Num((2//Nt) * sin(phase))
        end
    end
    for d in 1:M
        row_range = (d-1)*block_rows + 1 : d*block_rows
        col_range = (d-1)*Nt + 1 : d*Nt
        Γ_total[row_range, col_range] .= Γ_single
    end
    #ordering preserved: residual blocks are created per state in the same order
    rotated_system = Γ_total * [equation.lhs for equation in harmonic_system]
    return rotated_system
end

# 2D torus (hyper-time) projection: samples live on the product grid
# (θ1_i, θ2_j) = (2πi/Nt1, 2πj/Nt2), flattened with i outer / j inner to match the
# sampler. Orthogonality is exact per axis (Nt = 2·max|index|+1, both odd), for ANY
# tone ratio — ω1, ω2 never appear. Row order per state matches the 1D versions:
# [DC, cos₁, sin₁, ...] in basis-slot order.
function rotate_to_harmonic_frame(M, indices::Vector{Tuple{Int,Int}}, Nt1::Int, Nt2::Int, harmonic_system)
    block_rows = 2 * length(indices) + 1
    Nt = Nt1 * Nt2
    total_rows = M * block_rows
    total_cols = M * Nt
    Γ_total = zeros(Num, total_rows, total_cols)
    Γ_single = Matrix{Num}(undef, block_rows, Nt)
    for i in 0:(Nt1-1), j in 0:(Nt2-1)
        col = i * Nt2 + j + 1
        Γ_single[1, col] = Num(1//Nt)
        for (k, (m, n)) in enumerate(indices)
            phase = m * i * (2π / Nt1) + n * j * (2π / Nt2)
            Γ_single[2k, col]     = Num((2//Nt) * cos(phase))
            Γ_single[2k + 1, col] = Num((2//Nt) * sin(phase))
        end
    end
    for d in 1:M
        row_range = (d-1)*block_rows + 1 : d*block_rows
        col_range = (d-1)*Nt + 1 : d*Nt
        Γ_total[row_range, col_range] .= Γ_single
    end
    #ordering preserved: residual blocks are created per state in the same order
    rotated_system = Γ_total * [equation.lhs for equation in harmonic_system]
    return rotated_system
end

function diamond_truncation_indices(N::Int)
    indices = Tuple{Int, Int}[]
    for m in 0:N
        n_min = (m == 0) ? 0 : -N
        for n in n_min:N
            if m + abs(n) <= N
                push!(indices, (m, n))
            end
        end
    end
    return indices
end

function full_indices(Np::Int, Ni::Int; single_tone::Bool=false)
    indices = Tuple{Int, Int}[]
    max_bound = max(Np, Ni)
    for m in 0:max_bound
        # If single_tone is true, we only care about m, n=0
        range_n = single_tone ? (0:0) : (m == 0 ? 0 : -max_bound):max_bound
        
        for n in range_n
            is_pump = (m <= Np && n == 0) || (m == 0 && abs(n) <= Np)
            is_im = (m + abs(n) <= Ni)
            if is_pump || is_im
                push!(indices, (m, n))
            end
        end
    end
    return indices
end

function construct_fourier_basis(Np::Int, states::Vector{Num}; intermod_order=0, single_tone = false)
    K = length(states)
    @assert Np >= 1 "Need at least one pump harmonic"
    coeff_labels = 'A':'Z'
    intermod_order > Np ? indices = diamond_truncation_indices(intermod_order) : indices = full_indices(Np, intermod_order, single_tone = single_tone)
    #TODO: add optinal argument to remove DC term
    N_terms = length(indices) - 1
    dc_terms = @variables DC[1:K]
    d_dc_terms = @variables d_DC[1:K]
    d2_dc_terms = @variables d2_DC[1:K]
    #TODO: bigger system size
    fourier_basis_map = Dict{Num, FourierBasis}()
    @assert K < 13 "Have run out of letters in alphabet for state variables... will fix later"
    for k in 1:K
        cur_DC_term = DC[k]
        cos_labels, sin_labels = Symbol(coeff_labels[2*k-1]), Symbol(coeff_labels[2*k])
        cos_sym_arr = (@variables $cos_labels[1:N_terms])[1]
        sin_sym_arr = (@variables $sin_labels[1:N_terms])[1]

        cur_d_DC_term = d_DC[k]
        d_cos_labels, d_sin_labels = Symbol('d' * coeff_labels[2*k-1]), Symbol('d' * coeff_labels[2*k])
        d_cos_sym_arr = (@variables $d_cos_labels[1:N_terms])[1]
        d_sin_sym_arr = (@variables $d_sin_labels[1:N_terms])[1]

        cur_d2_DC_term = d2_DC[k]
        d2_cos_labels, d2_sin_labels = Symbol("d2" * coeff_labels[2*k-1]), Symbol("d2" * coeff_labels[2*k])
        d2_cos_sym_arr = (@variables $d2_cos_labels[1:N_terms])[1]
        d2_sin_sym_arr = (@variables $d2_sin_labels[1:N_terms])[1]

        #Create varaible map
        index_map = Dict{Tuple{Tuple{Int,Int},Symbol}, Num}()
        index_map[(indices[1], :DC)] = cur_DC_term
        for n in 1:N_terms
            index_map[(indices[n+1], :Cos)] = cos_sym_arr[n]
            index_map[(indices[n+1], :Sin)] = sin_sym_arr[n]
        end
        cur_basis = FourierBasis(cur_DC_term, cur_d_DC_term, cos_sym_arr, sin_sym_arr, d_cos_sym_arr, d_sin_sym_arr,
            cur_d2_DC_term, d2_cos_sym_arr, d2_sin_sym_arr, indices, index_map)
        fourier_basis_map[states[k]] = cur_basis
    end
    return fourier_basis_map
end
