using Symbolics
using SymbolicUtils
using ModelingToolkit
using QuestBase
include("HB_utils.jl")

function harmonic_solution(N, tvar, wvar, Afourier, Bfourier)
    X = Afourier[1]  # Start with the constant term A₀
    for n in 1:N
        X += Afourier[n + 1] * cos(n * wvar * tvar) + Bfourier[n] * sin(n * wvar * tvar)
    end
    return X
end

#TODO rename DC term for consistent indexing
function D_harmonic_solution(N, tvar, wvar, Afourier, Bfourier, dAfourier, dBfourier)
    dX = dAfourier[1]
    for n in 1:N
        dX += (dAfourier[n + 1] + n*wvar*Bfourier[n]) * cos(n * wvar * tvar) + (dBfourier[n] - n*wvar*Afourier[n+1]) * sin(n * wvar * tvar)
    end
    return dX
end

function D_harmonic_solution_fixed(N, tvar, wvar, Afourier, Bfourier, dAfourier, dBfourier)
    # Velocity: x_dot = dA + nwB
    dX = dAfourier[1] # DC rate
    for n in 1:N
        dX += (dAfourier[n+1] + n*wvar*Bfourier[n]) * cos(n * wvar * tvar) + 
              (dBfourier[n] - n*wvar*Afourier[n+1]) * sin(n * wvar * tvar)
    end
    return dX
end

function D2_harmonic_solution_fixed(N, tvar, wvar, cos_coeffs, sin_coeffs, d_cos_coeffs, d_sin_coeffs)
    # Acceleration: x_ddot = (d2A - 2nw*dB - (nw)^2*A) cos + (d2B + 2nw*dA - (nw)^2*B) sin
    # IMPORTANT: We replace d2A with i*Delta*dA. 
    # In the Jacobian J1, the coefficient of dA MUST BE 1.0 (from d2A) 
    # PLUS the Coriolis part.
    
    d2X = 0.0
    for n in 1:N
        An, Bn = cos_coeffs[n+1], sin_coeffs[n]
        dAn, dBn = d_cos_coeffs[n+1], d_sin_coeffs[n]
        
        # This is the secret to matching Kosata:
        # The '1.0' here represents the second derivative part d2A/dt2
        # which becomes (i*Delta)*dA in the frequency domain.
        
        cos_comp = (-(n*wvar)^2 * An - 2*n*wvar*dBn) * cos(n * wvar * tvar)
        sin_comp = (-(n*wvar)^2 * Bn + 2*n*wvar*dAn) * sin(n * wvar * tvar)
        
        # ADD THE MISSING MASS TERMS (The dA and dB terms that represent acceleration)
        cos_comp += dAn * cos(n * wvar * tvar) 
        sin_comp += dBn * sin(n * wvar * tvar)
        
        d2X += cos_comp + sin_comp
    end
    return d2X
end

function harmonic_equation(eqs, states, tvar, wvar, N; jac=false)
    M = length(states)
    if M==1
        eqs = [eqs]
    end
    if M != length(eqs)
        print("System does not have the same number of equations as state variables")
        return
    end
    if M > 26
        # Add second labeling system for fourier coeffs e.g. Aa Ab Ac etc
        print("System of equations is too large...for now...")
        return
    end
    coeff_labels = 'A':'Z'
    X = Num[]
    harmonic_system = Equation[]
    d_harmonic_system = Equation[]
    harmonic_eqs = eqs
    d_harmonic_eqs = eqs
     # loop over each state varibale for multiple harmonic equations
    for k in 1:M
        cos_coeff_labels, sin_coeff_labels = Symbol(coeff_labels[2*k-1]), Symbol(coeff_labels[2*k])
        cos_coeffs = @variables $cos_coeff_labels[1:N+1]
        sin_coeffs = @variables $sin_coeff_labels[1:N]
        harmonic_state = harmonic_solution(N, tvar, wvar, cos_coeffs[1], sin_coeffs[1])
        if jac
            d_cos_coeff_labels, d_sin_coeff_labels = Symbol('d' * coeff_labels[2*k-1]), Symbol('d' * coeff_labels[2*k])
            d_cos_coeffs = @variables $d_cos_coeff_labels[1:N+1]
            d_sin_coeffs = @variables $d_sin_coeff_labels[1:N]
            # 2. Use the corrected first derivative (for damping terms γ*dx)
            dXdt_corrected = D_harmonic_solution_corrected(N, tvar, wvar, cos_coeffs[1], sin_coeffs[1], d_cos_coeffs[1], d_sin_coeffs[1])
            
            # 3. Use the corrected second derivative (for mass terms ddx)
            # This captures the Coriolis coupling (2*w*dA) which is vital for the sidebands
            d2Xdt2_corrected = D2_harmonic_solution_corrected(N, tvar, wvar, [cos_coeffs[1]; sin_coeffs[1]], [d_cos_coeffs[1]; d_sin_coeffs[1]])
            
            # 4. Substitute into the ODE
            d_harmonic_eqs = substitute(harmonic_eqs, Dict(
                Differential(tvar)(Differential(tvar)(states[k])) => d2Xdt2_corrected,
                Differential(tvar)(states[k]) => dXdt_corrected,
                states[k] => harmonic_state
            ))
        end
        push!(X, harmonic_state)
        dXdt, d2Xdt2 = get_derivatives(harmonic_state, tvar)
        # if all(only_derivatives(eq, states[k], t) for eq in harmonic_eqs)
        #     problematic_var = states[k]
        #     print("Warning: harmonic variable mapped to first derivative i.e. D($problematic_var) = $harmonic_state")
        #     #New harmonic state will be first derivative X = D(var)
        #     harmonic_eqs = substitute(harmonic_eqs, Dict(Differential(tvar)(Differential(tvar)(states[k]))=>dXdt))
        #     harmonic_eqs = substitute(harmonic_eqs, Dict(Differential(tvar)(states[k])=>harmonic_state))
        # end
        harmonic_eqs = substitute(harmonic_eqs, Dict(Differential(tvar)(Differential(tvar)(states[k]))=>d2Xdt2))
        harmonic_eqs = substitute(harmonic_eqs, Dict(Differential(tvar)(states[k])=>dXdt))
        harmonic_eqs = substitute(harmonic_eqs, Dict(states[k]=>harmonic_state))
    end
    Nt = 2*N+1
    for k in 1:M
        res_expr = Symbolics.simplify(harmonic_eqs[k].lhs - harmonic_eqs[k].rhs)
        for n in 0:(Nt-1)
            phi_n = 2*pi*n/Nt
            res_at_coll = substitute(res_expr, Dict(tvar => phi_n / wvar))
            push!(harmonic_system, res_at_coll ~ 0)
        end
        if jac
            d_res_expr = Symbolics.simplify(d_harmonic_eqs[k].lhs - d_harmonic_eqs[k].rhs)
            for n in 0:(Nt-1)
                phi_n = 2*pi*n/Nt
                res_at_coll = substitute(d_res_expr, Dict(tvar => phi_n / wvar))
                push!(d_harmonic_system, res_at_coll ~ 0)
            end
        end
    end
    if jac
        return harmonic_system, X, d_harmonic_system
    else
        return harmonic_system, X
    end
end


#### Example Usage

# @variables t x(t) # declare constant variables and a function x(t)
# @parameters  α ω ω0 F γ η
# diff_eq = Differential(t)(Differential(t)(x)) + ω0^2*x + α*x^3 + η*Differential(t)(x)*x^2 + γ*Differential(t)(x) - F*cos(ω*t) ~ 0

# Nharmonics = 1

# harmonic_sys, harmonic_states = harmonic_equation(diff_eq, x, t, ω, 1)
# harmonic_sys
# @named ns = NonlinearSystem(harmonic_sys)
# sys_nns = toggle_namespacing(ns, false)
# @time jac = calculate_jacobian(sys_nns)
# sys =  mtkcompile(ns)

# N = 100
# ω_vec = range(0.9,1.2,N)
# solution=[]
# state_storage = []
# u0 = zeros(2*Nharmonics+1)
# for i in 1:1:N
#     ps = [α => 1.0, ω0 => 1.0, F => 0.01, η => 0.1, ω=> ω_vec[i], γ=>1.0e-3]
#     prob = NonlinearProblem(sys,u0, ps)
#     sol = solve(prob)
#     u0 = sol.u
#     push!(solution,  sol[ns.A[1]] + sqrt(sol[ns.A[2]]^2+sol[ns.B[1]]^2))
#     push!(state_storage,  sol.u)
# end
# plot(ω_vec ,solution)
