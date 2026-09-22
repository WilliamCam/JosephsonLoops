
#transient simulation of whole system
"""
    tsolve(model, u0, param_pairs, tspan; DAE=false, solver_opts=Rodas5(), guesses=nothing)

Solves the circuit in the time domain as an initial value problem. This is the alternative to
harmonic balance, useful for checking a steady state result and for transients that harmonic
balance cannot represent.

# Arguments
- `model`: the compiled model from [`build_circuit`](@ref).
- `u0`: initial state. Note that `build_circuit` returns an empty `u0`, so pass `guesses` here.
- `param_pairs`: parameter values.
- `tspan`: time span in NORMALISED time, so multiply a physical duration by `ωc`.

# Keywords
- `DAE::Bool = false`: build a `DAEProblem` instead of an `ODEProblem`.
- `solver_opts = Rodas5()`: the integrator.
- `guesses = nothing`: initialisation guesses for the unknowns.

# Returns
- the solution object. Index it with a model variable, for example `tsol[C1.i]`.

!!! note
    Timing instrumentation is built in and prints on every call. Time domain solves are far
    slower than harmonic balance for steady state problems.

# Example

tsol = tsolve(model, guesses, ps, (0.0, 1e-6) .* ωc; guesses = guesses)

"""
function tsolve(model, u0, param_pairs, tspan; DAE=false, solver_opts = Rodas5(), guesses = nothing, kwargs...)  
    if DAE
        prob = DAEProblem(model, merge(Dict(u0), Dict(param_pairs)), tspan; kwargs...)
    else
        if guesses !== nothing
            prob = ODEProblem(model, merge(Dict(u0), Dict(param_pairs)), tspan; guesses=guesses, kwargs...)
        else
            prob = ODEProblem(model, merge(Dict(u0), Dict(param_pairs)), tspan; kwargs...)
        end
    end
    sol = @time DifferentialEquations.solve(prob, solver_opts)
    return sol                                                  #Return the solved ODEProblem
end

#Plot a current or voltage of a component (resistor or capacitor)
function tplot(sol::ODESolution, c, model; units = "volts")
    if units == "amps"
        y = sol[c.i][2:end]
        ylabel = "Current (A)"
        label = string(c.i)
    elseif units == "volts"
        y = 1/(sol.t[2]-sol.t[1]) * Φ₀/(2.0*pi) * diff(sol[c.φ])
        ylabel = "Voltage  (V)"
        label = replace(string(c.φ), "φ" => "v")
    elseif units[1] == 'S'
        @assert length(units) == 3 "Error: Please state scattering parameter in form 'Sij'"
        port_i, i = units[2], parse(Int, units[2])
        port_j, j = units[3], parse(Int, units[3])
        S = get_HB_scattering_matrix(model, port_i, port_j)
        y = sol[S[i,j]]
        ylabel = units
        label = nothing
    end
    plot(sol.t[2:end], y, xlabel = "Time (s)", ylabel = ylabel, label = label)
end

# The scattering matrix helper lives in src/harmonic balance/utils.jl as
# get_HB_scattering_matrix. The copy that used to sit here read component fields
# (Rsrc.R, dθ) that no component in this package defines, so it threw on every model
# build_circuit produces.

#solve for the frequency response of some load component when subject to an AC source, by performing an ensemble of transient simulations
function ensemble_fsolve(
        model::ODESystem, u0, tspan, fspan, param_pairs, source,  load; 
        NPts = 1000, Ntraj = 100, alg = Rodas5(), units = "volts", kwargs...
    )
    tsaves = LinRange(tspan[1],tspan[2], NPts)
    ω_vec = 2*pi .* LinRange(fspan[1], fspan[2], Ntraj) 
    prob = ODEProblem(model, u0, tspan, param_pairs, saveat = tsaves; kwargs...)

    function RMS(x)
        return sqrt(mean((x .- mean(x)).^2))
    end

    function RMS_volts(sol,i)
        push!(logger, 1)
        println(string(Ntraj-length(logger)))
        (RMS(1/(sol.t[2]-sol.t[1])*Φ₀/(2*pi)*diff(sol[load.sys.θ])),false)
    end

    function RMS_amps(sol,i)
        push!(logger, 1)
        println(string(Ntraj-length(logger)))
        (RMS(sol[load.sys.i]),false)
    end
    if units == "volts"
        output_func = RMS_volts
    elseif units == "amps"
        output_func = RMS_amps
    end

    ω_index = findfirst(isequal(source.sys.ω), parameters(model))
    function prob_func(prob, i ,repeat)
        prob.p[ω_index] = ω_vec[i]
        prob
    end

    logger = []
    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=output_func)
    sol = DifferentialEquations.solve(ensemble_prob,alg, EnsembleSerial(), trajectories=Ntraj)
    return sol
end

function ensemble_parameter_sweep(
    model::ODESystem, u0, tspan, pspan, param_pairs, parameter,  load; 
    NPts = 1000, Ntraj = 100, alg = Rodas5(), units = "volts", Parallel = false, DAE = false, kwargs...
    )
    tsaves = LinRange(tspan[1],tspan[2], NPts)
    p_vec = LinRange(pspan[1], pspan[2], Ntraj)
    
    if Parallel == true
        method = EnsembleThreads()
    else
        method = EnsembleSerial()
    end

    if DAE == true
        dae_model = dae_index_lowering(model)
        prob = ODAEProblem(dae_model,  u0, tspan, param_pairs, saveat = tsaves; kwargs...)
    else
        prob = ODEProblem(model, u0, tspan, param_pairs, saveat = tsaves; kwargs...)
    end

    function RMS(x)
        return sqrt(mean((x .- mean(x)).^2))
    end

    function RMS_volts(sol,i)
        push!(logger, 1)
        println(string(Ntraj-length(logger)))
        (RMS(1/(sol.t[2]-sol.t[1])*Φ₀/(2*pi)*diff(sol[load.sys.θ])),false)
    end

    function RMS_amps(sol,i)
        push!(logger, 1)
        println(string(Ntraj-length(logger)))
        (RMS(sol[load.sys.i]),false)
    end
    if units == "volts"
        output_func = RMS_volts
    elseif units == "amps"
        output_func = RMS_amps
    end

    p_index = findfirst(isequal(parameter), parameters(model))
    function prob_func(prob, i ,repeat)
        prob.p[1][p_index] = p_vec[i]
        prob
    end

    logger = []
    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=output_func)
    sol = DifferentialEquations.solve(ensemble_prob,alg, method, trajectories=Ntraj)
    return sol
end






    

