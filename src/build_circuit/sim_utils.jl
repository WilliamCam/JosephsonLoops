using ModelingToolkit
using Statistics: mean
using SymbolicIndexingInterface: parameter_index

const FLUX_QUANTUM = 2.067833848e-15

function _component_system(c)
    try
        return getproperty(c, :sys)
    catch err
        if err isa ArgumentError
            return c
        else
            rethrow()
        end
    end
end

function _parameter_index(model, sym_or_idx)
    if sym_or_idx isa ModelingToolkit.ParameterIndex
        return sym_or_idx
    end
    idx = parameter_index(model, sym_or_idx)
    if idx === nothing
        error("Parameter $(sym_or_idx) is not part of the provided model.")
    end
    return idx
end

function _build_op_map(model, u0, param_pairs)
    state_map = Dict{Any,Any}(unknowns(model) .=> 0.0)
    if all(x -> x isa Pair, u0)
        for (var, val) in u0
            state_map[var] = val
        end
    elseif !isempty(u0)
        for (var, val) in zip(unknowns(model), u0)
            state_map[var] = val
        end
    end
    merge(state_map, Dict(param_pairs))
end

#transient simulation of whole system
function tsolve(model, u0, tspan, param_pairs; alg = Rodas5(), kwargs...)      
    prob = ODEProblem(model, _build_op_map(model, u0, param_pairs), tspan; kwargs...)   #Create an ODEProblem to solve for a specified time
    sol = solve(prob, alg)
    return sol                                                  #Return the solved ODEProblem
end

#Plot a current or voltage of a component (resistor or capacitor)
function tplot(sol::ODESolution, c; units = "volts")
    comp = _component_system(c)
    if units == "amps"
        y = sol[comp.i][2:end]
        ylabel = "Current (A)"
        label = string(comp.i)
    else
        y = 1/(sol.t[2]-sol.t[1]) * FLUX_QUANTUM/(2.0*pi) * diff(sol[comp.I,])
        ylabel = "Voltage  (V)"
        label = replace(string(comp.I,), "I," => "v")
    end
    plot(sol.t[2:end], y, xlabel = "Time (s)", ylabel = ylabel, label = label)
end

#solve for the frequency response of some load component when subject to an AC source, by performing an ensemble of transient simulations
function ensemble_fsolve(
        model::ODESystem, u0, tspan, fspan, param_pairs, source,  load; 
        NPts = 1000, Ntraj = 100, alg = Rodas5(), units = "volts", kwargs...
    )
    tsaves = LinRange(tspan[1],tspan[2], NPts)
    I_vec = 2*pi .* LinRange(fspan[1], fspan[2], Ntraj) 
    prob = ODEProblem(model, _build_op_map(model, u0, param_pairs), tspan, saveat = tsaves; kwargs...)

    load_comp = _component_system(load)
    source_comp = _component_system(source)
    source_param_idx = _parameter_index(model, source_comp.I)

    function RMS(x)
        return sqrt(mean((x .- mean(x)).^2))
    end

    function RMS_volts(sol,i)
        push!(logger, 1)
        println(string(Ntraj-length(logger)))
        (RMS(1/(sol.t[2]-sol.t[1])*FLUX_QUANTUM/(2*pi)*diff(sol[load_comp.I,])),false)
    end

    function RMS_amps(sol,i)
        push!(logger, 1)
        println(string(Ntraj-length(logger)))
        (RMS(sol[load_comp.i]),false)
    end
    if units == "volts"
        output_func = RMS_volts
    elseif units == "amps"
        output_func = RMS_amps
    end

    function prob_func(prob, i, repeat)
        new_p = copy(prob.p)
        new_p[source_param_idx] = I_vec[i]
        return remake(prob; p = new_p)
    end

    logger = []
    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=output_func)
    sol = solve(ensemble_prob,alg, EnsembleSerial(), trajectories=Ntraj)
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
        prob = ODAEProblem(dae_model,  _build_op_map(model, u0, param_pairs), tspan, saveat = tsaves; kwargs...)
    else
        prob = ODEProblem(model, _build_op_map(model, u0, param_pairs), tspan, saveat = tsaves; kwargs...)
    end

    load_comp = _component_system(load)

    function RMS(x)
        return sqrt(mean((x .- mean(x)).^2))
    end

    function RMS_volts(sol,i)
        push!(logger, 1)
        println(string(Ntraj-length(logger)))
        (RMS(1/(sol.t[2]-sol.t[1])*FLUX_QUANTUM/(2*pi)*diff(sol[load_comp.I,])),false)
    end

    function RMS_amps(sol,i)
        push!(logger, 1)
        println(string(Ntraj-length(logger)))
        (RMS(sol[load_comp.i]),false)
    end
    if units == "volts"
        output_func = RMS_volts
    elseif units == "amps"
        output_func = RMS_amps
    end

    parameter_idx = _parameter_index(model, parameter)
    function prob_func(prob, i, repeat)
        new_p = copy(prob.p)
        new_p[parameter_idx] = p_vec[i]
        return remake(prob; p = new_p)
    end

    logger = []
    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func, output_func=output_func)
    sol = solve(ensemble_prob,alg, method, trajectories=Ntraj)
    return sol
end
