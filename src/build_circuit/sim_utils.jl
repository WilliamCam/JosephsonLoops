function parameter_set(ps, param::Num, value::Float64)
    return push!(ps, param => value)
end

#Solve initial conditions
function solve_ini(model, old_u0, t_end, ps; alg = Rodas5(), kwargs...)
    tspan_ini = (0.0, t_end)                                #Create a timespan
    prob = ODEProblem(model, old_u0, tspan_ini, ps, save_everystep = false; kwargs...)  #Create an ODEProblem to solve for a specified time only saving the final component variable values
    sol = solve(prob, alg)                                       #Solve the ODEProblem
    new_u0 = sol[:,end]                                     #Set the new initial conditions to the 
    return new_u0                                           #return the new intial conditions
end

#transient simulation of whole system
function tsolve(model, u0, tspan, param_pairs; alg = Rodas5(), kwargs...)      
    prob = ODEProblem(model, u0, tspan, param_pairs; kwargs...)   #Create an ODEProblem to solve for a specified time
    sol = solve(prob, alg)
    return sol                                                  #Return the solved ODEProblem
end

#Plot a current or voltage of a component (resistor or capacitor)
function tplot(sol::ODESolution, c::Component; units = "volts")
    if units == "amps"
        y = sol[c.sys.i][2:end]
        ylabel = "Current (A)"
        label = string(c.sys.i)
    else
        y = 1/(sol.t[2]-sol.t[1]) * Φ₀/(2.0*pi) * diff(sol[c.sys.θ])
        ylabel = "Voltage  (V)"
        label = replace(string(c.sys.θ), "θ" => "v")
    end
    plot(sol.t[2:end], y, xlabel = "Time (s)", ylabel = ylabel, label = label)
end

#solve for the frequency response of some load component when subject to an AC source, by performing an ensemble of transient simulations
function ensemble_fsolve(
        model::ODESystem, u0, tspan, fspan, param_pairs, source,  load::Component; 
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
    sol = solve(ensemble_prob,alg, EnsembleSerial(), trajectories=Ntraj)
    return sol
end

function ensemble_parameter_sweep(
    model::ODESystem, u0, tspan, pspan, param_pairs, parameter,  load::Component; 
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
    sol = solve(ensemble_prob,alg, method, trajectories=Ntraj)
    return sol
end






    

