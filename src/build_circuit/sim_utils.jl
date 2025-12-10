
#transient simulation of whole system
function tsolve(model::System, u0::Vector{Pair{Num, Float64}}, 
    param_pairs::Vector{Pair{Num, Float64}}, tspan::Tuple{Float64,Float64}; DAE=false, solver_opts = Rodas5(),  kwargs...)  
    if DAE
        prob = DAEProblem(model, merge(Dict(u0), Dict(param_pairs)), tspan; kwargs...)
    else
        prob = ODEProblem(model, merge(Dict(u0), Dict(param_pairs)), tspan; kwargs...)
    end
    sol = @time solve(prob, solver_opts)
    return sol                                                  #Return the solved ODEProblem
end

#Plot a current or voltage of a component (resistor or capacitor)
function tplot(sol::ODESolution, c, model; units = "volts")
    if units == "amps"
        y = sol[c.i][2:end]
        ylabel = "Current (A)"
        label = string(c.i)
    elseif units == "volts"
        y = 1/(sol.t[2]-sol.t[1]) * Φ₀/(2.0*pi) * diff(sol[c.θ])
        ylabel = "Voltage  (V)"
        label = replace(string(c.θ), "θ" => "v")
    elseif units[1] == 'S'
        @assert length(units) == 3 "Error: Please state scattering parameter in form 'Sij'"
        port_i, i = units[2], parse(Int, units[2])
        port_j, j = units[3], parse(Int, units[3])
        S = get_scattering_matrix(model,port_i,port_j)
        y = sol[S[i,j]]
        ylabel = units
        label = nothing
    end
    plot(sol.t[2:end], y, xlabel = "Time (s)", ylabel = ylabel, label = label)
end

function get_scattering_matrix(model::System,i::Char,j::Char)
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
    return a / b
end


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






    

