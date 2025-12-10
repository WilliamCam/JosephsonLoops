#Netlist struct
"""
    process_netlist(loops::Vector{Vector{String}}; mutual_coupling::Vector{Tuple{Int64, Int64}} = [], ext_flux::Vector{Int64} = []) -> CircuitNetlist

Processes a circuit netlist defined by loops and creates a `CircuitNetlist` structure containing the components, mutual coupling information, and external flux values.

# Arguments
- `loops::Vector{Vector{String}}`: A vector of vectors, where each inner vector represents a loop containing component names as strings.
- `mutual_coupling::Vector{Tuple{Int64, Int64}}`: An optional vector of tuples representing mutual coupling between loops. Defaults to an empty vector.
- `ext_flux::Vector{Bool}`: An optional vector of external flux offsets associated with the loops. Defaults to an empty vector.

# Returns
- `CircuitNetlist`: A `CircuitNetlist` object that encapsulates the processed information, including loops, mutual coupling, external flux, component mappings, and identified branches.

# Example

loops = [["R1", "C1"], ["I1", "R1"], ["J1", "L1"]]
circuit = process_netlist(loops; mutual_coupling=[(1, 2)], ext_flux=[1])

"""
struct CircuitNetlist
    loops::Vector{Vector{String}}
    mutual_coupling::Vector{Tuple{Int64, Int64}}
    ext_flux::Vector{Bool}
    components::Dict{String, Vector{Int}}
    branches::Vector{String}
end

function find_components(loops::Vector{Vector{String}})
    numLoops = size(loops)[1]                      
    component_loop_mapping = Dict{String, Vector{Int}}()
    branches = Vector{String}()               
    for i in 1:numLoops                        
        for j in 1:length(loops[i])            
            component_loop_mapping[loops[i][j]]=push!(get(component_loop_mapping, loops[i][j], []), i-1) 
        end
    end

    for comp in keys(component_loop_mapping)
        if (comp[1] in ['I', 'R', 'C', 'J', 'L', 'P'])
            push!(branches, comp)
        end
    end
    return component_loop_mapping, branches
end

function process_netlist(
    loops::Vector{Vector{String}}; 
    mutual_coupling::Vector{Tuple{Int,Int}} = Vector{Tuple{Int,Int}}(), 
    ext_flux::Vector{Bool} = Vector{Bool}()
)
    numLoops = size(loops)[1]
    if ext_flux == []
        flux_vector = fill(false, numLoops)
    else
        flux_vector = ext_flux
    end
    if length(flux_vector) != length(loops)
        throw(ArgumentError("Mutual flux vector must be equal to number of circuit loops"))
    end
    component_loop_mapping, branches = find_components(loops)
    circuit = CircuitNetlist(loops, mutual_coupling, flux_vector, component_loop_mapping, branches)
    return circuit
end

"""
    build_circuit(circuit::CircuitNetlist) -> ODESystem

Constructs a ModelingToolkit ODESystem from a CircuitNetlist by creating and connecting circuit components.

# Arguments
- `circuit::CircuitNetlist`: A CircuitNetlist object containing:
  - `loops`: Vector of vectors defining circuit loops
  - `components`: Dictionary mapping component names to their loop indices
  - `mutual_coupling`: Vector of tuples defining mutual inductance coupling
  - `branches`: Vector of component names that form branches

# Returns
- A compiled ModelingToolkit ODESystem representing the circuit equations

# Details
- Creates ModelingToolkit components based on first letter of component name:
  - 'R': Resistor
  - 'C': Capacitor 
  - 'J': Josephson Junction
  - 'I': Current Source
  - 'L': Inductor
- Automatically handles mutual inductance coupling
- Creates ground node connections for components in single loops
- Composes final system from component equations and connections

# Example

loops = [["R1", "C1"], ["I1", "R1"]]
circuit = process_netlist(loops)
model = build_circuit(circuit)

# Throws
- `ArgumentError`: If an unrecognized component type is found in the netlist
"""

function build_circuit(circuit::CircuitNetlist)
    #Load circuit from netlist
    loops = circuit.loops
    numLoops = size(loops)[1]
    component_loop_mapping = circuit.components
    mutual_coupling = circuit.mutual_coupling
    branches = circuit.branches
    ext_flux = circuit.ext_flux

    #Use meta programming to construct MTK component objects parsing their given name.
    built_components = Dict()                               
    for j in branches                                      
        println(j)                                          
        if (j[1] == 'R')                                    
            new_c = "@named $j = Resistor()"
            new_c = Meta.parse(new_c)                      
            new_c = eval(new_c)
            built_components[j] = new_c                     
        elseif (j[1] == 'C')                                
            new_c = "@named $j = Capacitor()"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        elseif (j[1] == 'J')                                
            new_c = "@named $j = JosephsonJunction()"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        elseif (j[1] == 'I')                                
            new_c = "@named $j = CurrentSource()"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        elseif (j[1] == 'L')                                
            new_c = "@named $j = Inductor()"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        elseif (j[1] == 'P')                                
            new_c = "@named $j = Port()"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        else
            throw(ArgumentError("Error: Netlist component was not recognised.Check docs for supported cirucit components and naming conventions."))
        end
    end

    for n in mutual_coupling
        eval(Meta.parse("@named M" * string(n[1])*string(n[2]) * "= Inductor()"))
        built_components["M" *string(n[1])*string(n[2])] = eval(Meta.parse("M" * string(n[1])*string(n[2]))) 
    end
    for (i,k) in enumerate(ext_flux)
        if k
            eval(Meta.parse("@named Φₑ" * string(i) * "= ExternalFlux()"))
            built_components["Φₑ" * string(i)] = eval(Meta.parse("Φₑ" * string(i)))
        end
    end                                 
    already_in_loop = String[]
    eqs=Equation[]
    ground_loop_connectables = []
    for i in 1:length(loops)
        #TODO: assert that ports are always grounded
        connectables = []
        for component_name in keys(component_loop_mapping)
            if i-1 in component_loop_mapping[component_name]
                comp_system = built_components[component_name]
                if component_name in already_in_loop
                    push!(connectables, comp_system.out)
                else
                    push!(connectables, comp_system.in)
                    push!(already_in_loop, component_name)
                end
            end
        end
        for M in mutual_coupling
            if i in M
                component_name = "M$(M[1])$(M[2])"
                M_sys = built_components[component_name]
                if component_name in already_in_loop
                    push!(connectables, M_sys.out)
                else
                    push!(connectables, M_sys.in)
                    push!(already_in_loop, component_name)
                end
            end
        end
        if ext_flux[i]
            component_name = "Φₑ$(i)"
            k_sys = built_components[component_name]
            push!(connectables, k_sys.in)
            push!(ground_loop_connectables, k_sys.out)
        end
        push!(eqs, connect(connectables...))
    end

    
    for (comp, vals) in component_loop_mapping
        if length(vals) == 1
            comp_system = built_components[comp]
            push!(ground_loop_connectables, comp_system.out)
        end
    end

    eval(Meta.parse("@named ground" * "= GroundLoop()"))
    built_components["ground"] = eval(Meta.parse("ground"))
    g_sys = built_components["ground"]
    push!(ground_loop_connectables, g_sys.g)
    push!(eqs, connect(ground_loop_connectables...))
    
    sys = [component[2] for component in built_components]
    @named _model = System(eqs, t)
    @named model = compose(_model,sys)
    #TODO: Make calling Dict more readable i.e. c[2]
    print(built_components)


    new_model = mtkcompile(model)

    guesses = Pair{Num,Float64}[]
    u0 = Pair{Num,Float64}[]
    for state_var in unknowns(new_model)
        push!(guesses, state_var => 0.0)
        #TODO: Addition of ports creates initialization
        #push!(u0, D(state_var) => 0.0)
    end

    return new_model, u0, guesses                             
end





# loops = [["I1", "C1"],
# ["C1", "R2"]]
# ext_flux = [false, false]

# circuit = process_netlist(loops, ext_flux=ext_flux)

# #cirucit model ODAE system and initial condition vector are created.
# model = build_circuit(circuit)

