#Netlist struct
struct CircuitNetlist
    loops::Vector{Vector{String}}
    mutualInd::Vector{Tuple{Int64, Int64}}
    comps::Dict{String, Vector{Int}}
    params::Dict{String, Float64}
    junctions::Vector{String}
    σB::Matrix{Int64}
    phaseParity::Dict{String, Vector{Int64}}
end

#Find component parameters
function find_components(loops::Vector{Vector{String}})
    numLoops = size(loops)[1]                      #Number of loops in circuit
    componentLoopDict = Dict{String, Vector{Int}}()
    componentParamDict = Dict{String, Float64}()                    #Dictionary with components as keys and loops as values (used to find unique elements)
    junctions = Vector{String}()                              #Stores the names of the junctions
    for i in 1:numLoops                         #Iterate through all loops to find unique circuit components     
        for j in 1:length(loops[i])             #Iterate components in current loop
            componentLoopDict[loops[i][j]]=push!(get(componentLoopDict, loops[i][j], []), i-1) #Forms dict with unique circuit elements
        end
    end

    for comp in keys(componentLoopDict)         #Finds circiut component parameters
        if (comp in keys(componentParamDict))   #If component already has parameter skip
            if (comp[1] in ['I', 'R', 'C', 'J'])
                push!(junctions, comp)
            end
            continue
        elseif (comp[1] == 'I')
            push!(junctions, comp)                 #Gather data about current source
            componentParamDict[comp] = 1.0
        elseif (comp[1] == 'R')                 #Gather data about resistor
            push!(junctions, comp)
            componentParamDict[comp] = 1.0
        elseif (comp[1] == 'C')                 #Gather data about capacitor
            push!(junctions, comp)
            componentParamDict[comp] = 1.0
        elseif (comp[1] == 'L')                 #Gather data about inductor
            componentParamDict[comp] = 1.0
        elseif (comp[1] == 'J')                 #Gather data about josephson junction
            push!(junctions, comp)
            componentParamDict[comp] = 1.0
        end
    end
    return componentLoopDict, componentParamDict, junctions
end

#Use existing circuit data to form k, L, σA, σB and componentPhaseDirection
function process_netlist(loops::Vector{Vector{String}}; mutualInd = Vector[])

    numLoops = size(loops)[1]
    
    componentLoopDict, componentParamDict, junctions = find_components(loops)

    #Initialise σB matrices, and componentPhaseDirection dictionary            
    σB = zeros(Int64, 0, numLoops)           
    componentPhaseDirection = Dict{String, Vector{Int64}}() 

    ### Algorithm for finding σB & σA & componentPhaseDirection
    for i in 1:length(junctions)                                #Iterate through junctions
        current_row = Vector{Int64}()
        for j in 1:numLoops                                     #Iterate through loops
            if junctions[i] in loops[j]                         #If current junction is in current loop
                junc_loops = get(componentLoopDict, junctions[i], -1)   #Array containing the loops in which the current junction is present
                loop_count = 0
                for n in 1:length(junc_loops)
                    loop_count = loop_count + junc_loops[n]     #Sum the loop number of the loops in which the current junction is present
                end
                #If the sum the loop number of the loops in which the current junction is present is greater than the current loop number 
                #the direction of the phase is positive as the component must be on the RHS or bottom of the loop --- check readme.txt
                if (loop_count/length(junc_loops) >= (j-1))     
                    push!(current_row, 1)                       #Positive θ direction
                else
                    push!(current_row, -1)                      #Negative θ direction
                end
            else
                push!(current_row, 0)                           #No θ direction as this component does not exist in loop j
            end
        end
        componentPhaseDirection[junctions[i]] = current_row     #Push current_row to componentPhaseDirection dict
        σB = [σB; current_row']                                 #Push current_row to σB matrix
    end

    #Set matrices as transpose of existing matrices
    σA = transpose(σB)
    circuit = CircuitNetlist(loops, mutualInd, componentLoopDict, componentParamDict, junctions, σB, componentPhaseDirection)
    return circuit
end
#Build the circuit based on the open file
function build_circuit(circuit::CircuitNetlist)
    #Load circuit from netlist
    loops = circuit.loops
    numLoops = size(loops)[1]
    componentLoopDict = circuit.comps
    CPD = circuit.phaseParity
    mutualInd = circuit.mutualInd
    junctions = circuit.junctions
    σA = transpose(circuit.σB)

                                                         #Array to store loops that have been built using MTK
    built_components = OrderedDict()                        #Dictionary to store components that have been built with MTK
    for j in junctions                                      #Iterate through juncrions
        println(j)                                          #Display junction name
        if (j[1] == 'R')                                    #Built resistor case
            param = get(CPD, j, 0)
            new_c = "@named $j = Resistor()"
            new_c = Meta.parse(new_c)                       #Using metaprogramming to generate components with unique names and parameters
            new_c = eval(new_c)
            built_components[j] = new_c                     #Push component to built components dictionary
        elseif (j[1] == 'C')                                #Built capacitor case
            param = get(CPD, j, 0)
            new_c = "@named $j = Capacitor()"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        elseif (j[1] == 'J')                                #Built Josephson Junction case
            params = get(CPD, j, 0)
            new_c = "@named $j = JosephsonJunction()"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        elseif (j[1] == 'I')                                #Built Josephson Junction case
            params = get(CPD, j, 0)
            new_c = "@named $j = CurrentSource()"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        end
    end
    ### Algorithm for finding L
    L = zeros(Num, numLoops, numLoops)
    for j in 1:numLoops                                         #Iterate through all loops
        current_row = []
        for i in 1:numLoops                                     #Second iteration through all loops
            Lij = 0                                           #Float storing the value of the (j,i) position in matrix L
            #SELF COUPLING
            for n in loops[i]                                   #Iterate through components in loop i
                if ((n[1] == 'L'))
                    if (j-1 in get(componentLoopDict, n, -1))   #If component n is also in the loop j
                        if (n[1] == 'L')
                            eval(Meta.parse("@named " *n* " = Inductor()"))
                            built_components[n] = eval(Meta.parse(n))    
                            param = eval(Meta.parse(n*".L"))     #Inductor case for setting param
                        end
                        if (i == j)
                            Lij = Lij + param     #Adjust Lij by the value of the inductance of component n
                        else
                            Lij = Lij - param     #Adjust Lij by the value of the inductance of component n
                        end
                    end
                end
            end
            #MUTUAL COUPLING
            for n in mutualInd
                eval(Meta.parse("@named M" * string(n[1])*string(n[2]) * "= Inductor()"))
                built_components["M" *string(n[1])*string(n[2])] = eval(Meta.parse("M" * string(n[1])*string(n[2]))) 
                param = eval(Meta.parse("M" * string(n[1])*string(n[2]) * ".L"))
                if ((i != j) && (i in n) && (j in n)) #If the two currently observed loops are not the same loop and are stated as having mutual inductance
                    Lij = Lij - param              #Adjust Lij by the value of the mutual inductance
                end
            end
            push!(current_row, Lij)                      #Lij is pushed to current_row 
        end 
        L[j,:] = current_row'                                #current_row is pushed to the L matrix
    end
                                 
    already_in_loop = String[]
    eqs=Equation[]
    for i in 1:length(loops)
        connectables = []
        for component_name in keys(circuit.comps)
            if i-1 in circuit.comps[component_name]
                comp_system = built_components[component_name]
                if component_name in already_in_loop
                    push!(connectables, comp_system.out)
                else
                    push!(connectables, comp_system.in)
                    push!(already_in_loop, component_name)
                end
            end
        end
        for M in mutualInd
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
        push!(eqs, connect(connectables...))
    end

    ground_loop_connectables = []
    for (comp, vals) in circuit.comps
        if length(vals) == 1
            comp_system = built_components[comp]
            push!(ground_loop_connectables, comp_system.out)
        end
    end
    @named ground = GroundLoop()
    push!(ground_loop_connectables, ground.g)
    push!(eqs, connect(ground_loop_connectables...))

    built_components["ground"] = ground


    @named _model = System(eqs, t)
    @named model = compose(_model, built_components.vals)                    
    new_model = mtkcompile(model)
    return new_model                             
end


