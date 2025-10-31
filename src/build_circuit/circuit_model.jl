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
            if (comp[1] in ['V', 'R', 'C', 'J'])
                push!(junctions, comp)
            end
            continue
        elseif (comp[1] == 'V')                 #Gather data about voltage source 
            push!(junctions, comp)
            componentParamDict[comp] = 1.0
        elseif (comp[1] == 'I')                 #Gather data about current source
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
function process_netlist(loops::Vector{Vector{String}},
    mutualInd = [])

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

"""Example Usage
loops = [
["J1","L1"],
["L2","C1"],
["C1","R1"],
["R1","I1"]
]

coupling = [(1,2)]

circuit = process_netlist(loops, coupling)
"""



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


    eqs = Equation[]                                        #Array to store equations
    built_loops = []
                                                         #Array to store loops that have been built using MTK
    for i in 1:numLoops                                     #Iterate through all loops
        println("loop $(i)")                              #Display loop name
        current_loop = false
        c_source = ""
        for comp in loops[i]                                #Determine if a loop is a curernt source loop
            if startswith(comp, 'I')
                current_loop = true
                c_source = comp
            end
        end
        if (current_loop)                                   #Build a current source loop
            new_l = "@named loop$(i) = build_current_source_loop()"
            new_l = Meta.parse(new_l)
            new_l = eval(new_l)
        else                                                #Build a normal loop (no current source)
            new_l = "@named loop$(i) = build_loop()"
            new_l = Meta.parse(new_l)                       #Using metaprogramming to generate loops with unique names and parameters
            new_l = eval(new_l)
        end
        push!(built_loops, new_l)                           #Push built loop to built_loops array
    end

    built_components = OrderedDict()                               #Dictionary to store components that have been built with MTK
    for j in junctions                                      #Iterate through juncrions
        println(j)                                          #Display junction name
        if (j[1] == 'R')                                    #Built resistor case
            param = get(CPD, j, 0)
            new_c = "@named $j = build_resistor()"
            new_c = Meta.parse(new_c)                       #Using metaprogramming to generate components with unique names and parameters
            new_c = eval(new_c)
            built_components[j] = new_c                     #Push component to built components dictionary
        elseif (j[1] == 'C')                                #Built capacitor case
            param = get(CPD, j, 0)
            new_c = "@named $j = build_capacitor()"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        elseif (j[1] == 'V')                                #Built voltage source case
            param = get(CPD, j, 0)
            new_c = "@named $j = build_voltage_source()"
            new_c = Meta.parse(new_c)
            new_c = eval(new_c)
            built_components[j] = new_c
        elseif (j[1] == 'J')                                #Built Josephson Junction case
            params = get(CPD, j, 0)
            new_c = "@named $j = build_JJ()"
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
                if ((n[1] == 'J') || (n[1] == 'L'))
                    if (j-1 in get(componentLoopDict, n, -1))   #If component n is also in the loop j
                        if (n[1] == 'J')
                            param = eval(Meta.parse(n * ".sys.L"))    #JJ case for setting param
                        elseif (n[1] == 'L')
                            eval(Meta.parse("@named " *n* " = build_inductance()"))
                            #built_components[n] = eval(Meta.parse(n))    
                            param = eval(Meta.parse(n*".sys.L"))     #Inductor case for setting param
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
                eval(Meta.parse("@named M" * string(n[1])*string(n[2]) * "= build_inductance()"))
                built_components["M" *string(n[1])*string(n[2])] = eval(Meta.parse("M" * string(n[1])*string(n[2]))) 
                param = eval(Meta.parse("M" * string(n[1])*string(n[2]) * ".sys.L"))
                if ((i != j) && (i in n) && (j in n)) #If the two currently observed loops are not the same loop and are stated as having mutual inductance
                    Lij = Lij - param              #Adjust Lij by the value of the mutual inductance
                end
            end
            push!(current_row, Lij)                      #Lij is pushed to current_row 
        end 
        L[j,:] = current_row'                                #current_row is pushed to the L matrix
    end
                                 
    old_sys = []                                            #Array to store system states                                          
    u0 = Pair{Num, Any}[]                               #Array to store system initial condionts  (Set to 0)

    θcomponents = OrderedDict()                                        #Array to store components with phase differnece θ
    for comp in built_components                            #Iterate through components to find component  system states and intial conditons
        push!(old_sys, comp[2].sys)
        if  !(comp[1][1] in ['L', 'M'])
            push!(u0, comp[2].sys.θ=>nothing)                   #θ initialised to 0
            push!(u0, comp[2].sys.i=>0.0)                   #i initialised to 0
            θcomponents[comp[1]] = comp[2]
            if (uppercase(comp[1][1]) in ['C', 'J', 'V', 'R'])   
                push!(u0, D(comp[2].sys.θ)=>0.0)            #D(θ) initialised to 0 for capacitors, JJs and voltage sources
            end
        end
    end
    for loop in built_loops                                 #Iterate through components to find loop system states
        push!(old_sys, loop.sys)
        push!(u0, loop.sys.iₘ=>nothing)                 #iₘ initialised to 0
    end
    sys = Vector{ODESystem}(old_sys)                        #Convert system states array to an ODESystem vector form

    #Functions from model_builder.jl to form appropriate equations
    
    add_loops!(eqs, built_loops, σA, θcomponents, L)
    current_flow(eqs, CPD, built_loops, θcomponents)
    
    
    

    @named _model  =  ODESystem(eqs, t)                     #Create an ODESystem with the existing equations

    @named model = compose(_model, sys)                    
    new_model = mtkcompile(model)                #structural_simplify Algorithm to improve performance
    return new_model, u0
                                 #Return structuraly simplified model and initial conditions
end