using JLD2, FileIO, Symbolics
#print readme file in REPL
function readme()
    print("readme.txt")
end

#Find component parameters
function find_components(numLoops, loops, componentParamDict)
    componentLoopDict = Dict()                  #Dictionary with components as keys and loops as values (used to find unique elements)
    junctions = []                              #Stores the names of the junctions
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
function process_netlist(name)
    #Open file in read mode and gather all existing data then close
    file = jldopen("$name.jld2", "a+")
    numLoops = read(file, "editing/numLoops")
    componentLoopDict = read(file, "editing/componentLoopDict")
    componentParamDict = read(file, "editing/componentParamDict")
    mutualInd = read(file, "editing/mutualInd")
    junctions = read(file, "editing/junctions")
    loops = read(file, "editing/loops")
    close(file)

    #Open file in write mode, clearing existing data
    file = jldopen("$name.jld2", "w") #Open file to write new data

    #Initialise σB matrices, and componentPhaseDirection dictionary            
    σB = zeros(Num, 0, numLoops)           
    componentPhaseDirection = Dict() 

    ### Algorithm for finding σB & σA & componentPhaseDirection
    for i in 1:length(junctions)                                #Iterate through junctions
        current_row = []
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

    #Save data to file and close
    file["editing/loops"] = loops
    file["editing/componentParamDict"] = componentParamDict
    file["editing/componentLoopDict"] = componentLoopDict
    file["editing/junctions"] = junctions
    file["editing/numLoops"] = numLoops
    file["editing/mutualInd"] = mutualInd
    file["matrices/σA"] = σA
    file["matrices/σB"] = σB
    file["matrices/componentPhaseDirection"] = componentPhaseDirection
    close(file)
end

#Creates a new .jdl2 file storing data about a circuit
function new_netlist(name)
    file = jldopen("$name.jld2", "w")

    loops = []                                  #Stores components in each loop
    mutualInd = []                              #Stores data regarding which loops are mutally coupled
    extern_flux = []                            #Stores loops external flux data

    println("Enter the number of loops in the circuit:")
    numLoops = readline()
    numLoops = parse(Int8, numLoops)            #Store number of Loops as an int8 (up to 128 loops)


    println("WARNING CURRENT SOURCES MUST BE ON THE OUTSIDE OF AN EXTERNAL LOOP")
    for i in 1:numLoops                         #Asks about circuit elements
        push!(loops, [])
        println("Enter all components in Loop $(i) one by one\n(Enter '~' when all components are listed)")
        while true
            input = readline()
            if (input == "~")  
                break
            end
            push!(loops[i], input)              #Push component into current loop
        end
    end

    componentParamDict = Dict()
    componentLoopDict, componentParamDict, junctions = find_components(numLoops, loops, componentParamDict) #Find component parameters and store in dicts

    println("Enter any mutually coupled loops:\nE.g. if loop 1 and 2 are coupled with mutual inductance M12 enter\n'1,2'\n(Enter '~' when all are listed)")
    while true
        input = readline()                      #Mutual flux input
        if (input == "~")                       #Error handling
            break
        end
        try
            currentMutual = split(input, ',')
            if (length(currentMutual) != 2)
                throw(error)
            end
            mutualTuple = (parse(Int8, currentMutual[1]), parse(Int8, currentMutual[2]))
            push!(mutualInd, mutualTuple)
        catch e
            if isa(e, ArgumentError)
                println(" --- One or more variable types are incorrect, loops must be Int")
            elseif e == error
                println(" --- Incorrect input length, enter 2 values seperated by commas ---")
            end
        end
    end

    #Sava data to file
    file["editing/loops"] = loops
    file["editing/componentParamDict"] = componentParamDict
    file["editing/componentLoopDict"] = componentLoopDict
    file["editing/junctions"] = junctions
    file["editing/numLoops"] = numLoops
    file["editing/mutualInd"] = mutualInd
    close(file)

    process_netlist(name)
end

#Edit a .jdl2 file ### Unsure if this is usable after Symbolic Variable Updates ####
function edit_netlist(name)
    file = jldopen("$name.jld2", "r+") #Open file in read mode and retrieve data then close
    numLoops = read(file, "editing/numLoops")
    componentLoopDict = read(file, "editing/componentLoopDict")
    componentParamDict = read(file, "editing/componentParamDict")
    mutualInd = read(file, "editing/mutualInd")
    junctions = read(file, "editing/junctions")
    loops = read(file, "editing/loops")
    L = read(file, "matrices/L")
    σA = read(file, "matrices/σA")
    σB = read(file, "matrices/σB")
    componentPhaseDirection = read(file, "matrices/componentPhaseDirection")
    close(file)

    file = jldopen("$name.jld2", "w") #Open file to overwrite data
    
    while true
        println("What would you like to edit?")
        println(" --- Enter P to change component parameters ---")
        println(" --- Enter L to change components in loops ---")
        println(" --- Enter M to change mutally coupled loops ---")
        println(" --- Enter K to change external flux through loops ---")
        println(" --- Enter ~ when finished editing ---")
        input = readline()
        if (input == "~")                           #Exit the edit function
            break
        end
        if (uppercase(input) == "P")                #Change component parameters
            while true
                display(componentParamDict)         #Display the current componentParamDict
                println()
                println("Enter a component followed by '>' followed by it's new parameter. E.g. to change Ra from 2 Ohm to 3 Ohm enter Ra>3")
                println("For JJ's enter comma seperated parameters e.g. J1>10e-6,27,3.24,50e-3")
                println("Enter ~ when finished editing parameters")
                input = readline()
                if (input == "~")
                    break                           #Exit edit component parameters
                end
                comp_param = split(input, '>')
                if (comp_param[1] in keys(componentParamDict))  #Check if component exists
                    if startswith(comp_param[1], 'J')           #Josephson Junction case
                        params = split(comp_param[2], ',')
                        param_a = []
                        for p in params
                            try
                                input = Meta.parse(p)
                            catch e
                                if isa(e, MethodError)
                                    input = parse(Float64, p)
                                end
                            end
                            push!(param_a, input)
                        end
                        componentParamDict[comp_param[1]] = param_a
                    elseif (startswith(comp_param[1], 'I') || startswith(comp_param[1], 'V')) #Current/Voltage source case
                        params = split(comp_param[2], ',')
                        if (length(params) == 2)
                            param_a = []
                            for p in params
                                try
                                    input = Meta.parse(p)
                                catch e
                                    if isa(e, MethodError)
                                        input = parse(Float64, p)
                                    end
                                end
                                push!(param_a, input)
                            end
                            componentParamDict[comp_param[1]] = param_a
                        else
                            try
                                input = Meta.parse(comp_param[2])
                            catch e
                                if isa(e, MethodError)
                                    input = parse(Float64, input)
                                end
                            end
                            componentParamDict[comp_param[1]] = input
                        end
                    else                                        #All other components
                        try
                            input = Meta.parse(comp_param[2])
                        catch e
                            if isa(e, MethodError)
                                input = parse(Float64, input)
                            end
                        end
                        componentParamDict[comp_param[1]] = input
                    end
                else
                    println("Component does not exist, try agian")
                end
            end
        elseif (uppercase(input) == "L")            #Change Loop structures
            while true
                for i in 1:length(loops)            #Display each loop
                    println("Loop $(i-1): $(loops[i])")
                end
                println("Which loop would you like to edit?\nEnter ~ when finished editing loops")
                input = readline()
                if (input == "~")
                    break                           #Exit edit loops
                end
                loop_num = parse(Int8, input)+1
                while true
                    println(loops[loop_num])        #Display loop components
                    println("Enter a component in the loop to remove it, enter a new component to add it to the loop\nEnter ~ when finished editing loop $(loop_num-1)")
                    input = readline()
                    if (input == "~")               #Exit edit current loop
                        break
                    elseif (input in loops[loop_num])   #If input exists in the loop remove it
                        deleteat!(loops[loop_num], findall(x->x==input,loops[loop_num]))
                    else                                #If does not input exists in the loop add it
                        push!(loops[loop_num], input)
                    end
                end
            end
            componentLoopDict, componentParamDict, junctions = find_components(numLoops, loops, componentParamDict)
        elseif (uppercase(input) == "M")            #Change mutual inductances
            while true
                display(mutualInd)                  #Display current mutual inductances
                println()
                println("Enter existing mutually coupled loops to remove, enter new mutually coupled loops to add\nEnter ~ when finished editing mutally coupled loops")
                input = readline()
                if (input == "~")
                    break
                end
                currentMutual = split(input, ',')
                mutualTuple = (parse(Int8, currentMutual[1]), parse(Int8, currentMutual[2]))
                mutualTuple = (mutualTuple, parse(Float64, currentMutual[3]))
                if (mutualTuple in mutualInd)       #If the input exists in the mutual inductances remove it
                    deleteat!(mutualInd, findall(x->x==mutualTuple,mutualInd))
                else
                    push!(mutualInd, mutualTuple)   #If the input does not exist in the mutual inductances add it
                end
            end
        elseif (uppercase(input) == "K")            #Change external flux
            display(k)                              #Display old external flux matrix
            println()
            println("Enter the new external flux through each loop:\nE.g. if there are 3 loops (0, 1, 2) and 0.6 of the external flux passes through loop 1 and the remaining flux passes through loop 2 enter \n'0,0.6,0.4'")
            input = readline()
            k = []                                  #Clears external flux matrix
            flux = split(input, ',')
            for i in 1:length(flux)
                push!(k, parse(Float64, flux[i]))   #Writes to external flux matrix 
            end
            display(k)                              #Display new external flux matrix
            println()
        end
    end

    #Write data to file
    file["editing/loops"] = loops
    file["editing/componentParamDict"] = componentParamDict
    file["editing/componentLoopDict"] = componentLoopDict
    file["editing/junctions"] = junctions
    file["editing/numLoops"] = numLoops
    file["editing/mutualInd"] = mutualInd
    file["matrices/L"] = L
    file["matrices/σA"] = σA
    file["matrices/σB"] = σB
    file["matrices/componentPhaseDirection"] = componentPhaseDirection
    close(file)

    process_netlist(name)   #Reprocess data in case any changes have been made to circuit structure
end

#User input to enter edit or create functions
while true  
    println(" --- Enter 'E' to edit an existing netlist  --- ")
    println(" --- Enter 'N' to create new netlist  --- ")
    println(" --- Enter '~' to exit program  --- ")
    input = readline()
    if (uppercase(input) == "E")
        println(" --- Enter filename (excluding '.jld2')  --- ")
        input = readline()
        edit_netlist("$input")
    elseif (uppercase(input) == "N")
        println(" --- Enter filename (excluding '.jld2')  --- ")
        input = readline()
        new_netlist("$input")
    elseif (input == "~")
        exit()
    end
end