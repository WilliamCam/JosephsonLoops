
#Display function uses
function help()
    println()
    println("Use 'help()' to display this")
    println()
    println("Use 'open_file(name)' to open an existing .jdl2 file to extract circuit data")
    println()
    println("Use 'symbolic_assign()' to input numerical values for variables that are still in symbolic form")
    println("   - Symbolically assigned variables are temporary only, upon using net_edit() to change parameters symbolic variables return")
    println()
    println("Use 'net_edit(name)' to edit a netlist")
    println("   - After using net_edit(), open_file() must be called again to load correct data")
    println("   - Component parameters can also be edited by entering 'CPD[component] = new_parameter', (THIS IS TEMPORARY)")
    println()
    println("Use 'build(Φext)' to build the circuit based on the open file where Φext is the total external flux through the circuit")
    println()
    println("Use 'solve_init(old_u0, t_end)' to solve initial conditions for the model")
    println()
    println("Use 'sim_solve(u0, t_init, t_end)' to solve for a given timespan")
    println("   - A specific solver algorithm can be used for sim_solve() by adding solver=ROS3P() as a parameter for example")
    println()
    println("Use 'single_plot(comp, param)' to plot the voltage or current through a component")
    println("   - To plot voltage enter 'V' as the param")
    println("   - To plot current enter 'i' as the param")
end

#Open an existing .jdl2 file to extract circuit data
function open_file(name)
    f_name = name
    try
        if endswith(f_name, ".jld2")
            file = jldopen(f_name, "r")
        else 
            file = jldopen("$f_name.jld2", "r")
        end
    catch e
        if isa(e, SystemError)
            println("File does not exist")
        end
        #break?
    end
    if endswith(f_name, ".jld2")
        file = jldopen(f_name, "r")
    else 
        file = jldopen("$f_name.jld2", "r")
    end
    #file = jldopen("$f_name.jld2", "r")
    global numLoops = read(file, "editing/numLoops")
    global componentLoopDict = read(file, "editing/componentLoopDict")
    global CPD = read(file, "editing/componentParamDict")
    global mutualInd = read(file, "editing/mutualInd")
    global junctions = read(file, "editing/junctions")
    global loops = read(file, "editing/loops")
    global σA = read(file, "matrices/σA")
    #σB = read(file, "matrices/σB")
    global componentPhaseDirection = read(file, "matrices/componentPhaseDirection")
    close(file)
    print("File loaded \n circuit parameters:")
    for param in keys(CPD)
        print(string(param) * "\n")
    end

end

#Build the circuit based on the open file
function build_circuit(; dae_system = false)

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
    current_flow(eqs, componentPhaseDirection, built_loops, θcomponents)
    
    
    

    @named _model  =  ODESystem(eqs, t)                     #Create an ODESystem with the existing equations

    @named model = compose(_model, sys)                     #Compose the existing ODESystem with the system states vector
    if dae_system == true
        new_model = mtkcompile(dae_index_lowering(model))
    else
        new_model = mtkcompile(model)                #structural_simplify Algorithm to improve performance
    end

    return new_model, u0
                                 #Return structuraly simplified model and initial conditions
end