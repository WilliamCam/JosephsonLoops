using ModelingToolkit, Symbolics

const Φ₀ = 2.067833848e-15              #Flux quantum

@independent_variables t                            #Time variable
D = Differential(t)                     #First Differential operation
D2 = Differential(t)^2                  #Second Differential operation

function build_component(;name)         #General component with phase and current states
    sts = @variables θ(t)=1.0 i(t)=1.0  
    ODESystem(Equation[], t, sts, []; name=name)
end

struct Component
    sys::ODESystem
end

struct Inductor
    sys::ODESystem
end

function build_inductance(;name, L = 1.0)
    ps = @parameters L = L
    sys = ODESystem(Equation[], t, [], ps; name = name)
    Inductor(sys)
end 

function build_resistor(;name, R = 1.0) #Builds ODESystem for resistor using Component structure
    @named component = build_component()    #Create new component
    @unpack θ, i = component 
    ps = @parameters R=R               #Extract variables from component
    eqs = [
            D(θ)~i*(2*pi*R)/Φ₀              #Differential equation defining θ and R relationship
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end

function build_capacitor(;name, C = 1.0) #Builds ODESystem for capacitor using Component structure
    @named component = build_component()    #Create new component
    @unpack θ, i = component                #Extract variables from component
    ps = @parameters C=C
     
    eqs = [
            D2(θ)~i*2*pi/(Φ₀*C)             #Differential equation defining θ and C relationship
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end

function build_JJ(;name, I0 = 1.0, R = 1.0, C = 1.0, L = 1.0) #Builds ODESystem for Josephson Junction using Component structure
    @named component = build_component()    #Create new component
    @unpack θ, i = component                #Extract variables from component
    ps = @parameters R=R C=C L=L I0=I0
    eqs = [
            D2(θ) ~ (i - I0*sin(θ) - D(θ)*Φ₀/(2*pi*R))*(2*pi)/(Φ₀*C)    #Differential equation defining θ I0, R and C relationship
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end

function build_voltage_source(;name, V = 1.0, ω = 0.0) #Builds ODESystem for Voltage Source (AC or DC) using Component structure
    @named component = build_component()    #Create new component
    @unpack θ, i = component                #Extract variables from component
    ps = @parameters V=V ω=ω
    eqs = [
            D(θ) ~ V*cos(ω*t)*2*pi/Φ₀      #Differential equation defining θ, ω and V relationship
          ]
    sys = extend(ODESystem(eqs, t, [], ps; name=name), component)
    Component(sys)
end

struct Loop
    sys::ODESystem
end

function build_loop(;name, Φₑ = 0.0) #Builds ODESystem for general loop using Loop structure
    sts = @variables iₘ(t) = 0.0
    ps = @parameters Φₑ=Φₑ
    sys = ODESystem(Equation[], t, sts, ps, name=name)
    Loop(sys)
end

struct CurrentSourceLoop
    sys::ODESystem
end

function build_current_source_loop(;name, I = 1.0, ω = 0.0, Φₑ = 0.0) #Builds ODESystem for Current Source (AC or DC) loop using Loop structure
    @named loop = build_loop(Φₑ = Φₑ)              #Create new loop
    @unpack sys = loop                      #Extract variables from loop
    @unpack iₘ, Φₑ = sys                   
    ps = @parameters I=I ω=ω
    eqs = [
            0 ~ iₘ - I*cos(ω*t)             #Equation defining relationship between θ, ω and I
          ]
    sys = extend(ODESystem(eqs, t, [], ps, name=name), loop.sys)
    CurrentSourceLoop(sys)
end

function add_loops!(eqs::Vector{Equation}, loops, σ, cs, L)
    for i in 1:length(loops)
        if occursin("CurrentSourceLoop", string(loops[i])) == false                    #Check that loop is not a current source loop
            push!(eqs,0 ~ loops[i].sys.Φₑ - dot([Φ₀/(2*pi)*c.sys.θ for (n, c) in cs], σ[i,:])-  dot(L[:,i], [l.sys.iₘ for l in loops])) #Find Φₗ = L . iₘ  for each loop        
        end
    end
end


#Function to form  i = σB . iₘ for a loop
function current_flow(eqs::Vector{Equation}, componentPhaseDirection, built_loops, built_components)
    for (comp, phase_array) in componentPhaseDirection                  #Loop through each junction where "comp" is the junction name and "phase_array" is the row of σB corresponding to the current junction
        σB_im = 0                                                       #Temp variable to store σB * im for the current junction
        for i in 1:length(phase_array)
            σB_im = σB_im + phase_array[i]*built_loops[i].sys.iₘ        #Adds direction of im from each loop through the junction
        end
        push!(eqs, 0 ~ σB_im - get(built_components, comp, -1).sys.i)   #Find current of junction i = σB . iₘ (direction and strength of current from each loop through junction) for each junction
    end
end