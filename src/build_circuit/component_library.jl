using ModelingToolkit, Symbolics

const Φ₀ = 2.067833848e-15              #Flux quantum

@independent_variables t                            #Time variable
D = Differential(t)                     #First Differential operation
D2 = Differential(t)^2                  #Second Differential operation

@connector Loop begin
    Φ(t), [connect = Flow]
    iₘ(t)
end

@mtkmodel GroundLoop begin
    @components begin
        g = Loop()
    end
    @equations begin
        g.iₘ ~ 0
    end
end

@mtkmodel Branch begin
    @components begin
        in = Loop()
        out = Loop()
    end
    @variables begin
        θ(t)
        i(t)
    end
    @equations begin
        i ~ in.iₘ - out.iₘ #Define branch current as difference between loop currents
        0 ~ in.Φ + out.Φ #Flux flowing from left loop is equal to flux entering the right loop
        θ ~ 2*pi/Φ₀*(in.Φ)  #define branch phase as flux flowing from input loop (just sets the parity)
    end
end

@mtkmodel Inductor begin
    @components begin
        in = Loop()
        out = Loop()
    end
    @parameters begin
        L = 1.0
    end
        @variables begin
        i(t)
    end
    @equations begin
        in.Φ ~  L*i
        i ~ in.iₘ - out.iₘ 
        0 ~ in.Φ + out.Φ #Flux flowing from left loop is equal to flux entering the right loop
    end
end

@mtkmodel ExternalFlux begin
    @components begin
        in = Loop()
    end
    @parameters begin
        Φₑ = 1.0
    end
    @equations begin
        Φₑ ~ -in.Φ  
    end
end

@mtkmodel Resistor begin
    @extend Branch()
    @parameters begin
        R = 1.0 # Sets the default resistance
    end
    @equations begin
        D(θ)~i*(2*pi*R)/Φ₀  
    end
end

@mtkmodel Capacitor begin
    @extend Branch()
    @parameters begin
        C = 1.0
    end
    @equations begin
        D2(θ)~i*2*pi/(Φ₀*C)  
    end
end

@mtkmodel JosephsonJunction begin
    @extend Branch()
    @parameters begin
        C = 1.0
        R = 1.0
        I0 = 1.0
    end
    @equations begin
        D2(θ) ~ (i - I0*sin(θ) - D(θ)*Φ₀/(2*pi*R))*(2*pi)/(Φ₀*C)  
    end
end

@mtkmodel CurrentSource begin
    @extend Branch()
    @parameters begin
        I=1.0
        ω=1.0
    end
    @equations begin
        i ~ I*sin(ω*t)
    end
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
