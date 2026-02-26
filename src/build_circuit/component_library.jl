using ModelingToolkit, Symbolics, ModelingToolkitStandardLibrary.Blocks

const Φ₀ = 2.067833848e-15              #Flux quantum

@independent_variables t                #Time variable
D = Differential(t)                     #First Differential operation
D2 = Differential(t)^2                  #Second Differential operation

#TODO: stable and dev component tagging to pass to circuit_model.jl
struct ComponentAPI
    accepted::Vector{Pair{String, Symbol}}
    dev::Vector{Pair{String, Symbol}}
end

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
        θ ~ 2*pi/Φ₀*(out.Φ)  #define branch phase as flux flowing from input loop (just sets the parity)
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
        out = Loop()
    end
    @parameters begin
        Φₑ = 1.0
    end
    @equations begin
        Φₑ ~ -in.Φ
        0 ~ in.Φ + out.Φ #Flux flowing from left loop is equal to flux entering  
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
        ω=1.0, [tunable=true]
    end
    @equations begin
        i ~ I*sin(ω*t)
    end
end

@mtkmodel Port begin
    @components begin
        Rsrc = Resistor()
        Isrc = CurrentSource()
        in = Loop()
        out = Loop()
    end
    @variables begin
        i(t), [irreducible=true]
        dθ(t), [irreducible=true]
    end
    @equations begin
        [
            connect(Isrc.in, Rsrc.in)
            i ~ Rsrc.i 
            dθ ~ D(Rsrc.θ)
            connect(Rsrc.out, in)
            connect(Isrc.out, out)
        ]
    end
end


#Depreciated or not in current use ---------------

@mtkmodel InputCurrentSource begin
    @extend Branch()
    @components begin
        I=RealInput()
    end
    @equations begin
        i ~ I.u
    end
end

@mtkmodel PumpTone begin
    @components begin
        output = RealOutput()
    end
    @parameters begin
        _A = 1.0
        _ω = 0.0, [tunable=true]
    end
    @equations begin
        output.u ~ _A*sin(_ω*t)
    end
end

@mtkmodel OutputFluxSense begin
    @components begin
        in = Loop()
        out = Loop()
    end
    @variables begin
        dθ(t)
    end
    @equations begin
        in.iₘ ~ out.iₘ #ideal flux sensor
        dθ ~ 2*pi/Φ₀*(out.Φ)
        dθ ~ 2*pi/Φ₀*(-in.Φ)    
    end
end

@mtkmodel AnalysisPort begin
    #TODO: force port to ground with structural params ?
    @parameters begin
        I = 1.0
        ω = 0.0
    end
    @components begin
        in = Loop()
        out = Loop()
        Rport = Resistor()
        source = InputCurrentSource()
        sensor = OutputFluxSense()
        pump_tone = PumpTone(_A=I, _ω=ω)
    end
    @equations begin
        #cirucit structure
        [
            connect(source.in, Rport.in)
            connect(Rport.out, sensor.in)
            connect(sensor.out, in)
            connect(source.out, out)
        ]
        #analysis points
        [
            connect(pump_tone.output, :drive, source.I)
        ]
    end
end
