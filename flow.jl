@enum GasType OXYGEN METHANE CARBON_DIOXIDE AIR HYDROGEN
mutable struct FlowState{T<:Real}
    Temperature::T
    Pressure::T
    Velocity::T
    CrossArea::T
    Gas::GasType
    Density::T
end


function Specific_Heats(g::GasType)
    if g == OXYGEN
        return (918.0, 659.0)
    elseif g == METHANE
        return (2220.0, 1710.0)
    elseif g == CARBON_DIOXIDE
        return (846.0, 657.0)
    elseif g == AIR
        return (1005.0, 718.0)
    elseif g == HYDROGEN
        return (14300.0, 10100.0)
    else
        error("UNDEFINED GAS")
    end 
end

R(Cp, Cv) = Cp - Cv
γ(Cp, Cv) = Cp/Cv

function Flow_State(T, P, V, A, Gas::GasType)
    Cp, Cv = Specific_Heats(Gas)
    R_val = R(Cp, Cv)
    ρ = P/(R_val*T)
    return FlowState(T, P, V, A, Gas, ρ)
end

Mass_Flow_Rate(fs::FlowState) = fs.Density * fs.Velocity * fs.CrossArea

function Mach(fs::FlowState)
    Cp, Cv = Specific_Heats(fs.Gas)
    γ_val = γ(Cp, Cv)
    R_val = R(Cp, Cv)
    a = sqrt(γ_val*R_val*fs.Temperature)
    return fs.Velocity/a
end

function Stagnation(fs::FlowState)
    Cp, Cv = Specific_Heats(fs.Gas)
    γ_val = γ(Cp, Cv)
    M = Mach(fs)
    T0 = fs.Temperature * (1 + (γ_val-1)/2 * M^2)
    P0 = fs.Pressure * (1 + (γ_val-1)/2 * M^2)^(γ_val/(γ_val-1))
    return (T0, P0)
end


