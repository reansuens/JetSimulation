# TEMPERATURE IS IN KELVIN
# PRESSURE IS IN Pa
# Density is in Kg/m^3 
# Area is in m^2 
# Velocity is in m/s 
@enum GasType OXYGEN METHANE CARBON_DIOXIDE AIR HYDROGEN

const GAS_PROPERTIES = Dict{GasType, NamedTuple{(:Cp, :Cv, :R, :γ), NTuple{4,Float64}}}(
    OXYGEN         => (Cp =  918.0, Cv =  659.0, R =  259.0, γ = 1.393),
    METHANE        => (Cp = 2220.0, Cv = 1710.0, R =  510.0, γ = 1.298),
    CARBON_DIOXIDE => (Cp =  846.0, Cv =  657.0, R =  189.0, γ = 1.289),
    AIR            => (Cp = 1005.0, Cv =  718.0, R =  287.0, γ = 1.400),
    HYDROGEN       => (Cp =14300.0, Cv =10100.0, R = 4124.0, γ = 1.415)
)

struct FlowState{T<:AbstractFloat}
    Temperature::T
    Pressure::T
    Velocity::T
    CrossArea::T
    Gas::GasType
    Density::T
    γ::T
    R::T
    Cp::T
    Mach::T
    m_dot::T
    τ_o::T
    P_o::T
end

function FlowState(τ::T, P::T, V::T, A::T, gas::GasType) where {T<:AbstractFloat}
    properties = GAS_PROPERTIES[gas]
    Density = P / (properties.R * τ)
    a = sqrt(properties.γ * properties.R * τ)
    Mach = V / a
    m_dot = Density * V * A
    τ_o = τ * (1 + (properties.γ - 1)/2 * Mach^2)
    P_o = P * (1 + (properties.γ - 1)/2 * Mach^2)^(properties.γ/(properties.γ-1))

    fs = FlowState{T}(
        τ, P, V, A, gas,
        Density, properties.γ, properties.R, properties.Cp,
        Mach, m_dot, τ_o, P_o
    )
    return fs
end

# Helper functions
Mass_Flow_Rate(fs::FlowState) = fs.m_dot

function Mach(fs::FlowState)
    return fs.Mach
end

function Stagnation(fs::FlowState)
    return (fs.τ_o, fs.P_o)
end

function Entropy(fs::FlowState, P_ref::T, τ_ref::T) where T<:AbstractFloat
    s = fs.Cp * log(fs.Temperature / τ_ref) - fs.R  * log(fs.Pressure / P_ref)
    return s
end

