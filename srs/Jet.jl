module Jet

using LinearAlgebra

const P_AMBIENT = 101325.0
const T_AMBIENT = 288.15
const g0 = 9.80665
const R  = 287.05         
const GAMMA = 1.4
const CP    = 1004.0       

function ISA_Temperature(h::Real)
    h = Float64(h)
    lapse = -0.0065
    if h <= 11000
        return T_AMBIENT + lapse * h
    elseif h <= 20000
        return T_AMBIENT + lapse * 11000
    elseif h <= 32000
        lapse2 = 0.001
        return (T_AMBIENT + lapse * 11000) + lapse2 * (h - 20000)
    else
        error("Altitude out of bounds (model valid up to 32 km)")
    end
end

function ISA_Pressure(h::Real)
    h = Float64(h)
    lapse = -0.0065
    if h <= 11000
        T = ISA_Temperature(h)
        return P_AMBIENT * (T / T_AMBIENT)^(-g0 / (lapse * R))
    elseif h <= 20000
        T11 = ISA_Temperature(11000)
        P11 = P_AMBIENT * (T11 / T_AMBIENT)^(-g0 / (lapse * R))
        return P11 * exp(-g0 * (h - 11000) / (R * T11))
    elseif h <= 32000
        T11 = ISA_Temperature(11000)
        P11 = P_AMBIENT * (T11 / T_AMBIENT)^(-g0 / (lapse * R))
        P20 = P11 * exp(-g0 * (20000 - 11000) / (R * T11))
        T   = ISA_Temperature(h)
        lapse2 = 0.001
        return P20 * (T / T11)^(-g0 / (lapse2 * R))
    else
        error("Altitude out of bounds (model valid up to 32 km)")
    end
end

struct FlowState
    Tt::Float64     # total temperature [K]
    Pt::Float64     # total pressure [Pa]
    V::Float64      # velocity [m/s]
    A::Float64      # area 
    gamma::Float64
    cp::Float64
    mdot::Float64   # mass flow rate [kg/s]
end

const AIR = (gamma = GAMMA, cp = CP)

function FlowState(Tt, Pt, V, A, gas = AIR)
    T  = Tt - V^2 / (2 * gas.cp)
    a = sqrt(gas.gamma * R * T)
    M = V / a

    
    P = Pt / (1 + (gas.gamma - 1)/2 * M^2)^(gas.gamma/(gas.gamma - 1))
    ρ = P / (R * T)
    mdot = ρ * V * A
    FlowState(Tt, Pt, V, A, gas.gamma, gas.cp, mdot)
end
function static_from_stagnation(Tt, Pt, M, γ = GAMMA)
    τ = 1 + (γ-1)/2 * M^2
    T = Tt / τ
    P = Pt / τ^(γ/(γ-1))
    return T, P
end

function diffuser_exit_flow(fs0::FlowState, A2::Float64, η_d::Float64 = 0.98)
    γ   = fs0.gamma
    cp  = fs0.cp
    
    T0 = fs0.Tt - fs0.V^2 / (2 * cp)
    a0 = sqrt(γ * R * T0)
    M0 = fs0.V / a0
    P0 = fs0.Pt / (1 + (γ-1)/2 * M0^2)^(γ/(γ-1))
    ρ0 = P0 / (R * T0)
    
    Tt2 = fs0.Tt
    
    recovery = 1 + η_d * (γ-1)/2 * M0^2
    Pt2 = fs0.Pt * recovery^(γ/(γ-1))
    
    M2_target = 0.45
    
    function mass_balance(M2)
        T2, P2 = static_from_stagnation(Tt2, Pt2, M2, γ)
        ρ2 = P2 / (R * T2)
        V2 = M2 * sqrt(γ * R * T2)
        return ρ2 * V2 * A2 - ρ0 * fs0.V * fs0.A
    end

    M2_low, M2_high = 0.1, 1.0
    for _ in 1:20
        M2_mid = (M2_low + M2_high) / 2
        if mass_balance(M2_low) * mass_balance(M2_mid) <= 0
            M2_high = M2_mid
        else
            M2_low = M2_mid
        end
    end
    M2 = (M2_low + M2_high) / 2
    
    T2, P2 = static_from_stagnation(Tt2, Pt2, M2, γ)
    V2 = M2 * sqrt(γ * R * T2)

    ρ2 = P2 / (R * T2)
    mdot2 = ρ2 * V2 * A2
    
    FlowState(Tt2, Pt2, V2, A2, γ, cp, mdot2)
end

function compressor(fs::FlowState, PR::Float64, η_c::Float64)
    γ, cp = fs.gamma, fs.cp
    Tt_in = fs.Tt
    Tt3s  = Tt_in * PR^((γ-1)/γ)
    ΔT    = (Tt3s - Tt_in) / η_c
    Tt3   = Tt_in + ΔT
    Pt3   = fs.Pt * PR
    
    T3, P3 = static_from_stagnation(Tt3, Pt3, 0.2, γ)  
    V3 = 0.2 * sqrt(γ * R * T3)
    FlowState(Tt3, Pt3, V3, fs.A, γ, cp, fs.mdot)
end

function combustor(fs::FlowState, Tt4_target::Float64, LHV::Float64)
    cp = fs.cp
    f  = cp * (Tt4_target - fs.Tt) / LHV          
    Pt4 = 0.96 * fs.Pt                           
    mdot_new = fs.mdot * (1 + f)
    
    T4, P4 = static_from_stagnation(Tt4_target, Pt4, 0.15, fs.gamma)
    V4 = 0.15 * sqrt(fs.gamma * R * T4)
    FlowState(Tt4_target, Pt4, V4, fs.A, fs.gamma, cp, mdot_new), f
end

function turbine(fs::FlowState, work_required::Float64, η_t::Float64)
    γ, cp = fs.gamma, fs.cp
    Tt_in = fs.Tt
    ΔT_ideal = work_required / cp
    ΔT_real  = ΔT_ideal / η_t  # Note: turbine efficiency definition
    Tt5      = Tt_in - ΔT_real
    Tt5s     = Tt_in - ΔT_ideal
    PR_t     = (Tt5s / Tt_in)^(γ/(γ-1))
    Pt5      = fs.Pt * PR_t
    # Velocity increases through turbine
    T5, P5 = static_from_stagnation(Tt5, Pt5, 0.3, γ) 
    V5 = 0.3 * sqrt(γ * R * T5)
    FlowState(Tt5, Pt5, V5, fs.A, γ, cp, fs.mdot)
end

function nozzle(fs::FlowState, P_amb::Float64, A_e::Float64)
    γ  = fs.gamma
    cp = fs.cp
    Pt = fs.Pt
    Tt = fs.Tt
    
    
    Pcrit_ratio = (2/(γ+1))^(γ/(γ-1))   # ~0.528
    
    if Pt / P_amb > 1/Pcrit_ratio
        # Choked flow
        Me = 1.0
        Pe = Pt * Pcrit_ratio
        Te = Tt * (2/(γ+1))
    else

        Pe = P_amb
        Me = sqrt(( (Pt/Pe)^((γ-1)/γ) - 1 ) * 2/(γ-1))
        Te = Tt / (1 + (γ-1)/2 * Me^2)
    end
    
    Ve = Me * sqrt(γ * R * Te)
    

    ρe = Pe / (R * Te)
    mdot_e = ρe * Ve * A_e
    

    if abs(mdot_e - fs.mdot) / fs.mdot > 0.01
        A_e_adj = A_e * (fs.mdot / mdot_e)
        mdot_e = fs.mdot
    else
        A_e_adj = A_e
    end
    
    FlowState(Tt, Pt, Ve, A_e_adj, γ, cp, mdot_e)
end

# ---------- Entropy for T-s diagram ----------
function entropy(T::Real, P::Real; Tref = T_AMBIENT, Pref = P_AMBIENT)
    CP * log(T / Tref) - R * log(P / Pref)
end

function entropy_stations(states)
    s = Float64[]
    for fs in states
        T_static = fs.Tt - fs.V^2 / (2 * fs.cp)
        if T_static <= 0
            push!(s, 0.0)
            continue
        end
        a_loc = sqrt(fs.gamma * R * T_static)
        M_loc = fs.V / a_loc
        if M_loc < 0 || !isfinite(M_loc)
            push!(s, 0.0)
            continue
        end
        P_static = fs.Pt / (1 + (fs.gamma-1)/2 * M_loc^2)^(fs.gamma/(fs.gamma-1))
        s_val = entropy(T_static, P_static)
        push!(s, s_val)
    end
    s_min = minimum(s)
    return s .- s_min
end

export FlowState, ISA_Temperature, ISA_Pressure, diffuser_exit_flow, compressor, combustor,
       turbine, nozzle, entropy_stations, GAMMA, R, AIR, CP

end 
