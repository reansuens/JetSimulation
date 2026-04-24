module Jet

export FlowState, ISA_Temperature, ISA_Pressure, static_from_stagnation,
       diffuser_exit_flow, compressor, combustor, turbine, nozzle,
       entropy_stations, station_static, inlet_state, inlet_state_with_shock,
       normal_shock, GAMMA, R, AIR, CP

const P_AMBIENT = 101325.0
const T_AMBIENT = 288.15
const g0 = 9.80665
const R  = 287.05
const GAMMA = 1.4
const CP = 1004.5
const AIR = (gamma = GAMMA, cp = CP)

# ISA model: valid from sea level to 32 km.
function ISA_Temperature(h::Real)
    h = Float64(h)
    h < 0 && error("Altitude must be >= 0 m")
    lapse1 = -0.0065
    if h <= 11000
        return T_AMBIENT + lapse1*h
    elseif h <= 20000
        return T_AMBIENT + lapse1*11000
    elseif h <= 32000
        return (T_AMBIENT + lapse1*11000) + 0.001*(h - 20000)
    else
        error("ISA model is valid only up to 32 km")
    end
end

function ISA_Pressure(h::Real)
    h = Float64(h)
    h < 0 && error("Altitude must be >= 0 m")
    lapse1 = -0.0065
    if h <= 11000
        T = ISA_Temperature(h)
        return P_AMBIENT * (T/T_AMBIENT)^(-g0/(lapse1*R))
    elseif h <= 20000
        T11 = ISA_Temperature(11000)
        P11 = P_AMBIENT * (T11/T_AMBIENT)^(-g0/(lapse1*R))
        return P11 * exp(-g0*(h - 11000)/(R*T11))
    elseif h <= 32000
        T11 = ISA_Temperature(11000)
        P11 = P_AMBIENT * (T11/T_AMBIENT)^(-g0/(lapse1*R))
        P20 = P11 * exp(-g0*(20000 - 11000)/(R*T11))
        T20 = ISA_Temperature(20000)
        T = ISA_Temperature(h)
        lapse2 = 0.001
        return P20 * (T/T20)^(-g0/(lapse2*R))
    else
        error("ISA model is valid only up to 32 km")
    end
end

struct FlowState
    Tt::Float64
    Pt::Float64
    V::Float64
    A::Float64
    gamma::Float64
    cp::Float64
    mdot::Float64
end

function static_from_stagnation(Tt::Real, Pt::Real, M::Real, gamma::Real = GAMMA)
    tau = 1 + (gamma - 1)/2 * M^2
    T = Tt / tau
    P = Pt / tau^(gamma/(gamma - 1))
    return T, P
end

function mach_from_TtV(Tt::Real, V::Real, gamma::Real = GAMMA, cp::Real = CP)
    T = Tt - V^2/(2cp)
    T <= 1 && return NaN
    return V / sqrt(gamma*R*T)
end

function station_static(fs::FlowState)
    T = fs.Tt - fs.V^2/(2fs.cp)
    M = mach_from_TtV(fs.Tt, fs.V, fs.gamma, fs.cp)
    if !isfinite(M) || M < 0
        return (T = T, P = NaN, M = M, rho = NaN)
    end
    P = fs.Pt / (1 + (fs.gamma - 1)/2*M^2)^(fs.gamma/(fs.gamma - 1))
    rho = P/(R*T)
    return (T = T, P = P, M = M, rho = rho)
end

function flow_from_mach(Tt, Pt, M, A, mdot, gamma = GAMMA, cp = CP)
    T, P = static_from_stagnation(Tt, Pt, M, gamma)
    V = M * sqrt(gamma*R*T)
    return FlowState(Float64(Tt), Float64(Pt), V, Float64(A), Float64(gamma), Float64(cp), Float64(mdot))
end

function inlet_state(Tamb, Pamb, M0, A0, mdot)
    Tt0 = Tamb * (1 + (GAMMA - 1)/2*M0^2)
    Pt0 = Pamb * (1 + (GAMMA - 1)/2*M0^2)^(GAMMA/(GAMMA - 1))
    V0 = M0 * sqrt(GAMMA*R*Tamb)
    FlowState(Tt0, Pt0, V0, A0, GAMMA, CP, mdot)
end

function normal_shock(M1::Real, gamma::Real = GAMMA)
    M1 = Float64(M1)
    M1 <= 1 && return (M2 = M1, pt_ratio = 1.0, p_ratio = 1.0, t_ratio = 1.0)
    g = gamma
    M2sq = (1 + (g - 1)/2*M1^2)/(g*M1^2 - (g - 1)/2)
    p2p1 = 1 + 2g/(g + 1)*(M1^2 - 1)
    rho2rho1 = ((g + 1)*M1^2)/(2 + (g - 1)*M1^2)
    t2t1 = p2p1/rho2rho1
    # Total pressure ratio across a normal shock, Pt2/Pt1.
    pt2pt1 = p2p1 * ((1 + (g - 1)/2*M2sq)^(g/(g - 1))) / ((1 + (g - 1)/2*M1^2)^(g/(g - 1)))
    return (M2 = sqrt(M2sq), pt_ratio = pt2pt1, p_ratio = p2p1, t_ratio = t2t1)
end

function inlet_state_with_shock(Tamb, Pamb, M0, A0, mdot; shock_on::Bool = true)
    fs_free = inlet_state(Tamb, Pamb, M0, A0, mdot)
    if shock_on && M0 > 1.0
        sh = normal_shock(M0, GAMMA)
        # Normal shock is adiabatic, so total temperature is unchanged and total pressure drops.
        fs_post = flow_from_mach(fs_free.Tt, fs_free.Pt*sh.pt_ratio, sh.M2, A0, mdot, GAMMA, CP)
        return fs_free, fs_post, sh
    end
    return fs_free, fs_free, (M2 = M0, pt_ratio = 1.0, p_ratio = 1.0, t_ratio = 1.0)
end

function diffuser_exit_flow(fs0::FlowState, A2::Float64; eta_d::Float64 = 0.98, M2::Float64 = 0.45)
    eta_d = clamp(Float64(eta_d), 0.01, 1.0)
    Tt2 = fs0.Tt
    Pt2 = fs0.Pt * eta_d
    return flow_from_mach(Tt2, Pt2, M2, A2, fs0.mdot, fs0.gamma, fs0.cp)
end

function compressor(fs::FlowState, PR::Float64, eta_c::Float64; M3::Float64 = 0.25)
    gamma, cp = fs.gamma, fs.cp
    PR <= 1 && error("Compressor pressure ratio must be > 1")
    eta_c <= 0 && error("Compressor efficiency must be positive")
    Tt3s = fs.Tt * PR^((gamma - 1)/gamma)
    Tt3 = fs.Tt + (Tt3s - fs.Tt)/eta_c
    Pt3 = fs.Pt * PR
    return flow_from_mach(Tt3, Pt3, M3, fs.A, fs.mdot, gamma, cp)
end

function combustor(fs::FlowState, Tt4_target::Float64, LHV::Float64; eta_b::Float64 = 0.99, pressure_loss::Float64 = 0.05, M4::Float64 = 0.20)
    LHV <= 0 && error("Fuel LHV must be positive")
    eta_b <= 0 && error("Burner efficiency must be positive")
    if Tt4_target <= fs.Tt
        f = 0.0
        Tt4 = fs.Tt
    else
        denom = eta_b*LHV - fs.cp*Tt4_target
        denom <= 0 && error("Requested burner temperature is physically impossible for this LHV/cp model")
        f = fs.cp*(Tt4_target - fs.Tt)/denom
        Tt4 = Tt4_target
    end
    Pt4 = fs.Pt * (1 - pressure_loss)
    mdot4 = fs.mdot * (1 + f)
    return flow_from_mach(Tt4, Pt4, M4, fs.A, mdot4, fs.gamma, fs.cp), f
end

function turbine(fs::FlowState, compressor_specific_work::Float64, eta_t::Float64, f::Float64; M5::Float64 = 0.35)
    gamma, cp = fs.gamma, fs.cp
    eta_t <= 0 && error("Turbine efficiency must be positive")
    deltaT_actual = compressor_specific_work / ((1 + f)*cp)
    Tt5 = fs.Tt - deltaT_actual
    Tt5 <= 300 && error("Turbine exit temperature became nonphysical; reduce PR or increase burner temperature")
    deltaT_isentropic = deltaT_actual / eta_t
    Tt5s = fs.Tt - deltaT_isentropic
    Pt5 = fs.Pt * (Tt5s/fs.Tt)^(gamma/(gamma - 1))
    return flow_from_mach(Tt5, Pt5, M5, fs.A, fs.mdot, gamma, cp)
end

function nozzle(fs::FlowState, Pamb::Float64, Ae::Float64; eta_n::Float64 = 0.97)
    gamma, cp = fs.gamma, fs.cp
    eta_n = clamp(Float64(eta_n), 0.01, 1.0)
    Ae <= 0 && error("Nozzle area must be positive")
    fs.Pt <= Pamb && error("Nozzle total pressure is below/near ambient pressure; no useful expansion")
    Pcrit = fs.Pt * (2/(gamma + 1))^(gamma/(gamma - 1))
    if Pcrit > Pamb
        Me = 1.0
        Pe = Pcrit
    else
        Pe = Pamb
        Me = sqrt(max(0, 2/(gamma - 1)*((fs.Pt/Pe)^((gamma - 1)/gamma) - 1)))
    end
    Te = fs.Tt/(1 + (gamma - 1)/2*Me^2)
    Ve_ideal = Me * sqrt(gamma*R*Te)
    Ve = sqrt(eta_n) * Ve_ideal
    return FlowState(fs.Tt, fs.Pt, Ve, Ae, gamma, cp, fs.mdot), Pe
end

function entropy(T::Real, P::Real; Tref = T_AMBIENT, Pref = P_AMBIENT)
    CP*log(T/Tref) - R*log(P/Pref)
end

function entropy_stations(states)
    vals = map(states) do fs
        st = station_static(fs)
        if st.T <= 0 || !isfinite(st.P) || st.P <= 0
            return NaN
        end
        entropy(st.T, st.P)
    end
    s0 = first(filter(isfinite, vals))
    return vals .- s0
end

end
