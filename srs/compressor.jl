function compressor_exit_flow(fs_in::FlowState{T}, PR::T, η_c::T) where T<:AbstractFloat
    γ = fs_in.γ
    R = fs_in.R
    Cp = fs_in.Cp

  T_isentropic = fs_in.Temperature * PR^((γ_val - 1)/γ_val)
  
    T_exit = fs_in.Temperature + (T_isentropic - fs_in.Temperature)/η_c
    P_exit = fs_in.Pressure * PR
    ρ_exit = P_exit / (R_val * T_exit) # (ideal gas)
    V_exit = fs_in.Velocity
    Mach_exit = V_exit / sqrt(γ_val * R_val * T_exit)
    τ_o_exit = T_exit * (1 + (γ_val - 1)/2 * Mach_exit^2)
    P_o_exit = P_exit * (1 + (γ_val - 1)/2 * Mach_exit^2)^(γ_val/(γ_val - 1))
    m_dot_exit = ρ_exit * V_exit * fs_in.CrossArea

    fs_out = FlowState{T}(
        T_exit, P_exit, V_exit, fs_in.CrossArea, fs_in.Gas,
        ρ_exit, γ_val, R_val, Cp_val,
        Mach_exit, m_dot_exit, τ_o_exit, P_o_exit
    )

    return fs_out
end

#PR = Pressure Ratio
function polytropic_temperature(T_in::T, PR::T, γ::T, η_p::T) where T<:AbstractFloat
    return T_in * PR^((γ - 1)/(γ * η_p))
end
