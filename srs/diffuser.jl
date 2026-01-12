function diffuser_exit_flow(fs_in::FlowState{T}, A_exit::T) where T<:AbstractFloat
    # Mass flow conservation: m_dot = ρ * V * A
    ρ_exit = fs_in.m_dot / (fs_in.CrossArea * fs_in.Velocity)  # temporary exit density placeholder

    V_exit = fs_in.m_dot / (ρ_exit * A_exit)

    T_exit = fs_in.Temperature * (1 - (fs_in.Velocity^2 - V_exit^2)/(2 * fs_in.Cp))

    P_exit = ρ_exit * fs_in.R * T_exit

    Mach_exit = V_exit / sqrt(fs_in.γ * fs_in.R * T_exit)
    τ_o_exit = T_exit * (1 + (fs_in.γ - 1)/2 * Mach_exit^2)
    P_o_exit = P_exit * (1 + (fs_in.γ - 1)/2 * Mach_exit^2)^(fs_in.γ/(fs_in.γ - 1))
    m_dot_exit = ρ_exit * V_exit * A_exit

    fs_out = FlowState{T}(
        T_exit, P_exit, V_exit, A_exit, fs_in.Gas,
        ρ_exit, fs_in.γ, fs_in.R, fs_in.Cp,
        Mach_exit, m_dot_exit, τ_o_exit, P_o_exit
    )
    return fs_out
end
