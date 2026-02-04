using GLMakie
include("src/Jet.jl")
using .Jet
const T_MAX_TURBINE = 1600.0
const LHV_FUEL      = 43e6
const η_COMP        = 0.85
const η_TURB        = 0.90
const PR_COMP       = 12.0

const GAMMA = Jet.GAMMA
const R     = Jet.R
const AIR   = Jet.AIR
function compute_cycle(h_val, M_val, Ai_val, Ad_val, An_val)
    h  = Float64(h_val)
    M  = Float64(M_val)
    Ai = Float64(Ai_val)
    Ad = Float64(Ad_val)
    An = Float64(An_val)

    T_amb = Jet.ISA_Temperature(h)
    P_amb = Jet.ISA_Pressure(h)
    a_amb = sqrt(GAMMA * R * T_amb)
    V0    = M * a_amb

    Tt0 = T_amb * (1 + (GAMMA-1)/2 * M^2)
    Pt0 = P_amb * (1 + (GAMMA-1)/2 * M^2)^(GAMMA/(GAMMA-1))
    fs0 = Jet.FlowState(Tt0, Pt0, V0, 1.0, AIR)

    fs2 = Jet.diffuser_exit_flow(fs0, Ad)
    fs3 = Jet.compressor(fs2, PR_COMP, η_COMP)
    w_comp = fs3.cp * (fs3.Tt - fs2.Tt)
    fs4, f = Jet.combustor(fs3, T_MAX_TURBINE, LHV_FUEL)
    fs5 = Jet.turbine(fs4, w_comp, η_TURB)
    fs6 = Jet.nozzle(fs5, Jet.ISA_Pressure(h), An)

    stations = [fs0, fs2, fs3, fs4, fs5, fs6]
    (stations = stations, PR = PR_COMP, f = f)
end

fig = Figure(size=(960, 600), fontsize=11, backgroundcolor=:gray95)

# Main grid: 1 row, 2 columns
main_grid = fig[1, 1] = GridLayout(colgap=15)

# Left panel: controls + performance
left_panel = GridLayout(rowgap=8, colgap=6)
main_grid[1, 1] = left_panel
colsize!(main_grid, 1, Fixed(280))

# Right panel: plots
right_panel = GridLayout(rowgap=10, colgap=10)
main_grid[1, 2] = right_panel

Label(left_panel[1, 1:2], "CONTROLS", fontsize=13, font=:bold, color=:dodgerblue)

param_names  = ["Alt [m]", "Mach", "A_in [m²]", "A_dif [m²]", "A_noz [m²]"]
param_ranges = [0:500:32000, 0.3:0.025:2.5, 0.2:0.05:4.0, 0.4:0.05:2.5, 0.1:0.02:2.5]
param_starts = [0.0, 0.7, 0.5, 0.8, 0.3]

sliders = []
for i in 1:5
    Label(left_panel[i+1, 1], param_names[i], width=70, halign=:right, fontsize=10)
    sl = Slider(left_panel[i+1, 2], range=param_ranges[i], startvalue=param_starts[i], width=140)
    Label(left_panel[i+1, 3], lift(x -> string(round(x,digits=2)), sl.value), width=45, fontsize=10)
    push!(sliders, sl)
end
sl_h, sl_M, sl_Ai, sl_Ad, sl_An = sliders
Label(left_panel[7, 1:3], "PERFORMANCE", fontsize=13, font=:bold, color=:darkgreen)

metrics_layout = GridLayout(rowgap=6, colgap=6, tellheight=false)
left_panel[8:10, 1:3] = metrics_layout

cycle_obs = lift(sl_h.value, sl_M.value, sl_Ai.value, sl_Ad.value, sl_An.value) do h, m, ai, ad, an
    compute_cycle(h, m, ai, ad, an)
end

thrust_obs = lift(cycle_obs) do res
    fs0 = res.stations[1]
    fs6 = res.stations[end]
    round((fs6.mdot*fs6.V - fs0.mdot*fs0.V)/1000, digits=1)
end

fuel_obs = lift(cycle_obs) do res
    round(res.f*1000, digits=2)
end

mdot_obs = lift(cycle_obs) do res
    round(res.stations[1].mdot, digits=2)
end

Label(metrics_layout[1,1], "Thrust", width=50, fontsize=10)
Label(metrics_layout[1,2], lift(t -> "$t kN", thrust_obs), fontsize=11)

Label(metrics_layout[2,1], "f/a", width=50, fontsize=10)
Label(metrics_layout[2,2], lift(f -> "$f ‰", fuel_obs), fontsize=11)

Label(metrics_layout[3,1], "mdot", width=50, fontsize=10)
Label(metrics_layout[3,2], lift(m -> "$m kg/s", mdot_obs), fontsize=11)

x_stations = [0,2,3,4,5,6]
station_names = ["Free","Diff","Comp","Comb","Turb","Noz"]
station_colors = [:steelblue, :lightblue, :dodgerblue, :orange, :red, :darkred]

ax1 = Axis(right_panel[1,1], title="Total Pressure", xlabel="Station", ylabel="Pt [kPa]",
    xticks=(x_stations, station_names), xticklabelrotation=π/6,
    titlesize=12, xlabelsize=10, ylabelsize=10, xticklabelsize=8)
    
ax2 = Axis(right_panel[1,2], title="Static Temperature", xlabel="Station", ylabel="T [K]",
    xticks=(x_stations, station_names), xticklabelrotation=π/6,
    titlesize=12, xlabelsize=10, ylabelsize=10, xticklabelsize=8)
    
ax3 = Axis(right_panel[2,1], title="Velocity", xlabel="Station", ylabel="V [m/s]",
    xticks=(x_stations, station_names), xticklabelrotation=π/6,
    titlesize=12, xlabelsize=10, ylabelsize=10, xticklabelsize=8)
    
ax4 = Axis(right_panel[2,2], title="T-s Diagram", xlabel="s [J/kg·K]", ylabel="T [K]",
    titlesize=12, xlabelsize=10, ylabelsize=10)

function update_plots!(res)
    stations = res.stations
    Pt_vals = [s.Pt/1000 for s in stations]
    T_vals  = [s.Tt - s.V^2/(2*s.cp) for s in stations]
    V_vals  = [s.V for s in stations]
    s_vals  = Jet.entropy_stations(stations)

    empty!(ax1)
    empty!(ax2)
    empty!(ax3)
    empty!(ax4)

    # Total Pressure
    barplot!(ax1, x_stations, Pt_vals, color=station_colors, width=0.6)
    for (i,v) in enumerate(Pt_vals)
        text!(ax1, x_stations[i], v, text=string(round(v,digits=1)), 
              align=(:center,:bottom), fontsize=7, offset=(0,2))
    end

    # Static Temperature
    lines!(ax2, x_stations, T_vals, color=:orange, linewidth=2)
    scatter!(ax2, x_stations, T_vals, color=station_colors, markersize=10)
    for (i,v) in enumerate(T_vals)
        text!(ax2, x_stations[i], v, text=string(round(v,digits=0)),
              align=(:center,:bottom), fontsize=7, offset=(0,2))
    end

    # Velocity
    band!(ax3, x_stations, zeros(length(V_vals)), V_vals, color=(:limegreen, 0.3))
    lines!(ax3, x_stations, V_vals, color=:green, linewidth=2.5)
    scatter!(ax3, x_stations, V_vals, color=:darkgreen, markersize=8)
    for (i,v) in enumerate(V_vals)
        text!(ax3, x_stations[i], v, text=string(round(v,digits=0)),
              align=(:center,:bottom), fontsize=7, offset=(0,2))
    end

    # T-s Diagram
    lines!(ax4, s_vals, T_vals, color=:darkorange, linewidth=2)
    scatter!(ax4, s_vals, T_vals, color=station_colors, markersize=10)
    for i in eachindex(s_vals)
        text!(ax4, s_vals[i], T_vals[i], text=string(i-1),
              align=(:left,:bottom), fontsize=7, offset=(4,4))
    end

    # Force axis limits update
    autolimits!(ax1)
    autolimits!(ax2)
    autolimits!(ax3)
    autolimits!(ax4)
end

on(cycle_obs) do res
    update_plots!(res)
end


update_plots!(compute_cycle(0, 0.7, 0.5, 0.8, 0.3))

display(fig)

