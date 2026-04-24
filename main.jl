using GLMakie
include("srs/Jet.jl")
using .Jet

const GAMMA = Jet.GAMMA
const R = Jet.R

function compute_cycle(h_val, M_val, mdot_air_val, Ai_val, Tt4_val,
                       LHV_MJ_val, eta_diff_val, eta_comp_val, eta_burn_val,
                       eta_turb_val, eta_nozz_val, PR_comp_val, A_nozz_val,
                       shock_on_val)
    h = Float64(h_val)
    M0 = Float64(M_val)
    mdot_air = Float64(mdot_air_val)
    Ai = Float64(Ai_val)
    Tt4 = Float64(Tt4_val)
    LHV = Float64(LHV_MJ_val)*1e6
    eta_diff = Float64(eta_diff_val)
    eta_comp = Float64(eta_comp_val)
    eta_burn = Float64(eta_burn_val)
    eta_turb = Float64(eta_turb_val)
    eta_nozz = Float64(eta_nozz_val)
    PR_comp = Float64(PR_comp_val)
    A_nozz = Float64(A_nozz_val)
    shock_on = Bool(shock_on_val)

    Tamb = Jet.ISA_Temperature(h)
    Pamb = Jet.ISA_Pressure(h)
    V0 = M0*sqrt(GAMMA*R*Tamb)

    fs_free, fs_inlet, shock = Jet.inlet_state_with_shock(Tamb, Pamb, M0, Ai, mdot_air; shock_on = shock_on)
    fs2 = Jet.diffuser_exit_flow(fs_inlet, Ai; eta_d = eta_diff, M2 = 0.45)
    fs3 = Jet.compressor(fs2, PR_comp, eta_comp; M3 = 0.25)

    w_comp = fs3.cp*(fs3.Tt - fs2.Tt)
    fs4, f = Jet.combustor(fs3, Tt4, LHV; eta_b = eta_burn, pressure_loss = 0.05, M4 = 0.20)
    fs5 = Jet.turbine(fs4, w_comp, eta_turb, f; M5 = 0.35)
    fs6, Pe = Jet.nozzle(fs5, Pamb, A_nozz; eta_n = eta_nozz)

    F_net = fs6.mdot*fs6.V - fs_free.mdot*fs_free.V + (Pe - Pamb)*A_nozz
    mdot_fuel = f*mdot_air
    TSFC = F_net > 0 ? mdot_fuel/F_net*3600 : NaN

    jet_power = 0.5*(fs6.mdot*fs6.V^2 - fs_free.mdot*fs_free.V^2)
    propulsive_eff = jet_power > 0 ? (F_net*V0)/jet_power : NaN
    thermal_eff = mdot_fuel > 0 ? jet_power/(mdot_fuel*LHV) : NaN
    overall_eff = mdot_fuel > 0 ? (F_net*V0)/(mdot_fuel*LHV) : NaN

    return (
        stations = [fs_free, fs_inlet, fs2, fs3, fs4, fs5, fs6],
        f = f,
        F_net = F_net,
        TSFC = TSFC,
        eta_th = thermal_eff,
        eta_prop = propulsive_eff,
        eta_overall = overall_eff,
        mdot_fuel = mdot_fuel,
        Pamb = Pamb,
        Tamb = Tamb,
        Pe = Pe,
        V0 = V0,
        Tt4 = Tt4,
        LHV = LHV,
        A_nozz = A_nozz,
        shock = shock,
        shock_active = shock_on && M0 > 1.0,
        PR_comp = PR_comp
    )
end

safe_cycle(args...) = try
    compute_cycle(args...), "OK"
catch err
    nothing, sprint(showerror, err)
end

set_theme!(theme_dark())
fig = Figure(size = (1540, 930), fontsize = 14, backgroundcolor = RGBf(0.055, 0.065, 0.085))
root = fig[1, 1] = GridLayout(colgap = 18)
controls = root[1, 1] = GridLayout(rowgap = 8)
plots = root[1, 2] = GridLayout(rowgap = 12, colgap = 12)
colsize!(root, 1, Fixed(410))

Label(controls[1, 1:3], "TURBOJET CYCLE SIMULATOR", fontsize = 24, font = :bold, color = :cyan)
Label(controls[2, 1:3], "Flight, burner, component losses, and inlet shock controls", fontsize = 13, color = :gray80)

function add_slider!(grid, row, name, range, start; digits = 1, unit = "", label_width = 135)
    Label(grid[row, 1], name, width = label_width, halign = :right)
    sl = Slider(grid[row, 2], range = range, startvalue = start, width = 165)
    Label(grid[row, 3], lift(x -> string(round(x, digits = digits), unit), sl.value), width = 90, halign = :left)
    return sl
end

sl_h = add_slider!(controls, 4, "Altitude", 0:500:32000, 0; digits = 0, unit = " m")
sl_M = add_slider!(controls, 5, "Mach", 0.05:0.025:2.20, 0.70; digits = 2)
sl_mdot = add_slider!(controls, 6, "Air mass flow", 5:0.5:150, 30; digits = 1, unit = " kg/s")
sl_Tt4 = add_slider!(controls, 7, "Burner Tt4", 900:10:1900, 1500; digits = 0, unit = " K")

Label(controls[8, 1], "A inlet", width = 135, halign = :right)
inlet_options = [0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00, 1.20, 1.50, 2.00]
menu_Ai = Menu(controls[8, 2], options = zip(string.(inlet_options), inlet_options), default = "0.5", width = 165)
Label(controls[8, 3], " m^2", width = 90, halign = :left)

Label(controls[9, 1], "Shock model", width = 135, halign = :right)
menu_shock = Menu(controls[9, 2], options = zip(["On", "Off"], [true, false]), default = "On", width = 165)
Label(controls[9, 3], "normal shock", width = 90, halign = :left, color = :gray75)

Label(controls[11, 1:3], "ADVANCED ENGINE PARAMETERS", fontsize = 17, font = :bold, color = :lightskyblue)
sl_LHV = add_slider!(controls, 12, "Fuel LHV", 38.0:0.5:46.0, 43.0; digits = 1, unit = " MJ/kg")
sl_eta_d = add_slider!(controls, 13, "Diffuser eta", 0.80:0.005:1.00, 0.99; digits = 3)
sl_eta_c = add_slider!(controls, 14, "Compressor eta", 0.60:0.005:0.95, 0.85; digits = 3)
sl_eta_b = add_slider!(controls, 15, "Burner eta", 0.80:0.005:1.00, 0.99; digits = 3)
sl_eta_t = add_slider!(controls, 16, "Turbine eta", 0.60:0.005:0.95, 0.90; digits = 3)
sl_eta_n = add_slider!(controls, 17, "Nozzle eta", 0.70:0.005:1.00, 0.97; digits = 3)
sl_PR = add_slider!(controls, 18, "Compressor PR", 2.0:0.25:30.0, 12.0; digits = 2)
sl_A_noz = add_slider!(controls, 19, "Nozzle area", 0.05:0.01:1.00, 0.30; digits = 2, unit = " m^2")

status_label = Label(controls[21, 1:3], "OK", color = RGBAf(0.6, 1.0, 0.6, 1.0), fontsize = 11)

cycle_obs = lift(sl_h.value, sl_M.value, sl_mdot.value, menu_Ai.selection, sl_Tt4.value,
                 sl_LHV.value, sl_eta_d.value, sl_eta_c.value, sl_eta_b.value,
                 sl_eta_t.value, sl_eta_n.value, sl_PR.value, sl_A_noz.value,
                 menu_shock.selection) do h, m, md, ai, tt4, lhv, ed, ec, eb, et, en, pr, an, shock_on
    res, msg = safe_cycle(h, m, md, ai, tt4, lhv, ed, ec, eb, et, en, pr, an, shock_on)
    status_label.text[] = "Status: " * msg
    status_label.color[] = res === nothing ? RGBAf(1.0, 0.65, 0.0, 1.0) : RGBAf(0.6, 1.0, 0.6, 1.0)
    return res
end

metric_grid = controls[23:31, 1:3] = GridLayout(rowgap = 6)
Label(metric_grid[1, 1:2], "LIVE PERFORMANCE", fontsize = 16, font = :bold, color = :lightgreen)

fmt(f) = lift(cycle_obs) do res
    res === nothing ? nothing : f(res)
end

function metric!(row, title, obs, unit = "")
    Label(metric_grid[row, 1], title, width = 135, halign = :right, color = :gray80)
    Label(metric_grid[row, 2], lift(x -> x === nothing ? "--" : string(x, unit), obs), fontsize = 17, font = :bold, halign = :left)
end

metric!(2, "Net thrust", fmt(r -> round(r.F_net/1000, digits = 2)), " kN")
metric!(3, "TSFC", fmt(r -> isfinite(r.TSFC) ? round(r.TSFC, digits = 4) : NaN), " kg/N h")
metric!(4, "Fuel-air", fmt(r -> round(1000r.f, digits = 2)), " per mille")
metric!(5, "Fuel flow", fmt(r -> round(r.mdot_fuel, digits = 3)), " kg/s")
metric!(6, "Thermal eta", fmt(r -> round(100r.eta_th, digits = 1)), " %")
metric!(7, "Propulsive eta", fmt(r -> round(100r.eta_prop, digits = 1)), " %")
metric!(8, "Shock Pt loss", fmt(r -> round(100*(1 - r.shock.pt_ratio), digits = 1)), " %")
metric!(9, "Nozzle Pe/Pamb", fmt(r -> round(r.Pe/r.Pamb, digits = 3)), "")


x = [0, 1, 2, 3, 4, 5, 6]
names = ["Free", "Shock/In", "Diff", "Comp", "Burn", "Turb", "Nozz"]
colors = [:deepskyblue, :lightskyblue, :skyblue, :dodgerblue, :orange, :tomato, :red]

axP = Axis(plots[1, 1], title = "Total Pressure", xlabel = "Station", ylabel = "Pt [kPa]", xticks = (x, names), xticklabelrotation = pi/7)
axT = Axis(plots[1, 2], title = "Total and Static Temperature", xlabel = "Station", ylabel = "Temperature [K]", xticks = (x, names), xticklabelrotation = pi/7)
axV = Axis(plots[2, 1], title = "Velocity Distribution", xlabel = "Station", ylabel = "V [m/s]", xticks = (x, names), xticklabelrotation = pi/7)
axS = Axis(plots[2, 2], title = "T-s Diagram", xlabel = "Relative entropy [J/kg K]", ylabel = "Static T [K]")
axM = Axis(plots[3, 1], title = "Mach Number", xlabel = "Station", ylabel = "M [-]", xticks = (x, names), xticklabelrotation = pi/7)
axF = Axis(plots[3, 2], title = "Energy Summary", xlabel = "Quantity", ylabel = "Power [MW]", xticks = (1:3, ["Fuel", "Jet", "Useful"]))

function annotate_series!(ax, xs, ys; digits = 0)
    for (xi, yi) in zip(xs, ys)
        isfinite(yi) || continue
        text!(ax, xi, yi, text = string(round(yi, digits = digits)), align = (:center, :bottom), fontsize = 10)
    end
end

function update_plots!(res)
    res === nothing && return
    for ax in [axP, axT, axV, axS, axM, axF]
        empty!(ax)
    end

    st = res.stations
    static = Jet.station_static.(st)
    Pt = [s.Pt/1000 for s in st]
    Tt = [s.Tt for s in st]
    Ts = [q.T for q in static]
    V = [s.V for s in st]
    M = [q.M for q in static]
    srel = Jet.entropy_stations(st)

    barplot!(axP, x, Pt, color = colors, width = 0.65)
    annotate_series!(axP, x, Pt; digits = 0)

    lines!(axT, x, Tt, linewidth = 3, color = :orange, label = "Total T")
    scatter!(axT, x, Tt, markersize = 11, color = colors)
    lines!(axT, x, Ts, linewidth = 2, linestyle = :dash, color = :cyan, label = "Static T")
    axislegend(axT, position = :lt)

    band!(axV, x, zeros(length(V)), V, color = (:green, 0.18))
    lines!(axV, x, V, linewidth = 3, color = :lightgreen)
    scatter!(axV, x, V, markersize = 10, color = colors)
    annotate_series!(axV, x, V; digits = 0)

    lines!(axS, srel, Ts, linewidth = 3, color = :gold)
    scatter!(axS, srel, Ts, markersize = 12, color = colors)
    for i in eachindex(srel)
        isfinite(srel[i]) && isfinite(Ts[i]) || continue
        text!(axS, srel[i], Ts[i], text = names[i], fontsize = 10, align = (:left, :bottom))
    end

    lines!(axM, x, M, linewidth = 3, color = :violet)
    scatter!(axM, x, M, markersize = 10, color = colors)
    hlines!(axM, [1.0], color = (:white, 0.45), linestyle = :dash)
    annotate_series!(axM, x, M; digits = 2)

    fuel_power = res.mdot_fuel*res.LHV/1e6
    jet_power = 0.5*(st[end].mdot*st[end].V^2 - st[1].mdot*st[1].V^2)/1e6
    useful_power = res.F_net*res.V0/1e6
    vals = [fuel_power, jet_power, useful_power]
    barplot!(axF, 1:3, vals, color = [:orange, :cyan, :lightgreen], width = 0.65)
    annotate_series!(axF, 1:3, vals; digits = 2)

    foreach(autolimits!, [axP, axT, axV, axS, axM, axF])
end

on(cycle_obs) do res
    update_plots!(res)
end

initial, _ = safe_cycle(0, 0.70, 30.0, 0.50, 1500.0, 43.0, 0.99, 0.85, 0.99, 0.90, 0.97, 12.0, 0.30, true)
update_plots!(initial)

display(fig)
