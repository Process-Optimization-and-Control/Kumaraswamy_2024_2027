using LaTeXStrings

####################################################
# READ THE LEDAFLOW OUTPUT DATA
####################################################
full_output_csv_path = joinpath(case_dir, "full_lf_output.csv")
df = CSV.read(full_output_csv_path, DataFrame;
    header=10,            # column names start  from row 10  
    skipto=12,            # skip the units row (line 11)  
    delim = ',',
    normalizenames =false # keep names like "Pressure@Line 1 P 500m" 
    )

time = df[!, "Time"]

choke_valve_openings = Matrix(df[!, [
    "Opening fraction@WH Choke 1", 
    "Opening fraction@WH Choke 2", 
    "Opening fraction@WH Choke 3", 
    "Opening fraction@WH Choke 4"]])

pump_speed = df[!, "Speed@Pump"]

wellhead_flows = Matrix(df[!, [
        "MFR - total liquid@Wellbore1 P MFR", 
        "MFR - total liquid@Wellbore2 P MFR", 
        "MFR - total liquid@Wellbore3 P MFR", 
        "MFR - total liquid@Wellbore4 P MFR"]]).*(-1)



####################################################
# INITIALISATION AND SETPOINT DATA   
####################################################

# LedaFlow runs a 2000 s initialisation period at the start of the simulation
# during which all manipulated inputs (choke openings, pump speed) are held
# constant -- the MPC controller is not active yet. The setpoint trajectory is
# therefore shifted by `controller_on_time` so that the first setpoint sample
# lands at the instant the controller is switched on, not at t = 0.
controller_on_time = 2000   # s, end of LedaFlow init period = controller turned on

# Setpoints are sampled every DT seconds, with sample 1 at the controller-on time
sp_time = collect(0:DT:DT*(n_samples_setpoint-1)) .+ controller_on_time
well_flow_sps = hcat(well_1_flow_sp, well_2_flow_sp, well_3_flow_sp, well_4_flow_sp)

# Phase boundaries: times where any well's setpoint changes value.
# With :steppost the new value takes effect at sp_time[k], so that is
# where the boundary is drawn.
n_sp = size(well_flow_sps, 1)
change_ks = [k for k in 2:n_sp if any(well_flow_sps[k, :] .!= well_flow_sps[k-1, :])]
phase_bounds = sp_time[change_ks]
phase_edges = [controller_on_time; phase_bounds; time[end]]


function annotate_sp_phases(plt_annotate; fs = 16, y_offset=-0.1)
    # Plt - Initialization and Setpoint Phase Annotation
    # `fs` is the phase-label font size: annotations are not clipped to the
    # plot frame, so narrower panels (e.g. in multi-panel layouts) need a
    # smaller font to keep the text inside.

    # Mark where the controller is switched on (end of the LedaFlow init period).
    # A labelled, dotted vline gives a dotted sample in the legend, so the legend
    # entry itself is shown as a dotted line.
    vline!(plt_annotate, [controller_on_time],
    linestyle = :dot,
    color = :black,
    lw = 2,
    label = "Controller on")

    #vline!(plt9, phase_bounds, linestyle = :dash, color = :grey, lw = 1.5, label = "")

    # Label each phase midway between its boundaries, near the top of the axis
    ylo, yhi = ylims(plt_annotate)
    for i in 1:length(phase_edges)-1
        annotate!(plt_annotate, (phase_edges[i] + phase_edges[i+1]) / 2,
                yhi + y_offset * (yhi - ylo),
                # \mathbf renders in Computer Modern bold; \! removes the
                # math-mode space GR inserts between the letter and the digit.
                text(L"\mathbf{P\!%$i}", fs, "Computer Modern", :black))
    end

    for i in 2:2:length(phase_edges)-1
        vspan!(plt_annotate, [phase_edges[i], phase_edges[i+1]],
            color = :grey10, alpha = 0.15, label = "")
    end
end


####################################################
# PLOT DATA   
####################################################

# Default Plot Settings
tick_positions = 0:7200:time[end]  # ticks every 2 hours
tick_labels = string.(round.(Int, tick_positions ./ 3600))

default(
    legendfont=(16, "Computer Modern"), 
    tickfont=(16, "Computer Modern", :bold), 
    guidefont=(16, "Computer Modern", :bold),
    lw=3,
    yformatter=:plain,
    xformatter=:plain,
    legend=:outerright,
    size=(700, 500),
    grid=false,
    xticks=(tick_positions, tick_labels),
    xlabel="Time (hr)",
    dpi=750
)

####################################################
# MVs   
####################################################

# Choke Valve Openings Plot
plt1 = plot(time, choke_valve_openings, 
label = ["Well 1" "Well 2" "Well 3" "Well 4"], 
ylabel = "Valve Opening",
legend=:bottomright
)

# Pump Speed Plot 
plt2 = plot(time, pump_speed,
ylabel = "Pump Speed (rpm)",
legend = nothing,
)

plot!(plt2, [controller_on_time, time[end]], fill(speed_pump_max, 2),
                linestyle = :dash,
                color = :red,
                label = "High Limit")
plot!(plt2, [controller_on_time, time[end]], fill(speed_pump_min, 2),
                linestyle = :dash,
                color = :blue,
                label = "Low Limit")

annotate_sp_phases(plt1)  # add phase annotations to the choke valve plot
annotate_sp_phases(plt2)  # add phase annotations to the pump speed plot
####################################################
# CONSTRAINTS 
####################################################
bhp_meas_vec = Matrix(df[!, [
        "BHP - Zone 1@Well 1", 
        "BHP - Zone 1@Well 2", 
        "BHP - Zone 1@Well 3", 
        "BHP - Zone 1@Well 4"]])

pump_outlet_press = df[!, "Pressure@Mainline P 500m"]  

plt3 = plot(time, pump_outlet_press,
ylabel = "Pump Outlet Pressure (bar)",
label = ["Pump Outlet" "High Limit" "Low Limit"],
legend = :right
)
plot!(plt3, [controller_on_time, time[end]], fill(pipe_pressure_max[2], 2),
                linestyle = :dash,
                color = :red,
                label = "High Limit")
plot!(plt3, [controller_on_time, time[end]], fill(pipe_pressure_min[2], 2),
                linestyle = :dash,
                color = :blue,
                label = "Low Limit")
annotate_sp_phases(plt3)

                
plt4 = plot(time, bhp_meas_vec,
label = ["Well 1" "Well 2" "Well 3" "Well 4"],
ylabel = "Downhole Pressure (bar)",
legend =:right
)
plot!(plt4, [controller_on_time, time[end]], fill(pipe_pressure_max[52], 2),
                linestyle = :dash,
                color = :red,
                label = "High Limit")
#= plot!(plt4, [controller_on_time, time[end]], fill(pipe_pressure_min[52], 2),
                linestyle = :dash,
                color = :blue,
                label = "Low Limit") =#
annotate_sp_phases(plt4)

####################################################
# INLET PRESSURES TO AT MAINLINE, LINE1, LINE2   
####################################################

p_inlet_mainline = df[!, "Pressure@Mainline P 500m"]
p_inlet_line1 = df[!, "Pressure@Line 1 P 500m"]   
p_inlet_line2 = df[!, "Pressure@Line 2 P 500m"]   

plt5 = plot(time, hcat(p_inlet_mainline, p_inlet_line1, p_inlet_line2),
label = ["Mainline" "Line 1" "Line 2"],
ylabel = "Inlet Pressure (bar)"
)

####################################################
# WELLHEAD PRESSURES  
# Choke Valve Outlet Pressures (Left) and Inlet Pressures (Right)
####################################################
wellhead_pressures = Matrix(df[!,[
    "Pressure Left@WH Choke 1",
    "Pressure Left@WH Choke 2",
    "Pressure Left@WH Choke 3",
    "Pressure Left@WH Choke 4"]])

plt6 = plot(time, wellhead_pressures,
label = ["Well 1" "Well 2" "Well 3" "Well 4"],
ylabel = "Wellhead Pressure (bar)"
)

choke_vlv_inlet_pressures = Matrix(df[!,[
    "Pressure Right@WH Choke 1",
    "Pressure Right@WH Choke 2",
    "Pressure Right@WH Choke 3",
    "Pressure Right@WH Choke 4"]])

plt7 = plot(time, choke_vlv_inlet_pressures,
label = ["Well 1" "Well 2" "Well 3" "Well 4"],
ylabel = "Choke Valve Inlet Pressure (bar)"
)

####################################################
# INLET FLOW RATES AT MAINLINE, LINE1, LINE2   
####################################################

w_inlet_mainline = df[!, "MFR - total liquid@Mainline P 500m"].*(-1)
w_inlet_line1 = df[!, "MFR - total liquid@Line 1 P 500m"].*(-1)
w_inlet_line2 = df[!, "MFR - total liquid@Line 2 P 500m"].*(-1)

plt8 = plot(time, hcat(w_inlet_mainline, w_inlet_line1, w_inlet_line2),
label = ["Mainline" "Line 1" "Line 2"],
ylabel = "Inlet Flow Rate (kg/s)"
)



plt9 = plot(time, wellhead_flows,
label = ["Well 1" "Well 2" "Well 3" "Well 4"],
ylabel = "Wellhead Flow Rate (kg/s)",
legend =:right
)


plot!(plt9, sp_time, well_flow_sps,
seriestype = :steppost,        # setpoints are piecewise-constant
linestyle = :dot,
color = ["red" "red" "red" "blue"],             # match each setpoint to its well's color
label = ["Well 1 SP" "Well 2 SP" "Well 3 SP" "Well 4 SP"],
xlims = (0, time[end])         # clip the setpoint horizon to the sim length
)

annotate_sp_phases(plt9)  # add phase annotations to the wellhead flow plot

####################################################
# JOIN PLOTS 3 & 4 (Pump Outlet Pressure and BHP) AS ONE FIGURE  
####################################################

# Pump-outlet + downhole pressure as one 2-panel figure.
# Each panel is roughly half the width of a standalone figure, so the phase
plt34 = plot(plt3, plt4, layout = (1, 2), size = (1200, 500),
left_margin = 10Plots.mm, right_margin = 10Plots.mm,
top_margin = 10Plots.mm, bottom_margin = 10Plots.mm)

####################################################
# JOIN PLOTS 1 & 2 (Pump Speed and Choke Valves) AS ONE FIGURE  
####################################################

# Pump-outlet + downhole pressure as one 2-panel figure
plt12 = plot(plt1, plt2, layout = (1, 2), size = (1200, 500),
left_margin = 10Plots.mm, right_margin = 10Plots.mm,
top_margin = 10Plots.mm, bottom_margin = 10Plots.mm)



plots = [plt1 plt2 plt3 plt4 plt5 plt6 plt7 plt8 plt9 plt34 plt12]
plot_dir = joinpath(case_dir, "lf_plots")
mkpath(plot_dir)

for (i, plt) in enumerate(plots)
    display(plt)
    savefig(plt, joinpath(plot_dir, "plt_$i.png"))
end