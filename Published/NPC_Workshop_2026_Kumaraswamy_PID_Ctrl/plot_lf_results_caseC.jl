function plot_lf_case_results(case_name)

    ####################################################
    # READ THE LEDAFLOW OUTPUT DATA  
    ####################################################
    file_name = case_name * "_trends.csv"
    full_output_csv_path = joinpath(@__DIR__, file_name)

    df = CSV.read(joinpath(@__DIR__, full_output_csv_path), DataFrame;
    header=10,            # column names start  from row 10  
    skipto=12,            # skip the units row (line 11)  
    delim = ',',
    missingstring = ["--"],   # LedaFlow null marker at t=0
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
    # SETPOINT DATA
    ####################################################
    # NOTE: these overlays must mirror the SETPOINT SCHEDULE in ctrl_scheme_A.jl.
    # t_ss offsets the controller clock by the 100 s steady-state run that precedes it.
    t_ss = 100.0
    # Pump pressure setpoint (bar) - held at nominal in this case
    pc_pump_sp = [197.8 for x in time]
    # Well flow setpoints (kg/s) - all wells held at nominal in this case
    well_flow_sp = [125.0 for x in time]
    # Well pressure override setpoints (bar) - single step
    pc_override_sp = [x >= 500+t_ss ? 243.0 : 248.0 for x in time]

    ####################################################
    # PLOT DATA
    ####################################################

    # Default Plot Settings
    tick_positions = 0:1200:time[end] 
    tick_labels = string.(round.(Int, tick_positions))

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
        xlabel="Time (s)",
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

    hline!(plt2, [3900],
                    linestyle = :dash,
                    color = :red,
                    label = "High Limit")
    hline!(plt2, [2800],
                    linestyle = :dash,
                    color = :blue,
                    label = "Low Limit")


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
    label = "Pressure",
    legend = :topright
    )

    plot!(plt3, time, pc_pump_sp,
                    linestyle = :dash,
                    color = :red,
                    label = "SP")


    plt4 = plot(time, bhp_meas_vec,
    label = ["Well 1" "Well 2" "Well 3" "Well 4"],
    ylabel = "Downhole Pressure (bar)",
    legend =:topright
    )
    plot!(plt4, pc_override_sp,
                    linestyle = :dash,
                    color = :red,
                    label = "Override SP")

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

    w_inlet_mainline = df[!, "MFR - total liquid@Mainline P 19500m"].*(-1)

    # No minimum-flow override in this scheme, so there is no mainline flow
    # setpoint to overlay here - the trace is unconstrained.
    plt8 = plot(time, w_inlet_mainline,
    ylabel = "Mainline Flow Rate (kg/s)",
    label = "Flow Rate",
    )

    plt9 = plot(time, wellhead_flows[:, 2:4],
    label = ["Well 2" "Well 3" "Well 4"],
    ylabel = "Wellhead Flow Rate(s) (kg/s)",
    legend =:topright
    )
    plot!(plt9, time, well_flow_sp,
                    linestyle = :dash,
                    color = :red,
                    label = "SP: Wells 2-4")

    plt10 = plot(time, wellhead_flows[:, 1],
    ylabel = "Wellhead 1 Flow Rate (kg/s)",
    legend =:topright
    )


    ####################################################
    # JOIN PLOTS 1 & 2 (Pump Speed and Choke Valves) AS ONE FIGURE  
    ####################################################

    # Pump-outlet + downhole pressure as one 2-panel figure
    plt12 = plot(plt1, plt2, layout = (1, 2), size = (1200, 500),
    left_margin = 10Plots.mm, right_margin = 10Plots.mm,
    top_margin = 10Plots.mm, bottom_margin = 10Plots.mm)



    plots = [plt1 plt2 plt3 plt4 plt5 plt6 plt7 plt8 plt9 plt10 plt12]
    plot_dir = joinpath(@__DIR__, case_name)
    mkpath(plot_dir)

    for (i, plt) in enumerate(plots)
        display(plt)
        savefig(plt, joinpath(plot_dir, "plt_$i.png"))
    end

end
