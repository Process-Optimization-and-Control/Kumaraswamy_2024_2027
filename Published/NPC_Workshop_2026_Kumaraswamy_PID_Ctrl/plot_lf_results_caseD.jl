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
    bhp_meas_vec = Matrix(df[!, [
            "BHP - Zone 1@Well 1", 
            "BHP - Zone 1@Well 2", 
            "BHP - Zone 1@Well 3", 
            "BHP - Zone 1@Well 4"]])


    ####################################################
    # SETPOINT DATA   
    ####################################################
    t_ss = 100.0
    # Pump pressure setpoint (bar)
    pc_pump_sp = [x > 500+t_ss && x < 12500+t_ss ? 205.0 : x >= 12500+t_ss ? 190.0 : 197.8 for x in time]
    # Mainline min well flow override setpoint (kg/s)
    min_mainline_flow_sp = [490.0 for x in time]
    # Well flow setpoints (kg/s)    
    well_flow_sp = [125.0 for x in time]
    well4_flow_sp = [x >= 4500+t_ss && x < 8500+t_ss ? 135.0 : 125.0 for x in time]
    # Well pressure override setpoint (bar)
    pc_override_sp = [248.0 for x in time]

    ####################################################
    # PLOT DATA   
    ####################################################

    # Default Plot Settings
    tick_positions = 0:3600:time[end] 
    tick_labels = string.(round.(Int, tick_positions))
    # ---- Poster color scheme: Okabe-Ito (colorblind-safe, high-contrast) ----
    okabe_ito = [
        colorant"#0072B2",  # blue
        colorant"#D55E00",  # vermillion
        colorant"#009E73",  # green
        colorant"#CC79A7",  # reddish purple
        colorant"#E69F00",  # orange
        colorant"#56B4E9",  # sky blue
        colorant"#F0E442",  # yellow
        colorant"#000000",  # black
    ]
    # Default well colors 
    well_colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]

    default(
        palette = okabe_ito,
        legendfont=(18, "Computer Modern"),
        tickfont=(16, "Computer Modern", :bold),
        guidefont=(20, "Computer Modern", :bold),
        titlefont=(22, "Computer Modern", :bold),
        lw=4,
        yformatter=:plain,
        xformatter=:plain,
        legend=:right,
        size=(1000, 650),
        grid=false,
        framestyle=:box,
        left_margin=8mm,
        right_margin=8mm,
        top_margin=6mm,
        bottom_margin=8mm,
        xlabel="Time (s)",
        dpi=1500,
        xticks=(tick_positions, tick_labels)
    )    
    ####################################################
    # MVs   
    ####################################################

    # Choke Valve Openings Plot
    plt1 = plot(time, choke_valve_openings, 
    label = ["Well 1" "Well 2" "Well 3" "Well 4"], 
    ylabel = "Valve Opening",
    legend=:right
    )

    # Pump Speed Plot 
    plt2 = plot(time, pump_speed,
    ylabel = "Pump Speed (rpm)",
    label = "Speed",
    legend = :bottomright
    )


    hline!(plt2, [2800, 3900], color = :grey30, linestyle = :dot, lw = 5,
       linewidth = 2.5, label = "MV Limits")

    ####################################################
    # Controlled Variables (CVs) 
    ####################################################

    pump_outlet_press = df[!, "Pressure@Mainline P 500m"]  

    plt3 = plot(time, pump_outlet_press,
    ylabel = "Pump Outlet Pressure (bar)",
    label = "Pressure",
    legend = :right
    )

    plot!(plt3, time, pc_pump_sp,
                    linestyle = :dot, lw = 5,
                    color = :grey30,
                    label = "SP")

    plt4 = plot(time, bhp_meas_vec[:, 2:4],
    label = ["Well 2" "Well 3" "Well 4"],
    color = [okabe_ito[2] okabe_ito[3] okabe_ito[4]],
    ylabel = "Downhole Pressure (bar)",
    legend =:right
    )
    plot!(plt4, pc_override_sp,
                    linestyle = :dot, lw = 5,
                    color = :grey30,
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

    plt8 = plot(time, w_inlet_mainline,
    ylabel = "Mainline Flow Rate (kg/s)",
    label ="Flow Rate",
    )

    plot!(plt8, time, min_mainline_flow_sp,
                    linestyle = :dot, lw = 5,
                    color = :grey30,
                    label = "Override SP")


    plt9 = plot(time, wellhead_flows[:, 1:4],
    label = ["Well 1" "Well 2" "Well 3" "Well 4"],
    color  = [okabe_ito[1] okabe_ito[2] okabe_ito[3] okabe_ito[4]],
    ylabel = "Wellhead Flow Rate(s) (kg/s)",
    legend =:topright
    )
    plot!(plt9, time, well4_flow_sp,
                    linestyle = :dot, lw = 5,
                    color = okabe_ito[4],
                    label = "SP: Well 4")
    plot!(plt9, time, well_flow_sp,
                    linestyle = :dot, lw = 5,
                    color = :grey30,
                    label = "SP: Wells 2 & 3")


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


    # Keep this list identical to plot_lf_results_caseA/B.jl so that plt_N.svg means
    # the same panel in every case folder.
    plots = [plt1 plt2 plt3 plt4 plt5 plt6 plt7 plt8 plt9 plt10 plt12]
    plot_dir = joinpath(@__DIR__, case_name)
    mkpath(plot_dir)

    for (i, plt) in enumerate(plots)
        display(plt)
        savefig(plt, joinpath(plot_dir, "plt_$i.svg"))
    end

end
