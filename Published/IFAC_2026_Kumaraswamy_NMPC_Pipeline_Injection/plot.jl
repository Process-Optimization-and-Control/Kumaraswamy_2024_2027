####################################################
# PLOTTING AND SAVING RESULTS 
#################################################### 

function plot_save(results)
    
    # DEFAULT SETTINGS FOR PLOTS 

    default(legendfont=(16, "Computer Modern"), tickfont=(16, "Computer Modern", :bold), guidefont=(16, "Computer Modern", :bold))

    ####################################################
    # PLOTS FOR PRESSURE, FLOW, VALVE OPENING, CAVITATION INDEX, PRESSURE DROP ACROSS VALVE 
    ####################################################    
    data1 = results["p_in"][Name.(pipe_seg_inlet), :]'
    ymin, ymax = extrema(data1)
    plt1 = (
    legend=pipe_seg_names,
    var = data1,
    ylabel = "Inlet Pressure (bar)",
    xlabel = "Time (hr)",
    ylims = (0.9*ymin, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )
    data2 = results["p_in"][Name.(well_idx),:]'
    ymin, ymax = extrema(data2)
    plt2 = (
    legend=permutedims(["Wellhead $i" for i in 1:length(well_idx)]),
    var = data2,
    ylabel = "Pressure (bar)",
    xlabel = "Time (hr)",
    ylims = (0.9*ymin, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )

    data3 = results["w_in"][Name.(pipe_seg_inlet), :]'
    ymin, ymax = extrema(data3)
    plt3 = (
    legend=pipe_seg_names,
    var = data3,
    ylabel = "Flow (kg/s)",
    xlabel = "Time (hr)",
    ylims = (0.9*ymin, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )

    data4 = results["w_in"][Name.(well_idx),:]'
    ymin, ymax = extrema(data4)
    plt4 = (
    legend=permutedims(["Wellhead $i" for i in 1:length(well_idx)]),
    var = data4,
    ylabel = "Flow (kg/s)",
    xlabel = "Time (hr)",
    ylims = (0.98*ymin, 1.02*ymax),
    xlims =(-Inf, Inf) 
    )

    data5 = results["choke_vlv_op"][Name.(choke_vlv_idx),:]'
    ymin, ymax = extrema(data5)
    plt5 = (
    legend=permutedims(["Well $i" for i in 1:length(choke_vlv_idx)]),
    var = data5,
    ylabel = "Valve OP",
    xlabel = "Time (hr)",
    ylims = (0.97*ymin, 1.03*ymax),
    xlims =(-Inf, Inf) 
    )

    data6 = results["p_out"][Name.(well_idx),:]'
    ymin, ymax = extrema(data6)
    plt6 = (
    legend=permutedims(["Downhole $i" for i in 1:length(well_idx)]),
    var = data6,
    ylabel = "Pressure (bar)",
    xlabel = "Time (hr)",
    ylims = (0.97*ymin, 1.03*ymax),
    xlims =(-Inf, Inf) 
    )

    data7 = results["cavitation_idx"][Name.(choke_vlv_idx),:]'
    ymin, ymax = extrema(data7)
    plt7 = (
    legend=permutedims(["Valve $i" for i in 1:length(choke_vlv_idx)]),
    var = data7,
    ylabel = "Cavitation Index",
    xlabel = "Time (hr)",
    ylims = (0.9*ymin, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )

    data8 = results["deltap_wellhead_choke"][Name.(choke_vlv_idx),:]'
    ymin, ymax = extrema(data8)
    plt8 = (
    legend=permutedims(["Valve $i" for i in 1:length(choke_vlv_idx)]),
    var = data8,
    ylabel = "Pressure Drop Across Valve (bar)",
    xlabel = "Time (hr)",
    ylims = (0, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )
    
    data9 = results["w_in"][Name.(choke_vlv_idx),:]'
    ymin, ymax = extrema(data9)
    plt9 = (
    legend=permutedims(["Valve $i" for i in 1:length(choke_vlv_idx)]),
    var = data9,
    ylabel = "Flow through Valve",
    xlabel = "Time (hr)",
    ylims = (0, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )

    ####################################################
    # PLOTS FOR PUMPS IN THE SYSTEM 
    ##################################################
    pump_plots = []
    for i in pump_idx
        k = findfirst(==(i), pump_idx)

        pump_work_data = results["pump_work"][Name(i), :]
        ymin, ymax = extrema(pump_work_data)
        pump_work_plt = (
        legend = nothing,
        var = pump_work_data,
        ylabel = "Work Done by Pump $k (MWh)",
        xlabel = "Time (hr)",
        ylims = (0.9*ymin, 1.1*ymax),
        xlims =(-Inf, Inf) 
        )

        pump_pressure_data = [results["p_in"][Name(i), :], results["p_out"][Name(i), :]]
        ymin, ymax = extrema(vcat(pump_pressure_data...))
        pump_pressure_plt = (
        legend = ["Inlet" "Outlet"],
        var = pump_pressure_data,
        ylabel = "Pressure (bar)",
        xlabel = "Time (hr)",
        ylims = (0.9*ymin, 1.1*ymax),
        xlims =(-Inf, Inf) 
        )
    
        pump_speed_data = results["speed_pump"][Name(i), :]
        ymin, ymax = extrema(pump_speed_data)
        pump_speed_plt = (
        legend = nothing,
        var = pump_speed_data,
        ylabel = "Pump Speed, rpm",
        xlabel = "Time (hr)",
        ylims = (0.9*ymin, 1.1*ymax),
        xlims =(-Inf, Inf) 
        )

        head_data = results["head_meters"][Name(i), :]
        ymin, ymax = extrema(head_data)
        head_plt = (
        legend = nothing,
        var = head_data,
        ylabel = "Head $k, m",
        xlabel = "Time (hr)",
        ylims = (0.9*ymin, 1.1*ymax),
        xlims =(-Inf, Inf) 
        )
        

        push!(pump_plots, pump_work_plt)
        push!(pump_plots, pump_pressure_plt)
        push!(pump_plots, pump_speed_plt)
        push!(pump_plots, head_plt)

    end 

    # INDICATE WHICH PUMP PLOTS ARE BAR PLOTS AND WHICH ARE LINE PLOTS
    pump_plots_bar_indicator = [1, 0, 0, 0]

    # LIST OF ALL PLOTS TO BE GENERATED
    plots = [plt1, plt2, plt3, plt4, plt5, plt6, plt7, plt8, plt9]

    # MOVE INTO RESULTS DIRECTORY 
    cd(case_dir)

    function plot_data(x)
        tick_positions = 0:3600:NT  # ticks every hour
        tick_labels = string.(tick_positions ./ 3600)

        plt = plot(0:DT:NT, x.var,
        xlabel=x.xlabel, ylabel=x.ylabel, label = x.legend, lw=3.5, legend=:outerright, yformatter=:plain, xformatter =:plain, dpi=750, ylims=x.ylims, xlims=x.xlims, 
        xticks=(tick_positions, tick_labels),
        grid=false)
        display(plt)
        savefig(plt, CASE_NAME*"_plt"*string(plt_number)*".png")
        plt_number = plt_number + 1 
    end 

    function plot_bar(x)
        tick_positions = 0:3600:NT  # ticks every hour
        tick_labels = string.(tick_positions ./ 3600)

        # Convert NamedArray to regular Vector
        y_data = Vector(x.var)
        x_data = collect(0:DT:NT)

        plt = bar(x_data, y_data,
            xlabel=x.xlabel, ylabel=x.ylabel, label=x.legend, 
            legend=:outerright, yformatter=:plain, xformatter=:plain, 
            dpi=750, ylims=x.ylims, xlims=x.xlims, 
            xticks=(tick_positions, tick_labels),
            grid=false)

        display(plt)
        savefig(plt, CASE_NAME*"_plt"*string(plt_number)*".png")
        plt_number = plt_number + 1 
    end 

    
    plt_number = 1
    for x in plots 
        plot_data(x)
    end 

    for x in pump_plots 
        if pump_plots_bar_indicator[findfirst(==(x), pump_plots)] == 1
            plot_bar(x)
        else
            plot_data(x)
        end
    end 

    ####################################################
    # SAVE DATA IN JSON FORMAT
    ####################################################

    for varname in model_vars
        data = results[varname]
        open(CASE_NAME*"_"*varname*".json", "w") do io
            JSON.print(io, data)
        end 
    end 

    ####################################################
    # SWITCH TO ORIGINAL DIRECTORY
    ####################################################
    cd(original_dir)

end



