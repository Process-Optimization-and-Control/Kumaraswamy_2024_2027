####################################################
# ALGORITHM FOR PLOTTING AND SAVING RESULTS 
#################################################### 

export plot_save

function plot_save(results)

    default(legendfont=(16, "Computer Modern"), tickfont=(16, "Computer Modern", :bold), guidefont=(16, "Computer Modern", :bold))
    plot_size = (700, 500)
    pipe_seg_names = ["Mainline" "Line 1" "Line 2" "Well 1 " "Well 2" "Well 3" "Well 4"]

    ####################################################
    # PRESSURE PLOTS 
    #################################################### 
    
    dsel = last(axes(results["p_in"])[end-1])
    data1 = results["p_in"][Name.(pipe_seg_inlet), dsel, :]'
    ymin, ymax = extrema(data1)
    plt1 = (
    legend=pipe_seg_names,
    var = data1,
    ylabel = "Arc Inlet Pressure (bar)",
    xlabel = "Time (hr)",
    ylims = (0.9*ymin, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )
    
    dsel = last(axes(results["p_in"])[end-1])
    data2 = results["p_in"][Name.(well_idx), dsel, :]'
    ymin, ymax = extrema(data2)
    plt2 = (
    legend=permutedims(["Wellhead $i" for i in 1:length(well_idx)]),
    var = data2,
    ylabel = "Pressure (bar)",
    xlabel = "Time (hr)",
    ylims = (0.9*ymin, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )

    dsel = last(axes(results["w_in"])[end-1])
    data3 = results["w_in"][Name.(pipe_seg_inlet), dsel, :]'
    ymin, ymax = extrema(data3)
    plt3 = (
    legend=pipe_seg_names,
    var = data3,
    ylabel = "Flow (kg/s)",
    xlabel = "Time (hr)",
    ylims = (0.9*ymin, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )

    dsel = last(axes(results["w_in_corrected"])[end-1])
    data4 = results["w_in_corrected"][Name.(well_idx), dsel, :]'

    println(data4)
    ymin, ymax = extrema(data4)
    plt4 = (
    legend=permutedims(["Wellhead $i" for i in 1:length(well_idx)]),
    var = data4,
    ylabel = "Flow (kg/s)",
    xlabel = "Time (hr)",
    ylims = (0.98*ymin, 1.02*ymax),
    xlims =(-Inf, Inf) 
    )

    data5 = results["choke_vlv_op"][Name.(choke_vlv_idx), :]'
    println(data5)
    plt5 = (
    legend=permutedims(["Valve $i" for i in 1:length(choke_vlv_idx)]),
    var = data5,
    ylabel = "Valve OP",
    xlabel = "Time (hr)",
    ylims = (-Inf, Inf),
    xlims =(-Inf, Inf) 
    )

    dsel = last(axes(results["p_out"])[end-1])
    data6 = results["p_out"][Name.(well_idx), dsel, :]'
    ymin, ymax = extrema(data6)
    plt6 = (
    legend=permutedims(["Downhole $i" for i in 1:length(well_idx)]),
    var = data6,
    ylabel = "Pressure (bar)",
    xlabel = "Time (hr)",
    ylims = (0.97*ymin, 1.03*ymax),
    xlims =(-Inf, Inf) 
    )

    dsel = last(axes(results["cavitation_idx"])[end-1])
    data7 = results["cavitation_idx"][Name.(choke_vlv_idx), dsel, :]'
    ymin, ymax = extrema(data7)
    plt7 = (
    legend=permutedims(["Valve $i" for i in 1:length(choke_vlv_idx)]),
    var = data7,
    ylabel = "Cavitation Index",
    xlabel = "Time (hr)",
    ylims = (cavitation_min*0.9, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )

    dsel = last(axes(results["deltap_wellhead_choke"])[end-1])
    data9 = results["deltap_wellhead_choke"][Name.(choke_vlv_idx), dsel, :]'
    ymin, ymax = extrema(data9)
    plt9 = (
    legend=permutedims(["Valve $i" for i in 1:length(choke_vlv_idx)]),
    var = data9,
    ylabel = "Pressure Drop Across Valve (bar)",
    xlabel = "Time (hr)",
    ylims = (0, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )
    
    dsel = last(axes(results["w_in"])[end-1])
    data10 = results["w_in"][Name.(choke_vlv_idx), dsel, :]'
    ymin, ymax = extrema(data10)
    plt10 = (
    legend=permutedims(["Valve $i" for i in 1:length(choke_vlv_idx)]),
    var = data10,
    ylabel = "Flow through Valve",
    xlabel = "Time (hr)",
    ylims = (0, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )


    dsel = last(axes(results["p_reservoir"])[end-1])
    data11 = results["p_reservoir"][dsel, :]
    ymin, ymax = extrema(data11)
    plt11 = (
    legend=nothing,
    var = data11,
    ylabel = "Reservoir Pressure (bar)",
    xlabel = "Time (hr)",
    ylims = (0, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )

    dsel = last(axes(results["p_in"])[end-1])
    data12 = results["p_in"][Name.(choke_vlv_idx), dsel, :]'
    ymin, ymax = extrema(data12)
    plt12 = (
    legend=permutedims(["Valve $i" for i in 1:length(choke_vlv_idx)]),
    var = data12,
    ylabel = "Inlet Pressure (bar)",
    xlabel = "Time (hr)",
    ylims = (0.9*ymin, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )

    dsel = last(axes(results["p_out"])[end-1])
    data13 = results["p_out"][Name.(choke_vlv_idx), dsel, :]'
    ymin, ymax = extrema(data13)
    plt13 = (
    legend=permutedims(["Valve $i" for i in 1:length(choke_vlv_idx)]),
    var = data13,
    ylabel = "Outlet Pressure (bar)",
    xlabel = "Time (hr)",
    ylims = (0.9*ymin, 1.1*ymax),
    xlims =(-Inf, Inf) 
    )
    dsel = last(axes(results["p_node"])[end-1])
    node_pressure = Array(results["p_node"][Name.(pressure_con_idx), dsel, :])   # size: n_nodes x n_time

    # same backward difference used in the constraint
    pressure_roc = hcat(
        fill(NaN, length(pressure_con_idx), 1),
        diff(node_pressure; dims=2) ./ DT
    )'

    valid = pressure_roc[.!isnan.(pressure_roc)]
    ymin, ymax = extrema(valid)

    plt14 = (
    legend = permutedims(["Node $i" for i in pressure_con_idx]),
    var = pressure_roc,
    ylabel = "Pressure ROC (bar/s)",
    xlabel = "Time (hr)",
    ylims = (1.1*ymin, 1.1*ymax),
    xlims = (-Inf, Inf)
    )


    pump_plots = []
    for i in pump_idx
            # Which Tank? 
        k = findfirst(==(i), pump_idx)

        pump_work_data = vec(sum(results["pump_work"][Name(i), :, :], dims=1))   # total work per interval = sum of collocation-node shares
        ymin, ymax = extrema(pump_work_data)
        pump_work_plt = (
        legend = nothing,
        var = pump_work_data,
        ylabel = "Work Done by Pump $k (MWh)",
        xlabel = "Time (hr)",
        ylims = (0.9*ymin, 1.1*ymax),
        xlims =(-Inf, Inf) 
        )

        dsel = last(axes(results["p_in"])[end-1])
        pump_pressure_data = [results["p_in"][Name(i), dsel, :], results["p_out"][Name(i), dsel, :]]
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
        ylabel = "Pump $k Speed, rpm",
        xlabel = "Time (hr)",
        ylims = (0.9*ymin, 1.1*ymax),
        xlims =(-Inf, Inf) 
        )

        dsel = last(axes(results["head_meters"])[end-1])
        head_data = results["head_meters"][Name(i), dsel, :]
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
    pump_plots_bar_indicator = repeat([1, 0, 0, 0], length(pump_idx))
    plots = [plt1, plt2, plt3, plt4, plt5, plt6, plt7, plt9, plt10, plt11, plt12, plt13, plt14]
    # MOVE INTO RESULTS DIRECTORY
    mkpath(joinpath(@__DIR__, CASE_NAME))
    cd(case_dir)

    function plot_data(x)
        timesec = DT:DT:NP           # still in seconds
        tick_positions = DT:3600:NP  # ticks every hour
        tick_labels = string.(round.(Int, tick_positions ./ 3600))

        plt = plot(DT:DT:NP, x.var,
        xlabel=x.xlabel, ylabel=x.ylabel, label = x.legend, lw=3.5, legend=:outerright, yformatter=:plain, xformatter =:plain, dpi=750, ylims=x.ylims, xlims=x.xlims, 
        xticks=(tick_positions, tick_labels),
       size=plot_size,
        grid=false)
        if plt_number == 7     
            # add constraint lines per pipe
            hline!(plt, [cavitation_min],
                linestyle = :dash,
                color = :red,
                label = "Low Limit")

        end

        display(plt)
        savefig(plt, CASE_NAME*"_plt"*string(plt_number)*".png")
        plt_number = plt_number + 1

    end 
    
    function plot_bar(x)
        tick_positions = DT:3600:NP  # ticks every hour
        tick_labels = string.(tick_positions ./ 3600)

        # Convert NamedArray to regular Vector
        y_data = Vector(x.var)
        x_data = collect(DT:DT:NP)

        plt = bar(x_data, y_data,
            xlabel=x.xlabel, ylabel=x.ylabel, label=x.legend, 
            legend=:outerright, yformatter=:plain, xformatter=:plain, 
            dpi=750, ylims=x.ylims, xlims=x.xlims, 
            xticks=(tick_positions, tick_labels),
            size=plot_size,
            grid=false)

        display(plt)
        savefig(plt, CASE_NAME*"_plt"*string(plt_number)*".png")
        plt_number = plt_number + 1 
    end 

    plt_number = 1
    for x in plots 
        plot_data(x)
    end 

    for (i, x) in enumerate(pump_plots)
        if pump_plots_bar_indicator[i] == 1
            plot_bar(x)
        else
            plot_data(x)
        end
    end

    ####################################################
    # SWITCH TO ORIGINAL DIRECTORY
    ####################################################
    cd(original_dir)

end



