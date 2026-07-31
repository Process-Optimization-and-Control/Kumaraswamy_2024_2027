##############################
# FUNCTION TO PLOT AND SAVE DATA   
############################## 

function plot_save_results(fplt_ylim, hplt_ylim, Hprofile, Fprofile, alpha_tk, original_dir)

    ##############################
    # PLOT RESULTS
    ##############################

    flabels = [L"F_{%$i}" for i in 0:O-1]
    flabels = permutedims(flabels)

    hlabels = [L"h_{%$i}" for i in 1:NI]
    hlabels = permutedims(hlabels)


    plt_fmin, plt_fmax = fplt_ylim
    plt_hmin, plt_hmax = hplt_ylim

    default(size=(600,400), legendfont=(12, "times"), tickfont=(12, "times"), guidefont=(14, "times"))
    plt1 = plot(1:nt, [Fprofile[:,i] for i in 1:O],
        xlabel=L"Time\,(min)", ylabel=L"F\,(m^{3}/min)", lw=2.5, 
        label=flabels, ylimits = (plt_fmin, plt_fmax),
        yformatter=:plain, legend=:outerright, linetype=:steppre, dpi=750)
    plt2 = plot(1:nt+1, [Hprofile[:,i]/hcap[i] for i in 1:NI], 
        xlabel=L"Time\,(min)", ylabel=L"Level\,(\%)", lw=2.5, 
        label=hlabels, ylimits = (plt_hmin, plt_hmax),
        yformatter=:plain, legend=:outerright, linetype=:steppre, dpi=750)
    display(plt1)
    display(plt2)
    composite_figure = plot(plt1, plt2, layout=(1, 2), size=(2000, 500))
    display(composite_figure)
    ##############################
    # SAVE RESULTS
    ############################## 
    # CREATE FOLDER FOR RESULTS (In Original Directory)    
    cd(original_dir)
    mkpath(case_name*"_Results")

    # MOVE INTO CASE FOLDER 
    cd(original_dir*"/"*case_name*"_Results")

    # SAVE FIGURES (In Results Directory)
    savefig(plt1, case_name*"_Flow.png")
    savefig(plt2, case_name*"_Inventory.png")

    # SAVE INDIVIDUAL FLOW & HEIGHT PLOTS 
    for i in 1:O
        k = i-1
        plt = plot(1:nt, [Fprofile[:,i]],
        xlabel=L"Time\,(min)", ylabel=L"F\,(m^{3}/min)",   lw=2.5, 
        label=L"F_{%$k}", ylimits=(plt_fmin, plt_fmax),
        yformatter=:plain, legend=:outerright, linetype=:steppre, dpi=750)
        display(plt)
        savefig(plt, case_name*"_F"*string(k)*".png")
    end

    for i in 1:NI
        plt = plot(1:nt+1, [Hprofile[:,i]/hcap[i]],
        xlabel=L"Time\,(min)", ylabel=L"Level\,(\%)", lw=2.5,
        label=L"h_{%$i}", ylimits=(plt_hmin, plt_hmax),
        yformatter=:plain, legend=:outerright, linetype=:steppre, dpi=750)
        #display(plt)
        savefig(plt, case_name*"_h"*string(i)*".png")
    end

    # SAVE RESULTS (In Results Directory)
    writedlm(case_name*"_Inventory_profile.csv", Hprofile)
    writedlm(case_name*"_Flow_profile.csv", Fprofile)

    ##############################
    # SAVE PARAMETERS FOR THIS CASE 
    ##############################

    open(case_name*"_parameters.txt", "w") do file
        write(file, "NI = $(NI)\n")
        write(file, "O = $(O)\n")
        write(file, "area, m2 = $(a)\n")
        write(file, "height, m = $(hcap)\n")
        write(file, "Connectivity Matrix M:")
        write(file, "\n")
        writedlm(file, M, ',')
        write(file, "Producer Nodes (Connection to Flows) =")
        write(file, "\n")
        writedlm(file, P, ',')
        write(file, "Consumer Nodes (Connection to Flows) =")
        write(file, "\n")
        writedlm(file, C, ',')
        write(file, "Algorithm Weights")
        write(file, "Flow Rate Alpha (Optimization Weights) = $(alpha_f)\n")
        write(file, "Edge Weights = $(EI)\n")
        write(file, "Tank Inventory Weights = $(alpha_tk)\n")
        write(file, "MPC Parameters")
        write(file, "\n")
        write(file, "Prediction Horizon, mins = $(np)\n")
        write(file, "Control Horizon, mins = $(nm)\n")
        write(file, "Total Simulation Time, mins = $(nt)\n")
        write(file, "Penalty for Constraints = $(penalty)\n")
        write(file, "Penalty Multiplier (Height Constraints, Fmin Constraints, Deltah Constraints) = $(penalty_mult)\n")
        write(file, "Objective for MPC = $(obj_mpc)\n")
        write(file, "Model Constraints\n")
        write(file, "fmin Constraint Binary = $(fmin_binary)\n")
        write(file, "deltah Constraint Binary = $(deltah_binary)\n")
        write(file, "deltah Constraint  = $(deltah)\n")
        write(file, "hmax = $(hmax)\n")
        write(file, "hmin = $(hmin)\n")
        write(file, "fmin = $(fmin)\n")
        write(file, "Maximum Flow Constraint Vector\n")
        writedlm(file, fmax_vec', ',')
        write(file, "Height Disturbance Vector\n")
        writedlm(file, Bd_vec', ',')
        write(file, "Initial Conditions\n")
        write(file, "Initial Flow at f0 = $(f0)\n")
    end

    #######################################
    # IMPORTANT: RESTORE ORIGINAL DIRECTORY 
    #######################################

    cd(original_dir)

end 


