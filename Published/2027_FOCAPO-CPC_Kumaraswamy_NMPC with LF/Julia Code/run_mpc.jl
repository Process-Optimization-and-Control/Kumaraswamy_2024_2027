####################################################
# CONTROLLER STEP: SOLVE THE NMPC FOR ONE SAMPLE
#
# Called once per closed-loop sample, before run_simulator. The model itself was built once by
# mpc_model_build(); this function only injects everything that changes from sample to sample
# (initial condition, warm start, previously applied input, bias terms, parameter estimates,
# objective), solves, and returns the first input move plus the quantities run_simulator needs.
####################################################
function run_optimiser(initial_con, initial_guess, sim_t, model, prev_input, inj_est, VLV_2_CV_est, bias_well, bias_mainline, bias_pressure)

    ##############################
    # ASSIGN INITIAL GUESS
    # On the first sample this comes from the steady-state solution; afterwards it is the previous
    # solution shifted forward by one sample (see EXTRACT RESULTS below). 
    ##############################
    for varname in keys(initial_guess)
        var = model[Symbol(varname)]
        vals = JuMP.Containers.DenseAxisArray(initial_guess[varname], axes(var)...)
        set_start_value.(var, vals)
    end

    # Extend previous control inputs to full prediction horizon for setting initial guess
    choke_extended = repeat(Array(prev_input["choke_vlv_op"]), 1, T_mpc)
    speed_extended = repeat(Array(prev_input["speed_pump"]), 1, T_mpc)

    # Assign Initial Guess for Control Inputs 
    set_start_value.(model[:choke_vlv_op], choke_extended)
    set_start_value.(model[:speed_pump], speed_extended)

    ##############################
    # ASSIGN INITIAL CONDITIONS TO STATES
    # This is the state feedback: the start of the first collocation element is pinned to the value
    # reconstructed from the plant measurements in run_sim.jl. Every other element start is tied to
    # the previous element's end by a constraint in model_mpc.jl.
    ##############################
    for varname in state_vars
        var = model[Symbol(varname)]
        tuple_of_idx = axes(var)
        for idx in Iterators.product(tuple_of_idx...)
            # If variable is at t=1 and d=1, assign value from initial condition
            if last(idx) == 1 && idx[end-1] == 1
                fix.(var[idx...], initial_con[Name.((idx))...]; force=true)
            end
        end
    end

    ##############################
    # FIX PREVIOUSLY APPLIED INPUT FOR t=1 ROC CONSTRAINT
    # Without this the first move of each horizon would be unconstrained, and the closed loop could
    # step the inputs by more per sample than the rate limit allows.
    ##############################
    for k in choke_vlv_idx
        fix(model[:u_prev_choke][k], prev_input["choke_vlv_op"][Name(k), 1]; force=true)
    end
    for k in pump_idx
        fix(model[:u_prev_pump][k], prev_input["speed_pump"][Name(k), 1]; force=true)
    end

    ##############################
    # FIX INJECTIVITY, VALVE CV ESTIMATES, AND BIAS TERM
    # All of these arrive from run_sim.jl already resolved: if the corresponding flag in the case
    # file is off, the value passed in is the nominal parameter or a zero bias, so this block is
    # identical whether or not the offset-free layer is active.
    ##############################

    # Injectivity Values
    for k in well_idx 
        fix.(model[:injectivity][k, :, :], inj_est[Name(k)]; force=true)
    end 

    # Valve CV Estimates 
    for k in all_vlv_idx
        fix(model[:VLV_2_CV][k], VLV_2_CV_est[k]; force=true)
    end

    # Fix Bias Terms  
    # Bias term defined for wells 
    for k in well_idx 
        fix(model[:flow_bias][k], bias_well[Name(k)]; force=true)
    end 

    # Bias term defined for mainline flow measurement
    fix(model[:mainline_flow_bias], bias_mainline; force=true)

    # Bias term defined for pressure measurements at constrained nodes
    for k in pressure_con_idx
        fix(model[:p_bias][k], bias_pressure[Name(k)]; force=true)
    end

    ##############################
    # SPECIAL CONSTRAINTS
    # This rebuilds the objective with the setpoint slice for this sample time
    ##############################
    if special_constraints_dynamic !== nothing
        special_constraints_dynamic(model, sim_t)
    end

    ##############################
    # SOLVE MODEL
    ##############################
    set_optimizer_attribute(model, "print_level", 0)     # silence Ipopt; convergence is logged below instead
    set_attribute(model, "linear_solver", "ma27")        # HSL solver; requires a licensed HSL build of Ipopt
    set_optimizer_attribute(model, "mu_strategy", "adaptive")  
    optimize!(model)
    mpc_status = JuMP.termination_status(model)

    if mpc_status != MOI.LOCALLY_SOLVED && mpc_status != MOI.OPTIMAL
        println("MPC Model did not converge to global or local optima")
        println(mpc_status)
    else
        println("MPC Model Converged")
    end

    ##############################
    # LOG MPC STATUS TO TEXT FILE
    ##############################
    # Append the termination status for every MPC run to a single log file in case_dir.
    status_log_path = joinpath(case_dir, "mpc_status.txt")
    open(status_log_path, sim_t == DT ? "w" : "a") do io
        println(io, "t=$(sim_t)  status=$(mpc_status)")
    end

    ##############################
    # EXTRACT RESULTS
    ##############################
    new_initial_guess = Dict()
    predictions = Dict{String, Any}()

    # Bias-corrected outputs are saved for plotting but are not warm-started
    for varname in corrected_terms
        var = model[Symbol(varname)]
        data = Array(value.(var))
        predictions[varname] = NamedArray(data, axes(var))   # full-horizon prediction at this step
    end

    for varname in model_vars
        var = model[Symbol(varname)]
        n_dims = ndims(var)
        tuple_idx = axes(var)
        data = Array(value.(var))
        data_new = similar(data)
        predictions[varname] = NamedArray(data, axes(var))   # full-horizon prediction at this step

        # EXTRACT DATA FOR NEW INITIAL GUESS
        # Shift the solution one sample towards the present so it lines up with the next horizon:
        # sample t takes the value solved for t+1, and the final sample is duplicated to fill the
        # gap left at the end. 
        # The assumption is that time is always the last index.
        slicer_from = ntuple(d -> d == n_dims ? (2:T_mpc) : Colon(), n_dims)
        slicer_to = ntuple(d -> d == n_dims ? (1:(T_mpc-1)) : Colon(), n_dims)
        slicer_final = ntuple(d -> d == n_dims ? T_mpc : Colon(), n_dims)

        data_new[slicer_to...] = data[slicer_from...]
        data_new[slicer_final...] = data[slicer_final...]
        new_initial_guess[varname] = NamedArray(data_new, tuple_idx)
    end

    ##############################
    # SAVE PREDICTIONS TO FILE
    ##############################
    results_dir = joinpath(case_dir, "mpc_pred")
    jldsave(joinpath(results_dir, "pred_t$(sim_t).jld2"); predictions)

    ##############################
    # SAVE OPTIMAL MPC INPUT VALUES
    # Receding horizon: only the first move of the horizon is applied to the plant.
    ##############################
    optimal_input = Dict()

    for varname in u_vars
        var = model[Symbol(varname)]
        n_dims = ndims(var)
        data = Array(value.(var))
        opt = data[ntuple(d -> d == n_dims ? (1:1) : Colon(), n_dims)...]
        var_axes = axes(var)
        new_axes = (var_axes[1:end-1]..., 1:1)
        optimal_input[varname] = NamedArray(opt, new_axes)
    end 
    
    ##############################
    # SAVE DENSITY VALUES FOR FLOWS THROUGH CHOKE VALVES 
    ##############################
    # This is needed for online estimation of VLV_2_CV parameter in the simulator
    density_chokes = [value(model[:rho_out][nearest_pipe_idx[k], D, 1]) for k in choke_vlv_idx]
    
    ##############################
    # EXTRACT PREDICTED WELL OUTPUT FLOWS FOR CALCULATING BIAS
    ##############################
    # Taken at d = D, t = 1: the end of the first sample, which is the instant the plant is advanced
    # to next. run_sim.jl differences these against the measurements taken at that same instant.
    # Note these are the raw model predictions
    # bias update sees the full plant-model mismatch
    w_pred = NamedArray(Array(value.(model[:w_in][well_idx, D, 1])), (well_idx,))
    w_mainline_pred = value(model[:w_in][2, D, 1])
    pressure_pred = NamedArray(Array(value.(model[:p_node][pressure_con_idx, D, 1])), (pressure_con_idx,))

    return optimal_input, new_initial_guess, mpc_status, density_chokes, w_pred, w_mainline_pred, pressure_pred 
end



