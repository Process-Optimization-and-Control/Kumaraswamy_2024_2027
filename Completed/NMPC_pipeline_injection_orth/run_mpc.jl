####################################################
# NMPC DRIVER
# Solves the NMPC for one closed-loop sample: warm start the model, fix the initial condition,
# attach the case objective, solve, and hand back the first move together with a shifted warm
# start for the next sample.
####################################################

function run_optimiser(initial_con, initial_guess, sim_t, model, prev_input)

    ##############################
    # ASSIGN INITIAL GUESS FOR ALL MODEL VARIABLES
    ##############################
    # Everything except the inputs is warm started from the previous NMPC solution, shifted one
    # element forward.
    for varname in setdiff(model_vars, u_vars)
        var = model[Symbol(varname)]
        vals = JuMP.Containers.DenseAxisArray(initial_guess[varname], axes(var)...)
        set_start_value.(var, vals)
    end

    # The inputs have no collocation dimension, so they are warm started separately, from the whole
    # input trajectory the previous solve returned rather than from the single move that was applied.
    # The model object is reused across samples, so that trajectory is still on it and is read back
    # with value(); shifting it one element forward matches what is done to the states above
    # On the first sample there is no previous solve to read, so the applied move is held flat.
    if has_values(model)
        for varname in u_vars
            var = model[Symbol(varname)]
            prev_traj = Array(value.(var))
            shifted = similar(prev_traj)
            shifted[:, 1:(T_mpc-1)] = prev_traj[:, 2:T_mpc]
            shifted[:, T_mpc] = prev_traj[:, T_mpc]
            set_start_value.(var, shifted)
        end
    else
        choke_extended = repeat(Array(prev_input["choke_vlv_op"]), 1, T_mpc)
        speed_extended = repeat(Array(prev_input["speed_pump"]), 1, T_mpc)

        set_start_value.(model[:choke_vlv_op], choke_extended)
        set_start_value.(model[:speed_pump], speed_extended)
    end

    ##############################
    # ASSIGN INITIAL CONDITIONS TO RHO AND P_NODE 
    ##############################
    # With orthogonal collocation the initial condition is the d = 1 slice of the first finite
    # element: no model equation is written there, so whatever sits at (d=1, t=1) is fixed to the
    # plant state returned by the simulation model. In practice that is rho and p_node only - they
    # are the two variables declared over d = 1:D. Everything else is declared over d = 2:D, has no
    # d = 1 entry, and so falls through this loop untouched.
    for varname in keys(initial_con)
        var = model[Symbol(varname)]
        tuple_of_idx = axes(var)
        for idx in Iterators.product(tuple_of_idx...)
            # If variable is at t=1 and d=1, assign value from initial condition
            if last(idx) == 1 && idx[end-1] == 1
                fix.(var[idx...], initial_con[varname][Name.((idx))...]; force=true)
            end
        end
    end

    ##############################
    # FIX PREVIOUSLY APPLIED INPUT FOR t=1 ROC CONSTRAINT
    ##############################
    # The inputs carry no d = 1 initial condition, so the rate of change limit on the first move is
    # written against these two variables instead.
    for k in choke_vlv_idx
        fix(model[:u_prev_choke][k], prev_input["choke_vlv_op"][Name(k), 1]; force=true)
    end
    for k in pump_idx
        fix(model[:u_prev_pump][k], prev_input["speed_pump"][Name(k), 1]; force=true)
    end

    ##############################
    # SPECIAL CONSTRAINTS
    ##############################
    # The case file supplies the objective function (and any case specific constraints) here.
    if special_constraints_dynamic !== nothing
        special_constraints_dynamic(model, sim_t)
    end

    ##############################
    # SOLVE MODEL
    ##############################
    set_optimizer_attribute(model, "print_level", 0)
    set_attribute(model, "linear_solver", "ma27")
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
    # EXTRACT NEW INITIAL GUESS FOR NEXT MPC RUN
    ##############################
    new_initial_guess = Dict()

    for varname in setdiff(model_vars, u_vars)
        var = model[Symbol(varname)]
        n_dims = ndims(var)
        tuple_idx = axes(var)
        data = Array(value.(var))
        data_new = similar(data)

        # EXTRACT DATA FOR NEW INITIAL GUESS
        # ASSUMPTION: Time indexes are continuous
        slicer_from = ntuple(d -> d == n_dims ? (2:T_mpc) : Colon(), n_dims)
        slicer_to = ntuple(d -> d == n_dims ? (1:(T_mpc-1)) : Colon(), n_dims)
        slicer_final = ntuple(d -> d == n_dims ? T_mpc : Colon(), n_dims)

        data_new[slicer_to...] = data[slicer_from...]
        data_new[slicer_final...] = data[slicer_final...]
        new_initial_guess[varname] = NamedArray(data_new, tuple_idx)
    end

    ##############################
    # EXTRACT OPTIMAL INPUTS
    ##############################
    # Only the first move of the horizon is applied to the plant.
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

    return optimal_input, new_initial_guess, mpc_status
end
