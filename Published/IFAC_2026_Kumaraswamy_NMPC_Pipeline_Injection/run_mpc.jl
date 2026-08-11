
function run_optimiser(initial_con, initial_guess, sim_t, model)

    
    ##############################
    # ASSIGN INITIAL GUESS FOR ALL MODEL VARIABLES 
    ##############################
    for varname in model_vars
        var = model[Symbol(varname)]  
        vals = JuMP.Containers.DenseAxisArray(initial_guess[varname], axes(var)...)
        set_start_value.(var, vals)
    end     

    ##############################
    # ASSIGN INITIAL CONDITIONS TO STATE VARIABLES AND U_VARS 
    ##############################
    for varname in union(state_vars, u_vars)
        var = model[Symbol(varname)]
        n_dims = ndims(var)
        tuple_idx_init = ntuple(d -> d == n_dims ? 1 : Colon(), n_dims)
        vals = JuMP.Containers.DenseAxisArray(initial_con[varname], axes(var)...)
        fix.(var[tuple_idx_init...], vals[tuple_idx_init...]; force=true)
    end 

    ##############################
    # SPECIAL CONSTRAINTS 
    ##############################
    if special_constraints_dynamic !== nothing 
        special_constraints_dynamic(model, sim_t)
    end 

    ##############################
    # SOLVE MODEL 
    ##############################
    set_optimizer_attribute(model, "print_level", 2)
    set_attribute(model, "linear_solver", "ma57")
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
    # EXTRACT ALL PREDICTIONS FOR LOGGING 
    ##############################
    new_initial_guess = Dict()
    predictions = Dict{String, Any}()

    for varname in model_vars
        var = model[Symbol(varname)]
        n_dims = ndims(var)
        tuple_idx = axes(var)
        data = Array(value.(var))
        data_new = similar(data)
        predictions[varname] = NamedArray(data, tuple_idx)

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

    optimal_input = Dict()
    for varname in u_vars
        var = model[Symbol(varname)]
        n_dims = ndims(var)
        data = Array(value.(var))
        opt = data[ntuple(d -> d == n_dims ? (2:2) : Colon(), n_dims)...]
        new_axes = (axes(var)[1:end-1]..., 1:1)
        optimal_input[varname] = NamedArray(opt, new_axes)
    end


    ##############################
    # SAVE PREDICTIONS TO FILE
    ##############################
    mpc_pred_dir = joinpath(case_dir, "mpc_pred")
    jldsave(joinpath(mpc_pred_dir, "pred_t$(sim_t).jld2"); predictions)

    return optimal_input, new_initial_guess
end



