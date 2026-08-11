####################################################
# SIMULATION MODEL DRIVER
# Advances the plant model by one sample. The inputs are held at the move the NMPC just returned,
# the differential states are set to where the plant currently is, and the model is solved as a
# square feasibility problem for the state at the end of the interval.
####################################################

function run_simulator(initial_con, initial_guess, optimal_input, sim_t, model)

    ##############################
    # FIX U VARIABLES FOR t = 1 and 2
    # At t = 1, input comes from initial condition 
    # At t = 2, input comes from the optimal input returned by the NMPC
    ##############################
    # The NMPC returns a single move per input. It is applied at the start of the interval and
    # held across it (zero order hold), so the same value is fixed at both trapezoidal end points.
    for varname in u_vars
        var = model[Symbol(varname)]
        for k in axes(var)[1]
            fix(var[k, 1], initial_con[varname][Name(k), 1]; force=true)
            fix(var[k, 2], optimal_input[varname][Name(k), 1]; force=true)
        end
    end


    ##############################
    # ASSIGN INITIAL CONDITIONS TO STATES
    ##############################
    # Only the differential state (state_vars) are fixed. Every other variable at t=1 is recomputed by the
    # algebraic equations from these states and the applied inputs
    for varname in state_vars
        var = model[Symbol(varname)]
        for k in axes(var)[1]
            fix(var[k, 1], initial_con[varname][Name(k), 1]; force=true)
        end
    end

    ##############################
    # ASSIGN INITIAL GUESS FOR ALL MODEL VARIABLES
    ##############################
    # Every time step is warm started from the predictions of the previous simulation 
    for varname in keys(initial_guess)
        varname in u_vars && continue
        var = model[Symbol(varname)]
        n_dims = ndims(var)
        for idx in Iterators.product(axes(var)...)
            idx_guess = Base.setindex(idx, 1, n_dims) # same index, but taken at t = 1
            set_start_value(var[idx...], initial_guess[varname][Name.(idx_guess)...])
        end
    end

    ##############################
    # SOLVE MODEL
    ##############################
    set_optimizer_attribute(model, "print_level", 0)
    set_attribute(model, "linear_solver", "ma57")
    optimize!(model)
    sim_status = JuMP.termination_status(model)

    if sim_status != MOI.LOCALLY_SOLVED && sim_status != MOI.OPTIMAL
        println("Simulation Model did not converge to global or local optima")
        println(sim_status)
    else
        println("Simulation Model converged")
    end

    ##############################
    # LOG SIMULATION STATUS TO TEXT FILE
    ##############################
    # Append the termination status for every simulation run to a single log file in case_dir.
    status_log_path = joinpath(case_dir, "sim_status.txt")
    open(status_log_path, sim_t == DT ? "w" : "a") do io
        println(io, "t=$(sim_t)  status=$(sim_status)")
    end

    ##############################
    # EXTRACT NEW SIMULATION RESULTS
    ##############################
    # Every model variable at the end of the interval. This is the plant snapshot the notebook
    # appends to the results being plotted.
    sim_results = Dict()

    for varname in model_vars
        var = model[Symbol(varname)]
        n_dims = ndims(var)
        tuple_of_labels = axes(var)
        tuple_of_labels_snapshot = (tuple_of_labels[1:(end-1)]..., Base.OneTo(1))
        data = Array(value.(var))

        # EXTRACT DATA AT t = T_sim, the end of the interval
        slicer = ntuple(d -> d == n_dims ? (T_sim:T_sim) : Colon(), n_dims)
        sim_results[varname] = NamedArray(data[slicer...], tuple_of_labels_snapshot)
    end

    # Only the differential states carry over as the initial condition for the next sample.
    new_initial_con = Dict(varname => sim_results[varname] for varname in union(state_vars, u_vars))

    return new_initial_con, sim_results
end
