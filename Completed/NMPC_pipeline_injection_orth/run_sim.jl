####################################################
# SIMULATION MODEL DRIVER
# Advances the plant model by one sample, i.e. by one collocation finite element. The inputs are
# held at the move the NMPC just returned, the d = 1 slice is set to where the plant currently is,
# and the model is solved as a square feasibility problem for the state at the end of the element.
####################################################

function run_simulator(initial_con, optimal_input, sim_t, model)

    println(optimal_input)
    ##############################
    # FIX OPTIMAL INPUTS FROM MPC
    ##############################
    # The NMPC returns a single move per input. It is applied at the start of the element and held
    # across it (zero order hold), so one value per input is fixed for the whole element.
    for varname in u_vars
        var = model[Symbol(varname)]
        for k in axes(var)[1]
            fix(var[k, 1], optimal_input[varname][Name(k), 1]; force=true)
        end
    end

    ##############################
    # ASSIGN INITIAL CONDITIONS TO STATES AND ALGEBRAIC VARIABLES
    ##############################
    # The initial condition is the d = 1 slice, which after the move to d = 2:D declarations only
    # rho still has: it is fixed to the plant state at the end of the previous element. Every
    # other variable has no d = 1 entry, so it falls through to the else branch and is merely warm
    # started from that same state at each collocation point.
    for varname in keys(initial_con)
        var = model[Symbol(varname)]
        tuple_of_idx = axes(var)
        for idx in Iterators.product(tuple_of_idx...)
            # If variable is at d=1, assign value from initial condition
            if idx[end-1] == 1
                fix.(var[idx...], initial_con[varname][Name.((idx))...]; force=true)
            # For all other cases, only assign initial guess
            # Use initial condition as initial guess
            else
                # initial_con is a single snapshot: its collocation axis holds one label, 1,
                idx_modified = Base.setindex(idx, 1, length(idx)-1)
                set_start_value(var[idx...], initial_con[varname][Name.((idx_modified))...])
            end
        end
    end

    ##############################
    # SOLVE MODEL
    ##############################
    @objective(model, Min, 0.0)
    set_optimizer_attribute(model, "print_level", 0)
    set_attribute(model, "linear_solver", "ma27")
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
    # Every model variable at d = D of the element, i.e. at the end of the interval. This is both
    # the plant snapshot the notebook appends to the results and the initial condition for the
    # next sample. The inputs are excluded: they are carried separately, as optimal_input.
    # `sim_arr` is a plain Array, so the collocation dimension is indexed by position, not by
    # label: d = D is its LAST position whether the variable was declared over 1:D or 2:D.
    new_initial_con = Dict()

    for varname in setdiff(model_vars, u_vars)
        sim_arr = Array(value.(model[Symbol(varname)]))
        var_axes = axes(model[Symbol(varname)])
        n = ndims(sim_arr)
        d_end = size(sim_arr, n-1)
        ic_data = sim_arr[ntuple(i -> i == n-1 ? (d_end:d_end) : (i == n ? (1:1) : Colon()), n)...]
        new_axes = (var_axes[1:end-2]..., 1:1, 1:1)
        new_initial_con[varname] = NamedArray(ic_data, new_axes)
    end

    return new_initial_con, sim_status
end
