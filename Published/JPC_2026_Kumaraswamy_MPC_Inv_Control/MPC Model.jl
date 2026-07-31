using JuMP
using Ipopt
using LinearAlgebra
using HiGHS # Alternative Optimizer to IPOPT 


function create_mpc_model(alpha_tk)
    
    ##############################
    # CALCULATE A AND B MATRICES  
    ##############################
    A = I(NI)
    B = M./a

    ##############################
    # MPC PARAMETERS  
    ##############################
    # Prediction Horizon 
    T = np

    # Time Discount Factor 
    gamma = 0.5

    ##############################
    # INITIAL GUESS 
    ##############################
    # Use initial condition as initial guess 
    fguess = fill(f0, T)
    hguess = fill(h0, T)

    ##############################
    # CONSUMER INDEX  
    ############################## 
    # Indexes over which to multiply weights in the objective function
    nc = size(C)[1]
    nc_idx = zeros(nc)
    for k in 1:nc
        idx = findfirst(x -> x == 1, C[k, :])
        nc_idx[k] = idx
    end
    nc_idx = Int.(nc_idx)

    ##############################
    # BUILD MODEL 
    ##############################
    # Create a model & Select Optimizer 

    if optimizer == "Ipopt"
        model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "print_level", 0)  # Set verbosity to 0 to suppress output
        set_optimizer_attribute(model, "max_iter", 100000)  # Increase iteration limit  
    end

    if optimizer == "HiGHS"
        model = Model(HiGHS.Optimizer)
        set_optimizer_attribute(model, "log_to_console", false)
        set_optimizer_attribute(model, "output_flag", false)   
    end 

    if optimizer != "Ipopt" && optimizer != "HiGHS"
        println("Optimizer not assigned properly. Select either Ipopt or HiGHS")
    end 

    # Define variables 
    @variable(model, h[ni=1:NI, t=1:T], start=hguess[t][ni])
    @variable(model, f[o=1:O, t=1:T]>=0, start=fguess[t][o])
    @variable(model, hslack[ni=1:NI, t=1:T]>=0)
    @variable(model, hsurplus[ni=1:NI, t=1:T]>=0)
    @variable(model, fsurplus[o=1:O, t=1:T]>=0)
    @variable(model, deltahslack1[ni=1:NI, t=1:T]>=0)
    @variable(model, deltahslack2[ni=1:NI, t=1:T]>=0)

    # Disturbances 
    @variable(model, fmax[o=1:O]>=0)
    @variable(model, Bd[ni=1:NI]>=0)

    # Input Regularization 
    @variable(model, reg_slack1[o=1:O, t=1:T]>=0)
    @variable(model, reg_slack2[o=1:O, t=1:T]>=0)
    @variable(model, prev_f[o=1:O]>=0, start=f0[o])
    
    vars = (h = h, f = f)

    # Model Equations  
    @constraint(model, [t in 1:T-1], h[:, t+1] .- A*h[:, t] .- B*f[:, t] .- Bd[:] == 0)

    # Max and Min Inventory Constraints (Soft Constraints)
    @constraint(model, [t in 2:T], h[:, t] <= hmax + hslack[:, t])
    @constraint(model, [t in 2:T], h[ :,t] + hsurplus[:, t] >= hmin)

    # Max Flow Rate (Hard Constraint)
    @constraint(model, [t in 1:T], f[:,t] <= fmax)

    # Max Height Change Constraint
    if deltah_binary == 1
        @constraint(model, [i in 1:NI, t in 2:T], (h[i, t] - h[i, t-1]) <= deltah[i] + deltahslack1[i, t]) 
        @constraint(model, [i in 1:NI, t in 2:T], -(h[i, t] - h[i, t-1]) <= deltah[i] + deltahslack2[i, t]) 
    end

    # Min Flow Rate (Soft Constraint) 
    if fmin_binary == 1
        @constraint(model, [i in 1:O, t in 1:T], f[i, t] + fsurplus[i, t] >= fmin[i])
    end 

    # Input Usage Constraint Beyond Control Horizon 
    @constraint(model, [t in nm+1:T], f[:, t] .- f[:, t-1] == 0)

    # Additional / Special Case-Specific Constraints 
    if special_constraints !== nothing 
        special_constraints(model, vars, T)
    end 
    
    # EQUALITY CONSTRAINTS INCLUDED IN OBJECTIVE 
    if obj_mpc == 5
        eq_idx, equality_penalty = equality_con_param
        equality_objective = 0
        # Create indexed slack and surplus variables
        @variable(model, equality_slack[i in eq_idx, j in eq_idx, t=1:T] >= 0)
        @variable(model, equality_surplus[i in eq_idx, j in eq_idx, t=1:T] >= 0)
        for i in eq_idx
            for j in eq_idx
                if i == j 
                    @constraint(model, [t=1:T], equality_slack[i, j, t] == 0)
                    @constraint(model, [t=1:T], equality_surplus[i, j, t] == 0)
                end  
                if i != j 
                    @constraint(model, [t=1:T], (f[i, t] - f[j, t]) <= equality_slack[i, j, t])
                    @constraint(model, [t=1:T], (f[i, t] - f[j, t]) + equality_surplus[i, j, t] >= 0)
                    equality_objective += sum(equality_slack[i, j, t] + equality_surplus[i, j, t] for t in 1:T)*equality_penalty
                end
            end
        end 
    end

    # Input Regularization Constraint 
    @constraint(model, [o in 1:O, t in 2:T], (f[o, t] - f[o, t-1]) <= reg_slack1[o, t])
    @constraint(model, [o in 1:O, t in 2:T], -(f[o, t] - f[o, t-1]) <= reg_slack2[o, t])
    @constraint(model, [o in 1:O], (prev_f[o] - f[o, 1]) <= reg_slack1[o, 1])
    @constraint(model, [o in 1:O], -(prev_f[o] - f[o, 1]) <= reg_slack2[o, 1])

    ##############################
    # DEFINE OBJECTIVE FUNCTION 
    ##############################

    if obj_mpc == 0 
        # UNREACHABLE SETPOINT OBJECTIVE FUNCTION 
        @objective(model, Min, sum(sqrt(sum( (gamma^t)*(unreachable_fsp[o] - f[o, t])^2  for o in 1:O)) for t in 1:T )  + penalty_mult[1]*penalty*sum(hslack[ni,t] for ni in 1:NI for t in 1:T) + penalty_mult[1]*penalty*sum(hsurplus[ni,t] for ni in 1:NI for t in 1:T) + penalty_mult[2]*penalty*sum(fsurplus[o,t] for o in 1:O for t in 1:T) + penalty_mult[3]*penalty*sum(deltahslack1[ni,t] for ni in 1:NI for t in 1:T) + penalty_mult[3]*penalty*sum(deltahslack2[ni,t] for ni in 1:NI for t in 1:T) )
    end 

    if obj_mpc == 1
        # WEIGHTED OBJECTIVE FUNCTION   
        # AREA INCLUDED 
        # HEIGHT SOFT CONSTRAINTS: SLACK / SURPLUS VARIABLES 
        @objective(model, Min, -sum((sum(f[idx,t]*alpha_f[idx] for idx in nc_idx) + sum(a[ni]*h[ni, t]*alpha_tk[ni] for ni in 1:NI))*gamma^t for t in 1:T) + penalty_mult[1]*penalty*sum(hslack[ni,t] for ni in 1:NI for t in 1:T) + penalty_mult[1]*penalty*sum(hsurplus[ni,t] for ni in 1:NI for t in 1:T) + penalty_mult[2]*penalty*sum(fsurplus[o,t] for o in 1:O for t in 1:T) + penalty_mult[3]*penalty*sum(deltahslack1[ni,t] for ni in 1:NI for t in 1:T) + penalty_mult[3]*penalty*sum(deltahslack2[ni,t] for ni in 1:NI for t in 1:T) )
    end 
    
    ##############################
    # OBJ FUNCTIONS TO RESOLVE MULTIPLICITY OF SOLUTIONS
    ##############################
    
    if obj_mpc == 4
        # WEIGHTED OBJECTIVE FUNCTION   
        # AREA INCLUDED 
        # HEIGHT SOFT CONSTRAINTS: SLACK / SURPLUS VARIABLES 
        # INPUT REGULARIZATION: VERY SMALL PENALTY (1E-4) - IMPLEMENTED AS SOFT CONSTRAINTS
        @objective(model, Min, -sum((sum(f[idx,t]*alpha_f[idx] for idx in nc_idx) + sum(a[ni]*h[ni, t]*alpha_tk[ni] for ni in 1:NI))*gamma^t for t in 1:T) + penalty_mult[1]*penalty*sum(hslack[ni,t] for ni in 1:NI for t in 1:T) + penalty_mult[1]*penalty*sum(hsurplus[ni,t] for ni in 1:NI for t in 1:T) + penalty_mult[2]*penalty*sum(fsurplus[o,t] for o in 1:O for t in 1:T) + penalty_mult[3]*penalty*sum(deltahslack1[ni,t] for ni in 1:NI for t in 1:T) + penalty_mult[3]*penalty*sum(deltahslack2[ni,t] for ni in 1:NI for t in 1:T) + input_reg_penalty*sum((reg_slack1[o, t] + reg_slack2[o, t]) for o in 1:O for t in 1:T )  )
    end 

    if obj_mpc == 5
        # WEIGHTED OBJECTIVE FUNCTION   
        # AREA INCLUDED 
        # HEIGHT SOFT CONSTRAINTS: SLACK / SURPLUS VARIABLES 
        # EQUALITY CONSTRAINTS FOR FLOW INDEX 
        @objective(model, Min, -sum((sum(f[idx,t]*alpha_f[idx] for idx in nc_idx) + sum(a[ni]*h[ni, t]*alpha_tk[ni] for ni in 1:NI))*gamma^t for t in 1:T) + penalty_mult[1]*penalty*sum(hslack[ni,t] for ni in 1:NI for t in 1:T) + penalty_mult[1]*penalty*sum(hsurplus[ni,t] for ni in 1:NI for t in 1:T) + penalty_mult[2]*penalty*sum(fsurplus[o,t] for o in 1:O for t in 1:T) + penalty_mult[3]*penalty*sum(deltahslack1[ni,t] for ni in 1:NI for t in 1:T) + penalty_mult[3]*penalty*sum(deltahslack2[ni,t] for ni in 1:NI for t in 1:T) + equality_objective)
    end 

    if obj_mpc == 7
        # WEIGHTED OBJECTIVE FUNCTION   
        # AREA INCLUDED 
        # HEIGHT SOFT CONSTRAINTS: SLACK / SURPLUS VARIABLES 
        # HEIGHT SOFT CONSTRAINTS ONLY FOR TANKS AND NOT PROCESS NODES 
        storage_tk_idx = [1, 2, 3, 5, 7]
        @objective(model, Min, -sum((sum(f[idx,t]*alpha_f[idx] for idx in nc_idx) + sum(a[ni]*h[ni, t]*alpha_tk[ni] for ni in storage_tk_idx))*gamma^t for t in 1:T) + penalty_mult[1]*penalty*sum(hslack[ni,t] for ni in storage_tk_idx for t in 1:T) + penalty_mult[1]*penalty*sum(hsurplus[ni,t] for ni in storage_tk_idx for t in 1:T) + penalty_mult[2]*penalty*sum(fsurplus[o,t] for o in 1:O for t in 1:T) + penalty_mult[3]*penalty*sum(deltahslack1[ni,t] for ni in storage_tk_idx for t in 1:T) + penalty_mult[3]*penalty*sum(deltahslack2[ni,t] for ni in storage_tk_idx for t in 1:T) )
    end 
    


    return model
end 