using JuMP
using Ipopt
using LinearAlgebra
using HiGHS

function create_mpc_model(alpha_tk)

    #STATES
    # Voulme (5 tanks), Xb, Ss, So
    #Xb_tk, Ss_tk (4 tanks) 
    

    # CALCULATE A AND B MATRICES  
    A = I(NI)
    B = M  # No division by area for volume-based model
    

    # MPC PARAMETERS  
    T = npred  # Prediction horizon
    gamma = 0.98  # Time discount factor

    # INITIAL CONDITION AS INTIAL QUESS FOR OPTIMIZER
    fguess = fill(f0, T)
    Vguess = fill(V0, T)
    Xb_tk_guess = repeat(reshape(Xb_tk0, 1, NI, 1), 4, 1, T) 
    Ss_tk_guess = repeat(reshape(Ss_tk0, 1, NI, 1), 4, 1, T) 

    # INDEXES OVER WHICH TO MULTILY WEIGHTS IN THE OBJECTIVE FUNCTION
    nc = size(C)[1]
    nc_idx = zeros(nc)
    for k in 1:nc
        idx = findfirst(x -> x == 1, C[k, :])
        nc_idx[k] = idx
    end
    nc_idx = Int.(nc_idx)


    # BUILD MODEL 
    if optimizer == "Ipopt"
        model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "print_level", 0)
        set_optimizer_attribute(model, "max_iter", 100000)
    elseif optimizer == "HiGHS"
        model = Model(HiGHS.Optimizer)
        set_optimizer_attribute(model, "log_to_console", false)
        set_optimizer_attribute(model, "output_flag", false)
    else
        println("Optimizer not assigned properly. Select either Ipopt or HiGHS")
    end


    # DEFINE VARIABLES
    @variable(model, V[ni=1:NI, t=1:T], start=Vguess[t][ni])
    @variable(model, f[o=1:O, t=1:T]>=0, start=fguess[t][o])
   
    @variable(model, RAS[1:T] >= 0, start=RAS0) #guess value for MPC
    @variable(model, Kla[1:T]>= 0, start=Kla0) #guess value for MPC

    @variable(model, q_in[t=1:T] >= 0)
    @variable(model, q_out[t=1:T] >= 0)

    @variable(model, V_BRX[1:Col, 1:T])
    @variable(model, Xb[1:Col, 1:T])
    @variable(model, Ss[1:Col, 1:T])
    @variable(model, So[1:Col, 1:T])
    @variable(model, D[1:T])
    

    @variable(model, Vslack[ni=1:NI, t=1:T]>=0)
    @variable(model, Vsurplus[ni=1:NI, t=1:T]>=0)
    @variable(model, fsurplus[o=1:O, t=1:T]>=0) 
    @variable(model, Vslack_BRX[1:Col, 1:T] >= 0)
    @variable(model, Vsurplus_BRX[1:Col, 1:T] >= 0)
    @variable(model, Ssslack[1:Col, 1:T] >= 0)
    #@variable(model, RASslack[1:T] >= 0)
    #@variable(model, RASsurplus[1:T] >= 0)

    # COMPONENT VARIABLES
    @variable(model, Xb_tk[c=1:Col, ni=1:NI, t=1:T] >= 0) 
    @variable(model, Ss_tk[c=1:Col, ni=1:NI, t=1:T] >= 0)
    @variable(model, Xb_flow[c=1:Col, o=1:O, t=1:T] >= 0)
    @variable(model, Ss_flow[c=1:Col, o=1:O, t=1:T] >= 0)
    set_start_value.(Xb_tk, Xb_tk_guess)
    set_start_value.(Ss_tk, Ss_tk_guess)


    # DISTURBANCES
    @variable(model, fmax[o=1:O] >= 0)
    @variable(model, Bd[ni=1:NI]>=0)
    @variable(model, Ss_in_max >= 0)

    # INPUT REGULARIZATION
    @variable(model, reg_slack1[o=1:O, t=2:T] >= 0) 
    @variable(model, reg_slack2[o=1:O, t=2:T] >= 0) #for ras



    #Finding in and outflow of bioreactor, will not work with several inflows with different concentrations 
    idx_q_in = [o for o in 1:O if M[NI_BRX, o] > 0]
    idx_q_out = [o for o in 1:O if M[NI_BRX, o] < 0]
    @constraint(model, [t in 1:T], q_in[t] == sum(M[NI_BRX, o] * f[o, t] for o in idx_q_in))
    @constraint(model, [t in 1:T], q_out[t] == sum(-M[NI_BRX, o] * f[o, t] for o in idx_q_out))


    # MODEL EQUATIONS  
    tks = [i for i in 1:NI if i != NI_BRX]
    @constraint(model, [t in 1:T-1], V[tks, t+1] .- A[tks, tks]*V[tks, t] .- dt.*B[tks, :]*f[:, t] .- dt.*Bd[tks] .== 0) #Euler 

    @constraint(model, [t in 1:T], V[NI_BRX, t] == V_BRX[1, t])
    @constraint(model, [t=1:T], Xb_tk[:,NI_BRX,t] .== Xb[:,t])
    @constraint(model, [t=1:T], Ss_tk[:,NI_BRX,t] .== Ss[:,t])
    @constraint(model, [t=1:T], Xb_flow[:, 3, t] .== Xb_tk[:, NI_BRX, t])
    @constraint(model, [t=1:T], Ss_flow[:, 3, t] .== Ss_tk[:, NI_BRX, t])


    #Component equations for tank 1, 3, 4, 5 
    Moutlet = [findall(x -> x == -1, M[i, :]) for i in 1:size(M, 1)]

    tks_norm = [i for i in tks if i != 3] #tanks except bioreaactor and settling tank 3

    @constraint(model, [col in 1:Col, ni in tks_norm, t in 1:T], Xb_tk[col,ni,t] .== Xb_flow[col,Moutlet[ni],t])
    @constraint(model, [col in 1:Col, ni in tks_norm, t in 1:T], Ss_tk[col,ni,t] .== Ss_flow[col,Moutlet[ni],t])

    @constraint(model, [t in 1:T-1, ni in tks_norm], Adot'*Xb_tk[1:Col, ni, t] == dt*sum(M[ni, o]*f[o, t]*Xb_flow[2:Col, o, t] for o in 1:O))
    @constraint(model, [t in 1:T-1, ni in tks_norm], Adot'*Ss_tk[1:Col, ni, t] == dt*sum(M[ni, o]*f[o, t]*Ss_flow[2:Col, o, t] for o in 1:O))

    @constraint(model, [t=2:T], Xb_tk[1,:,t] .== Xb_tk[Col,:,t-1])
    @constraint(model, [t=2:T], Ss_tk[1,:,t] .== Ss_tk[Col,:,t-1])

    @constraint(model, [col in 1:Col, t in 1:T], Xb_tk[col, 3, t] .== Xb[col, t])
    @constraint(model, [col in 1:Col, t in 1:T], Ss_tk[col, 3, t] .== Ss[col, t])
    

    # Constraints flow concentrations
    @constraint(model, [t=1:T], Xb_flow[:,4,t] .== 0.10 * Xb_tk[:,3,t])
    @constraint(model, [t=1:T], Ss_flow[:,4,t] .== 0.10 * Ss_tk[:,3,t])

    @constraint(model, [t=1:T], Xb_flow[:,6,t] .== 0.90 * Xb_tk[:,3,t])
    @constraint(model, [t=1:T], Ss_flow[:,6,t] .== 0.90 * Ss_tk[:,3,t])

    # Ratio constraint RAS/inflow (RAS = MV) 
    @constraint(model, [t=1:T], RAS[t] * f[1, t] == f[recycle_idx, t])

    # Dilution rate
    @constraint(model, [t in 1:T], D[t] * V_BRX[1, t] == q_in[t])

    # Collocation model equation for states 
    @constraint(model, [t=1:T], Adot'*V_BRX[1:Col, t] .== dt * (q_in[t] .- q_out[t]))
   
    
    @constraint(model, [t=1:T], Adot'*Xb[1:Col, t] .== dt * (D[t] .* (Xb_flow[2:Col,2,t] .- Xb[2:Col, t]) .+ (mu_max .* Ss[2:Col, t] ./ (Ks .+ Ss[2:Col, t])) .* Xb[2:Col, t] .- b .* Xb[2:Col, t]))
    @constraint(model, [t=1:T], Adot'*Ss[1:Col, t] .== dt * (D[t] .* (Ss_flow[2:Col,2,t] .- Ss[2:Col, t]) .- (1 / Y) .* (mu_max .* Ss[2:Col, t] ./ (Ks .+ Ss[2:Col, t])) .* Xb[2:Col, t]))
    @constraint(model, [t=1:T], Adot'*So[1:Col, t] .== dt * (D[t] .* (So_in .- So[2:Col, t]) .+ Kla[t] .* (So_sat .- So[2:Col, t]) .- ((1 - Y)/Y) .* (mu_max .* Ss[2:Col, t] ./ (Ks .+ Ss[2:Col, t])) .* Xb[2:Col, t] .- b .* Xb[2:Col, t])) 
   

    # Collocation continuity: end of previous element = start of next
    @constraint(model, [t=2:T], V_BRX[1, t] == V_BRX[Col, t-1])
    @constraint(model, [t=2:T], Xb[1, t] == Xb[Col, t-1])
    @constraint(model, [t=2:T], Ss[1, t] == Ss[Col, t-1])
    @constraint(model, [t=2:T], So[1, t] == So[Col, t-1])

    # CONSTRAINTS
    @constraint(model, [t in 1:T], V[tks, t] .<= Vmax[tks] .+ Vslack[tks, t])
    @constraint(model, [t in 1:T], V[tks, t] .>= Vmin[tks] .- Vsurplus[tks, t])

    @constraint(model, [t in 1:T], V_BRX[:, t] .<= Vmax[NI_BRX] .+ Vslack_BRX[:, t])
    @constraint(model, [t in 1:T], V_BRX[:, t] .>= Vmin[NI_BRX] .- Vsurplus_BRX[:, t])

    @constraint(model, [t in 1:T], Ss[:, t] .<= Ssmax .+ Ssslack[:, t])
    @constraint(model, [o=1:O, t=1:T], f[o, t] <= fmax[o])
    @constraint(model, [t in 1:T], Ss_flow[:,1,t] .== Ss_in_max)

    @constraint(model, [t=1:T], Kla[t] <= Kla_max) 
    @constraint(model, [t=1:T], Kla[t] >= Kla_min) 

    #@constraint(model, [t=1:T], RAS[t] <= RAS_max + RASslack[t])
    #@constraint(model, [t=1:T], RAS[t] >= RAS_min - RASsurplus[t])

    @constraint(model, [t=1:T], f[4,t]==1/3*f[3,t])  #make tank three not drain
    @constraint(model, [t=1:T], f[6,t]==2/3*f[3,t])  #make tank three not drain


    # Min flow constraints (soft)
    if fmin_binary == 1
        @constraint(model, [o=1:O, t=2:T], f[o, t] + fsurplus[o, t] >= fmin[o])
    end

    # Input regularization constraints
    if fmin_binary == 2
        @constraint(model, [o=1:O, t=2:T], (f[o, t] - f[o, t-1]) <= reg_slack1[o, t])
        @constraint(model, [o=1:O, t=2:T], -(f[o, t] - f[o, t-1]) <= reg_slack2[o, t])
    end 
    ##############################
    # OBJECTIVE FUNCTION
    ##############################
    if obj_mpc == 9
        # Maximize bioreactor output flow, penalize volume deviations and substrate violations
        @objective(model, Min, 
            -sum((sum(f[idx, t] * alpha_f[idx] for idx in nc_idx) + sum(V[ni, t]*alpha_V[ni]./100 for ni in 1:NI)) * gamma^t for t in 1:T) +
            sum((So[col, t] .- So_sp)^2 for col in 1:Col for t in 1:T)*Q_so + 
            penalty_mult[3] * penalty * sum(Vslack[i, t] for i in tks for t in 1:T) +
            penalty_mult[3] * penalty * sum(Vsurplus[i, t] for i in tks for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vslack_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vsurplus_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[1] * penalty * sum(Ssslack[i, t] for i in 1:Col for t in 1:T) +
            sum((RAS[t]-RAS[t-1])^2 for t in 2:T)*Q_RAS
        )
    end


    if obj_mpc == 10
        # Maximize bioreactor output flow, penalize volume deviations and substrate violations
        @objective(model, Min, 
            -sum((sum(f[idx, t] * alpha_f[idx] for idx in nc_idx) + sum(V[ni, t]*alpha_V[ni]./100 for ni in 1:NI) + sum(Xb[col, t] for col in 1:Col for t in 1:T)) * gamma^t for t in 1:T) +
            sum((So[col, t] .- So_sp)^2 for col in 1:Col for t in 1:T)*Q_so + 
            penalty_mult[3] * penalty * sum(Vslack[i, t] for i in tks for t in 1:T) +
            penalty_mult[3] * penalty * sum(Vsurplus[i, t] for i in tks for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vslack_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vsurplus_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[1] * penalty * sum(Ssslack[i, t] for i in 1:Col for t in 1:T) +
            +sum((RAS[t]-RAS[t-1])^2 for t in 2:T)*Q_RAS
        )
    end

    if obj_mpc == 11
        # Maximize bioreactor output flow, penalize volume deviations and substrate violations
        #input regularizaation on flows and RAS
        @objective(model, Min, 
            -sum((sum(f[idx, t] * alpha_f[idx] for idx in nc_idx) + sum(V[ni, t]*alpha_V[ni]./100 for ni in 1:NI) + sum(Xb[col, t] for col in 1:Col for t in 1:T)) * gamma^t for t in 1:T) +
            sum((So[col, t] .- So_sp)^2 for col in 1:Col for t in 1:T)*Q_so + 
            penalty_mult[3] * penalty * sum(Vslack[i, t] for i in tks for t in 1:T) +
            penalty_mult[3] * penalty * sum(Vsurplus[i, t] for i in tks for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vslack_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vsurplus_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[1] * penalty * sum(Ssslack[i, t] for i in 1:Col for t in 1:T) +
            1*sum((reg_slack1[o, t] + reg_slack2[o, t]) for o in 1:O for t in 2:T) +
            sum((RAS[t]-RAS[t-1])^2 for t in 2:T)*Q_RAS
        )
    end

        if obj_mpc == 12
        # Maximize bioreactor output flow, penalize volume deviations and substrate violations
        @objective(model, Min, 
            -sum((sum(f[idx, t] * alpha_f[idx] for idx in nc_idx) + sum(V[ni, t]*alpha_V[ni]./100 for ni in 1:NI) + sum(Xb[col, t]*2 for col in 1:Col for t in 1:T)) * gamma^t for t in 1:T) +
            sum((So[col, t] .- So_sp)^2 for col in 1:Col for t in 1:T)*Q_so + 
            penalty_mult[3] * penalty * sum(Vslack[i, t] for i in tks for t in 1:T) +
            penalty_mult[3] * penalty * sum(Vsurplus[i, t] for i in tks for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vslack_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vsurplus_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[1] * penalty * sum(Ssslack[i, t] for i in 1:Col for t in 1:T) +
            +sum((RAS[t]-RAS[t-1])^2 for t in 2:T)*Q_RAS
        )
    end

    if obj_mpc == 13
        # Maximize  outflow, penalize volume deviations and substrate violations
        @objective(model, Min, 
            -sum((sum(f[idx, t] * alpha_f[idx] for idx in nc_idx) + sum(V[ni, t]*alpha_V[ni]./100 for ni in 1:NI)) * gamma^t for t in 1:T) +
            sum((So[col, t] .- So_sp)^2 for col in 1:Col for t in 1:T)*Q_so + 
            penalty_mult[3] * penalty * sum(Vslack[i, t] for i in tks for t in 1:T) +
            penalty_mult[3] * penalty * sum(Vsurplus[i, t] for i in tks for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vslack_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vsurplus_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[1] * penalty * sum(Ssslack[i, t] for i in 1:Col for t in 1:T)
        )
    end

    if obj_mpc == 13
        # Maximize  outflow, penalize volume deviations and substrate violations
        @objective(model, Min, 
            -sum((sum(f[idx, t] * alpha_f[idx] for idx in nc_idx) + sum(V[ni, t]*alpha_V[ni]./100 for ni in 1:NI)) * gamma^t for t in 1:T) +
            sum((So[col, t] .- So_sp)^2 for col in 1:Col for t in 1:T)*Q_so + 
            penalty_mult[3] * penalty * sum(Vslack[i, t] for i in tks for t in 1:T) +
            penalty_mult[3] * penalty * sum(Vsurplus[i, t] for i in tks for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vslack_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vsurplus_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[1] * penalty * sum(Ssslack[i, t] for i in 1:Col for t in 1:T)
        )
    end

    if obj_mpc == 14
        # Maximize  outflow, penalize volume deviations and substrate violations
        @objective(model, Min, 
            -sum((sum(f[idx, t] * alpha_f[idx] for idx in nc_idx) + sum(V[ni, t]*alpha_V[ni]./100 for ni in 1:NI)) * gamma^t for t in 1:T) +
            sum((So[col, t] .- So_sp)^2 for col in 1:Col for t in 1:T)*Q_so + 
            penalty_mult[3] * penalty * sum(Vslack[i, t] for i in tks for t in 1:T) +
            penalty_mult[3] * penalty * sum(Vsurplus[i, t] for i in tks for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vslack_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vsurplus_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[1] * penalty * sum(Ssslack[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[4] * penalty * sum(RASslack[t] for t in 1:T) +
            penalty_mult[4] * penalty * sum(RASsurplus[t] for t in 1:T)
        )
    end

    if obj_mpc == 15
        # Maximize  outflow, penalize volume deviations and substrate violations
        @objective(model, Min, 
            -sum((sum(f[idx, t] * alpha_f[idx] for idx in nc_idx) + sum(V[ni, t]*alpha_V[ni]./100 for ni in 1:NI)) * gamma^t for t in 1:T) +
            sum((So[col, t] .- So_sp)^2 for col in 1:Col for t in 1:T)*Q_so + 
            sum((RAS[t]-RAS_sp)^2 for t in 1:T)*Q_RAS +
            penalty_mult[3] * penalty * sum(Vslack[i, t] for i in tks for t in 1:T) +
            penalty_mult[3] * penalty * sum(Vsurplus[i, t] for i in tks for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vslack_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vsurplus_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[1] * penalty * sum(Ssslack[i, t] for i in 1:Col for t in 1:T) 
        )
    end


    return model

end