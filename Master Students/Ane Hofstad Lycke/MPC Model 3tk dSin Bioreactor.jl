using JuMP
using Ipopt
using LinearAlgebra
using HiGHS

function create_mpc_model(alpha_tk)
    

    # CALCULATE A AND B MATRICES  
    A = I(NI)
    B = M  # No division by area for volume-based model
    

    # MPC PARAMETERS  
    T = npred  # Prediction horizon
    gamma = 0.9  # Time discount factor


    # INITIAL CONDITION AS INTIAL QUESS FOR OPTIMIZER
    fguess = fill(f0, T)
    Vguess = fill(V0, T)


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
    @variable(model, Ss_in[t=1:T] >= 0)

    @variable(model, q_in[t=1:T] >= 0)
    @variable(model, q_out[t=1:T] >= 0)

    @variable(model, V_BRX[1:Col, 1:T])
    @variable(model, Xb[1:Col, 1:T])
    @variable(model, Ss[1:Col, 1:T])
    @variable(model, So[1:Col, 1:T])
    #@variable(model, Kla[1:Col, 1:T])
    @variable(model, Kla[t=1:T] >= 0)
    @variable(model, D[1:T])

    @variable(model, Vslack[ni=1:NI, t=1:T]>=0)
    @variable(model, Vsurplus[ni=1:NI, t=1:T]>=0)
    @variable(model, fsurplus[o=1:O, t=1:T]>=0) 
    @variable(model, Vslack_BRX[1:Col, 1:T] >= 0)
    @variable(model, Vsurplus_BRX[1:Col, 1:T] >= 0)
    @variable(model, Ssslack[1:Col, 1:T] >= 0)

#=
    Set V2 = Vbrx[1, t] 
    Set F2 = qin 
    Set F3 = qout 
=#
    @constraint(model, [t in 1:T], V[2, t] == V_BRX[1, t])
    @constraint(model, [t in 1:T], q_in[t] == f[2, t])
    @constraint(model, [t in 1:T], q_out[t] == f[3, t])


    # DISTURBANCES
    @variable(model, fmax[o=1:O] >= 0)
    @variable(model, Bd[ni=1:NI]>=0)
    @variable(model, Ss_in_max >= 0)


    # INPUT REGULARIZATION
    @variable(model, reg_slack1[o=1:O, t=2:T] >= 0)
    @variable(model, reg_slack2[o=1:O, t=2:T] >= 0)


    # MODEL EQUATIONS  
    tks = [i for i in 1:NI if i != NI_BRX]
    @constraint(model, [t in 1:T-1], V[tks, t+1] .- A[tks, tks]*V[tks, t] .- dt.*B[tks, :]*f[:, t] .- dt.*Bd[tks] .== 0)
    
    # Dilution rate
    @constraint(model, [t in 1:T], D[t] * V_BRX[1, t] == q_in[t])

    #Koh = 0.20 #*****************************************

    # Collocation model equation for states 
    @constraint(model, [t=1:T], Adot'*V_BRX[1:Col, t] .== dt * (q_in[t] .- q_out[t]))
    @constraint(model, [t=1:T], Adot'*Xb[1:Col, t] .== dt * (D[t] .* ((Xb_in) .- Xb[2:Col, t]) .+ (mu_max .* Ss[2:Col, t] ./ (Ks .+ Ss[2:Col, t]) ).* Xb[2:Col, t] .- b .* Xb[2:Col, t]))
    @constraint(model, [t=1:T], Adot'*Ss[1:Col, t] .== dt * (D[t] .* (Ss_in[t] .- Ss[2:Col, t]) .- (1 / Y) .* (mu_max .* Ss[2:Col, t] ./ (Ks .+ Ss[2:Col, t]) ) .* Xb[2:Col, t]))
    @constraint(model, [t=1:T], Adot'*So[1:Col, t] .== dt * (D[t] .* (So_in .- So[2:Col, t]) .+ Kla[t] .* (So_sat .- So[2:Col, t]) .- ((1 - Y)/Y) .* (mu_max .* Ss[2:Col, t] ./ (Ks .+ Ss[2:Col, t]) ) .* Xb[2:Col, t] .- b .* Xb[2:Col, t])) 
    #.* (So[2:Col, t])./ (Koh .+ So[2:Col, t]) 

    # Collocation continuity: end of previous element = start of next
    @constraint(model, [t=2:T], V_BRX[1, t] == V_BRX[Col, t-1])
    @constraint(model, [t=2:T], Xb[1, t] == Xb[Col, t-1])
    @constraint(model, [t=2:T], Ss[1, t] == Ss[Col, t-1])
    @constraint(model, [t=2:T], So[1, t] == So[Col, t-1])

    # CONSTRAINTS
    @constraint(model, [t in 2:T], V[tks, t] .<= Vmax[tks] .+ Vslack[tks, t])
    @constraint(model, [t in 2:T], V[tks, t] .>= Vmin[tks] .- Vsurplus[tks, t])

    @constraint(model, [t in 2:T], V_BRX[:, t] .<= Vmax[NI_BRX] .+ Vslack_BRX[:, t])
    @constraint(model, [t in 2:T], V_BRX[:, t] .>= Vmin[NI_BRX] .- Vsurplus_BRX[:, t])

    @constraint(model, [t in 2:T], Ss[:, t] .<= Ssmax .+ Ssslack[:, t])
    @constraint(model, [o=1:O, t=1:T], f[o, t] <= fmax[o])
    @constraint(model, [t in 1:T], Ss_in[t] == Ss_in_max)

    @constraint(model, [t=1:T], Kla[t] <= Kla_max) 
    @constraint(model, [t=1:T], Kla[t] >= Kla_min) 


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
            sum(sum((So[:, t] .- So_sp).^2) for t in 1:T)*Q_so + 
            penalty_mult[3] * penalty * sum(Vslack[i, t] for i in tks for t in 1:T) +
            penalty_mult[3] * penalty * sum(Vsurplus[i, t] for i in tks for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vslack_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[2] * penalty * sum(Vsurplus_BRX[i, t] for i in 1:Col for t in 1:T) +
            penalty_mult[4] * penalty * sum(fsurplus[o, t] for o in 1:O for t in 1:T) +
            penalty_mult[1] * penalty * sum(Ssslack[i, t] for i in 1:Col for t in 1:T)
        )
    end

    return model
end

                           
