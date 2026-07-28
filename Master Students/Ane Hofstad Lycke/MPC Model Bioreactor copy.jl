using JuMP
using Ipopt
using LinearAlgebra
using HiGHS # Alternative Optimizer to IPOPT 



function create_mpc_model(alpha_tk)
    
    ##############################
    # CALCULATE A AND B MATRICES  
    ##############################
    A = I(NI) #identity matrix of size NIxNI 
    B = M./a  # inflow outflow matrix divided by area of tanks

    ##############################
    # MPC PARAMETERS  
    ##############################
    # Prediction Horizon 
    T = npred

    # Time Discount Factor 
    gamma = 0.98


    # --- Parameters Bioreactor ---
    mu_max = 6.0  #day^-1
    Ks = 20.0 #g COD m^-3
    b = 0.62 #day^-1
    Y = 0.67 #g cell COD formed (g COD oxidizedx)^-1
    So_sat = 10.0 #g O2 m^-3

    # --- Influent ---
    Xb_in = 0.0
    Ss_in = 200.0
    So_in = 2.0

    Kla = 120.0 # d^-1
    Vmax = 1000
    Vmin = 970
    Ssmax = 10.0


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



    #Adot matrix for 4 collocation points (3rd order Gauss- Radau) from lagrange polynomial 
    #x(0) is known. K = 3. the start point in the next step = the interpolated solution of end point in previous time step
    Adot =
    [-4.139387691339813   1.739387691339811  -3.000000000000002
    3.224744871391587  -3.567840084690404   5.531972647421805
    1.167840084690405   0.775255128608409  -7.531972647421807
    -0.253197264742181   1.053197264742181   5.000000000000000]    
    # Collocation points BIOREACTOR
    Col = 4
    # State variables BIOREACTOR 
    @variable(model, V[1:Col, 1:T])
    @variable(model, Xb[1:Col, 1:T])
    @variable(model, Ss[1:Col, 1:T])
    @variable(model, So[1:Col, 1:T])

    @variable(model, Vslack[1:Col,1:T]>=0)
    @variable(model, Vsurplus[1:Col,1:T]>=0)
    @variable(model, Ssslack[1:Col,1:T]>=0)

    # Inputs BIOREACTOR
    @variable(model, q_in[1:T] >= 0)
    @variable(model, q_out[1:T] >= 0)
    # Dilution rate D = q_in/V 
    @variable(model, D[1:T])

    # Disturbances 
    @variable(model, qmaxin[1:T]>=0)
    @variable(model, qmaxout[1:T]>=0)

    # Input Regularization 
    #@variable(model, reg_slack1[o=1:O, t=2:T]>=0)
    #@variable(model, reg_slack2[o=1:O, t=2:T]>=0)

    # Model Equations  

    # Collocation Dynamic Constraints  BIOREACTOR
    @constraint(model, [t=1:T], Adot'*V[1:Col, t] .== dt .* (q_in[t] .- q_out[t])) 
    # Dilution rate BIOREACTOR
    @constraint(model, [t=1:T], D[t] * V[1, t] == q_in[t])

    @constraint(model, [t=1:T], Adot'*Xb[1:Col, t] .== dt * (D[t] .* (Xb_in .- Xb[2:Col, t]) .+ ( mu_max .* Ss[2:Col, t] ./ (Ks .+ Ss[2:Col, t]) ) .* Xb[2:Col, t] .- b .* Xb[2:Col, t]))
    @constraint(model, [t=1:T], Adot'*Ss[1:Col, t] .== dt * (D[t] .* (Ss_in .- Ss[2:Col, t]) .- (1 / Y) .* (mu_max .* Ss[2:Col, t] ./ (Ks .+ Ss[2:Col, t])) .* Xb[2:Col, t]))
    @constraint(model, [t=1:T], Adot'*So[1:Col, t] .== dt * (D[t] .* (So_in .- So[2:Col, t]) .+ Kla .* (So_sat .- So[2:Col, t]) .- ((1 - Y)/Y) .* (mu_max .* Ss[2:Col, t] ./ (Ks .+ Ss[2:Col, t])) .* Xb[2:Col, t] .- b .* Xb[2:Col, t])) 
    # first point in next finite element equals last step in previous time step BIOREACTOR
    @constraint(model, [t=2:T], V[1, t] == V[Col, t-1])
    @constraint(model, [t=2:T], Xb[1, t] == Xb[Col, t-1])
    @constraint(model, [t=2:T], Ss[1, t] == Ss[Col, t-1])
    @constraint(model, [t=2:T], So[1, t] == So[Col, t-1])

   
    @constraint(model, [t in 2:T], V[:, t] <= Vmax .+ Vslack[:, t])
    @constraint(model, [t in 2:T], V[:, t] >= Vmin .- Vsurplus[:, t])
    @constraint(model, [t in 2:T], Ss[:, t] <= Ssmax .+ Ssslack[:, t])

    # Max Flow Rate (Hard Constraint)
    @constraint(model, [t in 1:T], q_in[t] <= qmaxin[t])
    @constraint(model, [t in 1:T], q_out[t] <= qmaxout[t])

    # Max Height Change Constraint
    if deltah_binary == 1
        @constraint(model, [i in 1:NI, t in 2:T], (h[i, t] - h[i, t-1]) <= deltah[i] + deltahslack1[i, t]) 
        @constraint(model, [i in 1:NI, t in 2:T], -(h[i, t] - h[i, t-1]) <= deltah[i] + deltahslack2[i, t]) 
    end

    # Min Flow Rate (Soft Constraint) 
    if fmin_binary == 1
        @constraint(model, [i in 1:O, t in 2:T], f[i, t] + fsurplus[i, t] >= fmin[i])
    end 


    # Input Regularization  
   # @constraint(model, [o in 1:O, t in 2:T], (f[o, t] - f[o, t-1]) <= reg_slack1[o, t])
    #@constraint(model, [o in 1:O, t in 2:T], -(f[o, t] - f[o, t-1]) <= reg_slack2[o, t])

    ##############################
    # DEFINE OBJECTIVE FUNCTION 
    ##############################

    if obj_mpc == 9
        @objective(model, Min, -sum((q_out[t]*5 + V[1, t]/100)*gamma^t for t in 1:T) 
                            + penalty_mult[2]*penalty*sum(Vslack[1,t] for t in 1:T) 
                            + penalty_mult[2]*penalty*sum(Vsurplus[1,t] for t in 1:T) 
                            + penalty_mult[3]*penalty*sum(Ssslack[1,t] for t in 1:T)
        )
    end
    return model
end 