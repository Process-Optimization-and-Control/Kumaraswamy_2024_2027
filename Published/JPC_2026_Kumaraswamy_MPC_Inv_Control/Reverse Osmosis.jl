##############################
# NOTES
##############################
# 1. This is a Reverse Osmosis Case.
# 2. This script contain all the required parameters to run the MPC & dynamic model. 

##############################
# MODEL PARAMETERS  
##############################
# a - Tank Areas
# hcap - Full Tank Height 
# NI - Number of Tank Inventories 
# O - Number of Flows or "Edges" 
# unreachable_fsp - Unreachable Flow Setpoint 
# M - Connectivity Matrix
# P - Producer Nodes (Where are the producer nodes connected to?)
# C - Consumer Node Indexes (Where are the consumer nodes connected to?)
# alpha_f - Weights for Outflow to Consumers 
# EI - Edge Weights (Used to Calculate Distance from the Consumer)
# special_constraints - additional constraints to the model
# case_name - Name for Case


##############################
# DEFINE MODEL PARAMETERS 
##############################
a = [2.0, 2.0, 5.0, 0.1, 2.0, 0.1, 2.0]
hcap = [5.0; 5.0; 6.0; 0.1; 4.0; 0.1; 4.0]
NI = 7
O = 10 
unreachable_fsp = 100*ones(O)
M = [1 0 -1 0 0 0 0 0 0 0; 
     0 1 0 -1 0 0 0 0 0 0; 
     0 0 1 1 -1 0 0 0 0 1; 
     0 0 0 0 1 -1 0 0 0 0; 
     0 0 0 0 0 1 -1 0 0 0; 
     0 0 0 0 0 0 1 -1 -1 0; 
     0 0 0 0 0 0 0 0 1 -1] 
println(size(M))
P = zeros(2, O)
P[1, 1] = -1 
P[2, 2] = -1
C = [0 0 0 0 0 0 0 1 0 0]
alpha_f = ones(O) 
alpha_f[8] = 1.5 
EI = ones(O)
EI[10] = 1.8
EI[3] = 1.2

# Flow Ratio Constraints
function flow_ratio_constraints(model, vars, T)
    @constraint(model, [t in 1:T], vars.f[6, t] == 0.90*vars.f[5, t])
    @constraint(model, [t in 1:T], vars.f[9, t] == 0.40*vars.f[7, t])
    @constraint(model, [t in 1:T], vars.f[8, t] == 0.60*vars.f[7, t])

end
special_constraints = flow_ratio_constraints
case_name = "Reverse Osmosis"

##############################
# MPC PARAMETERS 
##############################
np = 10
nm = 3
nt = 80
penalty = 1000
penalty_mult = [1.0, 0.5, 0.3]
obj_mpc = 7
optimizer = "HiGHS"

##############################
# MPC Parameters for Objectives 
# Relevant only for specific objectives 
##############################
equality_con_param=([], 0)
input_reg_penalty=1e-3

##############################
# FIXED CONSTRAINTS / DISTURBANCE 
# NOT VARYING WITH TIME 
##############################
hmax = hcap*0.9
hmin = hcap*0.1
fmin = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
deltah = [0; 0; 0; 0; 0; 0; 0]
fmin_binary = 0
deltah_binary = 0
h0 = hmax
f0 = [1.8; 21.9; 1.8; 21.9; 37.033; 33.33; 33.33; 20.0; 13.33; 13.33]

##############################
# TIME-VARYING CONSTRAINTS / DISTURBANCE  
##############################
Bd1 = [0; 0; 0; 0; 0; 0; 0]
fmax1 = [1.8; 25.0; 1.8; 25.0; 45.0; 40.5; 40.5; 20.0; 16.2; 20.0]
fmax2 = [1.8; 20.0; 1.8; 25.0; 45.0; 40.5; 40.5; 20.0; 16.2; 20.0]

##############################
# FMAX AND BD MATRICES  
##############################
fmax_vec = zeros(O, nt)
Bd_vec = zeros(NI, nt)

for i in 1:nt
    if i <= 10
        fmax_vec[:, i] = fmax1;
    end
    
    if i > 10 && i <= 60
        fmax_vec[:, i] = fmax2;
    end

    if i > 60 
        fmax_vec[:, i] = fmax1;
    end

    Bd_vec[:, i] = Bd1;
    
end
