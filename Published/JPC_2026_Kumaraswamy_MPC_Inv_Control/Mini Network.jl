##############################
# NOTES
##############################
# 1. This is a Mini Network Case.
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
a = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
hcap = [2.3; 2.3; 2.3; 2.3; 2.3; 2.3; 2.3; 2.3]
NI = 8
O = 13
unreachable_fsp = 20*ones(O)
M = zeros(NI, O)
M[1, 1] = 1 
M[1, 3] = -1 
M[1, 13] = 1
M[2, 2] = 1 
M[2, 4] = -1
M[3, 3] = 1
M[3, 4] = 1
M[3, 5] = -1 
M[3, 6] = -1
M[3, 7] = -1
M[4, 5] = 1
M[4, 8] = -1 
M[5, 6] = 1
M[5, 9] = -1 
M[6, 7] = 1 
M[6, 10] = -1 
M[7, 8] = 1 
M[7, 9] = 1
M[7, 11] = -1
M[7, 12] = -1 
M[8, 12] = 1 
M[8, 13] = -1
P = zeros(2, O)
P[1, 1] = -1 
P[2, 2] = -1
C = zeros(2, O)
C[1, 10] = 1
C[2, 11] = 1
alpha_f = ones(O) 
alpha_f[10] = 1.5
alpha_f[11] = 1.2  
EI = ones(O)
EI[13] = 1.2
EI[3] = 1.2
EI[8] = 1.1
EI[9] = 1.2
special_constraints = nothing
case_name = "Mini Network Scenario"

##############################
# MPC PARAMETERS 
##############################
np = 10
nm = 3
nt = 120    
penalty = 10
penalty_mult = [1.0, 0.5, 0.3]
obj_mpc = 4
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
fmin = zeros(O)
fmin[3] = 0.35
deltah = 0.1*hmax
fmin_binary = 1
deltah_binary = 1
h0 = hmax
f0 = [0.4; 0.4; 0.4; 0.4; 0.25; 0.25; 0.3; 0.25; 0.25; 0.3; 0.5; 0.0; 0.0]

##############################
# TIME-VARYING CONSTRAINTS / DISTURBANCE  
##############################
Bd1 = zeros(NI)
fmax1 = [1.667; 1.667; 1.428; 1.428; 0.375; 0.375; 0.375; 0.3; 0.3; 0.3; 0.5; 1.0; 1.0]
fmax2 = [0.3; 0.3; 1.428; 1.428; 0.375; 0.375; 0.375; 0.3; 0.3; 0.3; 0.5; 1.0; 1.0]

##############################
# FMAX AND BD MATRICES  
##############################
fmax_vec = zeros(O, nt)
Bd_vec = zeros(NI, nt)

for i in 1:nt
    if i <= 10
        fmax_vec[:, i] = fmax1;
    end
    
    if i > 10 && i <= 100
        fmax_vec[:, i] = fmax2;
    end

    if i > 100 
        fmax_vec[:, i] = fmax1;
    end

    Bd_vec[:, i] = Bd1;
end
