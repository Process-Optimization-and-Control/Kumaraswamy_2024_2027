##############################
# NOTES
##############################
# 1. This is a 5 Tank Split Join Case.
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
a = [2.0, 1.0, 1.0, 1.0, 2.0]
hcap = [3.2; 2.3; 2.3; 2.3; 3.2]
NI = 5
O = 8 
unreachable_fsp = [20; 20; 20; 20; 20; 20; 20; 20]
M = [1 -1 -1 -1 0 0 0 0; 
    0 1 0 0 -1 0 0 0; 
    0 0 1 0 0 -1 0 0;
    0 0 0 1 0 0 -1 0;
    0 0 0 0 1 1 1 -1;]
P = [-1 0 0 0 0 0 0 0]
C = [0 0 0 0 0 0 0 1]
alpha_f = ones(O) 
alpha_f[8] = 1.5
EI = ones(O)
EI[5] = 1.1
EI[6] = 1.2
EI[7] = 1.3
special_constraints = nothing
case_name = "5Tk_Split_Join_Case2"

##############################
# MPC PARAMETERS 
##############################
np = 10
nm = 3
nt = 80
penalty = 10
penalty_mult = [1.0, 0.5, 0.3]
obj_mpc = 5
optimizer = "HiGHS"

##############################
# MPC Parameters for Objectives 
# Relevant only for specific objectives 
##############################
equality_con_param=([5,6,7], 2)
input_reg_penalty=1e-3

##############################
# FIXED CONSTRAINTS / DISTURBANCE 
# NOT VARYING WITH TIME 
##############################
hmax = hcap*0.9
hmin = hcap*0.1
fmin = [0; 0; 0; 0; 0; 0; 0; 0]
deltah = [0; 0; 0; 0; 0]
fmin_binary = 0
deltah_binary = 0
h0 = hmax
f0 = [0.75; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.75]

##############################
# TIME-VARYING CONSTRAINTS / DISTURBANCE  
##############################
Bd1 = [0; 0; 0; 0; 0]
fmax1 = [1.428; 0.375; 0.375; 0.375; 0.3; 0.3; 0.3; 0.75]
fmax2 = [0.5; 0.375; 0.375; 0.375; 0.3; 0.3; 0.3; 0.75]

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
