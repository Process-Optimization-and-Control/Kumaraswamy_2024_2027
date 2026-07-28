export a 
export hcap
export NI
export O
export unreachable_fsp
export M
export P
export C
export alpha_f
export EI
export special_constraints
export case_name            

export npred
export nm
export nt       
export obj_mpc
export penalty
export penalty_mult
export optimizer    

export hmax
export hmin 
export fmin
export fmin_binary
export deltah
export deltah_binary
export h0
export f0
export fmax_vec
export Bd_vec
export qmaxin_vec
export qmaxout_vec

export min_variance_parameters
export equality_con_param
export input_reg_penalty

using Distributions

##############################
# NOTES
##############################
# 1. This is a simple 3 Tank Case.
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


a = [1.0, 1.5, 2.0]
hcap = [2.3; 2.8; 3.2]
NI = 3
O = 4
unreachable_fsp = [20; 20; 20; 20]
M = [1 -1 0  0;0 1 -1 0; 0 0 1 -1]
P = [-1 0 0 0]
C = [0 0 0 1]
np = size(P)[1]
nc = size(C)[1]
alpha_f = ones(O) 
alpha_f[4] = 1.5 
EI = ones(O)
special_constraints = nothing 
case_name = "1Tk_Bioreactor_copy"

##############################
# MPC PARAMETERS 
##############################
nt = 10.0 #days
dt = 60/1440 #120/1440 #step time in days
npred = Int(1/dt)
nm = 3
penalty = 5000
penalty_mult = [1.0, 0.5, 0.3]
obj_mpc = 9
optimizer = "Ipopt"

##############################
# MPC Parameters for Objectives 
# Relevant only for specific objectives 
##############################
min_variance_parameters=(0, [], 0)
equality_con_param=([], 0)
input_reg_penalty=1e-3

##############################
# FIXED CONSTRAINTS / DISTURBANCE 
# NOT VARYING WITH TIME 
##############################
hmax = hcap*0.9
hmin = hcap*0.1
fmin = [0; 0; 0; 0; 0; 0]
deltah = [0; 0; 0; 0]
fmin_binary = 0
deltah_binary = 0
h0 = hmax
f0 = [1; 1; 1; 1]

##############################
# TIME-VARYING CONSTRAINTS / DISTURBANCE  
##############################
Bd1 = [0; 0; 0]
Bd2 = [0; 0; -0.05/a[3]]
fmax1 = [1.667; 1.428; 1.125; 1.0]
fmax2 = [0.833; 1.428; 1.125; 1.0]

q_inmax1 = 150.0
q_inmax2 = 5.0

##############################
# FMAX AND BD MATRICES  
##############################
qmaxin_vec = zeros(Int(nt/dt)+1)
qmaxout_vec = ones(Int(nt/dt)+1).*100.0
Bd_vec = zeros(NI, Int(nt))


for i in 1:(Int(nt/dt)+1)

    if i <= 70
        qmaxin_vec[i] = q_inmax1;
    end
    
    if i > 70 && i <= 90
        qmaxin_vec[i] = q_inmax2;
    end

    if i > 90
        qmaxin_vec[i] = q_inmax1;
    end
   
end