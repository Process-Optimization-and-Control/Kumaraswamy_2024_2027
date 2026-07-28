export a 
export vcap
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

export Vmax
export Vmin 
export fmin
export fmin_binary
export deltaV
export deltav_binary
export V0
export f0
export Xb_flow0, Ss_flow0, Xb_tk0, Ss_tk0
export fmax_vec
export Bd_vec
export Ss_in_vec
export qmaxin_vec
export qmaxout_vec

export min_variance_parameters
export equality_con_param
export input_reg_penalty

export recycle_idx, WAS_idx
export mu_max, Ks, b, Y, So_sat, So_sp, Q_so, Q_RAS
export Xb_in, Ss_in, So_in
export Kla0, RAS0
export Ssmax, Kla_min, Kla_max
export RAS_min, RAS_max
export dt, Col, Adot

export NI_BRX, RAS

using Distributions


##############################
# MODEL PARAMETERS - SIMPLE TANKS
##############################
# a - Tank Cross-sectional Areas
# hcap - Full Tank Heights (converted to volumes)
# vcap - Tank Volumes (hcap * a)
# NI - Number of Tank Inventories (3 tanks)
# O - Number of Flows (4 edges: inlet, tank1→2, tank2→3, outlet)
# M - Connectivity Matrix (for volume dynamics: dV/dt = M*f)
# P - Producer Nodes
# C - Consumer Node Indexes
# alpha_f - Weights for Outflow to Consumers
# case_name - Case identifier


NI = 5  # Number of inventory nodes (5 tanks)
NI_BRX = 2 # location of bioreactor
recycle_idx = 7 
WAS_idx = 8 
O = 8 
M = [1 -1 0 0 0 0 1 0; 
     0 1 -1 0 0 0 0 0; 
     0 0 1 -1 0 -1 0 0; 
     0 0 0 1 -1 0 0 0; 
     0 0 0 0 0 1 -1 -1]


unreachable_fsp = [20; 20; 20; 20]
P = [-1 0 0 0 0 0 0 0]
C = [0 0 0 0 1 0 0 1]
np = size(P)[1]
nc = size(C)[1]
alpha_f = ones(O) 
alpha_f[5] = 60
alpha_f[7] = 30
alpha_V = ones(NI)
alpha_V[NI_BRX] = 100
alpha_V[4] = 20
alpha_V[5] = 10
EI = ones(O)
special_constraints = nothing 
case_name = "3Tk_Recycle_WAS"


# MPC PARAMETERS 
nt = 10  # Total simulation time (days)
dt = 60/1440  # Time step in days (60 minutes)
npred = Int(1/dt)  # Prediction horizon
nm = 3 
penalty = 5000
penalty_mult = [1.0, 0.7, 0.5, 0.3]
obj_mpc = 15
optimizer = "Ipopt"


# BIOREACTOR PARAMETERS
mu_max = 6.0  # Maximum growth rate (1/day)
Ks = 20.0  # Half-saturation constant (g COD/m^3)
b = 0.62  # Decay rate (1/day)
Y = 0.67  # Yield coefficient (g cell COD / g COD oxidized)
So_sat = 10.0  # Oxygen saturation concentration (g O2/m^3)
So_sp = 3.0  # Setpoint for oxygen (g O2/m^3)
Q_so = 3.0  # Penalty for oxygen deviation

Ssmax = 15.0  # Maximum substrate (g COD/m^3)
Kla_min = 0.0  # Minimum oxygen transfer coefficient (1/day)
Kla_max = 360.0  # Maximum oxygen transfer coefficient (1/day)
#RAS_min = 1.25  # Minimum recycle ratio
#RAS_max = 1.25  # Maximum recycle ratio
RAS_sp = 1.25  # Setpoint for recycle ratio
Q_RAS = 2  # Penalty for RAS input regularization

# Collocation method parameters
Col = 4
Adot = [-4.139387691339813   1.739387691339811  -3.000000000000002
        3.224744871391587  -3.567840084690404   5.531972647421805
        1.167840084690405   0.775255128608409  -7.531972647421807
        -0.253197264742181   1.053197264742181   5.000000000000000]


# MPC Parameters for Objectives 
min_variance_parameters = (0, [], 0)
equality_con_param = ([], 0)
input_reg_penalty = 1e-3

Vcap = [600; 1200; 800; 400; 400]  # Tank volume capacity (m³)
Vmax = Vcap .* 0.9  
Vmin = Vcap .* 0.7  
fmin = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
deltaV = [0; 0; 0; 0; 0] 
fmin_binary = 0  
deltav_binary = 0  

# BIOREACTOR INLET
Xb_in = 0.0  # Inlet biomass (g COD/m^3)
So_in = 2.0  # Inlet oxygen (g O2/m^3)

#INITIAL STEADY STATE CONDITIONS
Xb0 = 8.82 #12.37 
So0 = 3.0
Ss0 = 2.6 #2.78 


# Initial conditions/Initial guesses for optimization
Kla0 = 1.23 #2
RAS0 = 1.25
V0 = Vmax  # Initial volumes = max volumes
f0 = [67; 150; 150; 50; 50; 100; 84; 17] #[75; 150; 150; 50; 50; 100; 75; 25] 
Xb_tk0 = [4.41; Xb0; Xb0; 0.88; 7.94] #[3.71; 12.37; 12.37; 1.23; 11.13] 
Ss_tk0 = [67.96; Ss0; Ss0; 0.26; 2.34] #[100; 2.78; 2.78; 0.28; 2.5] 
Xb_flow0 = [Xb_in; 4.41; Xb0*0.1; Xb0*0.1; Xb0*0.9; 7.94; 7.94; 7.94] #[0; 18.357; Xb0; Xb0*0.1; Xb0*0.1; Xb0*0.9; Xb0*0.9; Xb0*0.9] 
Ss_flow0 = [150; 67.96; Ss0; Ss0; Ss0; 2.34; 2.34; 2.34] #[150; 100; Ss0; Ss0*0.1; Ss0*0.1; Ss0*0.9; Ss0*0.9; Ss0*0.9]


# DISTURBANCE  
fmax1 = [200; 180; 180; 60; 50; 120; 300; 100] 
fmax2 = [2; 180; 180; 60; 50; 120; 300; 100]
Ss_in_max1 = 150.0
Ss_in_max2 = Ss_in_max1 #1500.0


# TIME-VARYING CONSTRAINT MATRICES  
fmax_vec = zeros(O, Int(nt/dt)+1)
Bd_vec = zeros(NI, Int(nt/dt)+1)
Ss_in_vec = zeros(Int(nt/dt)+1)

# Build time-varying constraint vectors
for i in 1:(Int(nt/dt)+1)

    if i <= 20
        fmax_vec[:, i] = fmax1
    elseif i <= 50
        fmax_vec[:, i] = fmax2
    else
        fmax_vec[:, i] = fmax1
    end

    if i <= 80
        Ss_in_vec[i] = Ss_in_max1
    elseif i <= 120
        Ss_in_vec[i] = Ss_in_max2
    else
        Ss_in_vec[i] = Ss_in_max1
    end
end