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
export deltav
export deltav_binary
export V0
export f0
export fmax_vec
export Bd_vec
export qmaxin_vec
export qmaxout_vec

export min_variance_parameters
export equality_con_param
export input_reg_penalty

export mu_max, Ks, Kla, b, Y, So_sat
export Xb_in, Ss_in, So_in
export Xb0, Ss0, So0
export Ssmax
export Col, Adot
export Vmax_bio, Vmin_bio, Ssmax_bio
export dt, Col, Adot

export NI_BRX

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


NI = 3  # Number of inventory nodes (3 tanks)
NI_BRX = 2 # location of bioreactor
O = 4   
M = [1 -1 0  0; 0 1 -1 0; 0 0 1 -1]
unreachable_fsp = [20; 20; 20; 20]
P = [-1 0 0 0]
C = [0 0 0 1]
np = size(P)[1]
nc = size(C)[1]
alpha_f = ones(O) 
alpha_f[4] = 30 
alpha_V = ones(NI)
alpha_V[NI_BRX] = 10
alpha_V[NI] = 5
EI = ones(O)
special_constraints = nothing 
case_name = "3Tk_Bioreactor"


# MPC PARAMETERS 
nt = 10.0  # Total simulation time (days)
dt = 60/1440  # Time step in days (60 minutes)
npred = Int(1/dt)  # Prediction horizon
nm = 3  
penalty = 5000
penalty_mult = [1.0, 0.7, 0.5, 0.3]
obj_mpc = 9
optimizer = "Ipopt"


# BIOREACTOR-SPECIFIC PARAMETERS
mu_max = 6.0  # Maximum growth rate (day⁻¹)
Ks = 20.0  # Half-saturation constant (g COD/m³)
Kla = 120.0  # Oxygen transfer coefficient (day⁻¹)
b = 0.62  # Decay rate (day⁻¹)
Y = 0.67  # Yield coefficient (g cell COD / g COD oxidized)
So_sat = 10.0  # Oxygen saturation concentration (g O₂/m³)

# BIOREACTOR INLET
Xb_in = 0.0  # Inlet biomass (g COD/m³)
Ss_in = 200.0  # Inlet substrate (g COD/m³)
So_in = 2.0  # Inlet oxygen (g O₂/m³)

Ssmax = 10.0  # Maximum substrate (g COD/m³)

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

Vcap = [600; 1200; 800]  # Tank volume capacity (m³)
Vmax = Vcap .* 0.9  
Vmin = Vcap .* 0.7  
fmin = [0; 0; 0; 0; 0; 0]
deltaV = [0; 0; 0; 0] 
fmin_binary = 0  
deltav_binary = 0  

# Initial conditions
V0 = Vmax  # Initial volumes = max volumes
f0 = [100; 100; 100; 100]  # Initial flows
Xb0 = 17.2
So0 = 9.865
Ss0 = 2.7

# DISTURBANCE  2
fmax1 = [200; 150; 120; 100]  
fmax2 = [2; 100; 100; 100] 


# TIME-VARYING CONSTRAINT MATRICES  
fmax_vec = zeros(O, Int(nt/dt)+1)
Bd_vec = zeros(NI, Int(nt/dt)+1)

# Build time-varying constraint vectors
for i in 1:(Int(nt/dt)+1)

    if i <= 5
        fmax_vec[:, i] = fmax1
    elseif i <= 50
        fmax_vec[:, i] = fmax2
    else
        fmax_vec[:, i] = fmax1
    end
end