##############################
# CASE DESCRIPTION   
##############################
# Setpoint Change in Well Flow Rates - Control 4 Wells  
# Setpoint Change in the Mainline Pump Flow Rate - Control Pump Speed 
# No Pump Work Done Minimisation
# No Offset Free Methods  
# Include Input Regularisation for both the pump speed and the 4 Choke Valves 

# Increased Pump Speed Limit 
# Speeded up first order filter for bias term 
# Increaed penalty on pressure constraint violation and cavitation constraint violation to avoid these violations in the presence of larger disturbances from the setpoint changes
# L1 Penalty for Input Regularisation 
# Same as case 5 with min. work done for pump in objective 

# NMPC Soft Constraints 
# 1. Pipe Pressure Constraints at Nodes 2, 52, 54, 65, 67 
# 2. Choke Valve Cavitation Index Constraint 

##############################
# MODEL PARAMETERS  
##############################
# original_dir: Directory for file storage  
# CASE_NAME: Case Title 
##############################
original_dir = pwd()
CASE_NAME = "case_final"
case_dir = original_dir*"/"*CASE_NAME
lf_case_id = "fdd84235-e7a2-4651-b946-d45cd4deb638"

##############################
# OFFSET-FREE CONTROL SETTINGS
##############################
well_flow_bias_bool = true # whether to include bias term for well flow rate measurements
well_flow_bias_filter = 0.7 # tune in (0,1); smaller = slower
mainline_flow_bias_bool = false # whether to include bias term for mainline flow rate measurements
mainline_flow_bias_filter = 0.7 # tune in (0,1); smaller = slower
pressure_bias_bool = true # whether to include bias term for pressure measurements at constrained nodes
pressure_bias_filter = 0.7 # tune in (0,1); smaller = slower
VLV_2_CV_est_bool = false # whether to estimate the valve CV
inj_est_bool = false # whether to estimate well injectivity

##############################
# NETWORK PROPERTIES
##############################
# 1KM DISCRETISATION MODEL 
# PIPE SEGMENTS INDEXES 
pipe_seg = [2:40, 41:49, 54:62, 51, 53, 64, 66]
pipe_seg_inlet = [2, 41, 54, 51, 53, 64, 66]
pipe_seg_outlet = [40, 49, 62, 51, 53, 64, 66]

# NODE AND SEGMENT TYPES & LOCATIONS 
# K: Total number of arc segments 
# J: Total number of nodes including nodes at the outlet of the last pipe segments (excluding nodes at pipe inlet) 
# choke_vlv_idx: Arc indexes of choke valves 
# vlv_idx: Arc indexes of other control valves 
# pipe_idx: Arc indexes of pipelines 
# inlet_idx: Arc indexes of flow supply to pipelines 
# well_idx: Arc indexes of pipes connected to the reservoir
# node_idx: Indexes of nodes for which mass balance equations need to be written (i.e. no accummulation at these nodes)
# nearest_pipe_idx: Index for pipe segments closest to choke valves (so that density can be referenced for in the choke valve pressure drop equation)

K = 66 # total number of arc segments 
J = 67 # total number of nodes including nodes at the outlet of the last pipe segments 
choke_vlv_idx = [50, 52, 63, 65]
vlv_idx = []
all_vlv_idx = vcat(choke_vlv_idx, vlv_idx)
pipe_idx = vcat(pipe_seg...)
well_idx = [51, 53, 64, 66]
pump_idx = [1]
node_idx = vcat(1:J)
inlet_node = [1] # inlet node(s)
outlet_node = [52, 54, 65, 67]
nearest_pipe_idx = Dict(
    50 => 49, 
    52 => 49, 
    63 => 62, 
    65 => 62)

##############################
# SEGMENT PROPERTIES 
##############################
# Theta: deg, Angle at which pipe segment is tilted 
# HEIGHT: m, Height 
# LENGTH: m, Length 
# DINNER: m, Inner Diameter 
# ACROSS: m2, Cross Sectional Area 
# V: m3, Volume
Theta = vcat(zeros(K, 1))
Theta[pipe_seg[1]] .= -0.44 # Slight downward slope in the mainline 
Theta[well_idx] .= -90
LENGTH = vcat(zeros(K, 1))
LENGTH[setdiff(pipe_idx, well_idx)] .= 1000
LENGTH[well_idx] .= 1200 # Length of wellbore segments 
DINNER = vcat(zeros(K, 1))
DINNER[pipe_seg[1]] .= 0.4572
DINNER[pipe_seg[2]] .= 0.2984
DINNER[pipe_seg[3]] .= 0.2984
DINNER[well_idx] .= 0.1571
ACROSS = 3.142*(DINNER./2).^2  
V = ACROSS.*LENGTH 

##############################
# NETWORK STRUCTURE 
##############################
# Dictionaries to map flows incoming and outgoing from a nod 
# Keys are nodes and values are the segments with inflows (for incoming dictionary) or outflows (for outgoing dictionary)
incoming = Dict{Int, Vector{Int}}()
for i in setdiff(node_idx, inlet_node)
    incoming[i] = [i-1]
end 

outgoing = Dict{Int, Vector{Int}}()
for i in setdiff(node_idx, outlet_node)
    outgoing[i] =  [i]
end 

outgoing[41] = [41, 54]
outgoing[50] = [50, 52]
outgoing[63] = [63, 65]

##############################
# DEFINE MPC TIME HORIZON PARAMETERS  
# NP - Prediction horizon in seconds 
# NM - Control horizon in number of intervals 
# NT - Total simulation time in seconds 
# DT - Time step 
# T_sim - Number of time steps for simulation model 
# T_mpc - Number of time steps for MPC 
# T_ss - Number of  time steps for steady-state model 
##############################
NP = 3600
NM = 4
NT = 3600*20
DT = 300
T_sim = 1 
T_mpc = Int(NP/DT)
T_ss = 1 

##############################
# OPERATIONAL PROPERTIES 
##############################
# CHOKE_VALVE_OP_SS: Steady State Valve Opening for all valves 
# VLV_2_CV: Control Valve CV (Convert Valve Opening to Flow Coefficient) 
# P_IN: Inlet Pressure to Pipeline 
# W_IN: Inlet flow to pipeline 
# P_RES: Reservoir Pressure (bar)
# PARAMETERS RELATED TO TANK & PUMP
# tank_capacity: m3/hr, tank for CO2 storage before further transport and injection 
# pump_efficiency: Efficiency of pump to calculate work down 

CHOKE_VALVE_OP_SS = vcat(zeros(K, 1))
CHOKE_VALVE_OP_SS[all_vlv_idx] .= 1.0
VLV_2_CV = vcat(zeros(K, 1))
VLV_2_CV[all_vlv_idx] = [0.0034, 0.0034, 0.0034, 0.0034]
P_IN = 100
W_IN = 500
P_RES = 120 
pump_efficiency = vcat(zeros(K, 1))
pump_efficiency[pump_idx] .= 0.85
well_injectivity = vcat(zeros(K, 1))
well_injectivity[well_idx] .= 1.0

##############################
# MPC VARIABLE DEFINITIONS   
##############################
# u_vars - Define input variables (MPC will provide optimal input values for these)
u_vars = ["choke_vlv_op", "speed_pump"]
state_vars = ["p_average"]

##############################
# MPC INPUT VARIABLE BOUNDS
##############################
choke_vlv_op_min = 0.0
choke_vlv_op_max = 1.0
# maximum rate of change for choke valve opening per second (used for input regularisation constraint)  
choke_vlv_op_roc = 1.0*(10.0^(-4)) 
speed_pump_min = 2800
speed_pump_max = 3900
# rpm/s (used for input regularisation constraint)
speed_pump_roc = 0.03 

##############################
# MPC CONTROLLED VARIABLES  
##############################
# pipe_pressure_min or pipe_pressure_max: bar, min or max allowable CO2 pressure in the pipeline 
# cavitation_min
pressure_con_idx = [2, 52, 54, 65, 67]
pipe_pressure_max = vcat(zeros(J, 1))
pipe_pressure_max[pressure_con_idx] .= [253, 253, 253, 253, 253]
pipe_pressure_min = vcat(zeros(J, 1))
pipe_pressure_min[pressure_con_idx] .= [180, 180, 180, 180, 180]
cavitation_min = 1.7

##############################
# MPC SETPOINTS 
##############################
n_samples_setpoint = length(DT:DT:(NT+NP)) + 1

# WELL FLOW SETPOINTS 
well_1_flow_sp = fill(125, n_samples_setpoint, 1)
well_2_flow_sp = fill(125, n_samples_setpoint, 1)
well_3_flow_sp = fill(125, n_samples_setpoint, 1)
well_4_flow_sp = fill(125, n_samples_setpoint, 1)

well_1_flow_sp[25:61] .= 130
well_1_flow_sp[62:101] .= 130
well_1_flow_sp[102:end] .= 120

well_2_flow_sp[25:61] .= 130
well_2_flow_sp[62:101] .= 130
well_2_flow_sp[102:end] .= 120

well_3_flow_sp[25:61] .= 130
well_3_flow_sp[62:101] .= 130
well_3_flow_sp[102:end] .= 120

well_4_flow_sp[25:61] .= 130
well_4_flow_sp[62:101] .= 140
well_4_flow_sp[102:end] .= 120

##############################
# Initial Guess for Steady State Model 
##############################
steady_state_initial_guess = 
    Dict(
    "p_in" => 120,
    "p_out" => 100,
    "w_in" => 500, 
    "w_out" => 500,
    "dpdz_in" => -0.01,
    "dpdz_out" => -0.01,
    "Z" => 0.15,
    "Z_in" => 0.15,
    "Z_out" => 0.15,
    "rho" => 900,
    "rho_in" => 900,
    "rho_out" => 900,
    "reynold_in" =>1e7,
    "reynold_out" => 1e7,
    "fdarcy_in" => 0.01, 
    "fdarcy_out" => 0.01, 
    "choke_vlv_op" => 1.0,
    "deltap_wellhead_choke" => 14,
    "p_node" => 100,
    "p_average" => 100,
    "q_tank_in" => 0.20,
    "q_tank_out" => 0.20,
    "tank_co2_volume" => 50,
    "rho_pump_inlet" => 922,
    "rho_pump_outlet" => 979,
    "pump_work" => 0.8,
    "cavitation_idx" => 0,
    "pump_cost" => 52,
    "electricity_price" =>100,
    "rho_tank" => 850, 
    "Z_tank" => 0.15,
    "speed_pump" => 3440,
    "head_meters" => 2000
    )  

##############################
# MPC CONSTRAINTS FOR OPEN LOOP / CLOSED LOOP SIMULATION 
##############################

function special_constraints_fn_dynamic(model, sim_t)
    
    ##############################
    # DISTURBANCES 
    ##############################
    start = Int(sim_t/DT) # Start from initial condition 

    # SET POINT TO MPC 
    well_1_flow_sp_mpc = well_1_flow_sp[start:(start + T_mpc - 1)]
    well_2_flow_sp_mpc = well_2_flow_sp[start:(start + T_mpc - 1)]
    well_3_flow_sp_mpc = well_3_flow_sp[start:(start + T_mpc - 1)]
    well_4_flow_sp_mpc = well_4_flow_sp[start:(start + T_mpc - 1)]

    ##############################
    # OBJECTIVE FUNCTION TO MPC 
    ##############################
    @objective(model, Min,
    100*sum((1.0 - model[:choke_vlv_op][k, t]) for k in choke_vlv_idx for t=1:T_mpc) # Choke Valve Opening Setpoint Tracking
    + 5/3*sum((model[:w_in_corrected][51, d, t] - well_1_flow_sp_mpc[t])^2 for d=2:D for t=1:T_mpc) # Well Flow Rate Setpoint 
    + 20/3*sum((model[:w_in_corrected][53, d, t] - well_2_flow_sp_mpc[t])^2 for d=2:D for t=1:T_mpc) # Well Flow Rate Setpoint 
    + 80/3*sum((model[:w_in_corrected][64, d, t] - well_3_flow_sp_mpc[t])^2 for d=2:D for t=1:T_mpc) # Well Flow Rate Setpoint 
    + 140/3*sum((model[:w_in_corrected][66, d, t] - well_4_flow_sp_mpc[t])^2 for d=2:D for t=1:T_mpc) # Well Flow Rate Setpoint 
    + 3500*sum(model[:pipe_pressure_slack][k, t] for k in pressure_con_idx for t=1:T_mpc)    # Pipe Pressure CV Constraint 
    + 3500*sum(model[:pipe_pressure_surplus][k, t] for k in pressure_con_idx for t=1:T_mpc)  # Pipe Pressure CV Constraint 
    + 2000*sum(model[:cavitation_slack][k, t] for k in choke_vlv_idx for t=1:T_mpc)    # Choke Valve Cavitation Index
    + 50*sum(model[:s_choke][k, t] for k in choke_vlv_idx for  t=1:T_mpc) # Input regularisation for moving the choke valve
    + 0.5*sum(model[:s_pump][k, t] for k in pump_idx for t=1:T_mpc) # Input regularisation for moving the pump speed
    )

end 

special_constraints_dynamic = special_constraints_fn_dynamic
