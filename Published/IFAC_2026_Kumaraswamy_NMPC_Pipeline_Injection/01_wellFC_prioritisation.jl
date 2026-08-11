##############################
# CASE DESCRIPTION   
##############################
# NMPC for a CO2 pipeline-injection system with 4 wells and a mainline pump 
# NMPC Objective: 
# 1. Minimise deviation of well flow rates from setpoints
# 2. Soft Constraints - Pipeline pressure, pipeline rate of change, and choke valve cavitation index
# 3. Input regularisation - Pump speed and Choke Valves 


##############################
# CASE IDENTIFIERS AND PATHS
##############################
# original_dir : working directory the notebook was launched from; all outputs are written
#                relative to it, and plot.jl cd()s back to it when it is done
# CASE_NAME    : case title, also the name of the results sub-folder
# case_dir     : results folder for this case
# mpc_pred_dir : sub-folder of case_dir holding the NMPC prediction dump for every sample
##############################
original_dir = pwd()
CASE_NAME = "01_wellFC_prioritisation"
case_dir = original_dir*"/"*CASE_NAME
mpc_pred_dir = joinpath(case_dir, "mpc_pred")

##############################
# NETWORK PROPERTIES
##############################
# The network is an arc-node graph. 
#
#   arc  1      pump
#   arcs 2:40   mainline, 39 km
#   node 41     mainline splits into Line 1 (arc 41) and Line 2 (arc 54)
#   arcs 41:49  Line 1, 9 km        
#   arcs 54:62  Line 2, 9 km
#   node 50     Line 1 splits to chokes 1 (arc 50) and 2 (arc 52)
#   node 63     Line 2 splits to chokes 3 (arc 63) and 4 (arc 65)
#   arcs 51, 53, 64, 66   wellbores 1-4, 1.2 km each, ending at the outlet nodes 52, 54, 65, 67
#
# PIPE SEGMENT INDEXES, GROUPED AS [mainline, Line 1, Line 2, wellbore 1, 2, 3, 4]
pipe_seg = [2:40, 41:49, 54:62, 51, 53, 64, 66]
pipe_seg_inlet = [2, 41, 54, 51, 53, 64, 66]   # first arc of each group (used for plotting)
pipe_seg_outlet = [40, 49, 62, 51, 53, 64, 66] # last arc of each group
pipe_seg_names = ["Mainline" "Line 1" "Line 2" "Well 1" "Well 2" "Well 3" "Well 4"]

# NODE AND SEGMENT TYPES & LOCATIONS
# K: Total number of arc segments
# J: Total number of nodes
# choke_vlv_idx: Arc indexes of choke valves
# vlv_idx: Arc indexes of other control valves (none in this network)
# all_vlv_idx: All valve arcs, i.e. choke_vlv_idx together with vlv_idx
# pipe_idx: Arc indexes of pipelines (mainline, both lines and all wellbores)
# pump_idx: Arc indexes of pumps
# well_idx: Arc indexes of the wellbores, i.e. the pipes connected to the reservoir
# node_idx: All node indexes.
# inlet_node: Node upstream of the pump, where the supply pressure P_IN is imposed
# outlet_node: Bottom-hole nodes, where the reservoir inflow relation closes the network
# nearest_pipe_idx: Index of the pipe segment immediately upstream of each choke valve, whose
#           outlet density is used in the choke valve pressure drop equation

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
# Theta: deg, Angle at which the pipe segment is tilted; negative means the segment runs downhill
# LENGTH: m, Length
# DINNER: m, Inner Diameter
# ACROSS: m2, Cross Sectional Area
# V: m3, Volume
# Entries are left at zero for arcs that are not pipes (the pump and the choke valves), which
# have no momentum or accumulation equations written for them.

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
# Dictionaries mapping the flows incoming to and outgoing from a node.
# Keys are nodes; values are the arcs feeding the node (incoming) or fed by it (outgoing).
# The default is a chain: arc i-1 enters node i and arc i leaves it. The three overrides below
# are the branch points, where one arc feeds two.
incoming = Dict{Int, Vector{Int}}()
for i in setdiff(node_idx, inlet_node)
    incoming[i] = [i-1]
end 

outgoing = Dict{Int, Vector{Int}}()
for i in setdiff(node_idx, outlet_node)
    outgoing[i] =  [i]
end 

outgoing[41] = [41, 54] # mainline splits into Line 1 and Line 2
outgoing[50] = [50, 52] # Line 1 splits into chokes 1 and 2
outgoing[63] = [63, 65] # Line 2 splits into chokes 3 and 4

##############################
# DEFINE MPC TIME HORIZON PARAMETERS  
# NP    - Prediction horizon in seconds
# NT    - Total closed-loop simulation time in seconds, excluding the steady-state initialisation
# DT    - s, sample time. Also the interval the simulation model is advanced by on each
#         closed-loop step
# T_sim - Number of time steps for simulation model 
# T_mpc - Number of time steps in the prediction horizon
# T_ss  - Number of time steps for the steady-state model (a single step, by definition)
##############################
NP = 3600
NT = 18000
DT = 300
T_sim = 2 # due to the trapezoidal integration scheme, the simulation model needs to be solved for 2 time steps to get the next time step solution
T_mpc = Int(NP/DT)
T_ss = 1 

##############################
# OPERATIONAL PROPERTIES 
##############################
# CHOKE_VALVE_OP_SS: Valve opening held fixed while the steady-state model is solved
# VLV_2_CV: Valve flow coefficient
# P_IN: bar, supply pressure at the pump inlet
# W_IN: kg/s, inlet flow to the pipeline. Only imposed in the steady-state model; in the NMPC
#       the flow follows from the pump curve and the network
# P_RES: bar, reservoir pressure
# pump_efficiency: Pump efficiency
# well_injectivity: kg/s/bar, nominal reservoir inflow coefficient, w = injectivity*(BHP - P_RES).

CHOKE_VALVE_OP_SS = vcat(zeros(K, 1)) # steady state valve opening for all valves
CHOKE_VALVE_OP_SS[all_vlv_idx] .= 1.0
VLV_2_CV = vcat(zeros(K, 1))
VLV_2_CV[all_vlv_idx] = [0.0034, 0.0034, 0.0034, 0.0034] # valve CV for all valves 
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
# u_vars     - Input variables; the MPC returns optimal values for these and they are the only
#              variables written to the simulation model
# state_vars - Differential states; these are the variables the initial condition is assigned to
#              at each MPC call, from the previous simulation model result
u_vars = ["choke_vlv_op", "speed_pump"]
state_vars = ["p_average"]

##############################
# MPC INPUT VARIABLE BOUNDS
##############################
# The *_roc values are per second and are multiplied by DT in model_mpc.jl to give the maximum
# movement allowed per sample. They are hard rate-of-change constraints; the separate L1 movement
# penalty in the objective is what smooths the inputs within those limits.
choke_vlv_op_min = 0.2
choke_vlv_op_max = 1.0
# maximum rate of change of the choke valve opening, per second
choke_vlv_op_roc = 1.0*(10.0^(-4))
speed_pump_min = 3000
speed_pump_max = 3700
# maximum rate of change of the pump speed, rpm/s
speed_pump_roc = 0.03

##############################
# MPC CONSTRAINED VARIABLES
##############################
# All of the constraints below are soft: they are relaxed by slack variables in model_mpc.jl and
# the violations are penalised in the objective function at the bottom of this file.
# pressure_con_idx: Nodes carrying a pressure constraint. Node 2 is the pump outlet; nodes
#           52, 54, 65 and 67 are the bottom-hole nodes of wells 1-4
# pipe_pressure_min / pipe_pressure_max: bar, min and max allowable CO2 pressure, entered per node
# cavitation_min: minimum allowable cavitation index across each choke valve
# pipe_pressure_roc: maximum allowable rate of change of the pipeline pressure, per second. The MPC model multiplies this by DT to get the maximum change allowed per sample.
# Applied only at pressure_con_idx nodes 
pressure_con_idx = [2, 52, 54, 65, 67]
pipe_pressure_max = vcat(zeros(J, 1))
pipe_pressure_max[pressure_con_idx] .= [250, 250, 250, 250, 250]
pipe_pressure_min = vcat(zeros(J, 1))
pipe_pressure_min[pressure_con_idx] .= [180, 180, 180, 180, 180]
cavitation_min = 1.7
pipe_pressure_roc = 0.005

##############################
# MPC SETPOINTS
##############################
# The setpoint trajectory has to cover the whole simulation plus one further prediction horizon,
# because the MPC at the last sample still looks NP seconds ahead.
n_samples_setpoint = length(DT:DT:(NT+NP)) + 1

# FLOW SETPOINTS 
well_1_flow_sp = fill(125, n_samples_setpoint, 1)
well_2_flow_sp = fill(125, n_samples_setpoint, 1)
well_3_flow_sp = fill(125, n_samples_setpoint, 1)
well_4_flow_sp = fill(125, n_samples_setpoint, 1)
    
well_4_flow_sp[25:end] .= 130
well_3_flow_sp[25:end] .= 130
well_2_flow_sp[25:end] .= 130
well_1_flow_sp[25:end] .= 130

##############################
# Initial Guess for Steady State Model
##############################
# One scalar per variable name; model_ss.jl uses it as the `start` value for every element of
# that variable. Keys with no matching variable in model_ss.jl are simply ignored.
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
    "rho_pump_inlet" => 922,
    "rho_pump_outlet" => 979,
    "pump_work" => 0.8,
    "cavitation_idx" => 0,
    "speed_pump" => 3440,
    "head_meters" => 2000,
    "q_pump" => 1950  # m3/h, i.e. W_IN/rho_pump_inlet*3600 at the nominal operating point
    )  
    
##############################
# MPC CONSTRAINTS FOR OPEN LOOP / CLOSED LOOP SIMULATION 
##############################

function special_constraints_fn_dynamic(model, sim_t)
    
    ##############################
    # SETPOINT SLICE FOR THIS SAMPLE
    ##############################
    # sim_t runs DT, 2*DT, ..., NT, so `start` runs 1, 2, ..., NT/DT and walks the setpoint
    # trajectory forward by one sample per closed-loop step.
    start = Int(sim_t/DT)

    # SET POINT TO MPC: one value per prediction step t = 1..T_mpc
    well_1_flow_sp_mpc = well_1_flow_sp[start:(start + T_mpc - 1)]
    well_2_flow_sp_mpc = well_2_flow_sp[start:(start + T_mpc - 1)]
    well_3_flow_sp_mpc = well_3_flow_sp[start:(start + T_mpc - 1)]
    well_4_flow_sp_mpc = well_4_flow_sp[start:(start + T_mpc - 1)]

    ##############################
    # OBJECTIVE FUNCTION TO MPC 
    ##############################
    @objective(model, Min, 
    5*sum((model[:w_in][51, t] - well_1_flow_sp_mpc[t])^2 for t=2:T_mpc) # Well Flow Rate Setpoint 
    + 20*sum((model[:w_in][53, t] - well_2_flow_sp_mpc[t])^2 for t=2:T_mpc) # Well Flow Rate Setpoint 
    + 80*sum((model[:w_in][64, t] - well_3_flow_sp_mpc[t])^2 for t=2:T_mpc) # Well Flow Rate Setpoint 
    + 200*sum((model[:w_in][66, t] - well_4_flow_sp_mpc[t])^2 for t=2:T_mpc) # Well Flow Rate Setpoint 
    + 1500*sum(model[:pipe_pressure_slack][k, t] for k in pressure_con_idx for t=2:T_mpc)    # Pipe Pressure CV Constraint 
    + 1500*sum(model[:pipe_pressure_surplus][k, t] for k in pressure_con_idx for t=2:T_mpc)  # Pipe Pressure CV Constraint 
    + 500*sum(model[:pipe_pressure_roc_slack][k, t] for k in pressure_con_idx for t=2:T_mpc)    # Pipe Pressure ROC Constraint 
    + 500*sum(model[:pipe_pressure_roc_surplus][k, t] for k in pressure_con_idx for t=2:T_mpc)  # Pipe Pressure ROC Constraint 
    + 800*sum(model[:cavitation_slack][k, t] for k in choke_vlv_idx for t=2:T_mpc)    # Choke Valve Cavitation Index  
    + 10*sum((model[:choke_vlv_op][50,t] - model[:choke_vlv_op][50,t-1])^2 for  t=2:T_mpc) # Input regularisation for moving the choke valve 
    + 10*sum((model[:choke_vlv_op][52,t] - model[:choke_vlv_op][52,t-1])^2 for t=2:T_mpc) # Input regularisation for moving the choke valve 
    + 10*sum((model[:choke_vlv_op][63,t] - model[:choke_vlv_op][63,t-1])^2 for  t=2:T_mpc) # Input regularisation for moving the choke valve 
    + 10*sum((model[:choke_vlv_op][65,t] - model[:choke_vlv_op][65,t-1])^2 for t=2:T_mpc) # Input regularisation for moving the choke valve 
    + sum(((model[:speed_pump][1, t] - model[:speed_pump][1, t-1])/100)^2 for t=2:T_mpc) # Input regularisation for moving the pump speed 
    )


end 

 
 special_constraints_dynamic = special_constraints_fn_dynamic
