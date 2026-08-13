####################################################
# NMPC INTERNAL MODEL
# First-principles dynamic model of the CO2 pipeline-injection network: mass and momentum
# balances per pipe segment, a pump curve, choke valve orifice equations and a linear reservoir
# inflow relation. Time is discretised by orthogonal collocation on finite elements.
#
# Every variable carries three indexes [k, d, t]:
#   k  arc or node index
#   d  collocation point within the finite element. The model equations are written at the
#      collocation points d = 2:D only, so that is the range every variable is declared over.
#      d = D is the end of the element and is the value reported for the sample.
#      The two exceptions are rho, the differential state, and p_node: both are declared over
#      d = 1:D because d = 1 (the start of the element) is referenced - by the collocation
#      derivative Adot'*rho and by the pressure rate-of-change constraint respectively.
#   t  finite element, 1:T_mpc, each of length DT
# The inputs (choke_vlv_op, speed_pump) have no collocation dimension: one value per finite
# element, held constant across it.
#
# The model is built ONCE and then reused for every closed-loop sample.
# All pressures are in bar and all flows in kg/s unless stated otherwise.
####################################################


function mpc_model_build()

    ##############################
    # DYNAMIC MODEL
    ##############################
    model = Model(Ipopt.Optimizer)

    ##############################
    # DEFINE VARIABLES
    ##############################
    # DEFINE GENERIC ARC VARIABLES
    @variable(model, p_in[k=1:K, d=2:D, t=1:T_mpc]>=0)
    @variable(model, p_out[k=1:K, d=2:D, t=1:T_mpc]>=0)
    @variable(model, w_in[k=1:K, d=2:D, t=1:T_mpc]>=0)
    @variable(model, w_out[k=1:K, d=2:D, t=1:T_mpc]>=0)

    # DEFINE PIPE SEGMENT VARIABLES AT INLET AND OUTLET
    @variable(model, reynold_in[k in pipe_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, reynold_out[k in pipe_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, fdarcy_in[k in pipe_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, fdarcy_out[k in pipe_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, dpdz_in[k in pipe_idx, d=2:D, t=1:T_mpc])
    @variable(model, dpdz_out[k in pipe_idx, d=2:D, t=1:T_mpc])

    # DEFINE LUMPED PIPE SEGMENT VARIABLES
    @variable(model, Z[k in pipe_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, Z_in[k in pipe_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, Z_out[k in pipe_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, rho[k in pipe_idx, d=1:D, t=1:T_mpc]>=500)
    @variable(model, rho_in[k in pipe_idx, d=2:D, t=1:T_mpc]>=500)
    @variable(model, rho_out[k in pipe_idx, d=2:D, t=1:T_mpc]>=500)
    @variable(model, p_average[k in pipe_idx, d=2:D, t=1:T_mpc]>=0)

    # DEFINE CHOKE VALVE VARIABLES
    @variable(model, choke_vlv_op[k in all_vlv_idx, t=1:T_mpc]>=0)
    @variable(model, deltap_wellhead_choke[k in all_vlv_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, cavitation_idx[k in choke_vlv_idx, d=2:D, t=1:T_mpc]>=0)

    # DEFINE NODE VARIABLES
    @variable(model, p_node[j in node_idx, d=1:D, t=1:T_mpc]>=0)

    # DEFINE PUMP VARIABLES
    @variable(model, pump_work[k in pump_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, rho_pump_inlet[k in pump_idx, d=2:D, t=1:T_mpc]>=200)
    @variable(model, rho_pump_outlet[k in pump_idx, d=2:D, t=1:T_mpc]>=200)
    @variable(model, head_meters[k in pump_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, q_pump[k in pump_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, speed_pump[k in pump_idx, t=1:T_mpc]>=0)

    # DEFINE RESERVOIR PARAMETERS
    # Well Injectivity & Reservoir Pressure.
    # These are declared as variables only so that they could be fixed to a new value at each
    # sample without rebuilding the model; they are never free for the optimiser.
    # p_reservoir and injectivity are currently fixed to constant values in this implementation,
    # but they could be made time-varying if desired.
    @variable(model, injectivity[k in well_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, p_reservoir[d=2:D, t=1:T_mpc]>=0)

    ##############################
    # COLLOCATION EQUATION FOR MASS BALANCE
    ##############################
    @constraint(model, [k in pipe_idx, t=1:T_mpc],
    Adot'*rho[k, 1:D, t]  .== (DT / V[k]) * (w_in[k, 2:D, t] - w_out[k, 2:D, t]))

    # Element continuity: each element starts where the previous one ended. Only the differential
    # state needs this 
    @constraint(model, [k in pipe_idx, t=2:T_mpc], rho[k, 1, t] == rho[k, D, t-1])

    # p_node is the one algebraic variable that also needs a d = 1 value, because the pipeline
    # pressure rate-of-change soft constraint below is written across the element as
    # p_node[k, D, t] - p_node[k, 1, t].
    @constraint(model, [j in node_idx, t=2:T_mpc], p_node[j, 1, t] == p_node[j, D, t-1])

    ##############################
    # PUMP EQUATIONS
    ##############################
    # r is the speed ratio against the 3000 rpm reference speed at which the pump curve was fitted.
    # The head/flow relation is applied in affinity-scaled form, head/r^2 as a quadratic in q/r, so
    # that one curve covers the whole speed range.
    @expression(model, r[k in pump_idx, t=1:T_mpc], speed_pump[k, t]/3000)
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_mpc], head_meters[k, d, t]/(r[k, t]^2) == -6*(10^-5)*(q_pump[k, d, t]/r[k, t])^2 - 1.0917*(q_pump[k, d, t]/r[k, t]) + 2859.6)
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_mpc], p_out[k, d, t] - p_in[k, d, t] == head_meters[k, d, t]*0.5*(rho_pump_inlet[k, d, t] + rho_pump_outlet[k, d, t])*G*(10^-5))
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_mpc], rho_pump_inlet[k, d, t] == p_in[k, d, t]*MW/R/poly2(p_in[k, d, t]*(10^5), fit_result_Z.param)/Tref)
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_mpc], rho_pump_outlet[k, d, t] == p_out[k, d, t]*MW/R/poly2(p_out[k, d, t]*(10^5), fit_result_Z.param)/Tref)
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_mpc], pump_work[k, d, t] == ( (p_out[k, d, t] - p_in[k, d, t])*(10^5)*w_in[k, d, t]/pump_efficiency[k]/(0.5*(rho_pump_inlet[k, d, t] + rho_pump_outlet[k, d, t])) )*DT/(3.6e9) ) # Pump work in MWh
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_mpc], w_in[k, d, t] == w_out[k, d, t]) # no accummulation across pump volumes
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_mpc], q_pump[k, d, t] == (w_in[k, d, t]/rho_pump_inlet[k, d, t])*3600) # m3/hour

    ##############################
    # MOMENTUM BALANCE EQUATIONS
    ##############################
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_mpc], dpdz_in[k, d, t] == -(fdarcy_in[k, d, t]*(w_in[k, d, t]^2)/(2*DINNER[k]*ACROSS[k]*ACROSS[k]*rho_in[k, d, t]))/1e5 - (rho_in[k, d, t]*G*sin(Theta[k]*pi/180))/1e5)
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_mpc], dpdz_out[k, d, t] == -(fdarcy_out[k, d, t]*(w_out[k, d, t]^2)/(2*DINNER[k]*ACROSS[k]*ACROSS[k]*rho_out[k, d, t]))/1e5 - (rho_out[k, d, t]*G*sin(Theta[k]*pi/180))/1e5)
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_mpc], p_out[k, d, t] == p_in[k, d, t] + (LENGTH[k]/2)*(dpdz_in[k, d, t] + dpdz_out[k, d, t]))

    ############################
    #  COMPRESSIBILITY FACTOR
    ############################
    # Z is taken from the quadratic CoolProp fit at the local pressure, and density follows from the
    # real gas law. Evaluated separately at the segment average, inlet and outlet.
    # NOTE: model_ss.jl additionally bounds Z <= 1; that bound is deliberately not imposed here.
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_mpc], Z[k, d, t] == poly2(p_average[k, d, t]*(10^5), fit_result_Z.param))
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_mpc], Z_in[k, d, t] == poly2(p_in[k, d, t]*(10^5), fit_result_Z.param))
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_mpc], Z_out[k, d, t] == poly2(p_out[k, d, t]*(10^5), fit_result_Z.param))
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_mpc], rho[k, d, t] == p_average[k, d, t]*MW/R/Z[k, d, t]/Tref)
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_mpc], rho_in[k, d, t] == p_in[k, d, t]*MW/R/Z_in[k, d, t]/Tref)
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_mpc], rho_out[k, d, t] == p_out[k, d, t]*MW/R/Z_out[k, d, t]/Tref)

    ############################
    # DEFINE FRICTION FACTOR & REYNOLDS NUMBER
    ############################
    # REYNOLDS NUMBER AND DARCY FRICTION FACTOR CALCULATION
    # Viscosity comes from the quadratic CoolProp fit.
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_mpc], reynold_in[k, d, t] ==  w_in[k, d, t]*DINNER[k]/(poly2(p_in[k, d, t]*(10^5), fit_result_mu.param)*ACROSS[k]))
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_mpc], reynold_out[k, d, t] ==  w_out[k, d, t]*DINNER[k]/(poly2(p_out[k, d, t]*(10^5), fit_result_mu.param)*ACROSS[k]))
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_mpc], fdarcy_in[k, d, t] == 1 / ( -1.8*log10(6.9/reynold_in[k, d, t] + (EPS/(DINNER[k]*3.7))^1.11) )^2)
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_mpc], fdarcy_out[k, d, t] == 1 / ( -1.8*log10(6.9/reynold_out[k, d, t] + (EPS/(DINNER[k]*3.7))^1.11) )^2)

    ############################
    # AVERAGE PRESSURE CALCULATION
    ############################
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_mpc], p_average[k, d, t] == (0.5*p_in[k, d, t] + (1/12)*dpdz_in[k, d, t]*LENGTH[k] + 0.5*p_out[k, d, t] - (1/12)*dpdz_out[k, d, t]*LENGTH[k])/1)

    ##############################
    # VALVE EQUATIONS
    ##############################
    # Valves are treated as zero-volume elements: the flow passes straight through and the only
    # effect is a pressure drop. The drop is forced to be non-negative so that flow cannot reverse
    # through a choke.
    # Orifice equation w = VLV_2_CV*opening*sqrt(rho*dP).
    # The density is taken from the outlet of the pipe segment immediately upstream of the valve.
    @constraint(model, [k in all_vlv_idx, d=2:D, t=1:T_mpc], (w_in[k, d, t]/(VLV_2_CV[k]*(choke_vlv_op[k, t])))^2/(rho_out[nearest_pipe_idx[k], d, t]) == deltap_wellhead_choke[k, d, t]*(10^5))
    @constraint(model, [k in all_vlv_idx, d=2:D, t=1:T_mpc], p_in[k, d, t] == p_out[k, d, t] + deltap_wellhead_choke[k, d, t])
    @constraint(model, [k in all_vlv_idx, d=2:D, t=1:T_mpc], w_in[k, d, t] == w_out[k, d, t])
    @constraint(model, [k in all_vlv_idx, d=2:D, t=1:T_mpc], p_in[k, d, t] >= p_out[k, d, t])

    # CAVITATION INDEX CALCULATION FOR CHOKE VALVES
    # Cavitation index = (p_in - p_vapour)/(p_in - p_out)
    @constraint(model, [k in choke_vlv_idx, d=2:D, t=1:T_mpc], cavitation_idx[k, d, t] == (p_in[k, d, t] - co2_vapor_pressure)/(p_in[k, d, t] - p_out[k, d, t]))

    ##############################
    # PIPELINE NETWORK - DEFINE JUNCTIONS / CONNECTIONS
    # MASS BALANCE
    ##############################
    # Interior nodes have no volume, so inflow equals outflow.
    @constraint(model, [j in setdiff(node_idx, inlet_node, outlet_node), d=2:D, t=1:T_mpc], sum(w_out[x, d, t] for x in incoming[j]) == sum(w_in[y, d, t] for y in outgoing[j]))

    ##############################
    # PIPELINE NETWORK
    # PRESSURE EQUALITY AT NODES
    ##############################
    # A node has a single pressure, shared by every arc meeting it. At a branch point this is what
    # couples the two downstream arcs together.
    @constraint(model, [j in setdiff(node_idx, outlet_node), k in outgoing[j], d=2:D, t=1:T_mpc], p_node[j, d, t] == p_in[k, d, t])
    @constraint(model, [j in setdiff(node_idx, inlet_node), k in incoming[j], d=2:D, t=1:T_mpc], p_node[j, d, t] == p_out[k, d, t])

    ##############################
    # BOUNDARY CONDITIONS FOR SYSTEM
    ##############################
    # Downstream: linear reservoir inflow, driven by the bottom-hole to reservoir pressure difference.
    @constraint(model, [k in well_idx, d=2:D, t=1:T_mpc], w_out[k, d, t] == (p_out[k, d, t] - p_reservoir[d, t])*injectivity[k, d, t])
    # Upstream: fixed supply pressure at the pump inlet.
    @constraint(model, p_in[1, 2:D, :] .== P_IN)

    # Reservoir pressure and injectivity held constant throughout; they are variables only so that
    # they could be updated per sample.
    for d in 2:D, k in well_idx, t in 1:T_mpc; fix(injectivity[k, d, t], well_injectivity[k]; force=true); end
    for d in 2:D, t in 1:T_mpc; fix(p_reservoir[d, t], P_RES; force=true); end

    ##############################
    # MV CONSTRAINTS
    # 1. Choke Valve Rate of Change Constraint
    # 2. Choke Valve Output Limits
    # 3. Pump Speed Output Limits
    # 4. Pump Speed Rate of Change Constraints
    ##############################
    @constraint(model, [k in choke_vlv_idx, t=2:T_mpc], choke_vlv_op[k, t] - choke_vlv_op[k, t-1] <= choke_vlv_op_roc*DT)
    @constraint(model, [k in choke_vlv_idx, t=2:T_mpc], choke_vlv_op[k, t-1] - choke_vlv_op[k, t] <= choke_vlv_op_roc*DT)
    @constraint(model, [k in choke_vlv_idx, t=1:T_mpc], choke_vlv_op[k, t] <= choke_vlv_op_max)
    @constraint(model, [k in choke_vlv_idx, t=1:T_mpc], choke_vlv_op[k, t] >= choke_vlv_op_min)
    @constraint(model, [k in pump_idx, t=1:T_mpc], speed_pump[k, t] <= speed_pump_max)
    @constraint(model, [k in pump_idx, t=1:T_mpc], speed_pump[k, t] >= speed_pump_min)
    @constraint(model, [k in pump_idx, t=2:T_mpc], speed_pump[k, t-1] - speed_pump[k, t] <= speed_pump_roc*DT)
    @constraint(model, [k in pump_idx, t=2:T_mpc], speed_pump[k, t] - speed_pump[k, t-1] <= speed_pump_roc*DT)

    # PREVIOUSLY APPLIED INPUT (set at each MPC call to initial_con value)
    @variable(model, u_prev_choke[k in choke_vlv_idx])
    @variable(model, u_prev_pump[k in pump_idx])

    # ROC AT t=1 (ties u[k,1] to previously applied input)
    @constraint(model, [k in choke_vlv_idx], choke_vlv_op[k, 1] - u_prev_choke[k] <= choke_vlv_op_roc*DT)
    @constraint(model, [k in choke_vlv_idx], u_prev_choke[k] - choke_vlv_op[k, 1] <= choke_vlv_op_roc*DT)
    @constraint(model, [k in pump_idx], speed_pump[k, 1] - u_prev_pump[k] <= speed_pump_roc*DT)
    @constraint(model, [k in pump_idx], u_prev_pump[k] - speed_pump[k, 1] <= speed_pump_roc*DT)

    ##############################
    # Soft Constraints
    # 1. Pipe Pressure
    # 2. Pipe Pressure Rate of Change
    # 3. Cavitation Index
    ##############################
    # All three are imposed at d = D, the end of each finite element.

    # PIPELINE PRESSURE CONSTRAINTS
    @variable(model, pipe_pressure_surplus[k in pressure_con_idx,  t=1:T_mpc]>=0)
    @variable(model, pipe_pressure_slack[k in pressure_con_idx,  t=1:T_mpc]>=0)
    @constraint(model, [k in pressure_con_idx,  t=1:T_mpc], p_node[k, D, t] + pipe_pressure_slack[k, t]>= pipe_pressure_min[k])
    @constraint(model, [k in pressure_con_idx,  t=1:T_mpc], p_node[k, D, t] <= pipe_pressure_max[k] + pipe_pressure_surplus[k, t])

    # PIPELINE PRESSURE ROC CONSTRAINTS (within each element: end minus start)
    # For t>=2, p_node[k, 1, t] == p_node[k, D, t-1] by algebraic continuity, so this matches the previous formulation.
    # For t=1, p_node[k, 1, 1] is fixed to the MPC's IC, giving (p_node[k, D, 1] - IC)/DT.
    @variable(model, pipe_pressure_roc_surplus[k in pressure_con_idx,  t=1:T_mpc]>=0)
    @variable(model, pipe_pressure_roc_slack[k in pressure_con_idx,  t=1:T_mpc]>=0)
    @constraint(model, [k in pressure_con_idx,  t=1:T_mpc], (p_node[k, D, t] - p_node[k, 1, t])/DT <= pipe_pressure_roc + pipe_pressure_roc_slack[k, t])
    @constraint(model, [k in pressure_con_idx,  t=1:T_mpc], -(p_node[k, D, t] - p_node[k, 1, t])/DT <= pipe_pressure_roc + pipe_pressure_roc_surplus[k, t])

    # CAVITATION INDEX CONSTRAINTS
    @variable(model, cavitation_slack[k in choke_vlv_idx,  t=1:T_mpc]>=0)
    @constraint(model, [k in choke_vlv_idx,  t=1:T_mpc], cavitation_idx[k, D, t] + cavitation_slack[k, t] >= cavitation_min)

    return model
end
