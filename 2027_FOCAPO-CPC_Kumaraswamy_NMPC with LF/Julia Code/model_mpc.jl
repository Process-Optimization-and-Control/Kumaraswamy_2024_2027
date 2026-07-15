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
    @variable(model, rho[k in pipe_idx, d=2:D, t=1:T_mpc]>=500)
    @variable(model, rho_in[k in pipe_idx, d=2:D, t=1:T_mpc]>=500)
    @variable(model, rho_out[k in pipe_idx, d=2:D, t=1:T_mpc]>=500)
    @variable(model, p_average[k in pipe_idx, d=1:D, t=1:T_mpc]>=0)

    # DEFINE CHOKE VALVE VARIABLES 
    @variable(model, choke_vlv_op[k in all_vlv_idx, t=1:T_mpc]>=0)
    @variable(model, deltap_wellhead_choke[k in all_vlv_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, cavitation_idx[k in choke_vlv_idx, d=2:D, t=1:T_mpc]>=0)

    # DEFINE NODE VARIABLES 
    @variable(model, p_node[j in node_idx, d=2:D, t=1:T_mpc]>=0)

    # DEFINE PUMP VARIABLES 
    @variable(model, pump_work[k in pump_idx, d=2:D, t=1:T_mpc]>=0) # Element work carried at d=D (endpoint); interior nodes fixed to 0 below
    @variable(model, rho_pump_inlet[k in pump_idx, d=2:D, t=1:T_mpc]>=200) 
    @variable(model, rho_pump_outlet[k in pump_idx, d=2:D, t=1:T_mpc]>=200) 
    @variable(model, head_meters[k in pump_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, q_pump[k in pump_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, speed_pump[k in pump_idx, t=1:T_mpc]>=0)

    # DEFINE RESERVOIR PARAMETERS 
    # Well Injectivity & Reservoir Pressure 
    @variable(model, injectivity[k in well_idx, d=2:D, t=1:T_mpc]>=0)
    @variable(model, p_reservoir[d=2:D, t=1:T_mpc]>=0)
    @variable(model, VLV_2_CV[k in all_vlv_idx]>=0)

    # OUTPUT MEASUREMENT BIAS TERM 
    @variable(model, flow_bias[k in well_idx])
    @variable(model, w_in_corrected[k in well_idx, d=2:D, t=1:T_mpc])
    @variable(model, mainline_flow_bias)
    @variable(model, w_in_corrected_mainline[d=2:D, t=1:T_mpc])
    @variable(model, p_bias[k in pressure_con_idx])
    @variable(model, p_node_corrected[k in pressure_con_idx, d=2:D, t=1:T_mpc])

    ##############################
    # COLLOCATION EQUATION FOR MASS BALANCE
    ##############################
    @constraint(model, [k in pipe_idx, t=1:T_mpc], 
    Adot'*p_average[k, 1:D, t]  .== DT / V[k] * (w_in[k, 2:D, t] - w_out[k, 2:D, t]) * R .* Z[k, 2:D, t] * Tref / MW )
    
    @constraint(model, [k in pipe_idx, t=2:T_mpc], p_average[k, 1, t] == p_average[k, D, t-1])


    ##############################
    # PUMP EQUATIONS    
    ##############################
    @expression(model, r[k in pump_idx, t=1:T_mpc], speed_pump[k, t]/3000)
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_mpc], head_meters[k, d, t]/(r[k, t]^2) == -6*(10^-5)*(q_pump[k, d, t]/r[k, t])^2 - 1.0917*(q_pump[k, d, t]/r[k, t]) + 2859.6)
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_mpc], p_out[k, d, t] - p_in[k, d, t] == head_meters[k, d, t]*0.5*(rho_pump_inlet[k, d, t] + rho_pump_outlet[k, d, t])*G*(10^-5))
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_mpc], rho_pump_inlet[k, d, t] == p_in[k, d, t]*MW/R/poly2(p_in[k, d, t]*(10^5), fit_result_Z.param)/Tref)
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_mpc], rho_pump_outlet[k, d, t] == p_out[k, d, t]*MW/R/poly2(p_out[k, d, t]*(10^5), fit_result_Z.param)/Tref)
    # Pump work in MWh over the element, evaluated at the endpoint node d=D (right-endpoint rule): power(d=D) * DT
    @constraint(model, [k in pump_idx, t=1:T_mpc], pump_work[k, D, t] == ( (p_out[k, D, t] - p_in[k, D, t])*(10^5)*w_in[k, D, t]/pump_efficiency[k]/(0.5*(rho_pump_inlet[k, D, t] + rho_pump_outlet[k, D, t])) )*DT/(3.6e9) )
    # Interior collocation nodes carry no work
    for k in pump_idx, d in 2:(D-1), t in 1:T_mpc
        fix(pump_work[k, d, t], 0.0; force=true)
    end
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
    # CONSTRAINTS FOR BULK MODULUS AND COMPRESSIBILITY FACTOR  
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
    @constraint(model, [k in all_vlv_idx, d=2:D, t=1:T_mpc], p_in[k, d, t] == p_out[k, d, t] + deltap_wellhead_choke[k, d, t])
    @constraint(model, [k in all_vlv_idx, d=2:D, t=1:T_mpc], w_in[k, d, t] == w_out[k, d, t])
    @constraint(model, [k in all_vlv_idx, d=2:D, t=1:T_mpc], p_in[k, d, t] >= p_out[k, d, t])

    @constraint(model, [k in all_vlv_idx, d=2:D, t=1:T_mpc], (w_in[k, d, t])^2 == ((VLV_2_CV[k]*choke_vlv_op[k, t])^2)*rho_out[nearest_pipe_idx[k], d, t]*deltap_wellhead_choke[k, d, t]*(10^5)) 

    # CAVITATION INDEX CALCULATION FOR CHOKE VALVES 
    @constraint(model, [k in choke_vlv_idx, d=2:D, t=1:T_mpc], cavitation_idx[k, d, t]*(p_in[k, d, t] - p_out[k, d, t]) == (p_in[k, d, t] - co2_vapor_pressure))

    ##############################
    # PIPELINE NETWORK - DEFINE JUNCTIONS / CONNECTIONS 
    # MASS BALANCE
    ##############################
    @constraint(model, [j in setdiff(node_idx, inlet_node, outlet_node), d=2:D, t=1:T_mpc], sum(w_out[x, d, t] for x in incoming[j]) == sum(w_in[y, d, t] for y in outgoing[j]))

    ##############################
    # PIPELINE NETWORK
    # PRESSURE EQUALITY AT NODES 
    ##############################
    @constraint(model, [j in setdiff(node_idx, outlet_node), k in outgoing[j], d=2:D, t=1:T_mpc], p_node[j, d, t] == p_in[k, d, t])
    @constraint(model, [j in setdiff(node_idx, inlet_node), k in incoming[j], d=2:D, t=1:T_mpc], p_node[j, d, t] == p_out[k, d, t])

    ##############################
    # BOUNDARY CONDITIONS FOR SYSTEM 
    ##############################
    @constraint(model, [k in well_idx, d=2:D, t=1:T_mpc], w_out[k, d, t] == (p_out[k, d, t] - p_reservoir[d, t])*injectivity[k, d, t])
    @constraint(model, p_in[1, 2:D, :] .== P_IN)

    # FIX RESERVOIR PRESSURE TO STEADY STATE VALUES FOR SIMULATION
    # INJECTIVITY ESTIMATED 
    for d in 2:D, t in 1:T_mpc; fix(p_reservoir[d, t], P_RES; force=true); end
    
    ##############################
    # BIAS CORRECTION FOR ALL MEASURED FLOWS
    ##############################
    # Flow Entering 4 Wells
    @constraint(model, [k in well_idx], w_in_corrected[k, :, :] .== flow_bias[k] .+ w_in[k, :, :])
    # Flow entering the mainline 
    @constraint(model, w_in_corrected_mainline[:, :] .== mainline_flow_bias .+ w_in[2, :, :])

    ##############################
    # BIAS CORRECTION FOR PRESSURES WITH CONSTRAINTS 
    ##############################
    @constraint(model, [k in pressure_con_idx], p_node_corrected[k, :, :] .== p_bias[k] .+ p_node[k, :, :])

    ##############################
    # MV CONSTRAINTS 
    # 1. Choke Valve Output Limits 
    # 2. Pump Speed Output Limits 
    # 3. Choke Valve Rate of Change Constraint 
    # 4. Pump Speed Rate of Change Constraints 
    # 5. Holding MV Constant After Control Horizon
    ##############################

    # MV Limits 
    @constraint(model, [k in choke_vlv_idx, t=1:NM], choke_vlv_op[k, t] <= choke_vlv_op_max)  
    @constraint(model, [k in choke_vlv_idx, t=1:NM], choke_vlv_op[k, t] >= choke_vlv_op_min)      
    @constraint(model, [k in pump_idx, t=1:NM], speed_pump[k, t] <= speed_pump_max)
    @constraint(model, [k in pump_idx, t=1:NM], speed_pump[k, t] >= speed_pump_min)

    # Rate of Change Constraints
    @constraint(model, [k in choke_vlv_idx, t=2:NM], choke_vlv_op[k, t] - choke_vlv_op[k, t-1] <= choke_vlv_op_roc*DT)
    @constraint(model, [k in choke_vlv_idx, t=2:NM], choke_vlv_op[k, t-1] - choke_vlv_op[k, t] <= choke_vlv_op_roc*DT)
    @constraint(model, [k in pump_idx, t=2:NM], speed_pump[k, t-1] - speed_pump[k, t] <= speed_pump_roc*DT)
    @constraint(model, [k in pump_idx, t=2:NM], speed_pump[k, t] - speed_pump[k, t-1] <= speed_pump_roc*DT)

    # Rate of Change Constraints at t=1 (ties first time step to previously applied input)
    # Variable for previously applied input (set at each MPC call to initial_con value)
    @variable(model, u_prev_choke[k in choke_vlv_idx])
    @variable(model, u_prev_pump[k in pump_idx])

    # ROC AT t=1 (ties u[k,1] to previously applied input)
    @constraint(model, [k in choke_vlv_idx], choke_vlv_op[k, 1] - u_prev_choke[k] <= choke_vlv_op_roc*DT)
    @constraint(model, [k in choke_vlv_idx], u_prev_choke[k] - choke_vlv_op[k, 1] <= choke_vlv_op_roc*DT)
    @constraint(model, [k in pump_idx], speed_pump[k, 1] - u_prev_pump[k] <= speed_pump_roc*DT)
    @constraint(model, [k in pump_idx], u_prev_pump[k] - speed_pump[k, 1] <= speed_pump_roc*DT)

    # Holding MV Constant After Control Horizon
    @constraint(model, [k in choke_vlv_idx, t=(NM+1):T_mpc], choke_vlv_op[k, t-1]  == choke_vlv_op[k, t])
    @constraint(model, [k in pump_idx, t=(NM+1):T_mpc], speed_pump[k, t-1]  == speed_pump[k, t])

    ##############################
    # Soft Constraints 
    # 1. Pipe Pressure 
    # 2. Cavitation Index  
    ##############################

    # PIPELINE PRESSURE CONSTRAINTS      
    @variable(model, pipe_pressure_surplus[k in pressure_con_idx,  t=1:T_mpc]>=0)
    @variable(model, pipe_pressure_slack[k in pressure_con_idx,  t=1:T_mpc]>=0)
    @constraint(model, [k in pressure_con_idx,  t=1:T_mpc], p_node_corrected[k, D, t] + pipe_pressure_slack[k, t]>= pipe_pressure_min[k])
    @constraint(model, [k in pressure_con_idx,  t=1:T_mpc], p_node_corrected[k, D, t] <= pipe_pressure_max[k] + pipe_pressure_surplus[k, t])

    # CAVITATION INDEX CONSTRAINTS
    @variable(model, cavitation_slack[k in choke_vlv_idx,  t=1:T_mpc]>=0)
    @constraint(model, [k in choke_vlv_idx,  t=1:T_mpc], cavitation_idx[k, D, t] + cavitation_slack[k, t] >= cavitation_min)

    # L1 Penalty on Input Movement (Input Regularisation)
    @variable(model, s_choke[k in choke_vlv_idx, t=1:T_mpc] >= 0)
    @variable(model, s_pump[k in pump_idx, t=1:T_mpc] >= 0)
    @constraint(model, [k in choke_vlv_idx, t=2:T_mpc], s_choke[k,t] >=  choke_vlv_op[k,t] - choke_vlv_op[k,t-1])
    @constraint(model, [k in choke_vlv_idx, t=2:T_mpc], s_choke[k,t] >= -(choke_vlv_op[k,t] - choke_vlv_op[k,t-1]))
    @constraint(model, [k in pump_idx, t=2:T_mpc], s_pump[k,t] >=  speed_pump[k,t] - speed_pump[k,t-1])
    @constraint(model, [k in pump_idx, t=2:T_mpc], s_pump[k,t] >= -(speed_pump[k,t] - speed_pump[k,t-1]))
    # L1 Penalty on Input Movement at t=1 (ties input movement to previously applied input)
    @constraint(model, [k in choke_vlv_idx], s_choke[k, 1] >=  choke_vlv_op[k, 1] - u_prev_choke[k])
    @constraint(model, [k in choke_vlv_idx], s_choke[k, 1] >= -(choke_vlv_op[k, 1] - u_prev_choke[k]))
    @constraint(model, [k in pump_idx], s_pump[k, 1] >=  speed_pump[k, 1] - u_prev_pump[k])
    @constraint(model, [k in pump_idx], s_pump[k, 1] >= -(speed_pump[k, 1] - u_prev_pump[k]))
    return model
end
 
