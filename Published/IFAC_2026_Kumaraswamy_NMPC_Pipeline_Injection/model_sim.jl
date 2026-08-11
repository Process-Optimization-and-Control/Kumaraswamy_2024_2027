
function sim_model_build()

    ##############################
    # SIMULATION MODEL  
    ##############################
    model = Model(Ipopt.Optimizer)

    ##############################
    # DEFINE VARIABLES   
    ##############################
    # DEFINE GENERIC ARC VARIABLES 
    @variable(model, p_in[k=1:K, t=1:T_sim]>=0)
    @variable(model, p_out[k=1:K, t=1:T_sim]>=0)
    @variable(model, w_in[k=1:K, t=1:T_sim]>=0)
    @variable(model, w_out[k=1:K, t=1:T_sim]>=0)

    # DEFINE PIPE SEGMENT VARIABLES AT INLET AND OUTLET 
    @variable(model, reynold_in[k in pipe_idx, t=1:T_sim]>=0)
    @variable(model, reynold_out[k in pipe_idx, t=1:T_sim]>=0)
    @variable(model, fdarcy_in[k in pipe_idx, t=1:T_sim]>=0)
    @variable(model, fdarcy_out[k in pipe_idx, t=1:T_sim]>=0)
    @variable(model, dpdz_in[k in pipe_idx, t=1:T_sim])
    @variable(model, dpdz_out[k in pipe_idx, t=1:T_sim])

    # DEFINE LUMPED PIPE SEGMENT VARIABLES 
    @variable(model, Z[k in pipe_idx, t=1:T_sim]>=0)
    @variable(model, Z_in[k in pipe_idx, t=1:T_sim]>=0)
    @variable(model, Z_out[k in pipe_idx, t=1:T_sim]>=0)
    @variable(model, rho[k in pipe_idx, t=1:T_sim]>=500)
    @variable(model, rho_in[k in pipe_idx, t=1:T_sim]>=500)
    @variable(model, rho_out[k in pipe_idx, t=1:T_sim]>=500)
    @variable(model, p_average[k in pipe_idx, t=1:T_sim]>=0)

    # DEFINE CHOKE VALVE VARIABLES 
    @variable(model, choke_vlv_op[k in choke_vlv_idx, t=1:T_sim]>=0)
    @variable(model, deltap_wellhead_choke[k in choke_vlv_idx, t=1:T_sim]>=0)
    @variable(model, cavitation_idx[k in choke_vlv_idx, t=1:T_sim]>=0)

    # DEFINE NODE VARIABLES 
    @variable(model, p_node[j in node_idx, t=1:T_sim]>=0)

    # DEFINE PUMP VARIABLES 
    @variable(model, pump_work[k in pump_idx, t=1:T_sim]>=0)
    @variable(model, rho_pump_inlet[k in pump_idx, t=1:T_sim]>=200) 
    @variable(model, rho_pump_outlet[k in pump_idx, t=1:T_sim]>=200) 
    @variable(model, head_meters[k in pump_idx, t=1:T_sim]>=0)
    @variable(model, q_pump[k in pump_idx, t=1:T_sim]>=0)
    @variable(model, speed_pump[k in pump_idx, t=1:T_sim]>=0)

    # DEFINE RESERVOIR PARAMETERS
    # Well Injectivity & Reservoir Pressure.
    # These two are declared as variables only so that they can be fixed to a new value at each
    # sample without rebuilding the model; they are never free for the optimiser.
    @variable(model, injectivity[k in well_idx, t=1:T_sim]>=0)
    @variable(model, p_reservoir[t=1:T_sim]>=0)

    ##############################
    # MASS BALANCE EQUATIONS 
    ##############################
    @constraint(model, [k in pipe_idx, t=1:T_sim-1], V[k]*(rho[k, t+1] - rho[k, t]) == 0.5*((w_in[k, t] - w_out[k, t]) + (w_in[k, t+1] - w_out[k, t+1]))*DT)

    ##############################
    # PUMP EQUATIONS    
    ##############################
    # r is the speed ratio against the 3000 rpm reference speed at which the pump curve was fitted.
    # The head/flow relation is applied in affinity-scaled form, head/r^2 as a quadratic in q/r, so
    # that one curve covers the whole speed range.
    @expression(model, r[k in pump_idx, t=1:T_sim], speed_pump[k, t]/3000)
    @constraint(model, [k in pump_idx, t=1:T_sim], head_meters[k, t]/(r[k, t]^2) == -6*(10^-5)*(q_pump[k, t]/r[k, t])^2 - 1.0917*(q_pump[k, t]/r[k, t]) + 2859.6)
    @constraint(model, [k in pump_idx, t=1:T_sim], p_out[k, t] - p_in[k, t] == head_meters[k, t]*0.5*(rho_pump_inlet[k, t] + rho_pump_outlet[k, t])*G*(10^-5))
    @constraint(model, [k in pump_idx, t=1:T_sim], rho_pump_inlet[k, t] == p_in[k, t]*MW/R/poly2(p_in[k, t]*(10^5), fit_result_Z.param)/Tref)
    @constraint(model, [k in pump_idx, t=1:T_sim], rho_pump_outlet[k, t] == p_out[k, t]*MW/R/poly2(p_out[k, t]*(10^5), fit_result_Z.param)/Tref)
    @constraint(model, [k in pump_idx, t=1:T_sim], pump_work[k, t] == ( (p_out[k, t] - p_in[k, t])*(10^5)*w_in[k, t]/pump_efficiency[k]/(0.5*(rho_pump_inlet[k, t] + rho_pump_outlet[k, t])) )*DT/(3.6e9) ) # Pump work in MWh
    @constraint(model, [k in pump_idx, t=1:T_sim], w_in[k, t] == w_out[k, t]) # no accummulation across pump volumes 
    @constraint(model, [k in pump_idx, t=1:T_sim], q_pump[k, t] == (w_in[k, t]/rho_pump_inlet[k, t])*3600)

    ##############################
    # MOMENTUM BALANCE EQUATIONS   
    ##############################
    @constraint(model, [k in pipe_idx, t=1:T_sim], dpdz_in[k, t] == -(fdarcy_in[k, t]*(w_in[k, t]^2)/(2*DINNER[k]*ACROSS[k]*ACROSS[k]*rho_in[k, t]))/1e5 - (rho_in[k, t]*G*sin(Theta[k]*pi/180))/1e5)
    @constraint(model, [k in pipe_idx, t=1:T_sim], dpdz_out[k, t] == -(fdarcy_out[k, t]*(w_out[k, t]^2)/(2*DINNER[k]*ACROSS[k]*ACROSS[k]*rho_out[k, t]))/1e5 - (rho_out[k, t]*G*sin(Theta[k]*pi/180))/1e5)
    @constraint(model, [k in pipe_idx, t=1:T_sim], p_out[k, t] == p_in[k, t] + (LENGTH[k]/2)*(dpdz_in[k, t] + dpdz_out[k, t]))
    
    ############################
    #  COMPRESSIBILITY FACTOR
    ############################
    # Z is taken from the quadratic CoolProp fit at the local pressure, and density follows from the
    # real gas law. Evaluated separately at the segment average, inlet and outlet.
    # NOTE: model_ss.jl additionally bounds Z <= 1; that bound is deliberately not imposed here.
    @constraint(model, [k in pipe_idx, t=1:T_sim], Z[k, t] == poly2(p_average[k, t]*(10^5), fit_result_Z.param))
    @constraint(model, [k in pipe_idx, t=1:T_sim], Z_in[k, t] == poly2(p_in[k, t]*(10^5), fit_result_Z.param))
    @constraint(model, [k in pipe_idx, t=1:T_sim], Z_out[k, t] == poly2(p_out[k, t]*(10^5), fit_result_Z.param))
    @constraint(model, [k in pipe_idx, t=1:T_sim], rho[k, t] == p_average[k, t]*MW/R/Z[k, t]/Tref)
    @constraint(model, [k in pipe_idx, t=1:T_sim], rho_in[k, t] == p_in[k, t]*MW/R/Z_in[k, t]/Tref)
    @constraint(model, [k in pipe_idx, t=1:T_sim], rho_out[k, t] == p_out[k, t]*MW/R/Z_out[k, t]/Tref)

    ############################
    # DEFINE FRICTION FACTOR & REYNOLDS NUMBER   
    ############################
    # REYNOLDS NUMBER AND DARCY FRICTION FACTOR CALCULATION
    # Viscosity comes from the quadratic CoolProp fit.
    @constraint(model, [k in pipe_idx, t=1:T_sim], reynold_in[k, t] ==  w_in[k, t]*DINNER[k]/(poly2(p_in[k, t]*(10^5), fit_result_mu.param)*ACROSS[k]))
    @constraint(model, [k in pipe_idx, t=1:T_sim], reynold_out[k, t] ==  w_out[k, t]*DINNER[k]/(poly2(p_out[k, t]*(10^5), fit_result_mu.param)*ACROSS[k]))
    @constraint(model, [k in pipe_idx, t=1:T_sim], fdarcy_in[k, t] == 1 / ( -1.8*log10(6.9/reynold_in[k, t] + (EPS/(DINNER[k]*3.7))^1.11) )^2)
    @constraint(model, [k in pipe_idx, t=1:T_sim], fdarcy_out[k, t] == 1 / ( -1.8*log10(6.9/reynold_out[k, t] + (EPS/(DINNER[k]*3.7))^1.11) )^2)

    ############################
    # AVERAGE PRESSURE CALCULATION    
    ############################
    @constraint(model, [k in pipe_idx, t=1:T_sim], p_average[k, t] == (0.5*p_in[k, t] + (1/12)*dpdz_in[k, t]*LENGTH[k] + 0.5*p_out[k, t] - (1/12)*dpdz_out[k, t]*LENGTH[k])/1)

    ##############################
    # VALVE EQUATIONS
    ##############################
    # Valves are treated as zero-volume elements: the flow passes straight through and the only
    # effect is a pressure drop. The drop is forced to be non-negative so that flow cannot reverse
    # through a choke.
    # Orifice equation w = VLV_2_CV*opening*sqrt(rho*dP).
    # The density is taken from the outlet of the pipe segment immediately upstream of the valve.
    @constraint(model, [k in all_vlv_idx, t=1:T_sim], (w_in[k, t]/(VLV_2_CV[k]*(choke_vlv_op[k, t])))^2/(rho_out[nearest_pipe_idx[k], t]) == deltap_wellhead_choke[k, t]*(10^5))

    @constraint(model, [k in all_vlv_idx, t=1:T_sim], p_in[k, t] == p_out[k, t] + deltap_wellhead_choke[k, t])
    @constraint(model, [k in all_vlv_idx, t=1:T_sim], w_in[k, t] == w_out[k, t])
    @constraint(model, [k in all_vlv_idx, t=1:T_sim], p_in[k, t] >= p_out[k, t])

    # CAVITATION INDEX CALCULATION FOR CHOKE VALVES
    # Cavitation index = (p_in - p_vapour)/(p_in - p_out)
    @constraint(model, [k in all_vlv_idx, t=1:T_sim], cavitation_idx[k, t] == (p_in[k, t] - co2_vapor_pressure)/(p_in[k, t] - p_out[k, t]))

    ##############################
    # PIPELINE NETWORK - DEFINE JUNCTIONS / CONNECTIONS 
    # MASS BALANCE
    ##############################
    # Interior nodes have no volume, so inflow equals outflow. 
    @constraint(model, [j in setdiff(node_idx, inlet_node, outlet_node), t=1:T_sim], sum(w_out[x, t] for x in incoming[j]) == sum(w_in[y, t] for y in outgoing[j]))

    ##############################
    # PIPELINE NETWORK
    # PRESSURE EQUALITY AT NODES 
    ##############################
    # A node has a single pressure, shared by every arc meeting it. At a branch point this is what
    # couples the two downstream arcs together.
    @constraint(model, [j in setdiff(node_idx, outlet_node), k in outgoing[j], t=1:T_sim], p_node[j, t] == p_in[k, t])
    @constraint(model, [j in setdiff(node_idx, inlet_node), k in incoming[j], t=1:T_sim], p_node[j, t] == p_out[k, t])

    ##############################
    # BOUNDARY CONDITIONS FOR SYSTEM
    ##############################
    # Downstream: linear reservoir inflow, driven by the bottom-hole to reservoir pressure difference.
    @constraint(model, [k in well_idx, t=1:T_sim], w_out[k, t] == (p_out[k, t] - p_reservoir[t])*injectivity[k, t])
    # Upstream: fixed supply pressure at the pump inlet. 
    @constraint(model, [t=1:T_sim], p_in[1, t] == P_IN)

    # FIX WELL INJECTIVITY AND RESERVOIR PRESSURE AT THEIR NOMINAL VALUES
    for k in well_idx, t in 1:T_sim; fix(injectivity[k, t], well_injectivity[k]; force=true); end
    for t in 1:T_sim; fix(p_reservoir[t], P_RES; force=true); end

    # OBJECTIVE FUNCTION 
    @objective(model, Min, 1)


    return model
end



