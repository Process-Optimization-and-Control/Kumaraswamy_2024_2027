####################################################
# STEADY-STATE MODEL
# The same first-principles network as model_mpc.jl with the accumulation terms removed, so there
# is no time dimension (t is a singleton, T_ss = 1, kept only so that the variable shapes match
# those of the NMPC and simulation models). It is solved once at the start of a run and serves
# two purposes:
#   1. its states and inputs are the initial condition the closed-loop run starts from
#   2. its full solution, repeated over the horizon, is the initial guess for the first NMPC and
#      simulation solves (case notebook)
# It is a feasibility problem, not an optimisation: the inputs are fixed and the objective is a
# constant, so Ipopt only has to find a consistent operating point.
####################################################

function steady_state()

    ##############################
    # STEADY STATE MODEL 
    ##############################
    model = Model(Ipopt.Optimizer)

    ##############################
    # DEFINE VARIABLES   
    ##############################
    # DEFINE GENERIC ARC VARIABLES 
    @variable(model, p_in[k=1:K, t=1:T_ss]>=0, start=steady_state_initial_guess["p_in"])
    @variable(model, p_out[k=1:K, t=1:T_ss]>=0, start=steady_state_initial_guess["p_out"])
    @variable(model, w_in[k=1:K, t=1:T_ss]>=0, start=steady_state_initial_guess["w_in"] )
    @variable(model, w_out[k=1:K, t=1:T_ss]>=0, start=steady_state_initial_guess["w_out"])

    # DEFINE PIPE SEGMENT VARIABLES AT INLET AND OUTLET 
    @variable(model, reynold_in[k in pipe_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["reynold_in"])
    @variable(model, reynold_out[k in pipe_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["reynold_out"])
    @variable(model, fdarcy_in[k in pipe_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["fdarcy_in"])
    @variable(model, fdarcy_out[k in pipe_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["fdarcy_out"])
    @variable(model, dpdz_in[k in pipe_idx, t=1:T_ss], start=steady_state_initial_guess["dpdz_in"])
    @variable(model, dpdz_out[k in pipe_idx, t=1:T_ss], start=steady_state_initial_guess["dpdz_out"])

    # DEFINE LUMPED PIPE SEGMENT VARIABLES 
    @variable(model, Z[k in pipe_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["Z"])
    @variable(model, Z_in[k in pipe_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["Z_in"])
    @variable(model, Z_out[k in pipe_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["Z_out"])
    @variable(model, rho[k in pipe_idx, t=1:T_ss]>=500, start=steady_state_initial_guess["rho"])
    @variable(model, rho_in[k in pipe_idx, t=1:T_ss]>=500, start=steady_state_initial_guess["rho_in"])
    @variable(model, rho_out[k in pipe_idx, t=1:T_ss]>=500, start=steady_state_initial_guess["rho_out"])
    @variable(model, p_average[k in pipe_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["p_average"])
    # DEFINE CHOKE VALVE VARIABLES 
    @variable(model, choke_vlv_op[k in choke_vlv_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["choke_vlv_op"])
    @variable(model, deltap_wellhead_choke[k in choke_vlv_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["deltap_wellhead_choke"])
    @variable(model, cavitation_idx[k in choke_vlv_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["cavitation_idx"])

    # DEFINE NODE VARIABLES 
    @variable(model, p_node[j in node_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["p_node"])

    # DEFINE PUMP VARIABLES 
    @variable(model, pump_work[k in pump_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["pump_work"])
    @variable(model, rho_pump_inlet[k in pump_idx, t=1:T_ss]>=200, start=steady_state_initial_guess["rho_pump_inlet"]) 
    @variable(model, rho_pump_outlet[k in pump_idx, t=1:T_ss]>=200, start=steady_state_initial_guess["rho_pump_outlet"]) 
    # speed_pump in particular must be started away from zero: the pump curve below is written in
    # affinity-scaled form with r = speed/3000 in the denominator, so a start of 0 makes head/r^2
    # blow up (~1e9) and Ipopt goes straight into restoration and converges to a locally
    # infeasible point. head_meters and q_pump are started near their expected values for the
    # same reason - to keep the first residual evaluation well scaled.
    @variable(model, head_meters[k in pump_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["head_meters"])
    @variable(model, q_pump[k in pump_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["q_pump"])
    @variable(model, speed_pump[k in pump_idx, t=1:T_ss]>=0, start=steady_state_initial_guess["speed_pump"])

    # DEFINE RESERVOIR PARAMETERS 
    # Well Injectivity & Reservoir Pressure 
    @variable(model, injectivity[k in well_idx, t=1:T_ss]>=0)
    @variable(model, p_reservoir[t=1:T_ss]>=0)

    ##############################
    # PUMP EQUATIONS    
    ##############################
    @expression(model, r[k in pump_idx, t=1:T_ss], speed_pump[k, t]/3000)
    @constraint(model, [k in pump_idx, t=1:T_ss], head_meters[k, t]/(r[k, t]^2) == -6*(10^-5)*(q_pump[k, t]/r[k, t])^2 - 1.0917*(q_pump[k, t]/r[k, t]) + 2859.6)
    @constraint(model, [k in pump_idx, t=1:T_ss], p_out[k, t] - p_in[k, t] == head_meters[k, t]*0.5*(rho_pump_inlet[k, t] + rho_pump_outlet[k, t])*G*(10^-5))
    @constraint(model, [k in pump_idx, t=1:T_ss], rho_pump_inlet[k, t] == p_in[k, t]*MW/R/poly2(p_in[k, t]*(10^5), fit_result_Z.param)/Tref)
    @constraint(model, [k in pump_idx, t=1:T_ss], rho_pump_outlet[k, t] == p_out[k, t]*MW/R/poly2(p_out[k, t]*(10^5), fit_result_Z.param)/Tref )
    @constraint(model, [k in pump_idx, t=1:T_ss], pump_work[k, t] == ( (p_out[k, t] - p_in[k, t])*(10^5)*w_in[k, t]/pump_efficiency[k]/(0.5*(rho_pump_inlet[k, t] + rho_pump_outlet[k, t])) )*DT/(3.6e9) ) # Pump work in MWh
    @constraint(model, [k in pump_idx, t=1:T_ss], w_in[k, t] == w_out[k, t]) # no accummulation across pump volumes 
    @constraint(model, [k in pump_idx, t=1:T_ss], q_pump[k, t] == (w_in[k, t]/rho_pump_inlet[k, t])*3600)
    #Note: q_pump is in m3/hour
    
    ##############################
    # MOMENTUM BALANCE EQUATIONS   
    ##############################
    @constraint(model, [k in pipe_idx, t=1:T_ss], dpdz_in[k, t] == -(fdarcy_in[k, t]*(w_in[k, t]^2)/(2*DINNER[k]*ACROSS[k]*ACROSS[k]*rho_in[k, t]))/1e5 - (rho_in[k, t]*G*sin(Theta[k]*pi/180))/1e5)
    @constraint(model, [k in pipe_idx, t=1:T_ss], dpdz_out[k, t] == -(fdarcy_out[k, t]*(w_out[k, t]^2)/(2*DINNER[k]*ACROSS[k]*ACROSS[k]*rho_out[k, t]))/1e5 - (rho_out[k, t]*G*sin(Theta[k]*pi/180))/1e5)
    @constraint(model, [k in pipe_idx, t=1:T_ss], p_out[k, t] == p_in[k, t] + (LENGTH[k]/2)*(dpdz_in[k, t] + dpdz_out[k, t]))

    ############################
    #  COMPRESSIBILITY FACTOR
    ############################
    # As in model_mpc.jl, but with the extra Z <= 1 bounds below to keep the feasibility solve from
    # wandering into a physically meaningless branch of the quadratic fit.
    @constraint(model, [k in pipe_idx, t=1:T_ss], Z[k, t] == poly2(p_average[k, t]*(10^5), fit_result_Z.param))
    @constraint(model, [k in pipe_idx, t=1:T_ss], Z_in[k, t] == poly2(p_in[k, t]*(10^5), fit_result_Z.param))
    @constraint(model, [k in pipe_idx, t=1:T_ss], Z_out[k, t] == poly2(p_out[k, t]*(10^5), fit_result_Z.param))
    @constraint(model, [k in pipe_idx, t=1:T_ss], rho[k, t] == p_average[k, t]*MW/R/Z[k, t]/Tref)
    @constraint(model, [k in pipe_idx, t=1:T_ss], rho_in[k, t] == p_in[k, t]*MW/R/Z_in[k, t]/Tref)
    @constraint(model, [k in pipe_idx, t=1:T_ss], rho_out[k, t] == p_out[k, t]*MW/R/Z_out[k, t]/Tref)
    @constraint(model, [k in pipe_idx, t=1:T_ss], Z[k, t] <= 1)
    @constraint(model, [k in pipe_idx, t=1:T_ss], Z_in[k, t] <= 1)
    @constraint(model, [k in pipe_idx, t=1:T_ss], Z_out[k, t] <= 1)

    ############################
    # DEFINE FRICTION FACTOR & REYNOLDS NUMBER   
    ############################
    # REYNOLDS NUMBER AND DARCY FRICTION FACTOR CALCULATION 
    @constraint(model, [k in pipe_idx, t=1:T_ss], reynold_in[k, t] ==  w_in[k, t]*DINNER[k]/(poly2(p_in[k, t]*(10^5), fit_result_mu.param)*ACROSS[k]))
    @constraint(model, [k in pipe_idx, t=1:T_ss], reynold_out[k, t] ==  w_out[k, t]*DINNER[k]/(poly2(p_out[k, t]*(10^5), fit_result_mu.param)*ACROSS[k]))
    @constraint(model, [k in pipe_idx, t=1:T_ss], fdarcy_in[k, t] == 1 / ( -1.8*log10(6.9/reynold_in[k, t] + (EPS/(DINNER[k]*3.7))^1.11) )^2)
    @constraint(model, [k in pipe_idx, t=1:T_ss], fdarcy_out[k, t] == 1 / ( -1.8*log10(6.9/reynold_out[k, t] + (EPS/(DINNER[k]*3.7))^1.11) )^2)

    ############################
    # AVERAGE PRESSURE CALCULATION    
    ############################
    @constraint(model, [k in pipe_idx, t=1:T_ss], p_average[k, t] == (0.5*p_in[k, t] + (1/12)*dpdz_in[k, t]*LENGTH[k] + 0.5*p_out[k, t] - (1/12)*dpdz_out[k, t]*LENGTH[k])/1)

    ##############################
    # VALVE EQUATIONS 
    ##############################
    @constraint(model, [k in all_vlv_idx, t=1:T_ss], (w_in[k, t]/(VLV_2_CV[k]*(choke_vlv_op[k, t])))^2/(rho_out[nearest_pipe_idx[k], t]) == deltap_wellhead_choke[k, t]*(10^5))
    @constraint(model, [k in all_vlv_idx, t=1:T_ss], p_in[k, t] == p_out[k, t] + deltap_wellhead_choke[k, t])
    @constraint(model, [k in all_vlv_idx, t=1:T_ss], w_in[k, t] == w_out[k, t])
    @constraint(model, [k in all_vlv_idx, t=1:T_ss], p_in[k, t] >= p_out[k, t])
    
    # CAVITATION INDEX CALCULATION FOR CHOKE VALVES
    # Reported for information only; unlike the NMPC there is no cavitation_min constraint here.
    @constraint(model, [k in all_vlv_idx, t=1:T_ss], cavitation_idx[k, t] == (p_in[k, t] - co2_vapor_pressure)/(p_in[k, t] - p_out[k, t]))

    ##############################
    # STEADY STATE EQUATIONS FOR PIPE SEGMENT
    ##############################
    # No accumulation, so what enters a pipe segment also leaves it.
    @constraint(model, [k in pipe_idx, t=1:T_ss], w_out[k, t] == w_in[k, t])

    ##############################
    # PIPELINE NETWORK - DEFINE JUNCTIONS / CONNECTIONS 
    # MASS BALANCE
    ##############################
    @constraint(model, [j in setdiff(node_idx, inlet_node, outlet_node), t=1:T_ss], sum(w_out[x, t] for x in incoming[j]) == sum(w_in[y, t] for y in outgoing[j]))

    ##############################
    # PIPELINE NETWORK
    # PRESSURE EQUALITY AT NODES 
    ##############################
    @constraint(model, [j in setdiff(node_idx, outlet_node), k in outgoing[j], t in 1:T_ss], p_node[j, t] == p_in[k, t])
    @constraint(model, [j in setdiff(node_idx, inlet_node), k in incoming[j], t in 1:T_ss], p_node[j, t] == p_out[k, t])

    ##############################
    # BOUNDARY CONDITIONS FOR SYSTEM 
    ##############################
    # Downstream: linear reservoir inflow, as in model_mpc.jl.
    @constraint(model, [k in well_idx, t=1:T_ss], w_out[k, t] == (p_out[k, t] - p_reservoir[t])*injectivity[k, t])
    # Upstream: BOTH the supply pressure and the total throughput are imposed here. Fixing the flow
    # is what pins down the operating point and hence the pump speed; the NMPC fixes only p_in and
    # lets the flow follow from the pump curve.
    for t in 1:T_ss; fix(p_in[1, t], P_IN; force=true); end
    for t in 1:T_ss; fix(w_in[1, t], W_IN; force=true); end

    ##############################
    # CHOKE VALVE OPENINGS AND PUMP SPEED LIMITS
    ##############################
    # The chokes are inputs here, not decisions: they are fixed at their nominal opening, leaving the
    # pump speed as the only degree of freedom the flow boundary condition can determine.
    for k in all_vlv_idx, t in 1:T_ss; fix(choke_vlv_op[k, t], CHOKE_VALVE_OP_SS[k]; force=true); end
    @constraint(model, [k in pump_idx, t=1:T_ss], speed_pump[k, t] <= speed_pump_max)
    @constraint(model, [k in pump_idx, t=1:T_ss], speed_pump[k, t] >= speed_pump_min)

    # FIX WELL INJECTIVITY AND RESERVOIR PRESSURE AT THEIR NOMINAL VALUES
    for k in well_idx, t in 1:T_ss; fix(injectivity[k, t], well_injectivity[k]; force=true); end
    for t in 1:T_ss; fix(p_reservoir[t], P_RES; force=true); end

    # Constant objective: this is a square feasibility problem, solved only for a consistent
    # operating point, so there is nothing to optimise.
    @objective(model, Min, 1)

    ##############################  
    # SOLVE MODEL 
    ############################## 
    set_optimizer_attribute(model, "print_level", 5)
    set_attribute(model, "linear_solver", "ma57")
        
    optimize!(model)    
    status = JuMP.termination_status(model)
    println(status)

    ##############################
    # PRINT ALL VARIABLE VALUES (comment out the loop below to silence)
    ##############################

    for var in all_variables(model)
        println("Variable $(name(var)) = ", value(var))
    end 

    return model

end



