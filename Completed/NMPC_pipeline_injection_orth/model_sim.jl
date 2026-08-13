####################################################
# SIMULATION MODEL
# The plant. It is the same first-principles network as model_mpc.jl, discretised by the same
# orthogonal collocation scheme, but solved over a single finite element (T_sim = 1) with the
# inputs fixed to the move the NMPC just returned. It is a square feasibility problem: there is
# nothing to optimise, only a consistent trajectory over the element to find.
####################################################

function sim_model_build()

    ##############################
    # SIMULATION MODEL
    ##############################
    model = Model(Ipopt.Optimizer)

    ##############################
    # DEFINE VARIABLES
    ##############################
    # DEFINE GENERIC ARC VARIABLES
    @variable(model, p_in[k=1:K, d=2:D, t=1:T_sim]>=0)
    @variable(model, p_out[k=1:K, d=2:D, t=1:T_sim]>=0)
    @variable(model, w_in[k=1:K, d=2:D, t=1:T_sim]>=0)
    @variable(model, w_out[k=1:K, d=2:D, t=1:T_sim]>=0)

    # DEFINE PIPE SEGMENT VARIABLES AT INLET AND OUTLET
    @variable(model, reynold_in[k in pipe_idx, d=2:D, t=1:T_sim]>=0)
    @variable(model, reynold_out[k in pipe_idx, d=2:D, t=1:T_sim]>=0)
    @variable(model, fdarcy_in[k in pipe_idx, d=2:D, t=1:T_sim]>=0)
    @variable(model, fdarcy_out[k in pipe_idx, d=2:D, t=1:T_sim]>=0)
    @variable(model, dpdz_in[k in pipe_idx, d=2:D, t=1:T_sim])
    @variable(model, dpdz_out[k in pipe_idx, d=2:D, t=1:T_sim])

    # DEFINE LUMPED PIPE SEGMENT VARIABLES
    @variable(model, Z[k in pipe_idx, d=2:D, t=1:T_sim]>=0)
    @variable(model, Z_in[k in pipe_idx, d=2:D, t=1:T_sim]>=0)
    @variable(model, Z_out[k in pipe_idx, d=2:D, t=1:T_sim]>=0)
    @variable(model, rho[k in pipe_idx, d=1:D, t=1:T_sim]>=500)
    @variable(model, rho_in[k in pipe_idx, d=2:D, t=1:T_sim]>=500)
    @variable(model, rho_out[k in pipe_idx, d=2:D, t=1:T_sim]>=500)
    @variable(model, p_average[k in pipe_idx, d=2:D, t=1:T_sim]>=0)

    # DEFINE CHOKE VALVE VARIABLES
    @variable(model, deltap_wellhead_choke[k in all_vlv_idx, d=2:D, t=1:T_sim]>=0)
    @variable(model, cavitation_idx[k in choke_vlv_idx, d=2:D, t=1:T_sim]>=0)

    # DEFINE NODE VARIABLES
    @variable(model, p_node[j in node_idx, d=2:D, t=1:T_sim]>=0)

    # DEFINE PUMP VARIABLES
    @variable(model, pump_work[k in pump_idx, d=2:D, t=1:T_sim]>=0)
    @variable(model, rho_pump_inlet[k in pump_idx, d=2:D, t=1:T_sim]>=200)
    @variable(model, rho_pump_outlet[k in pump_idx, d=2:D, t=1:T_sim]>=200)
    @variable(model, head_meters[k in pump_idx, d=2:D, t=1:T_sim]>=0)
    @variable(model, q_pump[k in pump_idx, d=2:D, t=1:T_sim]>=0)

    # DEFINE RESERVOIR PARAMETERS
    # Well Injectivity & Reservoir Pressure.
    # These two are declared as variables only so that they can be fixed to a new value at each
    # sample without rebuilding the model; they are never free for the optimiser.
    @variable(model, injectivity[k in well_idx, d=2:D, t=1:T_sim]>=0)
    @variable(model, p_reservoir[d=2:D, t=1:T_sim]>=0)

    # DEFINE MPC INPUTS - no collocation variables
    # These are fixed by run_sim.jl to the move the NMPC returned and held across the element.
    @variable(model, choke_vlv_op[k in all_vlv_idx, t=1:T_sim]>=0)
    @variable(model, speed_pump[k in pump_idx, t=1:T_sim]>=0)

    ##############################
    # COLLOCATION EQUATION FOR MASS BALANCE
    ##############################
    @constraint(model, [k in pipe_idx, t=1:T_sim],
    Adot'*rho[k, 1:D, t]  .== (DT / V[k]) * (w_in[k, 2:D, t] - w_out[k, 2:D, t]))

    ##############################
    # PUMP EQUATIONS
    ##############################
    # r is the speed ratio against the 3000 rpm reference speed at which the pump curve was fitted.
    # The head/flow relation is applied in affinity-scaled form, head/r^2 as a quadratic in q/r, so
    # that one curve covers the whole speed range.
    @expression(model, r[k in pump_idx, t=1:T_sim], speed_pump[k, t]/3000)
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_sim], head_meters[k, d, t]/(r[k, t]^2) == -6*(10^-5)*(q_pump[k, d, t]/r[k, t])^2 - 1.0917*(q_pump[k, d, t]/r[k, t]) + 2859.6)
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_sim], p_out[k, d, t] - p_in[k, d, t] == head_meters[k, d, t]*0.5*(rho_pump_inlet[k, d, t] + rho_pump_outlet[k, d, t])*G*(10^-5))
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_sim], rho_pump_inlet[k, d, t] == p_in[k, d, t]*MW/R/poly2(p_in[k, d, t]*(10^5), fit_result_Z.param)/Tref)
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_sim], rho_pump_outlet[k, d, t] == p_out[k, d, t]*MW/R/poly2(p_out[k, d, t]*(10^5), fit_result_Z.param)/Tref)
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_sim], pump_work[k, d, t] == ( (p_out[k, d, t] - p_in[k, d, t])*(10^5)*w_in[k, d, t]/pump_efficiency[k]/(0.5*(rho_pump_inlet[k, d, t] + rho_pump_outlet[k, d, t])) )*DT/(3.6e9) ) # Pump work in MWh
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_sim], w_in[k, d, t] == w_out[k, d, t]) # no accummulation across pump volumes
    @constraint(model, [k in pump_idx, d=2:D, t=1:T_sim], q_pump[k, d, t] == (w_in[k, d, t]/rho_pump_inlet[k, d, t])*3600) # m3/hour

    ##############################
    # MOMENTUM BALANCE EQUATIONS
    ##############################
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_sim], dpdz_in[k, d, t] == -(fdarcy_in[k, d, t]*(w_in[k, d, t]^2)/(2*DINNER[k]*ACROSS[k]*ACROSS[k]*rho_in[k, d, t]))/1e5 - (rho_in[k, d, t]*G*sin(Theta[k]*pi/180))/1e5)
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_sim], dpdz_out[k, d, t] == -(fdarcy_out[k, d, t]*(w_out[k, d, t]^2)/(2*DINNER[k]*ACROSS[k]*ACROSS[k]*rho_out[k, d, t]))/1e5 - (rho_out[k, d, t]*G*sin(Theta[k]*pi/180))/1e5)
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_sim], p_out[k, d, t] == p_in[k, d, t] + (LENGTH[k]/2)*(dpdz_in[k, d, t] + dpdz_out[k, d, t]))

    ############################
    #  COMPRESSIBILITY FACTOR
    ############################
    # Z is taken from the quadratic CoolProp fit at the local pressure, and density follows from the
    # real gas law. Evaluated separately at the segment average, inlet and outlet.
    # NOTE: model_ss.jl additionally bounds Z <= 1; that bound is deliberately not imposed here.
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_sim], Z[k, d, t] == poly2(p_average[k, d, t]*(10^5), fit_result_Z.param))
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_sim], Z_in[k, d, t] == poly2(p_in[k, d, t]*(10^5), fit_result_Z.param))
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_sim], Z_out[k, d, t] == poly2(p_out[k, d, t]*(10^5), fit_result_Z.param))
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_sim], rho[k, d, t] == p_average[k, d, t]*MW/R/Z[k, d, t]/Tref)
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_sim], rho_in[k, d, t] == p_in[k, d, t]*MW/R/Z_in[k, d, t]/Tref)
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_sim], rho_out[k, d, t] == p_out[k, d, t]*MW/R/Z_out[k, d, t]/Tref)

    ############################
    # DEFINE FRICTION FACTOR & REYNOLDS NUMBER
    ############################
    # REYNOLDS NUMBER AND DARCY FRICTION FACTOR CALCULATION
    # Viscosity comes from the quadratic CoolProp fit.
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_sim], reynold_in[k, d, t] ==  w_in[k, d, t]*DINNER[k]/(poly2(p_in[k, d, t]*(10^5), fit_result_mu.param)*ACROSS[k]))
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_sim], reynold_out[k, d, t] ==  w_out[k, d, t]*DINNER[k]/(poly2(p_out[k, d, t]*(10^5), fit_result_mu.param)*ACROSS[k]))
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_sim], fdarcy_in[k, d, t] == 1 / ( -1.8*log10(6.9/reynold_in[k, d, t] + (EPS/(DINNER[k]*3.7))^1.11) )^2)
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_sim], fdarcy_out[k, d, t] == 1 / ( -1.8*log10(6.9/reynold_out[k, d, t] + (EPS/(DINNER[k]*3.7))^1.11) )^2)

    ############################
    # AVERAGE PRESSURE CALCULATION
    ############################
    @constraint(model, [k in pipe_idx, d=2:D, t=1:T_sim], p_average[k, d, t] == (0.5*p_in[k, d, t] + (1/12)*dpdz_in[k, d, t]*LENGTH[k] + 0.5*p_out[k, d, t] - (1/12)*dpdz_out[k, d, t]*LENGTH[k])/1)

    ##############################
    # VALVE EQUATIONS
    ##############################
    # Valves are treated as zero-volume elements: the flow passes straight through and the only
    # effect is a pressure drop. The drop is forced to be non-negative so that flow cannot reverse
    # through a choke.
    # Orifice equation w = VLV_2_CV*opening*sqrt(rho*dP).
    # The density is taken from the outlet of the pipe segment immediately upstream of the valve.
    @constraint(model, [k in all_vlv_idx, d=2:D, t=1:T_sim], (w_in[k, d, t]/(VLV_2_CV[k]*(choke_vlv_op[k, t])))^2/(rho_out[nearest_pipe_idx[k], d, t]) == deltap_wellhead_choke[k, d, t]*(10^5))
    @constraint(model, [k in all_vlv_idx, d=2:D, t=1:T_sim], p_in[k, d, t] == p_out[k, d, t] + deltap_wellhead_choke[k, d, t])
    @constraint(model, [k in all_vlv_idx, d=2:D, t=1:T_sim], w_in[k, d, t] == w_out[k, d, t])
    @constraint(model, [k in all_vlv_idx, d=2:D, t=1:T_sim], p_in[k, d, t] >= p_out[k, d, t])

    # CAVITATION INDEX CALCULATION FOR CHOKE VALVES
    # Cavitation index = (p_in - p_vapour)/(p_in - p_out)
    @constraint(model, [k in choke_vlv_idx, d=2:D, t=1:T_sim], cavitation_idx[k, d, t] == (p_in[k, d, t] - co2_vapor_pressure)/(p_in[k, d, t] - p_out[k, d, t]))

    ##############################
    # PIPELINE NETWORK - DEFINE JUNCTIONS / CONNECTIONS
    # MASS BALANCE
    ##############################
    # Interior nodes have no volume, so inflow equals outflow.
    @constraint(model, [j in setdiff(node_idx, inlet_node, outlet_node), d=2:D, t=1:T_sim], sum(w_out[x, d, t] for x in incoming[j]) == sum(w_in[y, d, t] for y in outgoing[j]))

    ##############################
    # PIPELINE NETWORK
    # PRESSURE EQUALITY AT NODES
    ##############################
    # A node has a single pressure, shared by every arc meeting it. At a branch point this is what
    # couples the two downstream arcs together.
    @constraint(model, [j in setdiff(node_idx, outlet_node), k in outgoing[j], d=2:D, t=1:T_sim], p_node[j, d, t] == p_in[k, d, t])
    @constraint(model, [j in setdiff(node_idx, inlet_node), k in incoming[j], d=2:D, t=1:T_sim], p_node[j, d, t] == p_out[k, d, t])

    ##############################
    # BOUNDARY CONDITIONS FOR SYSTEM
    ##############################
    # Downstream: linear reservoir inflow, driven by the bottom-hole to reservoir pressure difference.
    @constraint(model, [k in well_idx, d=2:D, t=1:T_sim], w_out[k, d, t] == (p_out[k, d, t] - p_reservoir[d, t])*injectivity[k, d, t])
    # Upstream: fixed supply pressure at the pump inlet.
    @constraint(model, p_in[1, 2:D, :] .== P_IN)

    # FIX WELL INJECTIVITY AND RESERVOIR PRESSURE AT THEIR NOMINAL VALUES
    for d in 2:D, k in well_idx, t in 1:T_sim; fix(injectivity[k, d, t], well_injectivity[k]; force=true); end
    for d in 2:D, t in 1:T_sim; fix(p_reservoir[d, t], P_RES; force=true); end

    return model
end
