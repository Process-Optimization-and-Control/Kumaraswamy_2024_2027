####################################################
# DESCRIPTION
####################################################
# PID control structure (Scheme A - baseline, NO overrides):
#  - Wells 2-4: wellhead flow control via choke valve opening (velocity-form PI).
#  - Well 1:    choke on manual (fixed opening = 1.0).
#  - Mainline:  pressure control at 500 m via pump speed (velocity-form PI).
# Controllers use the proportional-on-measurement velocity form; outputs are
# clamped to valve/speed limits, and the velocity form is based on the previous
# CLAMPED output, which provides anti-windup.
# NOTE: this is the plain regulatory layer. Neither the per-well high-BHP override
# (min-select on choke) nor the mainline min-flow override (max-select on pump
# speed) used in Schemes B/C/D is implemented here, so well BHP and mainline flow
# are measured but unused. Tuning matches Scheme D: wells tau_c = 120 s,
# pump pressure tau_c = 500 s.

function run_pid_controller(new_measurements, old_measurements, old_outputs, old_clamped_outputs, dt, i)

    mainline_flow, mainline_pressure, well1_flow, well1_press, well1_bhp, well2_flow, well2_press, well2_bhp, well3_flow, well3_press, well3_bhp, well4_flow, well4_press, well4_bhp = new_measurements 
    mainline_flow_prev, mainline_pressure_prev, well1_flow_prev, well1_press_prev, well1_bhp_prev, well2_flow_prev, well2_press_prev, well2_bhp_prev, well3_flow_prev, well3_press_prev, well3_bhp_prev, well4_flow_prev, well4_press_prev, well4_bhp_prev = old_measurements 
    
    ##############################
    # OLD CLAMPED OUTPUTS
    # The velocity form integrates from the previous CLAMPED output (anti-windup),
    # so `old_outputs` is not used here - it is kept in the signature because the
    # run_pid_scheme_*.ipynb driver passes it.
    ##############################
    choke_vlv_op_2_old_clamped = old_clamped_outputs[2]
    choke_vlv_op_3_old_clamped = old_clamped_outputs[3]
    choke_vlv_op_4_old_clamped = old_clamped_outputs[4]
    pump_speed_old_clamped = old_clamped_outputs[5]

    ##############################
    # MVs on MANUAL  
    ##############################   

    # WELL 1 ON MANUAL
    choke_vlv_op_1     = 1.0

    ##############################
    # SETPOINT SCHEDULE
    #   Matches Phase 0 -> Phase 1 of Scheme D so that A is a like-for-like
    #   "no overrides" comparison against D on the same disturbance.
    #     Phase 0 [0   -  500): baseline, pump 197.8, all wells 125
    #     Phase 1 [500 - 7200): pump pressure 197.8 -> 205  (settle, pump tau_c 500)
    #   Run with nt = 7200 (dt = 30 -> 241 steps).
    ##############################

    # Pump pressure setpoint (bar)
    if i >= 500
        PC_pump_SP = 205.0
    else
        PC_pump_SP = 197.8        # Phase 0 baseline
    end

    # Well flow setpoints (kg/s) - all wells held at nominal
    well_flow_sp = 125.0

    ##############################
    # CONTROLLER CALCULATION
    ##############################

    ##############################
    # TUNING PARAMETERS FOR PID CONTROLLERS 
    # Kc_well - Proportional gain for flow control at wellhead
    # tau_I_well - Integral time constant for flow control at wellhead (seconds)
    # Kc_speed - Proportional gain for pump speed control
    # tau_I_speed - Integral time constant for pump speed control (seconds)
    # (No override parameters in this scheme - see Schemes B/C/D.)
    # All gains are SIMC: Kc = (tau_1/k)/(tau_c + theta), tau_I = min(tau_1, 4*(tau_c + theta))
    ##############################

    # tau_c = 120s
    Kc_well_2 = 4.27079e-5
    tau_I_well_2 = 0.1
    FC2_SP     = well_flow_sp

    # tau_c = 120s
    Kc_well_3 = 4.27079e-5
    tau_I_well_3 = 0.1
    FC3_SP     = well_flow_sp

    # tau_c = 120s
    Kc_well_4 = 4.27079e-5
    tau_I_well_4 = 0.1
    FC4_SP     = well_flow_sp

    # tau_c = 500s
    Kc_speed = 29.5426
    tau_I_speed = 1039.0

    # (pump pressure setpoint is set in the SETPOINT SCHEDULE block above)

    ##############################
    # CONTROLLER ERRORS    
    ##############################
    FC2_error = FC2_SP - well2_flow
    FC3_error = FC3_SP - well3_flow
    FC4_error = FC4_SP - well4_flow
    PC_pump_error = PC_pump_SP - mainline_pressure
    println("FC2_error: $FC2_error, FC3_error: $FC3_error, FC4_error: $FC4_error, PC_pump_error: $PC_pump_error")
    
    ##############################
    # FLOW CONTROLLER (PID) - VELOCITY FORM
    # Proportional term on measurement
    ##############################
    FC2_OP = -Kc_well_2*(well2_flow - well2_flow_prev) + (Kc_well_2/tau_I_well_2)*dt*FC2_error + choke_vlv_op_2_old_clamped
    FC3_OP = -Kc_well_3*(well3_flow - well3_flow_prev) + (Kc_well_3/tau_I_well_3)*dt*FC3_error + choke_vlv_op_3_old_clamped
    FC4_OP = -Kc_well_4*(well4_flow - well4_flow_prev) + (Kc_well_4/tau_I_well_4)*dt*FC4_error + choke_vlv_op_4_old_clamped        
    PC_pump_OP = -Kc_speed*(mainline_pressure - mainline_pressure_prev) + (Kc_speed/tau_I_speed)*dt*PC_pump_error + pump_speed_old_clamped

    println("FC2_OP: $FC2_OP, FC3_OP: $FC3_OP, FC4_OP: $FC4_OP, PC_pump_OP: $PC_pump_OP")

    ##############################
    # Clamp valve and speed outputs to min/max values
    ##############################
    speed_min = 2800.0
    speed_max = 3900.0

    valve_min = 0.2
    valve_max = 1.0

    FC2_OP_clamped = clamp(FC2_OP, valve_min, valve_max)
    FC3_OP_clamped = clamp(FC3_OP, valve_min, valve_max)
    FC4_OP_clamped = clamp(FC4_OP, valve_min, valve_max)
    PC_pump_OP_clamped = clamp(PC_pump_OP, speed_min, speed_max)

    outputs = [choke_vlv_op_1, FC2_OP, FC3_OP, FC4_OP, PC_pump_OP]
    clamped_outputs = [choke_vlv_op_1, FC2_OP_clamped, FC3_OP_clamped, FC4_OP_clamped, PC_pump_OP_clamped]

    println("outputs: $outputs, clamped_outputs: $clamped_outputs")
    return outputs, clamped_outputs 
end 