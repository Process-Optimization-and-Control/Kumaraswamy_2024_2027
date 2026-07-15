
function run_simulator(optimal_input, prev_input, sim_t, last_time_point, density_chokes, w_pred, w_mainline_pred, pressure_pred, bias_well, bias_mainline, bias_pressure)

    ##############################
    # EXTRACT OPTIMAL AND PREVIOUS INPUTS FROM MPC 
    ##############################  

    choke_valve_op_1 = optimal_input["choke_vlv_op"][Name.(50), 1]
    choke_valve_op_2 = optimal_input["choke_vlv_op"][Name.(52), 1]
    choke_valve_op_3 = optimal_input["choke_vlv_op"][Name.(63), 1]
    choke_valve_op_4 = optimal_input["choke_vlv_op"][Name.(65), 1]
    pump_speed = optimal_input["speed_pump"][Name.(1), 1]

    choke_valve_op_1_prev = prev_input["choke_vlv_op"][Name.(50), 1]
    choke_valve_op_2_prev = prev_input["choke_vlv_op"][Name.(52), 1]
    choke_valve_op_3_prev = prev_input["choke_vlv_op"][Name.(63), 1]
    choke_valve_op_4_prev = prev_input["choke_vlv_op"][Name.(65), 1]
    pump_speed_prev = prev_input["speed_pump"][Name.(1), 1]

    ##############################
    # CREATE RAMP PROFILE
    ##############################

    ramp_duration = DT        # seconds to ramp from prev -> optimal
    dt_ramp       = 1          # seconds, supply an opening every second

    # times start at last_time_point and advance 300 s (inclusive endpoint)
    times = collect(last_time_point:dt_ramp:(last_time_point + ramp_duration))

    # linear ramp from prev to optimal over the interval, capped to [0, 1]
    ramp(prev, opt) = clamp.(prev .+ (opt - prev) .* (times .- last_time_point) ./ ramp_duration, 0.0, 1.0)

    openings_1 = ramp(choke_valve_op_1_prev, choke_valve_op_1)
    openings_2 = ramp(choke_valve_op_2_prev, choke_valve_op_2)
    openings_3 = ramp(choke_valve_op_3_prev, choke_valve_op_3)
    openings_4 = ramp(choke_valve_op_4_prev, choke_valve_op_4)

    # format Julia vectors as JS array literals
    js_arr(v) = "[" * join(v, ", ") * "]"
    times_js      = js_arr(times)
    openings_1_js = js_arr(openings_1)
    openings_2_js = js_arr(openings_2)
    openings_3_js = js_arr(openings_3)
    openings_4_js = js_arr(openings_4)

    ##############################
    # WRITE JAVASCRIPT TO RUN LEDAFLOW SIMULATION
    ##############################
    trends_csv_path = joinpath(case_dir, "trends.csv")

    js = """
    var calcModule = ledaModules.CALCULATE();
    var fullRModule = ledaModules.FULLRESULTS();

    var caseID = "$lf_case_id";  
    var cc = new Compound(caseID);
    var Valve = cc.relation("child", "Leda1DnPhValve", true)
    var Pump = cc.relation("child", "Leda1DnPhPump", true)

    Valve[0].setStime($times_js);    Valve[0].setOpenFrac($openings_1_js)
    Valve[1].setStime($times_js);    Valve[1].setOpenFrac($openings_2_js)
    Valve[2].setStime($times_js);    Valve[2].setOpenFrac($openings_3_js)
    Valve[3].setStime($times_js);    Valve[3].setOpenFrac($openings_4_js)
    Pump[0].setSpeed([$pump_speed * 0.10472])

    //////// Running dynamic simulation of one time step //////////

    var currentCase = new LedaGeneralCase(caseID);
    var timeAdvance = $DT;
    currentCase.mycase.setTimeAdvance(timeAdvance);

    calcModule.setCaseId(caseID);
    calcModule.calculate();

    //////// Exctracting pressure results //////////

    fullRModule.setUuid(caseID);
    fullRModule.trendLoggersToSvFile("$trends_csv_path", ["Mainline P 500m", 
               "Mainline P 38500m", 
               "Line 1 P 500m",
               "Line 1 P 8500m",
               "Line 2 P 500m",
               "Line 2 P 8500m",
               "Wellbore1 P 600m",
               "Wellbore2 P 600m",
               "Wellbore3 P 600m",
               "Wellbore4 P 600m",
               "Wellbore1 P MFR",
               "Wellbore2 P MFR",
               "Wellbore3 P MFR",
               "Wellbore4 P MFR",
               "Well 1", 
               "Well 2", 
               "Well 3", 
               "Well 4",
               "WH Choke 1",
               "WH Choke 2",
               "WH Choke 3",
               "WH Choke 4"
]);
    """
    ##############################
    # RUN LEDAFLOW SIMULATION 
    ##############################
    softsh = resolve_softsh_path()
    script = joinpath(@__DIR__, "ledaflow_write.js")
    write(script, js)
    logfile = joinpath(case_dir, "leda.log")
    run(pipeline(`$softsh $script`, stdout=logfile, stderr=logfile, append=true))

    ############################## 
    # STATE ESTIMATION 
    # EXTRACT P AVERAGE LEDAFLOW RESULTS AND ASSIGN TO NEW INITIAL CONDITIONS
    ##############################
    df = CSV.read(trends_csv_path, DataFrame;
    header=10,            # column names start  from row 10  
    skipto=12,            # skip the units row (line 11)  
    delim = ',',
    normalizenames =false # keep names like "Pressure@Line 1 P 500m" 
    )

    p_average_mainline = Vector(df[end, ["Pressure@Mainline P 500m", "Pressure@Mainline P 38500m"]])   # size: n_time
    p_average_line1 = Vector(df[end, ["Pressure@Line 1 P 500m", "Pressure@Line 1 P 8500m"]])   # size: n_time
    p_average_line2 = Vector(df[end, ["Pressure@Line 2 P 500m", "Pressure@Line 2 P 8500m"]])   # size: n_time

    interp_mainline = linear_interpolation([2, 40], p_average_mainline)
    interp_line1 = linear_interpolation([41, 49], p_average_line1)
    interp_line2 = linear_interpolation([54, 62], p_average_line2)

    new_initial_con = NamedArray(fill(0.0, length(pipe_idx), 1, 1), (pipe_idx, 1:1, 1:1))
    new_initial_con[Name.(2:40), 1, 1] = interp_mainline(2:40)
    new_initial_con[Name.(41:49), 1, 1] = interp_line1(41:49)
    new_initial_con[Name.(54:62), 1, 1] = interp_line2(54:62)
    new_initial_con[Name.(well_idx), 1, 1] = Vector(df[end, ["Pressure@Wellbore1 P 600m", "Pressure@Wellbore2 P 600m", "Pressure@Wellbore3 P 600m", "Pressure@Wellbore4 P 600m"]])
    

    ############################## 
    # BIAS ESTIMATION  
    ##############################

    # Measurements of Flow going out of the Choke Valve 
    w_meas_vec = Vector(df[end, [
        "MFR - total liquid@Wellbore1 P MFR", 
        "MFR - total liquid@Wellbore2 P MFR", 
        "MFR - total liquid@Wellbore3 P MFR", 
        "MFR - total liquid@Wellbore4 P MFR"]]).*(-1)
    
    w_meas = NamedArray(w_meas_vec, (well_idx,))

    # Calculate Bias for Flow Rate through Choke Valves 
    L = well_flow_bias_filter
    # First Order Filter on the Bias Term 
    bias_well = (1 - L).*parent(bias_well) .+ L.*(parent(w_meas) .- parent(w_pred))
    bias_well = NamedArray(bias_well, (well_idx,)) # convert back to NamedArray with proper indexing

    # Measurements of Flow going into the Mainline 
    w_mainline_meas = df[end, "MFR - total liquid@Mainline P 500m"]
    # Calculate Bias for Mainline Flow Measurement
    L_mainline = mainline_flow_bias_filter
    bias_mainline =(1 - L_mainline)*bias_mainline + L_mainline*(w_mainline_meas - w_mainline_pred)

    # Calculate Bias for Pressure Measurements at Constrained Nodes
    # Bottom hole pressure measurements for each well
    bhp_meas_vec = Vector(df[end, [
        "BHP - Zone 1@Well 1", 
        "BHP - Zone 1@Well 2", 
        "BHP - Zone 1@Well 3", 
        "BHP - Zone 1@Well 4"]])

    pump_outlet_press = p_average_mainline[1] # use pressure at 500m as proxy for pump outlet pressure
    pressure_con_meas_vec = vcat(pump_outlet_press, bhp_meas_vec) # size: n_pressure_constrained_nodes
    L_pressure = pressure_bias_filter
    bias_pressure = (1 - L_pressure).*parent(bias_pressure) .+ L_pressure.*(parent(pressure_con_meas_vec) .- parent(pressure_pred))
    bias_pressure = NamedArray(bias_pressure, (pressure_con_idx,)) # convert back to NamedArray with proper indexing

    ############################## 
    # VALVE CV ESTIMATION   
    ##############################
    # Pressure Drop Measurements 
    p_drop_meas_vec = Vector(df[end, [
        "Pressure drop@WH Choke 1", 
        "Pressure drop@WH Choke 2", 
        "Pressure drop@WH Choke 3", 
        "Pressure drop@WH Choke 4"]]).*(-1)

    # Valve Openings Applied in the Simulation 
    choke_openings = [choke_valve_op_1, choke_valve_op_2, choke_valve_op_3, choke_valve_op_4]

    # Estimate Valve CVs based on the orifice equation: w = CV * opening * sqrt(pressure_drop * density)
    VLV_2_CV_est = vcat(zeros(K, 1))
    VLV_2_CV_est[choke_vlv_idx] = w_meas_vec./(choke_openings .* sqrt.(p_drop_meas_vec.*density_chokes.*1e5))

    ############################## 
    # INJECTIVITY ESTIMATION   
    ##############################
    # Measure flow rate injected into the reservoir
    w_inj_meas_vec = Vector(df[end, [
        "Total MFR - total@Well 1", 
        "Total MFR - total@Well 2", 
        "Total MFR - total@Well 3", 
        "Total MFR - total@Well 4"]])


    # Estimate Injectivity 
    inj_est = w_inj_meas_vec .* (-1)  ./  (bhp_meas_vec .- P_RES)
    inj_est = NamedArray(inj_est, (well_idx,))

    ############################## 
    # OVERRIDE BIAS AND PARAMETER ESTIMATES DEPENDING ON SETTINGS
    # Boolean values are obtained from case_directories 
    ##############################

    if inj_est_bool == false
        inj_est = well_injectivity[well_idx]
        inj_est = NamedArray(inj_est, (well_idx,))
    end

    if VLV_2_CV_est_bool == false
        VLV_2_CV_est = VLV_2_CV
    end

    if well_flow_bias_bool == false
        bias_well = zeros(length(well_idx))
        bias_well = NamedArray(bias_well, (well_idx,))
    end 

    if mainline_flow_bias_bool == false
        bias_mainline = 0.0
    end

    if pressure_bias_bool == false
        bias_pressure = zeros(length(pressure_con_idx))
        bias_pressure = NamedArray(bias_pressure, (pressure_con_idx,))
    end

    ############################## 
    # UPDATE LAST TIME POINT FOR NEXT SIMULATION RUN
    ############################## 
    last_time_point = df[end, "Time"]


    return new_initial_con, inj_est, VLV_2_CV_est, bias_well, bias_mainline, bias_pressure, last_time_point 
end



