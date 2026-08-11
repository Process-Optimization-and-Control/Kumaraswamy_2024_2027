####################################################
# IMPLEMENT PID CONTROLLERS EVERY TIME STEP AND EXTRACT RESULTS TO CSV FILE   
####################################################


function run_ledaflow_sim(clamped_outputs, old_clamped_outputs, dt, i)

    ##############################
    # UNPACK INPUTS 
    ##############################
    choke_vlv_op_1_prev = old_clamped_outputs[1]
    choke_vlv_op_2_prev = old_clamped_outputs[2]
    choke_vlv_op_3_prev = old_clamped_outputs[3]
    choke_vlv_op_4_prev = old_clamped_outputs[4]
        
    choke_vlv_op_1 = clamped_outputs[1]
    choke_vlv_op_2 = clamped_outputs[2]
    choke_vlv_op_3 = clamped_outputs[3]
    choke_vlv_op_4 = clamped_outputs[4]
    pump_speed = clamped_outputs[5]

    println("choke_vlv_op_1: $choke_vlv_op_1, choke_vlv_op_2: $choke_vlv_op_2, choke_vlv_op_3: $choke_vlv_op_3, choke_vlv_op_4: $choke_vlv_op_4, pump_speed: $pump_speed")
    println("choke_vlv_op_1_prev: $choke_vlv_op_1_prev, choke_vlv_op_2_prev: $choke_vlv_op_2_prev, choke_vlv_op_3_prev: $choke_vlv_op_3_prev, choke_vlv_op_4_prev: $choke_vlv_op_4_prev")
    output_csv_path = joinpath(@__DIR__, "lf_output.csv")

    js = """
        var calcModule = ledaModules.CALCULATE();
        var fullRModule = ledaModules.FULLRESULTS();

        var caseID = "$lf_case_id";  
        var cc = new Compound(caseID);
        var Valve = cc.relation("child", "Leda1DnPhValve", true)
        var Pump = cc.relation("child", "Leda1DnPhPump", true)
        
        Valve[0].setOpenFrac([$choke_vlv_op_1])
        Valve[1].setOpenFrac([$choke_vlv_op_2])
        Valve[2].setOpenFrac([$choke_vlv_op_3])
        Valve[3].setOpenFrac([$choke_vlv_op_4])
        Pump[0].setSpeed([$pump_speed * 0.10472])

        //////// Running dynamic simulation of one time step //////////

        var currentCase = new LedaGeneralCase(caseID);
        var timeAdvance = $dt;
        currentCase.mycase.setTimeAdvance(timeAdvance);

        calcModule.setCaseId(caseID);
        calcModule.calculate();

        //////// Exctracting results //////////

        fullRModule.setUuid(caseID);
        fullRModule.trendLoggersToSvFile("$output_csv_path");
    """

    softsh = resolve_softsh_path()
    script = joinpath(@__DIR__, "ledaflow_write.js")
    write(script, js)
    logfile = joinpath(@__DIR__, "ledaflow_write.log")
    run(pipeline(`$softsh $script`, stdout=logfile, stderr=logfile, append=true))

    df = CSV.read(output_csv_path, DataFrame;
        header=10,            # column names start  from row 10  
        skipto=12,            # skip the units row (line 11)  
        delim = ',',          # LedaFlow export is comma-delimited
        missingstring = ["--"], # LedaFlow writes "--" for no-data cells; treat as missing
        normalizenames =false, # keep names like "Pressure@Line 1 P 500m"
        )

    mainline_pressure = df[end, "Pressure@Mainline P 500m"]  
    mainline_flow = df[end, "MFR - total@Mainline P 19500m"]*(-1)
    well1_flow = df[end, "MFR - total liquid@Wellbore1 P MFR"]*(-1)
    well2_flow = df[end, "MFR - total liquid@Wellbore2 P MFR"]*(-1)
    well3_flow = df[end, "MFR - total liquid@Wellbore3 P MFR"]*(-1)
    well4_flow = df[end, "MFR - total liquid@Wellbore4 P MFR"]*(-1)
    well1_press = df[end, "Pressure@Wellbore1 P MFR"]
    well2_press = df[end, "Pressure@Wellbore2 P MFR"]
    well3_press = df[end, "Pressure@Wellbore3 P MFR"]
    well4_press = df[end, "Pressure@Wellbore4 P MFR"]
    well1_bhp = df[end, "BHP - Zone 1@Well 1"]
    well2_bhp = df[end, "BHP - Zone 1@Well 2"]
    well3_bhp = df[end, "BHP - Zone 1@Well 3"]
    well4_bhp = df[end, "BHP - Zone 1@Well 4"]

    new_measurements = [mainline_flow, mainline_pressure, well1_flow, well1_press, well1_bhp, well2_flow, well2_press,well2_bhp, well3_flow, well3_press, well3_bhp, well4_flow, well4_press, well4_bhp]

    return new_measurements
end 