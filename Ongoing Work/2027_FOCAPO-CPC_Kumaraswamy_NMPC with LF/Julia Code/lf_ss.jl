function run_ledaflow_ss(ss_model)

    # Extract Steady-State Pump and Choke Settings from the Steady State Model
    ss_pump_speed   = value.(ss_model[:speed_pump][1, 1])
    ss_choke_vlv_op_1 = value.(ss_model[:choke_vlv_op][choke_vlv_idx[1], 1])  
    ss_choke_vlv_op_2 = value.(ss_model[:choke_vlv_op][choke_vlv_idx[2], 1])  
    ss_choke_vlv_op_3 = value.(ss_model[:choke_vlv_op][choke_vlv_idx[3], 1])  
    ss_choke_vlv_op_4 = value.(ss_model[:choke_vlv_op][choke_vlv_idx[4], 1])  

    ss_csv_path = joinpath(case_dir, "steady_state_trends.csv")

    js = """
    var calcModule = ledaModules.CALCULATE();
    var fullRModule = ledaModules.FULLRESULTS();

    var caseID = "$lf_case_id";  

    //////// Reset case to initial conditions ////////

    var cc = new Compound(caseID);
    var Valve = cc.relation("child", "Leda1DnPhValve", true)
    var Pump = cc.relation("child", "Leda1DnPhPump", true)

    Valve[0].setOpenFrac($ss_choke_vlv_op_1)
    Valve[1].setOpenFrac($ss_choke_vlv_op_2)
    Valve[2].setOpenFrac($ss_choke_vlv_op_3)
    Valve[3].setOpenFrac($ss_choke_vlv_op_4)

    Pump[0].setSpeed([$ss_pump_speed * 0.10472])  // convert from rpm to rad/s

    //////// Running steady state simulation //////////

    // Purge results keeping only the first global sample without resetting the time
    calcModule.purge("KeepFirst", 1, caseID);

    var currentCase = new LedaGeneralCase(caseID);
    var timeAdvance = 2000;
    currentCase.mycase.setTimeAdvance(timeAdvance);

    calcModule.setCaseId(caseID);
    calcModule.calculate();


    //////// Exctracting results //////////

    fullRModule.setUuid(caseID);
    fullRModule.trendLoggersToSvFile("$ss_csv_path");

    """

    softsh = resolve_softsh_path()
    script = joinpath(@__DIR__, "ledaflow_ss.js")
    write(script, js)
    logfile = joinpath(case_dir, "leda_ss.log")
    run(pipeline(`$softsh $script`, stdout=logfile, stderr=logfile, append=true))

end 