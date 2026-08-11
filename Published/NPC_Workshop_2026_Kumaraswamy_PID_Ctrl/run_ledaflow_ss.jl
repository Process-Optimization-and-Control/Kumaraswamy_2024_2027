##############################
# RUN LEDAFLOW SIMULATION AT STEADY-STATE
##############################

ss_output_csv_path = joinpath(@__DIR__, "ss_lf_output.csv")

js = """
var calcModule = ledaModules.CALCULATE();
var fullRModule = ledaModules.FULLRESULTS();
var caseID = "$lf_case_id";  

//////// Reset case to initial conditions ////////

var cc = new Compound(caseID);
var Valve = cc.relation("child", "Leda1DnPhValve", true)
var Pump = cc.relation("child", "Leda1DnPhPump", true)

Valve[0].setOpenFrac([1.0])
Valve[1].setOpenFrac([1.0])
Valve[2].setOpenFrac([1.0])
Valve[3].setOpenFrac([1.0])

Pump[0].setSpeed([3458 * 0.10472])  // convert from rpm to rad/s

//////// Running steady state simulation //////////

// Purge results keeping only the first global sample without resetting the time
calcModule.purge("KeepFirst", 1, caseID);

var currentCase = new LedaGeneralCase(caseID);
var timeAdvance = 100;
currentCase.mycase.setTimeAdvance(timeAdvance);

calcModule.setCaseId(caseID);
calcModule.calculate();


//////// Exctracting results //////////

fullRModule.setUuid(caseID);
fullRModule.trendLoggersToSvFile("$ss_output_csv_path");

"""

softsh = resolve_softsh_path()
script = joinpath(@__DIR__, "ledaflow_ss.js")
write(script, js)
logfile = joinpath(@__DIR__, "ledaflow_ss.log")
run(pipeline(`$softsh $script`, stdout=logfile, stderr=logfile, append=true))

