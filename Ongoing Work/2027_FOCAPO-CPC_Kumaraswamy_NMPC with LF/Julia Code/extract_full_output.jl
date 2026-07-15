####################################################
# EXTRACT ALL LOGGER RESULTS TO A CSV FILE   
####################################################
full_output_csv_path = joinpath(case_dir, "full_lf_output.csv")

js = """
var fullRModule = ledaModules.FULLRESULTS();
var caseID = "$lf_case_id";  

//////// Exctracting results //////////

fullRModule.setUuid(caseID);
fullRModule.trendLoggersToSvFile("$full_output_csv_path");

"""

softsh = resolve_softsh_path()
script = joinpath(@__DIR__, "ledaflow_extract.js")
write(script, js)
logfile = joinpath(case_dir, "leda_extract.log")
run(pipeline(`$softsh $script`, stdout=logfile, stderr=logfile, append=true))
