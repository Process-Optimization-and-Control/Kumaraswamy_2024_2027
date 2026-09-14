####################################################
# EXTRACT ALL LOGGER RESULTS TO A CSV FILE
#
# Run once after a closed-loop simulation has finished, to dump every LedaFlow trend logger over
# the whole run to full_lf_output.csv. This is the file plot_lf_results.jl reads, so it has to be
# produced before the results of a case can be plotted. The per-sample trends.csv written during
# the run holds only the handful of loggers the controller needs and is overwritten each sample.
#
# This is a script, not a function: `include` it and it runs. It needs case_dir and lf_case_id, so
# a case file must be included first. LF_full_output_headers.txt lists the columns it produces.
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
