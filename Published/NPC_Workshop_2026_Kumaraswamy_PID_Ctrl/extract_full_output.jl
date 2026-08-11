####################################################
# EXTRACT ALL LOGGER RESULTS TO A CSV FILE   
####################################################

function extract_full_output(lf_case_id::String, csv_path::String)
    ##############
    # EXTRACT RESULTS TO CSV FILE
    ##############
    full_output_csv_path = joinpath(@__DIR__, csv_path)

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
    logfile = joinpath(@__DIR__, "ledaflow_extract.log")
    run(pipeline(`$softsh $script`, stdout=logfile, stderr=logfile, append=true))
end 