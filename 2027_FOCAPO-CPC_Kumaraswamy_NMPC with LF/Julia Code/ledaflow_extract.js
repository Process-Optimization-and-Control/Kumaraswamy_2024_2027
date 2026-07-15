var fullRModule = ledaModules.FULLRESULTS();
var caseID = "fdd84235-e7a2-4651-b946-d45cd4deb638";  

//////// Exctracting results //////////

fullRModule.setUuid(caseID);
fullRModule.trendLoggersToSvFile("/home/archanak/projects/00 Code Templates/2027_FOCAPO_CPC/model_mismatch/case_final/full_lf_output.csv");

