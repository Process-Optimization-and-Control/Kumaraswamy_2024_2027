var calcModule = ledaModules.CALCULATE();
var fullRModule = ledaModules.FULLRESULTS();

var caseID = "a8dd37ee-50ac-4210-991f-87d393811fd8";  

//////// Reset case to initial conditions ////////

var cc = new Compound(caseID);
var Valve = cc.relation("child", "Leda1DnPhValve", true)
var Pump = cc.relation("child", "Leda1DnPhPump", true)

Valve[0].setOpenFrac(1.0)
Valve[1].setOpenFrac(1.0)
Valve[2].setOpenFrac(1.0)
Valve[3].setOpenFrac(1.0)

Pump[0].setSpeed([3440.0201802713636 * 0.10472])  // convert from rpm to rad/s

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
fullRModule.trendLoggersToSvFile("/home/archanak/projects/00 Code Templates/2027_FOCAPO_CPC/model_mismatch/case_no_state_est/steady_state_trends.csv");

