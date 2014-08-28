source("Workspace/cheatgrass-abundance-model.R")
argv <- commandArgs(trailingOnly=T)
Rsw_loadClimateDb("/Users/ktaylora/dbWeatherData_GTD.sqlite3")
Rsw_CAM_run(sites=argv, Scenario="Current", years=c(1979,1983), debug=T, greppable=T)
