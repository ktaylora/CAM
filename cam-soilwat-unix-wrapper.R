source("cheatgrass-abundance-model.R")
argv <- commandArgs(trailingOnly=T)
Rsw_loadClimateDb("dbWeatherData_GTD.sqlite3")
Rsw_CAM_run(sites=argv, Scenario="Current", years=c(1979,1985), debug=T, greppable=T)
