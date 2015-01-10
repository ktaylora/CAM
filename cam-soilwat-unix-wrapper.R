source("cheatgrass-abundance-model.R")
argv <- commandArgs(trailingOnly=T)

Rsw_setConnection(argv[2])
o <- dbW_getWeatherData(Site_id=1);
  years <- c(o[[1]]@year, o[[length(o)]]@year);
    rm(o);
    
Rsw_CAM_run(sites=argv[1], Scenario="Current", years=years, debug=T, greppable=T)
