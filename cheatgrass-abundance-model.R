#
# Cheatgrass Abundance Model
#
# A mechanistic framework for describing cheatgrass growth (biomass and seedbank dynamics) over time.  This
# demographic framework is a frontend for SOILWAT that uses daily time-step soil hydrology information and
# a mass of information published in the scientific literature to simulate growth for bromus tectorum.
#
# Original Author: Kyle Taylor (kyle.a.taylor@gmail.com) [2014]
#
#

require(Rsoilwat)
require(snow)

firstLine <<- T        # used in producing greppable output
seedsContributed <<- F # indicates if our population has contributed seeds of its own volition. Used to address the "population/seedbank crash" bug.

#
# seedProduction()
# Our model opperates on a 1 m^2 scale.  Assuming this scale, 10,750 tillers/m^2 is a "high" density plot (Piemeisel, 1938).
# This assumes a plot is bare ground.  Theoretically, this density will be lower if there is other vegetation on the plot.  I
# am applying a linear correction factor to the max plot density that reflects % bare ground at a site.  Higher bare
# ground assumes the lower overall plot density, and thus a greater the amount of cheatgrass seed produced.

CAM_seedProduction <- function(session, bareGround=0.5){

  # Curve fit using a polynomial interpolation and observations from (Piemeisel, 1938; Young et al.,1987)
  # Assuming seeds produced f(population density), when limiting factors are removed.  At high densities, at least 10 viable seeds are produced per individual,  
  # while at low densities as many as 5000 seeds can be produced, such that:
  meanSeedsProduced<-function(x){ x<-10*(x^(-1.6385558412)); x[x>5000]<-5000; return(x) }

  # estimate corrected plot density
  plotDensity <- nrow(session)/(10750*bareGround) # (Piemeisel, 38; Young, 87)
  readyToSeed <- (session$lifestage == "established" & session$age > 20)
    plantAges <- session$age[readyToSeed]
    plantAGB  <- session$agBiomass[readyToSeed]
    
  if(sum(readyToSeed) > 0){
	# change the phenological signature of the plants that contributed seed
    session$lifestage[readyToSeed] <- "senescent"
    # simulate seed production using a truncated normal distribution
	nSeed <- round(meanSeedsProduced(plotDensity))
	  onSeed <- nSeed <- round(rep(nSeed,sum(readyToSeed))) 
		nSeed[nSeed < 1] <- 0 
		  nSeed <- round(nSeed * (plantAGB/(plantAges*0.00055749))) # treat nSeed as the maximum, and correct according to biomass accumulation
          #r<-round(runif(n=length(nSeed)))
          #a<-(sum(nSeed)/100)/sum(r/max(r)) 
          #nSeed <- round((r/max(r))*a)
		cat(" -- correction factor: (", median(onSeed), ")*(",(median(plantAGB)/(median(plantAges)*0.00055749)),")=",median(nSeed),"\n",sep="")
		cat(" -- seed production event: nSeed=",sum(nSeed),", nIndiv=", sum(readyToSeed), "\n", sep="")
		  seedsContributed <<- T
  } else { nSeed <- 0 } # if nothing is ready to seed, set to 0 and return whole population table back to user

  #return seeds and the session data back to the user
  return(list(session,sum(nSeed)))
}

#
# CAM_germination()
# Given a pool of some number of seeds (n), define the
# number of individuals that will successfully germinate.
#

# bareground is a vector of length 2, indicating number of current individuals in plot and the assumption of bareground (via LAI or prior knowledge)
# associated with the area

CAM_germination <- function(seeds=NULL, swp=0, sTemp=5, snowcover=0, bareGround=c(1,1), rainDays=0) {

  # sanity check input data
  if(is.null(seeds) | (nrow(seeds) < 1)){ cat(" -- error: germination request on empty seedbank.\n"); return(0) }

  #
  # Assuming that germination under suitable conditions is a function of seedbank age (Klemmedson and Smith (1964))
  # in which the optimal rate is achieved when seeds are ~ 4 months old (Klemmedson, 64).  This is treated as a back-of-the-envelope
  # calculation for after-ripening that affects CG over-winter.  As it stands, this could use a little more detailing.  Such as providing
  # a mechanism to relate snowcover and soil temperature to the after-ripening process.
  #

  logRate <- function(t0=1, t=179) {
   # n<-n*(log(t0/t)-log(1/t))*0.1966888
   # n<-n*(exp(t0/t)*0.3740752)
   # curve-fit polynomial
    # max <- = -2E-18*t^6 + 1E-14*t^5 - 4E-11*t^4 + 6E-08*t^3 - 4E-05*t^2 + 0.0107*t - 0.0877
    max <-(-0.0028*t^2+t-30)
      # this is a curve fit to a parabola maximum at 120 days (Klemmedson)
      # o <-(-0.0042*t0^2+t0+1)

      # this is a curve fit to a parabola that assumes no growth potential when the seed is < 1 month old, and
      # relaxes the assumption that the peak occurs at 120 days, and instead occurs ~170 days.  Allows the seed
      # to mature upto ~ 1 full year, which has been reported in the literature as maximum seed age for germination
      # (in the wild).

      o<-(-0.0028*t0^2+t0-30)
        o[o<0]<-0
          o<- o/(max/0.96) # scale our polynomial so that vertex of the porabola is y=0.96, assumes a maximum of ~96% germination rate (Mack and Pyke)

    return(o)
  }

  # if SWP is greater than -1.5 MPa and sTemp is greater than 0 celcius (Harris, 67; Pyke and Novak; Bradford et al., 2006)
  # Goodwin et al., (1996) note germination rates higher when SWP approaches 0 (1-2 day rate @ 0 MPa, 2-5 @ -1 MPa).
  # - if sTemp is greater than or equal to 29(c), germination fails (McCarlie, 2001)
  # - germinant growth / establishment stops at soil temperatures > 15(c) (Young,2000)
  # - germination itself is inhibited at > 30(c) (Harris, 1976)
  # - germination can occur at soil temperatures just above 0 (Evans, 1972)

  if(swp > -1.5 && swp < 0 &&
     sTemp >= 0 && sTemp <= 30 &&
     rainDays > 2 &&
     snowcover < 150)
  {
	# test: maximally germinate seeds older than 1 year (D. Schlaepfer)
	seeds$age[seeds$age > 178] <- 178
    # calculate the probability of germination for each seed, based on its age.
    r<-logRate(t0=seeds$age)
	# test: calculate the plot density coefficient to modify number of germinants
    cat("  -- conditions suitable for germination. \n")
	plotDensityBeta <- 1-(bareGround[1]/(10750*bareGround[2])) # High density in monocultures (Young et al.,78)
      if(plotDensityBeta < 0) { plotDensityBeta <- 0 }
	# assume a flat rate of germination for this event.  Find the unique values of germination % in the age groups within the seedbank,
	# take the median probability of germination amoung those groups, and germinate that median %, preferrentially taking the highest probability
	# seeds out of the stack first.
	nToGerminate<-round(median(unique(r))*nrow(seeds))
	  nToGerminate <- round(nToGerminate*plotDensityBeta) # test: bare ground correction
	     ## D. Schlaepfer hated this ... Re-implementing a self-thinning algorithm under mortality() to account
       ## for outrageous germination events.
       ##
       ## bug fix: don't germinate a ridiculous number of seeds just because the seedbank is large.  
	     ## max densities observed by Young et al. (78) were ~10000/m2.  Steward and Hall (46) report 15000/m2 in large monocultures.  
	     ## let's treat 10,000/m2 as a soft boundary and take a guess around 10000 by sampling the normal distribution
	     ## correctedMax <- round(10000*bareGround[2])
	     ## if(nToGerminate > correctedMax) { 
       ##   cat("  -- outlandish germination fix. Setting n ~=",correctedMax,"\n",sep=""); 
       ##   nToGerminate <- round(rnorm(n=1,mean=correctedMax,sd=1000)); 
       ## }
    cat("  -- germinants = (bareGroundBeta*nSuitableSeeds) = (", plotDensityBeta,")*(",round(median(unique(r))*nrow(seeds)),") = ", nToGerminate, "\n",sep="")
    if(nToGerminate > 0){
	  probs <- seq(from=1,to=0.1,by=-0.1)
	  sample <- NA
	  for(p in probs){
		  sample <- try(sample(x=which(r>p),size=nToGerminate,replace=F), silent=T)
		  if(class(sample) != "try-error") break
  	  }
  	  if(length(sample)<nToGerminate){ 
  	    cat("  -- not enough viable seeds in SB to satisfy n=",nToGerminate,"\n") 
  	    return(list(0, seeds))
  	  } else {
	    keep <- which(!(1:length(seeds$age)) %in% sample)
	      seeds <- seeds[keep,] # take out the germinated seeds from the bank
		    if(class(seeds) == "numeric") { seeds <- data.frame(age=seeds) }
	    return(list(length(sample), seeds)) # return the number of germinants and seedbank to user
	  }
	}
  }

  return(list(0,seeds)) # no germination?  return 0 and the seedbank in-tact
}

#
# CAM_phenology()
# initiate lifestage transitions
#

CAM_phenology <- function(n, establishmentSignal){

  # Bradford and Harris note atleast 15 days of suitable (non-drought) conditions needed
  # to "establish".  This isn't an physiologically meaningful number.  Find a way to relate
  # age, biomass/rootlength accumulation, and environmental signal to lifestage transition

  if(establishmentSignal >= 15){
	e <- which(n$age >= 20 & n$rootLength > 5 & n$lifestage == "seedling")
    if(length(e)>0) {
		n$lifestage[e] <- "established"
		#cat(" -- establishment event, n=",length(e),"\n",sep="")
        #cat(" -- establishment positions:",e,"\n")
	}
  }

  if(establishmentSignal > 1) { # relax our assumption about estab. signal for older plants
	e <- which(n$age >= 45 & n$rootLength > 3 & n$lifestage == "seedling")
	if(length(e)>0){
      n$lifestage[e] <- "established"
      #cat(" -- establishment event, n=",length(e),"\n",sep="")
      #cat(" -- establishment positions:",e,"\n")
    }
  }

  return(n)
}

#
# CAM_mortality()
#

CAM_mortality <- function(n, sTemp=0, droughtSignal=0){

  # sheley function treating carrying capacity as a function of median biomass.  
  # deriving this took some effort and is based on Sheley, 2004.  talk to k. taylor for notes
  sheley_k_upperLimit <- function(x) { 
    x<-x*1000; # the regression accepts biomass in grams
    exp((-0.1715*x)+9.4089); 
  }

  # thinning coefficient that relates our current population size to carrying capacity (defined above)
  # via an exponential function [x=(n/k)]. As the population approaches 20 times carrying capacity, 
  # out thinning factor (beta) goes to 1
  beta_t <- function(x) { (x^2)/(20^2) }

  # assume total population mortality if soil temperatures are <= -23 deg C (Lloyd, 1955)
  if(sTemp <= -23){
	cat("-- mass population freezing mortality event.\n")
    #n$lifestage <- "dead"; 
    n<-data.frame()
  }
  #
  # test: make the assumption that more than 20 days of sustained drought will kill successively larger portions of the established population 
  #
  else if(droughtSignal > 20){
    survivors <- (droughtSignal-20)/10 # normalized to a 0-1 scale
      survivors <- 1-(survivors)
        survivors <- sample(1:nrow(n), size=round(nrow(n)*survivors))
    cat(" -- drought mortality event:",nrow(n)-length(survivors)," established plants lost.\n")
    n<-n[survivors,]
  }
  # cull seedlings that have been exposed to 10 or more days of drought (Frasier, 1994)
  # adults are drought hardy, but must have accumulated some root depth to capitalize on water resources
  else if(droughtSignal > 10){
    # cull seedlings
    cull <- n$lifestage == "seedling"
    if(sum(cull) > 0){
      cat(" -- drought mortality event:",sum(cull)," seedlings lost.\n")
      #n$lifestage[which(cull)] <- "dead"
      n<-n[which(!cull),]
    }
  #
  # Individuals that have gone to see typically die after sudden drought exposure
  #
  } else if(droughtSignal > 5) { # derived from (Billings, 1952; Hall, 1949) 
    cull <- n$lifestage == "senescent"
    if(sum(cull) > 0){
      cat(" -- drought mortality event:",sum(cull)," senescents lost.\n")
      #n$lifestage[which(cull)] <- "dead"
      n<-n[which(!cull),]
    }    
  } 
  
  #
  # Self-thinning algorithm implemented to account for outrageous population growth.
  #
  k <- sheley_k_upperLimit(median(n$agBiomass)) # actual k
    k <- beta_t(nrow(n)/k) # solve for an appropriate kill coefficient
      k <- round(k*nrow(n)) # number to kill in population

  n <- n[order(n$agBiomass, decreasing=F),] # sort our table from smallest-to-largest
    n <- n[1:k,] # kill the small ones first

  #
  # Assume mortality for any individual older than 280 days (Hulbert, 55; Spence, 37; Harris, 67)
  #
  cull <- (n$age > 279)
  if(sum(cull)>0){
	  #n$lifestage[cull] <- "dead"
	  n<-n[which(!cull),]
    cat(" -- mortality event: EOL reached for", sum(cull), "individuals.\n")
  }
  #
  # test: remove n plants individuals due to disease, herbivory, etc
  #
  if(nrow(n)>1){
    nToExpire <- sample(1:nrow(n),size=round(nrow(n)*(0.1/365)))
      keep <- which(!(1:nrow(n) %in% nToExpire))
        n<-n[keep,]
  }
  
  return(n)
}

#
# CAM_rootGrowth()
#

CAM_rootGrowth <- function(n,sTemp,swp) {

  # f(swp) = beta parameter to apply to our empirical root growth rate maximum that accounts for SWP.
  # This parabolic function assumes roots grow more when water stress is highest
  # where beta max occurs @ 1.3 -SWP.  When SWP increases beyond this, I assume the plant is going to
  # close its stomata to prevent drought mortality, limiting net photosynthesis/production and slowing
  # the growth rate of its roots.  This function assumes rGr stops completely @ SWP <= -2.7
  swpBeta <- function(x) { x<-(-x); return(-(0.3049*(x^2)-(0.7926*x)-0.05152)) }

  # max CG root growth depth is ~130 cm for a single growing season when nitrogen
  # is not limiting, and ~ 90 cm under "normal" nitrogen conditions (Hulbert, 55; Spence, 37; Harris, 67)
  rGr <- 130/280 # daily root growth rate (cm) [110 cm / mean lifespan for species]

  # assuming root growth can happen at temps greater than or equal to 3 deg celcius (Harris, 67)
  if(sTemp > -3){
    notDead_bool <- (n$lifestage != "scenescent")
	growth <- swpBeta(swp)
      growth <- (growth/max(growth))*rGr
	    if(growth < 0){ growth <- 0 }
    n$rootLength[notDead_bool] <- n$rootLength[notDead_bool] + growth
  }

  return(n)
}

#
# Aboveground Growth
#

CAM_shootGrowth <- function(n,sTemp,snowcover=0){

  # above ground growth halted under snow cover
  if(snowcover < 150){
    # fit a parabolic function with the assumptions :
    # assume growth rate optimum @ ~16 degrees celcius, with a minimum at ~ 5 degrees celcius and inhibition around 29 degrees celcius
    # (McCarlie, 2001; Harris, 1976)
    gr<-function(x){return(-0.0167*(x^2)+(0.533*x)-2.0012633)}

    g<-gr(sTemp)
      if(g<0) g<-0; # our function is undefined when y<0; assume no growth beyond our defined range
        g<-g/2.251537 # add a link step to scale our output [0/1] relative to the curve maximum (x=16 degrees; y=2.251537)

    # apply our rate x max mass of observed cheatgrass growth / day in the field (3.22*10^-5 g/plant) [From: Harris,67; Klemmedson,64]
    # g<-g*0.00000322222
      g<-ifelse(g<1,g,1)*0.00055749 # back-of-the-envelope [From Hulbert, 1955]
      #g<-ifelse(g<1,g,1)*2.9 # (Klemmedson, 64)
      n$agBiomass[n$lifestage != "scenescent"] <- n$agBiomass[n$lifestage != "scenescent"] + g
  }
  return(n)
}

#
# cam :: MAIN
#

CAM_run <- function(n=1, session=NULL, maxSeedbankLife=(365*3), debug=F, greppable=F, hobble=0){

  stopifnot(!is.null(session))

  #
  # session data
  #

              seedbank <- data.frame(age=rep(1,n))
			population <- data.frame()
             cohortAge <- 0
         droughtSignal <- 0                   # number of consecutive days of ecological drought conditions
   establishmentSignal <- 0                   # number of consecutive days with SWP between -0.3 and -0.7 (15 Req, Harris)
        coldSnapSignal <- 0
  winterDormancySignal <- 0                   # number of consecutive days of winter dormancy
            rainSignal <- 0                   # number of consecutive days of rain
	                 doy <- 1:nrow(session)


  # functions
  calcDoy <- function(x) { x<-x+180; return(round((((x/365)-floor(x/365))*365)+1)) } # note: climate observations clipped to start at mid summer (doy=180) of first year
  calcYear <- function(x) { x<-x+180; return(floor(x/365)+1) }
  
  #
  # hypothetical structure of an "individual" in our population, after germination :
  # age (days), agBiomass=1, rootLength=0.1 (mm), lifestage=c("seedling","emergent","mature","senescent")
  #

  for(i in 1:nrow(session)){ 
  
   ##
   ## Check environmental signals (1 mm daily precipitation requirement [Pyke, 94])
   ##
  
   wet_bool <- (session[i,]$swp > -1.5 & session[i,]$swp < 0);      # textbook SWP levels for plant Available Water Content: -1.5 <-> -0.03
                                                                    # wet signaling will adjust our definition of drought and be used to satisfy establishment    
  
   # check for drought signal and Bradford/Harris Establishment Signal (Bradford & Lauenroth, 2006)
   # note, original HB signal: swp >= -0.7 & swp <= -0.3
     
   if(wet_bool) 
     { rainSignal <- rainSignal + 1; establishmentSignal <- establishmentSignal + 1; droughtSignal <- 0; }
   else
     { rainSignal <- establishmentSignal <- 0; droughtSignal <- droughtSignal+1; }
   
   ##
   ## If we have a standing crop, let's subject them to growth, phenology, and mortality
   ##
   
   if(nrow(population)>0) {

     if(sum(notDead_bool)>0){ 
       #
       # check for mortality conditions for existing individuals in the population pool.  Kill if appropriate
       #
       population <- CAM_mortality(population, droughtSignal=droughtSignal)

       #
       # check for dormancy release / initialize a growth phase for existing individuals in the population pool
       #
       population <- CAM_shootGrowth(population, sTemp=session[i,]$sTemp, snowcover=session[i,]$snowcover)
       population <- CAM_rootGrowth(population, sTemp=session[i,]$sTemp, swp=session[i,]$swp)

       #
       # check to see if today's environmental conditions trigger a lifestage transition
       #
       population <-CAM_phenology(population, establishmentSignal)
     }

   } 

   notDead_bool <- (population$lifestage != "scenescent") # update our vector of the undead
   
   ##
   ## simulate potential seedling germination for the day, accounting for seedbank size
   ##

   if(nrow(seedbank) > 0){
     g<-nrow(population[notDead_bool,])
       g <-CAM_germination(seeds=seedbank, swp=session[i,]$swp,
                      sTemp=session[i,]$sTemp, snowcover=session[i,]$snowcover,
                      bareGround=c(ifelse(g>0, g, 1), 0.5), rainDays=rainSignal)

     seedbank <- g[[2]] # update the seedbank
     rainSignal <- (g[[1]] == 0) * rainSignal # test: reset our "rain duration" signal to address rapid consecutive germination and seedbank exhaustion bug
   }

   #
   # add today's germinants to the population pool
   #

   if(g[[1]] > 0){
    population <- rbind(population, data.frame(age=rep(1,g[[1]]), germinationDOY=doy[rep(i,g[[1]])], agBiomass=rep(0.00001,g[[1]]),
                                                rootLength=rep(0.00001,g[[1]]), lifestage=rep("seedling",g[[1]])))
	  g[[1]] <- 0
   }

   notDead_bool <- (population$lifestage != "scenescent") # update our vector of the undead
   
   #
   # simulate potential seed production for the day for established individuals that are ripened
   #

   if(nrow(population)>0){
     s<-CAM_seedProduction(population,bareGround=0.5)
       population <- s[[1]]
         seedbank <- rbind(seedbank, data.frame(age=rep(1,s[[2]])))
	       rm(s)
   }

   ##
   ## adjust population age (days)
   ##

   # for individuals
   if(nrow(population)>0){
     population$lifestage <- as.vector(population$lifestage);
     population$age[notDead_bool] <- population$age[notDead_bool] + 1
   }

   # for the seedbank
   if(nrow(seedbank)>0){
	   # kill off those seeds in the seedbank who are set to expire
	   seedbank <- data.frame(age=seedbank$age[seedbank$age < maxSeedbankLife])
     # test: assume a daily percentage of the seedbank is lost to seed predation, mold/fungi, and other mortality factors.  
     # this should add-up to ~10% of yearly seedbank size over the course of 365 days, if we assume the seedbank size is static.
     #
     # "Pathogen infection imposes a race for survival between the seed and the pathogen: seeds that germinate slowly and become infected are killed, 
     # but those that germinate quickly can survive despite infection." (Beckstead et al. 2007) [Mordecai, 2013]
     agePrefMortalityRate<-function(age=1){ o<--0.00137*age^2+age-2; }
     max <- agePrefMortalityRate(365)
     a <- seedbank$age
       a[a>365] <- 365
         a <- agePrefMortalityRate(seedbank$age)/max
     nSeedToExpire <- round(length(a)*(0.1/365))
     # preferrentially kill older seeds first
     probs <- seq(from=1,to=0.1,by=-0.1)
     sample <- NA
     for(p in probs){
       sample <- try(sample(x=which(a>p),size=nSeedToExpire,replace=F), silent=T)
       if(class(sample) != "try-error") break
     }
     keep <- which(!which(seedbank$age == seedbank$age) %in% sample)
       seedbank <- data.frame(age=seedbank$age[keep]) # take out the affected seeds from the bank
     # add a day to the age of remaining seeds
     seedbank$age <- seedbank$age+1
   }
   
   ## debug
   if(debug) { 
     if(!greppable){
       cat("year(",calcYear(i),") doy(",calcDoy(i),") rain.sig(",rainSignal,") drt.sig(",droughtSignal,") p.size(", 
       ifelse(nrow(population[notDead_bool,]) > 0, nrow(population[notDead_bool,]), 0), ") ag.bio(",ifelse(nrow(population)>0, 
       mean(population$agBiomass),0),") p.MnAge(",ifelse(nrow(population[notDead_bool,]) > 0, median(population[notDead_bool,]$age, na.rm=T), 0),
       ") snow(",session[i,]$snowcover,") precip(", session[i,]$precip,") swp(",session[i,]$swp, ") sTemp(",session[i,]$sTemp,
       ") p.germnts(",g[[1]],")", " estb.sig.(", establishmentSignal, ")", " sb.size(",nrow(seedbank),") p.MnRtLen(", 
       ifelse(nrow(population[notDead_bool,]) > 0, mean(population[notDead_bool,]$rootLength, na.rm=T), 0), ")", "\n",sep="")
     } else {
       if(firstLine){ cat("year,doy,rain.sig,drt.sig,p.size,ag.mdBiomass,p.meanAge,snow,precip,swp,s.temp,p.germinants,estab.signal,sb.size,p.mnRtLen","\n",sep="");firstLine<<-F; }
       cat(calcYear(i),",",calcDoy(i),",",rainSignal,",",droughtSignal,",", 
       ifelse(nrow(population[notDead_bool,]) > 0, nrow(population[notDead_bool,]), 0), ",",ifelse(nrow(population[notDead_bool,])>0, 
       median(population$agBiomass[notDead_bool]),0),",",ifelse(nrow(population[notDead_bool,]) > 0, median(population[notDead_bool,]$age, na.rm=T), 0),
       ",",session[i,]$snowcover,",", session[i,]$precip,",",session[i,]$swp, ",",session[i,]$sTemp,
       ",",g[[1]],",", establishmentSignal, ",",nrow(seedbank),",", 
       ifelse(nrow(population[notDead_bool,]) > 0, mean(population[notDead_bool,]$rootLength, na.rm=T), 0), "\n",sep="")       
     }
     
     Sys.sleep(hobble);
   }

    # test: if there are no individuals and no seeds in seed bank, reintroduce a marginal number of seeds
    if(nrow(population[notDead_bool,]) < 1 & nrow(seedbank) < 1){
      if(!seedsContributed) { cat(" -- population crash: no seeds contributed to seedbank.\n"); break; } # if our population has never contributed seed, do not make the assumption of persistence in the seedbank
      cat(" -- population crash: seedbank and population are at zero. Reintroducing a n=",n," seeds...\n", sep="")
      seedbank <- data.frame(age=rep(1,n))
	  population <- data.frame()
    }
  }

  return(population)
}

#
# Rsw_loadClimateDb()
#

Rsw_loadClimateDb <- function(filePath){ if(file.exists(filePath)) Rsoilwat::dbW_setConnection(filePath) }

#
# Rsoilwat/cheatgrass abundance module interface
# Wrapper function for initiating Rsoilwat and then parsing output data for the CAM interface.
#

Rsw_CAM_run <- function(extent=NULL, sites=NULL, years=NULL, Scenario="Current", initialCG_N=500, maxSeedbankLife=(365*3), debug=F, greppable=T, hobble=F){

  # import default libraries
  require(rgdal)
  require(raster)

  # local variables
  if(is.null(sites)){
    sites <- vector() # if we weren't passed any sites explicitly, assume we will run all sites implemented in the climate db
  }

  #
  # extractSites()
  # optionally extract sites within the extent of a user-specified buffer-region, if asked
  #

  extractSites <- function(){
    # attempt to read-in sites from dbW
    sites<-try(dbW_getSiteTable())
    # convert to a spatial point data and extract accordingly
      if(class(sites) == "try-error") stop("error in dbW_getSiteTable().  Have you initialized the Rsoilwat weather database for this session?\n");
    sites <- SpatialPointsDataFrame(coords=data.frame(x=sites$Longitude,y=sites$Latitude),
                                    data=data.frame(Site_ID=sites$Site_id),
                                    proj4string=crs("+init=epsg:4326")) # make an assumption that climate data are always provided as lat/long
    sites <- raster::crop(sites, extent)
    return(as.vector(sites$Site_ID))
  }

  #
  # postProcessRswOutput()
  # post-process the output returned from Rsoilwat::sw_exec() into something that makes sense to CAM
  #
  
  postProcessRswOutput <- function(t){
    # extract information from Rsoilwat's output relavent to the CAM
    out_swp <- t$swp
      out_swp <- out_swp$dy[,3] # always assume the 3rd col. corresponds to the 1st sim. layer; this is the only soil layer we care about for invasive bromes.
        out_swp <- out_swp[180:length(out_swp)]  # crop the first 179 off of our sample, so we begin our sample in the Fall
    out_sTemp <- t$soil_temp
      out_sTemp <- out_sTemp$dy[,3]
        out_sTemp <- out_sTemp[180:length(out_sTemp)]
    out_precip <- t$precip
      out_precip <- out_precip$dy[,3]
        out_precip <- out_precip[180:length(out_precip)]
    out_snowcover <- t$snowpack
      out_snowcover <- out_snowcover$dy[,3]
        out_snowcover <- out_snowcover[180:length(out_snowcover)]
        
    # convert to units that are meaningful to CAM
    out<- data.frame(swp=-(out_swp/10),          # Bar -> -MPa
                     sTemp=out_sTemp,
                     precip=out_precip*10,       # cm -> mm
                     snowcover=out_snowcover*10) # cm -> mm
    return(out)
  }
  
  #
  # Rsw_call :: MAIN
  #

  # parse out sites within the extent of our research area, if askedout
  if(!is.null(extent)){
    if(class(extent) != "Extent" & length(sites) == 0){
      extent <- try(extent(extent))
        if(class(extent) == "try-error") stop("fatal: unable to determine extent of object passed to Rsw_call()\n");
    }
    # extract sites within extent
    sites <- extractSites()
  } else {
    # if the user didn't ask to crop our site selection and didn't specify sites explicitly at run-tim, attempt to read-in
    # sites directly from an active dbW connection
    if(length(sites) == 0){
      sites<-try(dbW_getSiteTable())
        sites <- as.vector(sites$Site_ID)
    }
  }

  # process our soilwat runs iteratively for each requested site
  for(s in sites){

    # fetch current site weather data and define run settings for Rsoilwat
    focal_rData <- sw_inputData() # dummy template
      swWeather_FirstYearHistorical(focal_rData) <- years[1]
    focal_wData <- dbW_getWeatherData(Site_id=s,Scenario=Scenario,startYear=years[1],endYear=years[2])

    # assign values to current run settings
    swYears_StartYear(focal_rData) <- years[1]
    swYears_EndYear(focal_rData) <- years[2]

	  # execute the soilwat run
    focal_outData <- Rsoilwat::sw_exec(data=focal_rData, weatherList=focal_wData, colNames=T)

	  # post-process our data for the CAM
    focal_outData <- postProcessRswOutput(focal_outData)

    # execute the CAM
    focal_outData <- CAM_run(n=initialCG_N, session=focal_outData, debug=debug, greppable=greppable, hobble=hobble)
    return() # debug: only running one site at a time right now
  }
}
