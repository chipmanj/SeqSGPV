# sgpvAM.R
# J Chipman

library(plyr)


sgpvAM <- function(){
  
  
  # 1 collect array of simulated data
  mcmcMonitoring <- replicate(n = 300, expr = {
      sgpvAMdata(rnorm, dataGenArgs = list(n=700), effectGeneration = 0.5,
                 deltaL2 = -0.5, deltaL1 = -0.15, deltaG1 = 0.15, deltaG2 = 0.5,
                 monitoringIntervalLevel = 0.05)
  })
  
  
  # 2 Make sure all simulations will continue until completion
  #   Look for stability of sgpv for last set of maxAlert patients
  
  getMore <- max(aaply(mcmcMonitoringFixed, .margins = 3, .fun = mcmcMonitoringEnoughCheck))
  while(getMore > 0){
    print("Generating more data to ensure simulations go to completion.")
    print(paste("Current observations = ", nrow(mcmcMonitoring[,,1]), " adding another ", getMore))
  
    
    aaply(mcmcMonitoring, .margins = 3, .fun = sgpvAMdata,
          rnorm, dataGenArgs = list(n=getMore), effectGeneration = 0.5,
          deltaL2 = -0.5, deltaL1 = -0.15, deltaG1 = 0.15, deltaG2 = 0.5,
          monitoringIntervalLevel = 0.05, exisitingY =)
      sgpvAMdata(rnorm, dataGenArgs = list(n=getMore), effectGeneration = 0.5,
                 deltaL2 = -0.5, deltaL1 = -0.15, deltaG1 = 0.15, deltaG2 = 0.5,
                 monitoringIntervalLevel = 0.05, exisitingY = )
    })
    
      
    
  }
  
  
  
  mcmcTest <- replicate(n = 3, expr = {
    sgpvAMdata(rnorm, dataGenArgs = list(n=750), effectGeneration = 0.5,
               deltaL2 = -0.5, deltaL1 = -0.15, deltaG1 = 0.15, deltaG2 = 0.5,
               monitoringIntervalLevel = 0.05)
  })
  
  
  waitWidths <- c(0.15, 0.20, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60)
  for (w in 1:length(waitWidths)){

    ww <- waithWidths[w]
    
    # 2 adaptively monitor simulated data across multiple burn ins
    mcmcEndOfStudy <- aaply(mcmcMonitoringFixed, .margins = 3, .fun = function(x) {
                              t(sgpvAMrules(data = x, waitWidth = ww,
                                            lookSteps = 1, maxAlertSteps = 100,
                                            monitoringIntervalLevel = 0.05))
    })
    
    
    
    # 3 Check that all simulations ran to completion. If not supplement data
    
    
    
    
    
    # 3 aggregate simulated data
    mcmcEndOfStudyAve <- t(aaply(mcmcEndOfStudy, .margins = c(2,3), .fun = mean))
    
  }
  

}







# 2 adaptively monitor simulated data across multiple burn ins
a <- aaply(mcmcMonitoringFixed, .margins = 3, .fun = function(x) {
  t(sgpvAMrules(data = x, waitWidth = 0.3,
                lookSteps = 1, maxAlertSteps = 100,
                monitoringIntervalLevel = 0.05))
})


# 3 aggregate simulated data
aa <- t(aaply(a, .margins = c(2,3), .fun = mean))

