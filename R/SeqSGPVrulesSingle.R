#' @title SeqSGPVrules
#'
#' @description Applies monitoring frequencies to PRISM monitoring and returns the final observed observation. SeqSGPVrules is called within SeqSGPV.
#'
#' @param data Single mcmc replicate
#' @param monitoringFrequency Matrix of monitoring frequency combinations
#' @param randomize TRUE if length(allocations) > 1
#'
#' @export
SeqSGPVrulesSingle <- function(data, monitoringFrequency, randomize){



  # 3 Affirm alert and report end of study (eos) operating characteristics (oc)
  #   - To be applied across various k (affirmation steps)
  affirmEndOfStudy <- function(i, monitoringFrequency, effectX){

    W <- monitoringFrequency[i,"W"]   # wait time
    S <- monitoringFrequency[i,"S"]   # monitoring steps
    A <- monitoringFrequency[i,"A"]   # required affirmation steps
    L <- monitoringFrequency[i,"L"]   # lag time to outcome
    N <- monitoringFrequency[i,"N"]   # maximum sample size


    # 1 Establish observations that surpass wait time and occur at lookSteps
    looks     <- W + (0:nrow(data)) * S
    monitor   <- data[,"n"][data[,"n"] %in% looks]


    # 2 Raise Alert
    alertNotROPE  <- which(data[,"sgpvROPE"] == 0)
    alertNotROME  <- which(data[,"sgpvROME"] == 0)

    # Stop times for being Not ROPE and being Not ROME
    stopNotROPE  <- alertNotROPE[(alertNotROPE - A) %in% alertNotROPE  & (alertNotROPE %in% monitor) ]
    stopNotROME  <- alertNotROME[(alertNotROME - A) %in% alertNotROME  & (alertNotROME %in% monitor) ]


    # Stop times with and without Lag
    if(N==Inf){
      stop    <- min(stopNotROPE, stopNotROME, na.rm = TRUE)
      stopLag <- stop + L
    } else {
      stop    <- min(stopNotROPE, stopNotROME, N, na.rm = TRUE)
      stopLag <- min(stop+L,N)
    }

    # End of study (immediate outcomes)
    eos                      <- data[data[,"n"]==stop,]
    eos["stopNotROPE"]       <- as.numeric(eos["sgpvROPE"]==0)
    eos["stopNotROME"]       <- as.numeric(eos["sgpvROME"]==0)
    eos["stopInconclusive"]  <- as.numeric(eos["stopNotROPE"]==0 & eos["stopNotROME"]==0)


    # Lag outcomes
    eosLag                     <- data[data[,"n"]==stopLag,]
    eosLag["stopNotROPE"]      <- as.numeric(eosLag["sgpvROPE"]==0)
    eosLag["stopNotROME"]      <- as.numeric(eosLag["sgpvROME"]==0)
    eosLag["stopInconclusive"] <- as.numeric(eosLag["stopNotROPE"]==0 &
                                             eosLag["stopNotROME"]==0)
    eosLag["stopInconsistent"] <-
      as.numeric( ( eos["stopNotROPE"] == 1 & eosLag["stopNotROPE"] != 1 )  |
                  ( eos["stopNotROME"] == 1 & eosLag["stopNotROME"] != 1 )  )

    eosLag["stopRejH0_YN"] <- as.numeric( ( eos["rejH0"] == 1 & eosLag["rejH0"] != 1 ) )
    eosLag["stopRejH0_NY"] <- as.numeric( ( eos["rejH0"] != 1 & eosLag["rejH0"] == 1 ) )



    # Keep theta and operating characteristics (drop y, trt, and sgpvs)
    keepStats <- c("n","est","bias","rejH0","cover", "stopNotROPE","stopNotROME","stopInconclusive")


    if(!is.na(stop)){
      oc <- c(data[1,"theta0"],
              data[1,effectX],
              eos[keepStats],
              lag = eosLag[c(keepStats,"stopInconsistent","stopRejH0_YN","stopRejH0_NY")])

      return(c(oc,wait=W,steps=S,affirm=A,lag=L,N=N))
    } else {
      return(NULL)
    }
  }

  if(randomize==FALSE){
    effectX <- "effect0"
  } else {
    effectX <- "effect1"
  }

  affirmedEnd <- t(sapply(1:nrow(monitoringFrequency),affirmEndOfStudy,monitoringFrequency=monitoringFrequency, effectX=effectX))


  return(affirmedEnd)

}
