#' @export
sgpvAMrulesSingle <- function(data, designLooks, getUnrestricted){



  # 3 Affirm alert and report end of study (eos) operating characteristics (oc)
  #   - To be applied across various k (affirmation steps)
  affirmEndOfStudy <- function(i, designLooks){

    W <- designLooks[i,"W"]   # wait time
    S <- designLooks[i,"S"]   # monitoring steps
    A <- designLooks[i,"A"]   # required affirmation steps
    L <- designLooks[i,"L"]   # lag time to outcome
    N <- designLooks[i,"N"]   # maximum sample size


    # 1 Establish observations that surpass wait time and occur at lookSteps
    looks     <- W + (0:nrow(data)) * S
    monitor   <- data[,"n"] %in% looks


    # 2 Raise Alert
    alertNotROPE  <- which(data[,"sgpvROPE"] == 0 & monitor==TRUE)
    alertNotROME  <- which(data[,"sgpvROME"] == 0 & monitor==TRUE)

    alertNotROPEAny  <- length(alertNotROPE) > 0
    alertNotROMEAny  <- length(alertNotROME) > 0



    # Stop times for being Not ROPE and being Not ROME
    stopNotROPE  <- alertNotROPE[alertNotROPE %in% (alertNotROPE + A)]
    stopNotROME  <- alertNotROME[alertNotROME %in% (alertNotROME + A)]

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
              data[1,"effect1"],
              eos[keepStats],
              lag = eosLag[c(keepStats,"stopInconsistent","stopRejH0_YN","stopRejH0_NY")])

      return(c(oc,wait=W,steps=S,affirm=A,lag=L,maxN=N))
    } else {
      return(NULL)
    }
  }

  affirmedEnd <- t(sapply(1:nrow(designLooks),affirmEndOfStudy,designLooks=designLooks))


  return(affirmedEnd)

}
