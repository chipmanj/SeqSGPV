#' @export
sgpvAMrulesSingle <- function(data,
                              waitJ,
                              lookS,
                              affirmK,
                              getUnrestricted,  maxN, lagN=0){


  # 1 Establish observations that surpass wait time and occur at lookSteps
  looks     <- waitJ + (0:nrow(data)) * lookS
  monitor   <- data[,"n"] %in% looks



  # 2 Raise Alert
  alertNotROPE  <- which(data[,"sgpvROPE"] == 0 & monitor==TRUE)
  alertNotROME  <- which(data[,"sgpvROME"] == 0 & monitor==TRUE)

  alertNotROPEAny  <- length(alertNotROPE) > 0
  alertNotROMEAny  <- length(alertNotROME) > 0



  # 3 Affirm alert and report end of study (eos) operating characteristics (oc)
  #   - To be applied across various k (affirmation steps)
  affirmEndOfStudy <- function(K){


    # Stop times for being Not ROPE and being Not ROME
    stopNotROPE  <- alertNotROPE[alertNotROPE %in% (alertNotROPE + K)]
    stopNotROME  <- alertNotROME[alertNotROME %in% (alertNotROME + K)]

    if(getUnrestricted==TRUE){

      stop <- min(stopNotROPE, stopNotROME, na.rm = TRUE)

      # Unrestricted sample size (immediate outcomes)
      eos                      <- data[data[,"n"]==stop,]
      eos["stopNotROPE"]       <- as.numeric(eos["sgpvROPE"]==0)
      eos["stopNotROME"]       <- as.numeric(eos["sgpvROME"]==0)
      eos["stopInconclusive"]  <- as.numeric(eos["stopNotROPE"]==0 & eos["stopNotROME"]==0)

      if(lagN > 0){
        # Unrestricted sample size (stopping and then observing lagged outcomes)
        eosLag                     <- data[data[,"n"]==stop + lagN,]
        eosLag["stopNotROPE"]      <- as.numeric(eosLag["sgpvROPE"]==0)
        eosLag["stopNotROME"]      <- as.numeric(eosLag["sgpvROME"]==0)
        eosLag["stopInconclusive"] <- as.numeric(eosLag["stopNotROPE"]==0 &
                                                 eosLag["stopNotROME"]==0)
        eosLag["stopInconsistent"] <-
          as.numeric( ( eos["stopNotROPE"] == 1 & eos["stopNotROME"] != 1 & eosLag["stopNotROPE"] != 1 )  |
                      ( eos["stopNotROPE"] != 1 & eos["stopNotROME"] == 1 & eosLag["stopNotROME"] != 1 )  )


      } else {
        eosLag <- NULL
      }

    } else eos <- eosLag <- NULL

    # MaxN stop
    if(!is.null(maxN)){


      stop <- min(stopNotROPE, stopNotROME, maxN, na.rm = TRUE)

      # Maximum sample size of maxN (immediate outcomes)
      eosMaxN                     <- data[data[,"n"]==stop,]
      eosMaxN["stopNotROPE"]      <- as.numeric(eosMaxN["sgpvROPE"]==0)
      eosMaxN["stopNotROME"]      <- as.numeric(eosMaxN["sgpvROME"]==0)
      eosMaxN["stopInconclusive"] <- as.numeric(eosMaxN["stopNotROPE"]==0 &
                                                eosMaxN["stopNotROME"]==0)

      if(lagN > 0){
        # Maximum sample size of maxN (stopping and then observing lagged outcomes)
        eosLagMaxN                     <- data[data[,"n"]==min(stop+lagN,maxN),]
        eosLagMaxN["stopNotROPE"]      <- as.numeric(eosLagMaxN["sgpvROPE"]==0)
        eosLagMaxN["stopNotROME"]      <- as.numeric(eosLagMaxN["sgpvROME"]==0)
        eosLagMaxN["stopInconclusive"] <- as.numeric(eosLagMaxN["stopNotROPE"]==0 &
                                                     eosLagMaxN["stopNotROME"]==0)
        if(data[,"n"]==stop+lagN){
          eosLagMaxN["stopInconsistent"] <-
            as.numeric( ( eosMaxN["stopNotROPE"] == 1 & eosMaxN["stopNotROME"] != 1 & eosLagMaxN["stopNotROPE"] != 1 )  |
                        ( eosMaxN["stopNotROPE"] != 1 & eosMaxN["stopNotROME"] == 1 & eosLagMaxN["stopNotROME"] != 1 )  )

        } else {
          eosLagMaxN["stopInconsistent"] <- 0
        }

      } else {
        eosLagMaxN <- NULL
      }

    } else eosMaxN <- eosLagMaxN <- NULL


    # Keep theta and operating characteristics (drop y, trt, and sgpvs)
    keepStats <- c("n","est","bias","rejH0","cover", "stopNotROPE","stopNotROME","stopInconclusive")


    if(!is.na(stop)){
      oc <- c(data[1,"theta0"],
              data[1,"effect1"],
              eos[keepStats],
              lag     = eosLag[    c(keepStats,"stopInconsistent")],
              maxN    = eosMaxN[     keepStats],
              lagMaxN = eosLagMaxN[c(keepStats,"stopInconsistent")])

      return(c(oc,affirmK=K))
    } else {
      return(NULL)
    }
  }

  affirmedEnd <- t(sapply(affirmK,affirmEndOfStudy))


  return(affirmedEnd)

}
