power2sampleT2E <- function(piecewiseSurvivalTime, med1, med2, kappa = 1, t0, 
                            accrualTime, accrualIntensity, dropout, dropoutTime = 12,
                            d, alpha = 0.05, M = 10 ^ 3, print.progress, 
                            seed = NA, wrho = 0, wgamma = 1, 
                            rho = c(0, 0, 1), gamma = c(0, 1, 0)){
    
  # piecewiseSurvivalTime:  See ?rpact::getSimulationSurvival.
  # med1:                   Median for experimental arm.
  # med0:                   Median of historical control.
  # kappa:                  Shape parameter of Weibull distribution.
  # t0:                     Milestone timepoint of interest.
  # accrualTime:            See ?rpact::getSimulationSurvival.
  # accrualIntensity:       See ?rpact::getSimulationSurvival.
  # dropout:                Drop-out rate, on same time scale as med1.
  # dropoutTime:            See ?rpact::getSimulationSurvival.
  # d:                      #events the final censored data should have 
  # alpha:                  Significance level (1-sided) of all the tests.
  # M:                      Number of simulation runs.
  # print.progress          Print every print.progess'th interation to indicate progress.
  #                         Set to M + epsilon if no progress to be shown.
  # seed:                   Set if you want to chose seed. 
  # wrho:   Fleming-Harrington weight function for weighted logrank (real number).
  # wgamma: Fleming-Harrington weight function for weighted logrank (real number).
  # rho:    Fleming-Harrington weight functions for maxcombo (can be vector).
  # gamma:  Fleming-Harrington weight functions for maxcombo (can be vector).
  
  require(rpact)

  if (is.na(seed) == FALSE){set.seed(seed)}
  
  # names of tests
  nams <- c("logrank", "wlogrank", "maxcombo", "milestone", "median", "RMST")
  
  # set up matrices to collect results
  pvals <- matrix(NA, nrow = M, ncol = length(nams))
  colnames(pvals) <- nams
    
  # ----------------------------------
  # simulate M clinical trials using rpact
  # ----------------------------------
  trial1 <- getSimulationSurvival(lambda1 = log(2) / med1, lambda2 = log(2) / med2, kappa = kappa, 
                                    piecewiseSurvivalTime = piecewiseSurvivalTime,
                                    dropoutRate1 = dropout[1], dropoutRate2 = dropout[2], 
                                    dropoutTime = dropoutTime, accrualTime = accrualTime, 
                                    accrualIntensity = accrualIntensity, plannedEvents = d, 
                                    maxNumberOfIterations = M,
                                    maxNumberOfRawDatasetsPerStage = M)
  trial2 <- getRawData(trial1)
    
  for (j in 1:M){
      
      time <- trial2[trial2$iterationNumber == j, "timeUnderObservation"]
      event <- as.numeric(trial2[trial2$iterationNumber == j, "event"])
      arm <- as.numeric(trial2[trial2$iterationNumber == j, "treatmentGroup"] == 1)    
      # 0 = control, 1 = treatment
      
      # ----------------------------------
      # remove patients that arrived after cutoff
      # ----------------------------------
      rem <- (time >= 0)
      time <- time[rem]
      event <- event[rem]
      arm <- arm[rem]
      
      # ----------------------------------
      # compute all tests
      # ----------------------------------
      pvals[j, ] <- pval2sampleT2E(time = time, event = event, arm = arm, t0 = t0, 
                                   wrho = wrho, wgamma = wgamma, rho = rho, gamma = gamma)
      
      if (j / print.progress == round(j / print.progress)){print(paste("run ", j, " of ", M, " done", sep = ""))}
  }
    
  # power to reject H0
  power <- apply(pvals <= alpha, 2, mean)
  return(power)
}
  





