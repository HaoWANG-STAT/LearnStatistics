power1sampleT2E <- function(med1, med0, kappa = 1, t0, accrualTime, accrualIntensity, dropout, 
                            dropoutTime = 12, d, alpha = 0.05, lambda0 = 1, M = 10 ^ 3,
                            print.progress = TRUE, seed = NA){
  
  # med1:             Median for experimental arm.
  # med0:             Median of historical control.
  # kappa:            Shape parameter of Weibull distribution.
  # t0:               Milestone timepoint of interest.
  # accrualTime:      See ?rpact::getSimulationSurvival.
  # accrualIntensity: See ?rpact::getSimulationSurvival.
  # dropout:          Drop-out rate, on same time scale as med1.
  # dropoutTime:      See ?rpact::getSimulationSurvival.
  # d:                #events the final censored data should have 
  # alpha:            Significance level (1-sided) of all the tests.
  # lambda0:          Null hypothesis hazard ratio.
  # M:                Number of simulation runs.
  # print.progress:   Print every print.progess'th interation to indicate progress.
  #                   Set to M + epsilon if no progress to be shown.
  # seed:             Set if you want to chose seed. 
  
  require(rpact)
  require(survival)
  require(bpcp)
  
  reject <- matrix(NA, nrow = M, ncol = 5)
  colnames(reject) <- c("milestone (cloglog)", "milestone (exact)", "median", "logrank", "exponential")

  if (is.na(seed) == FALSE){set.seed(seed)}
  
  S0t0 <- 1 - pexp(t0, rate = log(2) / med0)
  
  # -------------------
  # generate data
  # 2-arm trial in rpact with same median in both groups
  # -------------------
  trial1 <- getSimulationSurvival(median1 = med1, median2 = med1, kappa = kappa, 
                                  dropoutRate1 = dropout, dropoutRate2 = dropout, 
                                  dropoutTime = dropoutTime, accrualTime = accrualTime, 
                                  accrualIntensity = accrualIntensity,
                                  plannedEvents = d, maxNumberOfIterations = M,
                                  maxNumberOfRawDatasetsPerStage = M)
  trial2 <- getRawData(trial1)
  
  for (i in 1:M){
    
    time <- trial2[trial2$iterationNumber == i, "timeUnderObservation"]
    event <- as.numeric(trial2[trial2$iterationNumber == i, "event"])
    
    # ----------------------------------
    # remove patients that arrived after cutoff
    # ----------------------------------
    rem <- (time >= 0)
    time <- time[rem]
    event <- event[rem]

    # -------------------
    # test for milestone test based on log-minus-log transformation
    # see Klein et al (2007) and Nagashima et al (2021)
    # check if null hypothesis milestone is outside of 1 - 2 * alpha CI (1 - sided test)
    # alternative approach: bpcp::bpcp
    # -------------------
    n <- length(time)
    tmp <- (Surv(time, event) ~ rep(1, n))
    obj <- survfit(tmp, conf.int = 1 - 2 * alpha, conf.type = "log-log", type = "kaplan", 
                   error = "greenwood")
    # dat <- cbind(obj$time, obj$surv, obj$lower, obj$upper)
    # dat <- dat[dat[, 1] >= t0, ]
    # ci1 <- dat[1, 3:4] 
    # reject[i, "milestone (cloglog)"] <- (ci1[1] > S0t0)
    reject[i, "milestone (cloglog)"] <- (summary(obj, time = t0, extend = TRUE)$lower > S0t0)
    
    # exact 
    fitexact1 <- bpcp::bpcp(time = time, status = event, alpha = 2 * alpha, stype = "km", midp = FALSE)
    ci2 <- as.numeric(StCI(fitexact1, tstar = t0))[3:4]
    reject[i, "milestone (exact)"] <- (ci2[1] > S0t0)

    # -------------------
    # test for median based on mid-p value
    # check if null hypothesis median is outside of 1 - 2 * alpha CI (1 - sided test)
    # -------------------
    #ci4 <- as.numeric(median(fitexact1)[1, ][3:4])
    ci4 <- c(quantile(obj, probs = 0.5)$lower, quantile(obj, probs = 0.5)$upper)
    reject[i, "median"] <- (ci4[1] > med0)
    
    # -------------------
    # one-sample logrank test
    # historical control is weibull --> see ?pweibull for cumulative hazard
    # -------------------
    obse <- sum(event)
      
    # cumulative hazard evaluated at all event times
    expe <- sum(-pweibull(time, shape = kappa, scale = med0 / log(2), lower = FALSE, log = TRUE))
    teststat <- (obse - lambda0 * expe) / sqrt(lambda0 * expe)
    reject[i, "logrank"] <- (teststat < -qnorm(1 - alpha))

    # -------------------
    # one-sample parametric test based on exponential
    # Wald test statistic Cook & deMets, p. 206
    # -------------------
    null_rate <- log(2) / med0
    teststat2 <- (obse - null_rate * sum(time)) ^ 2 / obse
    reject[i, "exponential"] <- (teststat2 > qchisq(1 - alpha, df = 1))
    # -------------------

    if (i / print.progress == round(i / print.progress)){print(paste("run ", i, " of ", M, " done", sep = ""))}
  }
  
  # power to reject H0
  res <- apply(reject, 2, mean)
  return(res)
}