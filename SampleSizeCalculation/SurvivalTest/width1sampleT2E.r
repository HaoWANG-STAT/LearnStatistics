width1sampleT2E <- function(med1, kappa = 1, t0, accrualTime, accrualIntensity, 
                            dropout, dropoutTime = 12, d, conf.level = 0.95, M = 10 ^ 3, 
                            print.progress = TRUE, seed = NA, 
                            which = c("milestone", "nonparametric", "exponential")){
  
  # med1:             Median for single-arm arm.
  # kappa:            Shape parameter of Weibull distribution.
  # t0:               Milestone timepoint of interest.
  # accrualTime:      See ?rpact::getSimulationSurvival.
  # accrualIntensity: See ?rpact::getSimulationSurvival.
  # dropout:          Drop-out rate, on same time scale as med0.
  # dropoutTime:      See ?rpact::getSimulationSurvival.
  # d:                #events the final censored data should have 
  # conf.level:       2-sided confidence level of confidence intervals.
  # M:                Number of simulation runs.
  # print.progress    Print every print.progess'th interation to indicate progress.
  #                   Set to M + epsilon if no progress to be shown.
  # seed:             Set if you want to chose seed. 
  
  require(rpact)
  require(survival)
  require(bpcp)
  
  ci_milestone <- matrix(NA, ncol = 2, nrow = M)
  ci_median <- ci_milestone
  ci_median_exp <- ci_milestone
  milestone <- rep(NA, M)
  median <- rep(NA, M)
  median_exp <- rep(NA, M)
  
  if (is.na(seed) == FALSE){set.seed(seed)}
  
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

    n <- length(time)
    tmp <- (Surv(time, event) ~ rep(1, n))
    obj <- survfit(tmp, conf.int = conf.level, conf.type = "log-log", type = "kaplan", 
                   error = "greenwood")

    # -------------------
    # confidence interval based on cloglog trafo
    # see Klein et al (2007) and Nagashima et al (2021)
    # -------------------
    # dat <- data.frame(cbind(obj$time, obj$surv, obj$lower, obj$upper))
    # if (max(dat[, 1]) >= t0){
    #   dat <- dat[dat[, 1] >= t0, ]
    #   ci_milestone[i, ] <- as.numeric(dat[1, 3:4])
    #   milestone[i] <- dat[1, 2]
    # }
    fitsum <- summary(obj, time = t0, extend = TRUE) 
    ci_milestone[i, ] <-  c(fitsum$lower,fitsum$upper) 
    milestone[i] <- fitsum$surv
    
    # -------------------
    # nonparametric confidence interval for median
    # -------------------
    #fitBP <- bpcp(time = time, status = event, alpha = 1 - conf.level)
    #ci <- median(fitBP)[1, ][3:4]
    ci <- c(quantile(obj, probs = 0.5)$lower, quantile(obj, probs = 0.5)$upper)
    ci_median[i, ] <- ci
    median[i] <- quantile(obj, probs = 0.5)$quantile
    
    # -------------------
    # parametric confidence interval for median based on exponential assumption
    # -------------------
    r <- sum(event)
    mle_rate <- sum(event) / sum(time)
    med_exp <- log(2) / mle_rate
    ci_median_exp[i, ] <- med_exp + c(-1, 1) * qnorm((1 + conf.level) / 2) * med_exp / sqrt(r)
    median_exp[i] <- med_exp

    if (i / print.progress == round(i / print.progress)){print(paste("run ", i, " of ", M, " done", sep = ""))}
  }
  
  # all CIs and widths
  width <- cbind(apply(ci_milestone, 1, diff), apply(ci_median, 1, diff), apply(ci_median_exp, 1, diff))
  colnames(width) <- c("milestone", "median nonparametric", "median exponential")
  res <- list("milestone" = milestone, "ci_milestone" = ci_milestone, "median" = median, 
              "ci_median" = ci_median, "median_exp" = median_exp, "ci_median_exp" = ci_median_exp, 
              "width" = width)
  return(res)
}





