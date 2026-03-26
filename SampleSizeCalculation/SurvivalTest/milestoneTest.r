milestoneTest <- function(time, event, tmt, t0){

  require(survival)
  
  indA <- (tmt == names(table(tmt))[1])
  tmp <- Surv(time[indA], event[indA]) ~ tmt[indA]
  fit1 <- survfit(tmp, conf.type = 'log-log')
  dat.km1 <- cbind(fit1$time, fit1$surv, fit1$std.err)
  dat.km1 <- dat.km1[dat.km1[, 1] < t0, ]
  dat.km1 <- dat.km1[nrow(dat.km1), ]
  St1 <- dat.km1[2]
  se1 <- dat.km1[3]
  # fit1sum <- summary(fit1, time = t0, extend = TRUE)
  # St1 <-  fit1sum$surv
  # se1 <- fit1sum$std.err
  
  indB <- !indA
  tmp <- Surv(time[indB], event[indB]) ~ tmt[indB]
  fit2 <- survfit(tmp, conf.type = 'log-log')
  dat.km2 <- cbind(fit2$time, fit2$surv, fit2$std.err)
  dat.km2 <- dat.km2[dat.km2[, 1] < t0, ]
  dat.km2 <- dat.km2[nrow(dat.km2), ]
  St2 <- dat.km2[2]
  se2 <- dat.km2[3]
  # fit2sum <- summary(fit2, time = t0, extend = TRUE)
  # St2 <-  fit2sum$surv
  # se2 <- fit2sum$std.err
  
  # chi^2 test statistic
  # from Klein et al (2007), Stat. Med.
  #X1 <- (St1 - St2) ^ 2 / (St1 ^ 2 * se1 ^ 2 + St2 ^ 2 * se2 ^ 2)
  X3 <- (log(-log(St1)) - log(-log(St2))) ^ 2 / (se1 ^ 2 / (log(St1) ^ 2) + se2 ^ 2 / (log(St2) ^ 2))

  # p-value  
  # same as ComparisonSurv::Fixpoint.test(time = time, status = event, group = arm, t0 = t0)
  res <- 1 - pchisq(X3, df = 1)
  
  # if milestone is not reached in either group put the $p$-value to 1
  if (is.na(res)){res <- 1}
  
  return(res)
}
  