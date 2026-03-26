testMedianTwoGroupsT2E <- function(time, event, group){
  
  require(survival)
  
  # test for equality of medians for T2E endpoint according to Chen and Zhang (2016)
  n <- length(time)
  
  # compute overall median
  s_all <- Surv(time, event)
  km_all <- survfit(s_all ~ rep(1, n))
  med_all <- quantile(km_all, probs = 0.5)$quantile
  
  # objects to save results
  sigma <- rep(NA, 2)
  eta <- sigma
  medians <- sigma
  
  for (i in 1:2){
    
    ind.i <- (group == names(table(group))[i])
    ni <- sum(ind.i)
    time1 <- time[ind.i]
    event1 <- event[ind.i]
    # evaluate KM at overall median in each group
    s1 <- Surv(time1, event1)
    km1 <- survfit(s1 ~ rep(1, ni), conf.type = "plain")
    medians[i] <- quantile(km1, probs = 0.5)$quantile
    skm1 <- summary(km1)
    eta[i] <- tail(skm1$surv[skm1$time < med_all], 1)
    
    # Greenwood variance at the overall median
    v11 <- (tail(skm1$std.err[skm1$time < med_all], 1)) ^ 2
    
    # second summand of variance in Cheng and Zhang
    # event time closest to med_1
    t11 <- time1[event[ind.i] == 1]   
    t11 <- t11[t11 != medians[i]]
    t11 <- t11[which(abs(t11 - medians[i]) == min(abs(t11 - medians[i])))]
    
    # evaluate F1 at the two values
    F1 <- min(skm1$surv[skm1$time < medians[i]])
    F11 <- min(skm1$surv[skm1$time < t11])
    
    # now compute variance
    v12 <- (F1 - F11) ^ 2 / 2
    sigma[i] <- sqrt(v11 + v12)
  }
  
  # compute test statistic
  w <- 1 / (sigma ^ 2)
  h <- w / sum(w)
  C <- sum(w * (eta - sum(h * eta)) ^ 2)
  
  pval <- 1 - pchisq(C, df = 1)
  delta <- diff(medians)
  
  # if median is not reached in at least one group put the $p$-value to 1
  if (is.na(pval)){pval <- 1}
  
  res <- list("pval" = pval, "delta median" = delta)
  return(res)
}
