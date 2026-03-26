pval2sampleT2E <- function(time, event, arm, t0, wrho = 0, wgamma = 1, rho = c(0, 0, 1), gamma = c(0, 1, 0)){
    
  # time:   Time-to-event.
  # event:  Event indicator, 1 = event, 0 = censored.
  # arm:    Treatment arm.
  # t0:     Milestone timepoint of interest.
  # wrho:   Fleming-Harrington weight function for weighted logrank (real number).
  # wgamma: Fleming-Harrington weight function for weighted logrank (real number).
  # rho:    Fleming-Harrington weight functions for maxcombo (can be vector).
  # gamma:  Fleming-Harrington weight functions for maxcombo (can be vector).
  
    require(survival)
    require(nph)
    require(survRM2)

    nams <- c("logrank", "wlogrank", "maxcombo", "milestone", "median", "RMST")
    pvals <- rep(NA, length(nams))
    names(pvals) <- nams
    
    # -------------------
    # logrank p-value (score test from Cox regression)
    # -------------------
    logrank_pval <- summary(coxph(Surv(time, event) ~ arm))$sctest["pvalue"]
    pvals["logrank"] <- logrank_pval
  
    # -------------------
    # weighted logrank p-value
    # (rho, gamma) = 
    # (0, 0): standard logrank
    # (0, 1): more weight on late events
    # (1, 0): more weight on early
    # (1, 1): more weight on intermediate events
    # alternatively: survMisc::comp, 
    # -------------------
    wlogrank_pval <- nph::logrank.test(time = time, event = event, group = arm, rho = wrho, gamma = wgamma)
    pvals["wlogrank"] <- wlogrank_pval$test["p"][[1]]
  
    # -------------------
    # max-combo test
    # -------------------
    maxcombo_pval <- nph::logrank.maxtest(time = time, event = event, group = arm, rho = rho, gamma = gamma)
    pvals["maxcombo"] <- maxcombo_pval$pmult
  
    # -------------------
    # milestone, using cloglog trafo
    # -------------------
    # milestone_pval <- ComparisonSurv::Fixpoint.test(time = time, status = event, group = arm, t0 = t0)$test[, "pvalue"]
    # pvals["milestone1"] <- milestone_pval1[3]
    # milestone_pval <- bpcp::fixtdiff(time = time, status = event, group = arm, testtime = t0, trans = "cloglog")$p2  # not used b/c of warning in description file
    # pvals[i, "milestone"] <- as.numeric(milestone_pval)
    suppressWarnings(milestone_pval <- milestoneTest(time = time, event = event, tmt = arm, t0 = t0))
    pvals["milestone"] <- milestone_pval
  
    # -------------------
    # test for median difference
    # -------------------
    suppressWarnings(median_pval <- testMedianTwoGroupsT2E(time = time, event = event, group = arm))
    pvals["median"] <- median_pval$pval
  
    # -------------------
    # RMST
    # -------------------
    rmst_pval <- survRM2::rmst2(time = time, status = event, arm = arm, tau = NULL)  
    pvals["RMST"] <- rmst_pval$unadjusted.result[1, "p"]
    
    return(pvals)
}




