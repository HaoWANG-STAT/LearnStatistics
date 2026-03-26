# The code is used to calculate the DDCP from a event driven study design angle
# assuming the proportion of number of cases in the active arm follows a binomial distribution conditional on total number of cases 

library(extraDistr)
library(rpact)

# new full power study design
beta <- 0.2
riskratio <- 0.2
rate_pbo <- 0.10      # Placebo rate
rate_bxm <- rate_pbo * riskratio          # BXM infection rate
pi1 <- rate_bxm/(rate_bxm+rate_pbo)



# the total number of cases get from getSampleSizeRates does not use a conservative approach 
# (ie, there is no larger number of cases in the input n where power is less than 1-beta.)
#------eventDriven <- rpact::getSampleSizeRates(pi1 = pi1, groups = 1,
#------                                        normalApproximation = FALSE, thetaH0 = 0.5,
#------                                        sided = 1, alpha = 0.025, beta = beta) 
#------summary(eventDriven)
#------
#------
#------# get the total number of cases and minimum number of cases in active arm for a positive study
#------ncases <- eventDriven$numberOfSubjects
#------case_bxm_cutoff <- eventDriven$criticalValuesEffectScale * ncases


# nBinomial1Sample function provides an option that allows to use a conservative approach
ncases <- gsDesign::nBinomial1Sample(p0 = pi1, p1 = 0.5, alpha = 0.025, beta = beta, n = 10:100, conservative = TRUE, outtype = 2)$n
ncases
# find the maximum number of cases in treatment group that providing a positive result p<0.025 (can also get MDD)
for (i in ncases:1) {
  p.value <- binom.test(x = i, 
           n = ncases, 
           p = 0.5, 
           alternative = "less", 
           conf.level = 0.95)$p.value
  # if xx <0.025 breat and output i
  if (p.value < 0.025) {
    case_bxm_cutoff <- i
    break
  }
}

# MDD for RR
case_bxm_cutoff/(ncases-case_bxm_cutoff)
(case_bxm_cutoff+1)/(ncases-case_bxm_cutoff-1)

# Prior before study: beta(a, b) 
# we take here uniform priors to be uninformative
a <- 1
b <- 1

# Data observed in previous study BLOCKSTONE subgroup: BMX treated IP
e1 <- 5    # Number of cases in BXM arm in previous Phase 3 study
e2 <- 21     # Number of cases in placebo arm in previous Phase 3 study

# Posterior after previous study
alpha_post <- 1 + e1
beta_post <- 1 + e2

# DDCP is the probability of observing x cases or less in the new study
# while the predictive posterior distribution is a beta-binomial distribution
ddcp <- pbbinom(q = case_bxm_cutoff, size = ncases, alpha = alpha_post, beta = beta_post)
ddcp


# compute the DDCP for bridging study (the target total number of cases is 6)
ncases_bridge <- 6

ddcp_bridge1 <- pbbinom(ncases_bridge*0.5, size = ncases_bridge, alpha = alpha_post, beta = beta_post)

rr_crit <- 1 - (1-0.14)*0.5 # the risk ratio threshold preserving half observed effect
cases_bxm_crit <- rr_crit/(1+rr_crit) * ncases_bridge # the number of cases in the treatment group at the critical value
ddcp_bridge2 <- pbbinom(cases_bxm_crit, size = ncases_bridge, alpha = alpha_post, beta = beta_post)
c(ddcp, ddcp_bridge1, ddcp_bridge2)
