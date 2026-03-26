# load rpact
library(rpact)

# get the number of cases (infection) for a full power study
# pi1 is the proportion of cases occur in the treatment group among all cases across groups 
# the null hypothesis is that the treatment has no effect which is equivalent to pi1=0.5 (50% of cases in the treatment group)

riskratio <- 0.2
rate_pbo <- 0.10      # Placebo transmission rate
rate_bxm <- rate_pbo * riskratio          # BXM transmission rate
pi1 <- rate_bxm/(rate_bxm+rate_pbo)
totalCases <- rpact::getSampleSizeRates(pi1 = pi1, groups = 1,
                                               normalApproximation = FALSE, thetaH0 = 0.5,
                                               sided = 1, alpha = 0.025, beta = 0.2) 
summary(totalCases)

# change x to find the maximum number of cases in treatment group that providing a positive result (can also get MDD)
binom.test(x = ceiling(totalCases$maxNumberOfSubjects*pi1), 
           n = totalCases$maxNumberOfSubjects, 
           p = 0.5, 
           alternative = "less", 
           conf.level = 0.95)

# get sample size according to the total number of cases
# [factor] is a value to adjust the placebo infection rate 
factor <- seq(1, 0.5, -0.1)
N_bxm <- ceiling(totalCases$maxNumberOfSubjects*pi1/rate_bxm/factor)
N_pbo <- ceiling(totalCases$maxNumberOfSubjects*(1-pi1)/rate_pbo/factor)
N_tot <- N_bxm + N_pbo

# sample size under different placebo rates
data.frame(
  rate_pbo = factor*rate_pbo,
  N_tot = N_tot)


# get the number of cases for a bridging study ----
# get the initial number by setting alpha to 0.49999
totalCases_bridge <- rpact::getSampleSizeRates(pi1 = pi1, groups = 1,
                                        normalApproximation = FALSE, thetaH0 = 0.5,
                                        sided = 1, alpha = 0.49999, beta = 0.1) 
summary(totalCases_bridge)

N_cases_bridge <- totalCases_bridge$maxNumberOfSubjects

# get the sample size
# factor is a value to adjust the placebo rate 
factor <- seq(1, 0.5, -0.1)
N_bxm <- ceiling(totalCases_bridge$maxNumberOfSubjects*pi1/rate_bxm/factor)
N_pbo <- ceiling(totalCases_bridge$maxNumberOfSubjects*(1-pi1)/rate_pbo/factor)
N_tot <- N_bxm + N_pbo

# sample size under different placebo rates
data.frame(
  rate_pbo = factor*rate_pbo,
  N_tot = N_tot)


# get the probability of bridging ----
# change the following two parameters to test the robustness of the results
N_cases_bridge <- seq(3, 20, 1) # change total number of cases from 6 to 8
pi1 <- 0.2/(1+0.2)

rr_crit <- 1 - (1-0.14)*0.5 # the risk ratio threshold preserving half observed effect
cases_bxm_crit <- rr_crit/(1+rr_crit) * N_cases_bridge # the number of cases in the treatment group at the critical value
prob_bridge2 <- pbinom(q = cases_bxm_crit, size = N_cases_bridge, prob = pi1, lower.tail = TRUE)

prob_bridge1 <- pbinom(q = N_cases_bridge*0.5, size = N_cases_bridge, prob = pi1, lower.tail = TRUE)

# prob_bridge2 is not monotonic, however, it is always > 85% starting from 6 cases, thus could be a choice for a bridging study
data.frame(N_cases_bridge, prob_bridge1, prob_bridge2)
