# Calculate DDCP for the full power study and briding study
## Hao Wang, 2025-July-7

# Load bpp
library(bpp)

# load rpact
library(rpact)

# get the sample size for a full power study (power is set to 0.83 to approximate the sample size of 300 evaluable patients)
sampleSizeResult1 <- rpact::getSampleSizeRates(pi2 = 0.1, pi1 = 0.02, riskRatio = FALSE,
                                               sided = 1, alpha = 0.025, beta = 0.166) 
summary(sampleSizeResult1)

n1 <- ceiling(sampleSizeResult1$numberOfSubjects1)
n2 <- ceiling(sampleSizeResult1$numberOfSubjects2)

mdd <- sampleSizeResult1$criticalValuesEffectScale
mdd


# DDCP calculation for full power study ----
# refers to https://pages.github.roche.com/successR/successR/binary_bpp_design.html#ddcp-at-the-beginning-of-the-trial

set.seed(2021)

# number of simulations to approximate integration
M <- 10 ^ 5

# priors before Phase 2: beta(a1, b1) refers to intervention,
# beta(a2, b2) to control
# we take here uniform priors to be uninformative and match
# the approximate analysis of above
a1 <- 1
b1 <- 1
a2 <- 1
b2 <- 1

# cases and number of patients in previous Phase 3 BLOCKSTONE subgroup: BMX treated IP
# these frequencies correspond to pi10 and pi20 above
r21 <- 5
n21 <- 195 
r22 <- 21 
n22 <- 197

#r21 <- 7
#n21 <- 374
#r22 <- 51
#n22 <- 375

# posterior after previous Phase 3 - simply beta-updating
pi21 <- rbeta(M, a1 + r21, b1 + n21 - r21)
pi22 <- rbeta(M, a2 + r22, b2 + n22 - r22)

# compute power for every simulated pair and average
ddcp_m <- rep(NA, M)
i <- 0:n1

# in the last inequality above we have a "<" sign. Subtracting a small number in the first 
# argument of pbinom achieves this since pbinom rounds down the first argument 
# when it is not an integer
q <- n1 * i / n2 - n1 * abs(mdd) - 1e-10
for (m in 1:M){
  ddcp_m[m] <- 1 - sum((1 - pbinom(q = q, size = n2, p = pi21[m])) * 
                         dbinom(x = i, size = n2, p = pi22[m]))
}

ddcp <- mean(ddcp_m)
ddcp


# The following code calculates the DDCP for bridging study - showing same trend ----
n1 <- n2 <- 80
mdd_1 <- 0

# compute power for every simulated pair and average
ddcp_m <- rep(NA, M)
i <- 0:n1

# in the last inequality above we have a "<" sign. Subtracting a small number in the first 
# argument of pbinom achieves this since pbinom rounds down the first argument 
# when it is not an integer
q <- n1 * i / n2 - n1 * abs(mdd_1) - 1e-10
for (m in 1:M){
  ddcp_m[m] <- 1 - sum((1 - pbinom(q = q, size = n2, p = pi21[m])) * 
                         dbinom(x = i, size = n2, p = pi22[m]))
}

ddcp_1 <- mean(ddcp_m)
ddcp_1

# The following code calculates the DDCP for bridging study - preserving half effect -----
mdd_2 <- 0.043 # (1-0.14)*0.5*0.1
# compute power for every simulated pair and average
ddcp_m <- rep(NA, M)
i <- 0:n1

# in the last inequality above we have a "<" sign. Subtracting a small number in the first 
# argument of pbinom achieves this since pbinom rounds down the first argument 
# when it is not an integer
q <- n1 * i / n2 - n1 * abs(mdd_2) - 1e-10
for (m in 1:M){
  ddcp_m[m] <- 1 - sum((1 - pbinom(q = q, size = n2, p = pi21[m])) * 
                         dbinom(x = i, size = n2, p = pi22[m]))
}

ddcp_2 <- mean(ddcp_m)
ddcp_2