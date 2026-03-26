library(tidyverse)
library(bpp)
library(rpact)

# specifications of the pivotal trial
# set up a group-sequential design with interim at 63% and 84% of information
# significance levels at final analysis
design <- getDesignGroupSequential(sided = 2, alpha = 0.049, beta = 0.15,
                                   informationRates = c(.63, .84, 1),
                                   typeOfDesign = "asOF")

# number of events needed
sampleSize <- getSampleSizeSurvival(design = design, hazardRatio = 0.7)
nevents <- ceiling(sampleSize$maxNumberOfEvents)
nevents

# MDD at final
hrMDD <- sampleSize$criticalValuesEffectScaleLower[3, 1]
hrMDD

# Calculate DDCP at the beginning of the trial ----
thetas <- seq(0.3, 2.7, by = 0.01)

bpp0 <- bpp::bpp_t2e(prior = "normal", successHR = hrMDD, d = nevents,
                     priorHR = 0.7, priorsigma = sqrt(4/50))
bpp0


# validate DDCP using simulation
nsim <- 100000

power <- pnorm(log(hrMDD),log(0.7),sqrt(4/nevents))
power

prior <- rnorm(nsim, log(0.7), sqrt(4/50)) #Normal prior corresponding to information of 50 events in 1:1 randomized trial
post <- rep(0, nsim)
for(i in 1:nsim) {
  post[i] <- rnorm(1, prior[i], sqrt(4/nevents))
}

data.frame(prior, post) %>% ggplot()+
  geom_freqpoly(aes(prior), color="red", alpha=0.2)+
  geom_freqpoly(aes(post), color="green", alpha=0.3)

cpts <- sum(exp(post)<=hrMDD)/nsim
cpts


# Update prior and DDCP according to new external information ----
# effect estimate and SE
hr_ext <- 0.81
se_ext <- sqrt(4 / 245)

up1 <- NormalNormalPosterior(datamean = log(hr_ext), sigma = se_ext, n = 1, nu = log(0.7), 
                             tau = sqrt(4/50))

bpp1 <- bpp_t2e(prior = "normal", successHR = hrMDD, d = nevents,
                priorHR = exp(up1$postmean), priorsigma = up1$postsigma)

# Normal prior
up1


# Compute bpp after not stopping at interim ----
# specifications for the interim analysis
# local significance levels based on rpact output
alphas <- 2 * design$stageLevels
alphas
# number of events
nevents <- ceiling(as.vector(sampleSize$eventsPerStage))
nevents
# MDDs at interim and final based on rpact output
hrMDD <- as.vector(sampleSize$criticalValuesEffectScaleLower)
hrMDD

# efficacy boundary --> MDD at interim analysis
effi <- hrMDD[2]

# futility boundary --> chosen informally
futi <- 1


# assuming we do not stop at a blinded interim analysis for efficacy only:
bpp3 <- bpp_1interim_t2e(prior = "normal", successHR = hrMDD[3], d = nevents[c(2,3)],
                                   IntEffBoundary = effi, IntFutBoundary = 1, IntFixHR = 0.81, 
                                   priorHR = exp(log(0.7)), propA = 0.5, thetas = thetas, 
                                   priorsigma = sqrt(4/50)) # initial prior
bpp3$"BPP after not stopping at interim exact"


bpp4 <- bpp_1interim_t2e(prior = "normal", successHR = hrMDD[3], d = nevents[c(2,3)],
                                   IntEffBoundary = effi, IntFutBoundary = 1, IntFixHR = 0.81, 
                                   priorHR = exp(log(0.81)), propA = 0.5, thetas = thetas, 
                                   priorsigma = sqrt(4/245)) # prior based on interim data
bpp4$"BPP after not stopping at interim exact"

bpp5 <- bpp_1interim_t2e(prior = "normal", successHR = hrMDD[3], d = nevents[c(2,3)],
                         IntEffBoundary = effi, IntFutBoundary = 1, IntFixHR = 0.81, 
                         priorHR = exp(up1$postmean), propA = 0.5, thetas = thetas, 
                         priorsigma = up1$postsigma)  # updated prior according to interim data
bpp5$"BPP after not stopping at interim exact"





