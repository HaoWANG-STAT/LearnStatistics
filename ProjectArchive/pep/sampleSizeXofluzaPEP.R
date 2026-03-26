# Sample size calculation for Xofluza PEP study
# Author: Hao Wang
# Date: 14-Aug-2023

alpha <- 0.05
beta <- 0.1
p0 <- c(0.1, 0.1362, 0.1, 0.05)  #incidence of control
p1 <- c(0.04, 0.019, 0.02, 0.01) #incidence of treatment
rd <- p0-p1 #risk difference
rr <- p1/p0 #risk ratio
k <- 1 #randomization ratio

# for comparing risk difference using Z test, refer to sample size calculation in clinical research book----
n0.rd <- (qnorm(1-alpha/2)+qnorm(1-beta))^2/rd^2*(p1*(1-p1)/k + p0*(1-p0))
N.rd <- n0.rd*(k+1)
N.rd

# for comparing risk ratio, refer to page 458 of sample size calculation in clinical research book----
n0.rr <- (qnorm(1-alpha/2)+qnorm(1-beta))^2/log(rr)^2*((1-p1)/p1 + (1-p0)/p0)
N.rr <- n0.rr*(k+1)
N.rr

# calculation, refer to cancer epidemiology book page 123 ----
m0 <- 100 # initial sample size
m1 <- m0*k
a <- m1*p1
b <- m0*p0

se.rd.initial <- sqrt(p1*(1-p1)/m1 + p0*(1-p0)/m0) #standard deviation for risk difference
se.rd.require <- abs(rd)/(qnorm(1-alpha/2)+qnorm(1-beta))
scale.factor.rd <- (se.rd.initial/se.rd.require)^2
scale.factor.rd
N.rd.2 <- (m0+m1)*scale.factor.rd
N.rd.2

se.log.rr.initial <- sqrt(1/a+1/b-1/m1-1/m0) #standard deviation for log(rr)
se.log.rr.require <- abs(log(rr))/(qnorm(1-alpha/2)+qnorm(1-beta))
scale.factor.rr <- (se.log.rr.initial/se.log.rr.require)^2
scale.factor.rr
N.rr.2 <- (m0+m1)*scale.factor.rr
N.rr.2


# validation using function in rpact ----
## .getSampleSizeFixedRates is used in rpact to calculate sample size based on Farrington-Manning method
## .getEffectScaleBoundaryDataRates() is used to calcuate MDD

sampleSizeResult1 <- rpact::getSampleSizeRates(pi2 = 0.1, pi1 = 0.04, riskRatio = TRUE,
                                               sided = 1, alpha = 0.025, beta = 0.1) 
summary(sampleSizeResult1)


sampleSizeResult2 <- rpact::getSampleSizeRates(pi2 = 0.136, pi1 = 0.019, riskRatio = TRUE,
                                               sided = 1, alpha = 0.025, beta = 0.1)
summary(sampleSizeResult2)


sampleSizeResult3 <- rpact::getSampleSizeRates(pi2 = 0.05, pi1 = 0.01, riskRatio = TRUE,
                                               sided = 1, alpha = 0.025, beta = 0.2) 
summary(sampleSizeResult3)




# Non-inferiority
sampleSizeResult3 <- rpact::getSampleSizeRates(pi2 = 0.05, pi1 = 0.05, riskRatio = TRUE, thetaH0 = 1.5,
                                               sided = 1, alpha = 0.025, beta = 0.1)
summary(sampleSizeResult3)


# validation using function in gsDesign ----
gsDesign::nBinomial(p1 = .04, p2 = .1, delta0 = 0, alpha = .025, sided = 1, beta = .1, scale = "RR")


# validation using function power.prop.test ----
power.prop.test(p1 = .04, p2 = .1, power = .90)

# Repeat the power for CENTERSTONE study
power.prop.test(n = 2030/2, p1 = .14, p2 = .20)      ## => power = 0.95
power.prop.test(n = 2030/2, p1 = .14*.65, p2 = .20*.65)      ## => power = 0.80


# plot for power change with RR ----
library(tidyverse)

data.frame(rpact::getPowerRates(pi2 = 0.136, pi1 = 0.136*rr, riskRatio = TRUE, directionUpper = FALSE,
                                sided = 1, alpha = 0.025, maxNumberOfSubjects=350))$overallReject


parameters <- data.frame(rate.ratio=seq(0.15, 0.3, 0.05)) %>%
  cross_join(data.frame(placebo.rate=seq(0.07, 0.2, 0.01)))
parameters

rslt <- parameters %>% 
  group_by(placebo.rate, rate.ratio) %>%
  do(mutate(., power=data.frame(rpact::getPowerRates(pi2 = placebo.rate, pi1 = 0.1*rate.ratio, riskRatio = TRUE, 
                                                     directionUpper = FALSE, sided = 1, alpha = 0.025, 
                                                     maxNumberOfSubjects=325))$overallReject))


rslt %>% ggplot(aes(x=placebo.rate, y=power, color=factor(rate.ratio)))+
  geom_line()+
  theme_bw()+
  geom_hline(yintercept = 0.8, linetype=4, color="red")+
  labs(x = 'Placebo rate', 
       y = paste0('Power'),
       color = "Risk ratio"
       # title = 'Powered treatment effect under different sample size'
       )+
  scale_x_continuous(breaks = seq(0.07, 0.2, 0.01))


# plot for sample size change with RR ----
rpact::getSampleSizeRates(pi2 = 0.1, pi1 = 0.02, riskRatio = TRUE, thetaH0 = 0.6, 
                          sided = 1, alpha = 0.49999, beta = 0.2)


parameters2 <- data.frame(rate.ratio=seq(0.15, 0.3, 0.05)) %>%
  cross_join(data.frame(control.rate=seq(0.05, 0.2, 0.05)))
parameters2

rslt2 <- parameters2 %>% 
  group_by(control.rate, rate.ratio) %>%
  do(mutate(., sample.size=data.frame(rpact::getSampleSizeRates(pi2 = control.rate, pi1 = control.rate*rate.ratio, riskRatio = TRUE, 
                                                     sided = 1, alpha = 0.025, beta=0.2))$nFixed))

rslt2 %>% ggplot(aes(x=rate.ratio, y=sample.size, color=factor(control.rate)))+
  geom_line()+
  theme_bw()+
  #geom_hline(yintercept = 200, linetype=4, color="red")+
  labs(x = 'risk ratio', 
       y = paste0('Powered sample size to detect \n', 'corresponding risk ratio'),
       color = "Control rate"
       # title = 'Powered treatment effect under different sample size'
  )

# find MDD under fixed sample size (if control rate is lower than expected) ----
parameters3 <- data.frame(control.rate=seq(0.05, 0.1, 0.01)) %>%
  cross_join(data.frame(rate.ratio=seq(0.1, 0.4, 0.01)))
parameters3

rslt3 <- parameters3 %>% 
  group_by(control.rate, rate.ratio) %>%
  do(mutate(., power=data.frame(rpact::getPowerRates(pi2 = control.rate, pi1 = control.rate*rate.ratio, riskRatio = TRUE, 
                                                     directionUpper = FALSE, sided = 1, alpha = 0.025, 
                                                     maxNumberOfSubjects=274))$overallReject))


# plot for sample size change with placebo rate (fix BXM rate) ----

parameters4 <- data.frame(bxm.rate=0.025) %>%
  cross_join(data.frame(control.rate=seq(0.07, 0.15, 0.01)))
parameters4

rslt4 <- parameters4 %>% 
  group_by(bxm.rate, control.rate) %>%
  do(mutate(., sample.size=data.frame(rpact::getSampleSizeRates(pi2 = control.rate, pi1 = bxm.rate, riskRatio = TRUE, 
                                                                sided = 1, alpha = 0.025, beta=0.2))$nFixed))

rslt4 %>% ggplot(aes(x=control.rate, y=sample.size, color=factor(bxm.rate)))+
  geom_line()+
  theme_bw()+
  #geom_hline(yintercept = 200, linetype=4, color="red")+
  labs(x = 'Placebo rate', 
       y = paste0('Sample size'),
       color = "BXM rate"
       # title = 'Powered treatment effect under different sample size'
  )


# plot for sample size change with placebo rate (fix RR) ----
parameters5 <- data.frame(rate.ratio=c(0.15, 0.2, 0.25, 0.3)) %>%
  cross_join(data.frame(control.rate=seq(0.07, 0.2, 0.01)))
parameters5

rslt5 <- parameters5 %>% 
  group_by(rate.ratio, control.rate) %>%
  do(mutate(., sample.size=data.frame(rpact::getSampleSizeRates(pi2 = control.rate, pi1 = control.rate * rate.ratio, riskRatio = TRUE, 
                                                                sided = 1, alpha = 0.025, beta=0.2))$nFixed))

rslt5 %>% 
  ggplot(aes(x=control.rate, y=sample.size/0.9, color=factor(rate.ratio)))+ # divide by 0.9 to adjust for 10% drop out
  geom_line()+
  theme_bw()+
  #geom_hline(yintercept = 200, linetype=4, color="red")+
  labs(x = 'Placebo rate', 
       y = paste0('Sample size'),
       color = "Risk ratio"
       # title = 'Powered treatment effect under different sample size'
  ) +
  scale_x_continuous(breaks = seq(0.07, 0.2, 0.01)) +
  scale_y_continuous(breaks = seq(0, 900, 50))
