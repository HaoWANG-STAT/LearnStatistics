library(foreach)
library(tidyverse)
library(here)
require(sandwich)
require(lmtest)


# Generate data from binomial distribution

# Input: assumptions for simulation
# Output: 
#    - Estimated RR for the simulation
#    - p value for the simulation
#    - MDD for the simulation



simBinomial = function(null.rate = 0.1,
                       rate.ratio = 0.4,
                       N.total = 200,
                       N.arms = 2,
                       alpha = 0.05,
                       dropout = 0.1){
  
  N.perArm = floor((1-dropout)*N.total/N.arms)
  
  # simulated data
  d = data.frame(
    arm = 0,
    event = rbinom(N.perArm,size=1,prob=null.rate)
  ) %>% rbind(
    data.frame(
      arm = 1,
      event = rbinom(N.perArm,size=1,prob=null.rate*rate.ratio) 
    )
  )
  
  # Poisson regression
  mod = glm(event ~ arm, data = d, family=poisson(link=log))
  
  sd   = sqrt(summary(mod)$cov.scaled[2,2])
  coef = summary(mod)$coefficients[2]
  p = 2*(1 - pt(abs(coef/sd),df=2*N.total-2))
  mdd = exp(qt(alpha/2,mod$df.null, lower.tail=TRUE)*sd)
  
  
  # Modified Poisson regression
  p2 = coeftest(mod, vcov = sandwich)[2, 4] # robust P value
  sd2 = coeftest(mod, vcov = sandwich)[2, 2] # robust std.err
  mdd2 = exp(qt(alpha/2,mod$df.null, lower.tail=TRUE)*sd2)
  
  # summary results
  data.frame(
    rr = exp(coef),
    p = p,
    p2 = p2,
    mdd = mdd,
    mdd2 = mdd2
  )
}



# -> Power Simulation -----------------------------------------
dropout = 0
alpha = 0.05
placebo.rate = c(0.1)
exp.rate = seq(0.01, 0.02, 0.01)
rate.ratio = round(exp.rate/placebo.rate, 2)

nsim = 10000
sample.size = seq(200, 300, 20)


# set parameters for all scenarios
parameters = data.frame(placebo.rate, rate.ratio) %>%
  cross_join(data.frame(sample.size))

# run simulations
sim.data = parameters %>% 
  group_by(rate.ratio, sample.size) %>% 
  do(
    foreach(1:nsim, .combine='rbind') %do% simBinomial(null.rate = .$placebo.rate,
                                                       rate.ratio = .$rate.ratio, 
                                                       N.total = .$sample.size,
                                                       N.arms = 2, 
                                                       alpha = alpha, 
                                                       dropout = dropout)
  ) 

sim.rslt = sim.data %>% 
  summarise(power = mean(p<0.05)
            , power2=mean(p2<0.05)
            , mdd=median(mdd)
            , mdd2=median(mdd2)
            , .groups="keep")


save(sim.data, file = here("simdata.Rdata")
save(sim.rslt, file = here("simrslt.Rdata")

# check the distribution of simulated MDD
sim.data %>% filter(rate.ratio==0.2 & sample.size==200) %>%
  ggplot(aes(mdd2, fill=factor(rate.ratio*sample.size)))+ 
  geom_histogram(bins=10)



# run simulations for a specific assumption

sim.data.2 = data.frame(sample.size=seq(190, 220, 10)) %>% 
  group_by(sample.size) %>% 
  do(
    foreach(1:nsim, .combine='rbind') %do% simBinomial(null.rate = 0.136,
                                                       rate.ratio = 0.019/0.136, # observed data
                                                       N.total = .$sample.size,
                                                       N.arms = 2, 
                                                       alpha = alpha, 
                                                       dropout = dropout)
  ) 

sim.rslt.2 = sim.data.2 %>% 
  summarise(power = mean(p<0.05)
            , power2=mean(p2<0.05)
            , mdd=median(mdd)
            , mdd2=median(mdd2)
            , .groups="keep")

save(sim.data.2, file = here("simdata2.Rdata")
save(sim.rslt.2, file = here("simrslt2.Rdata")

sim.data.2 %>% filter(sample.size==200) %>%
  ggplot(aes(log(mdd2), fill=factor(sample.size)))+ 
  geom_histogram(bins=10)


sim.data.2 %>% filter(sample.size==200) 
