sim<-50000

# brain mets
# use 7 out of 11 as prior
p <- rbeta(sim, 7, 4)
N <- 55

# binom.test(30,55, p = 0.4, alternative = "two.sided") 

r.success <- 30 # 30 out of 55's 95 CI lower bound is gt 0.4

r.n <- rep(NA, sim)

for (i in 1:sim)
{
  r.n[i] <- rbinom(1, N, p[i])
}
cpts.bm1 <- sum(r.n>=r.success)/sim

# use 13 out of 21 as prior
p <- rbeta(sim, 13, 8)
N <- 55
r.success <- 30 # 30 out of 55's 95 CI lower bound is 0.4
r.n <- rep(NA, sim)
for (i in 1:sim)
{
  r.n[i] <- rbinom(1, N, p[i])
}
cpts.bm2 <- sum(r.n>=r.success)/sim

# allcomer (brain mets + no brain mets)
# use 13 out of 21 as prior
p <- rbeta(sim, 13, 8)
N <- 145
# binom.test(82,145, p = 0.48, alternative = "two.sided") 
r.success <- 82 # 82 out of 145's 95 CI lower bound is gt 0.48
r.n <- rep(NA, sim)
for (i in 1:sim)
{
  r.n[i] <- rbinom(1, N, p[i])
}
cpts.allcomer <- sum(r.n>=r.success)/sim

c(cpts.bm1, cpts.bm2, cpts.allcomer)



# ROS1+2L
# use 14 out of 40 as prior
p <- rbeta(sim, 14, 26)
N <- 66

# binom.test(17,66, p = 0.15, alternative = "two.sided") 

r.success <- 17 # 17 out of 66's 95 CI lower bound is gt 0.15

r.n <- rep(NA, sim)

for (i in 1:sim)
{
  r.n[i] <- rbinom(1, N, p[i])
}
cpts.2L <- sum(r.n>=r.success)/sim
cpts.2L
