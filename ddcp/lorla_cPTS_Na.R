library(LearnBayes)
library(extraDistr)
#library(Bolstad)

# Success Criteria
# at least 16 responders need to be observed out of 66 to meet the success criteria: lower of 95% CI excluding 15%
binom.test(17,66, p = 0.15, alternative = "two.sided") 
sum(dbinom(17:66,66,0.3))

# prior: observed 14 responders in 40 pts 
N <- 66
R.Target <- 17
N.observed <- 40  #Number of patients from observed data
R.observed <- 14   #Number of responders from observed data
pi=binom.test(R.observed, N.observed, p = 0.3, alternative = "two.sided")
cil=pi$conf.int[1]
cih=pi$conf.int[2]

quantile1=list(p=0.025, x=cil)
quantile2=list(p=0.975, x=cih)
priorbeta <- beta.select(quantile1,quantile2)
priorbeta  

prob1 <- pbbinom(R.Target-1,N,alpha=priorbeta[1], beta=priorbeta[2],lower.tail=F,log.p=F)

prob1



# Success Criteria
# at least 30 responders need to be observed out of 55 to meet the success criteria: lower of 95% CI excluding 40%
binom.test(30,55, p = 0.4, alternative = "two.sided") 
sum(dbinom(30:55,55,0.6))

# prior: observed 7 responders in 11 pts 
N.observed <- 11  #Number of patients from observed data
R.observed <- 7   #Number of responders from observed data
pi=binom.test(R.observed, N.observed, p = 0.4, alternative = "two.sided")
cil=pi$conf.int[1]
cih=pi$conf.int[2]

quantile1=list(p=0.025, x=cil)
quantile2=list(p=0.975, x=cih)
priorbeta <- beta.select(quantile1,quantile2)

N <- 55
R.Target <- 30
prob2 <- pbbinom(R.Target-1,N,alpha=priorbeta[1], beta=priorbeta[2],lower.tail=F,log.p=F)

prob2


# prior: observed 13 responders in 21 pts 
N.observed <- 21  #Number of patients from observed data
R.observed <- 13  #Number of responders from observed data
pi=binom.test(R.observed, N.observed, p = 0.5, alternative = "two.sided")
cil=pi$conf.int[1]
cih=pi$conf.int[2]

quantile1=list(p=0.025, x=cil)
quantile2=list(p=0.975, x=cih)
priorbeta <- beta.select(quantile1,quantile2)

priorbeta[1]/(priorbeta[1]+priorbeta[2])

N <- 145
R.Target <- 82
prob3 <- pbbinom(R.Target-1,N,alpha=priorbeta[1], beta=priorbeta[2],lower.tail=F,log.p=F)

prob3



#calculate analytical solutions using beta-binomial distribution
pbbinom(16,66,alpha=14, beta=26,lower.tail=F,log.p=F)

pbbinom(29,55,alpha=7, beta=4,lower.tail=F,log.p=F)
# sum(dbbinom(30:55,55,alpha=13, beta=8))

pbbinom(81,145,alpha=13, beta=8,lower.tail=F,log.p=F)


