## ---- include = FALSE---------------------------------------------------------
now <- as.character(as.POSIXlt(Sys.time()))
today <- as.Date(substr(now, 1, 10))
now <- paste(today, " at ", substr(now, 12, 19), sep = "")

library(knitr)
 
# bpcp: package for Brookmeyer-Crowley confidence interval for median
# rpact: use for simulation of single-arm trial

## set some knitr options
knitr::opts_chunk$set(echo = TRUE)

## read in functions
colFmt <- function(x, color = 'red'){

  require(knitr)
  
  # http://stackoverflow.com/questions/29067541/rmarkdown-how-to-change-the-font-color
  x <- as.character(x)
  outputFormat <- opts_knit$get("rmarkdown.pandoc.to")
  if(outputFormat == 'latex')
    paste("\\textcolor{", color, "}{", x, "}", sep = "")
  else if(outputFormat == 'html')
    paste("<font color='", color, "'>", x, "</font>", sep = "")
  else
    x
}

## specify paths
path      <- ""

## delete cached stuff --> otherwise, simulations are often not re-run for a different M
folders <- c("cache/html", "files/figure-html")
nam <- "survival1sample_"
for (i in 1:length(folders)){
  all <- list.files(path = paste(path, nam, folders[i], "/", sep = ""))
  if (length(all) > 0){file.remove(paste(path, nam, folders[i], "/", all, sep = ""))}
}


## ---- include = FALSE---------------------------------------------------------
# make the code available for download
purl("survival1sample.Rmd", quiet = TRUE)


## ---- include=TRUE, echo=TRUE, message = FALSE--------------------------------
## load packages
packs.html <- c("knitr", "survival", "rpact", "bpcp", "reporttools")
for (i in 1:length(packs.html)){library(packs.html[i], character.only = TRUE)}


## ---- include=TRUE, echo=TRUE-------------------------------------------------
# assumed effects
lambda0 <- 1       # null hypothesis for hazard  ratio
lambda1 <- 0.75    # planning alternative hypothesis for sample size computation
R <- 0.75          # true hazard ratio

# historical control rate and median
m_ctl <- 5
lam_ctl <- log(2) / m_ctl

# experimental rate and median
lam_exp <- lam_ctl * R
m_exp <- log(2) / lam_exp


## ---- echo = TRUE, results = 'asis', message = FALSE, fig.cap = "", fig.align = "center", fig.width = 7, fig.height = 5.5----
par(las = 1, mar = c(4.5, 4.5, 3, 1))
xs <- seq(0, 100, by = 0.01)
S_ctl <- 1 - pweibull(xs, shape = 1, scale = 1 / lam_ctl)
S_exp <- 1 - pweibull(xs, shape = 1, scale = 1 / lam_exp)

plot(0, 0, type = "n", xlim = c(0, 30), ylim = c(0, 1), main = "assumed survival functions", 
     xlab = "time", ylab = "even-free probability")
lines(xs, S_ctl, col = 1, lwd = 3)
lines(xs, S_exp, col = 2, lwd = 3)


## ---- include=TRUE, echo=TRUE-------------------------------------------------
alpha <- 0.025
beta <- 0.2
den <- qnorm(1 - alpha) + sqrt(lambda1) * qnorm(1 - beta)
nom <- 1 - lambda1
fink <- ceiling(lambda1 * (den / nom) ^ 2)
fink


## ---- include=TRUE, echo=FALSE------------------------------------------------
source("code/power1sampleT2E.r")


## ---- include=TRUE, echo=TRUE-------------------------------------------------
# n: number of patients, as in Table above
n <- 124

# recruitment per month
permonth <- 15

# define objects that can be submitted to rpact's "getSimulationSurvival" function
# within power1sampleT2E
accrualTime <- c(0, n / permonth)
accrualIntensity <- permonth

# apply drop-out (2.5% annually)
dropout.year <- 0.025
dropout.month <- -log(1 - dropout.year) / 12

# milestone time
t0 <- 6

# number of simulations
M <- 10 ^ 4


## ---- include=TRUE, echo=TRUE-------------------------------------------------
# compute power for all three approaches
# assess T1E --> chose data generating mechanism for experimental same as historical control
power1sampleT2E(med1 = m_ctl, med0 = m_ctl, kappa = 1, t0 = t0, accrualTime = accrualTime, 
                accrualIntensity = accrualIntensity, dropout = dropout.month, d = fink, 
                alpha = 0.05, lambda0 = 1, M = M, print.progress = M + 1, seed = 2022)


## ---- include=TRUE, echo=TRUE-------------------------------------------------
# compute power for all three approaches
power1sampleT2E(med1 = m_exp, med0 = m_ctl, kappa = 1, t0 = t0, accrualTime = accrualTime, 
                accrualIntensity = accrualIntensity, dropout = dropout.month, d = fink, 
                alpha = 0.05, lambda0 = 1, M = M, print.progress = M + 1, seed = 2022)


## ---- include=TRUE, echo=TRUE-------------------------------------------------
# parameter theta = 1 / rate
theta <- m_exp / log(2)

# maximal width of confidence interval, same scale as m (e.g. months)
delta <- 4   # so we want confidence interval [8, 12]

# confidence level 
conf.level <- 0.95

# compute necessary number of events
r <- 4 * (log(2) * qnorm((1 + conf.level) / 2) * theta / delta) ^ 2
r


## ---- include=TRUE, echo=FALSE------------------------------------------------
source("code/width1sampleT2E.r")


## ---- include=TRUE, echo=TRUE-------------------------------------------------
# compute width of confidence intervals
conf.level <- 0.95
d <- 80
out <- width1sampleT2E(med1 = m_exp, kappa = 1, t0 = t0, accrualTime = accrualTime, 
                accrualIntensity = accrualIntensity, dropout = dropout.month, d = d, 
                conf.level = conf.level, M = M, print.progress = M + 1, seed = 2022)

# the resulting object contains all CIs for median (nonparametric and exponential) and milestone, 
# as well as their widths
# we compute the median width of all M CIs
ws <- apply(out$width, 2, median)
ws


## ---- include=TRUE, echo=TRUE-------------------------------------------------
# recruitment interval 
t1 <- 12

# period of follow-up
t2 <- 12

# rate computed from theta above
lambda <- 1 / theta

# probability of an event
psi <- 1 - (exp(- lambda * t2) - exp(- lambda * (t1 + t2))) / (lambda * t1)

# number of patients
n <- r / psi
n


## ---- include = FALSE---------------------------------------------------------
## this chunk is here to delete large cached objects (likely some simulation stuff that is stored) 
## that are not needed for the html file to render properly, but which cause github 
## to reject the commit
path      <- ""

## delete cached stuff
folders <- c("cache/html")
nam <- "survival1sample_"
for (i in 1:length(folders)){
  
  all <- list.files(path = paste(path, nam, folders[i], "/", sep = ""))
  if (length(all) > 0){file.remove(paste(path, nam, folders[i], "/", all, sep = ""))}
}

