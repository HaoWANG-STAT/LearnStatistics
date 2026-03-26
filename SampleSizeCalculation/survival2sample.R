## ---- include = FALSE---------------------------------------------------------
now <- as.character(as.POSIXlt(Sys.time()))
today <- as.Date(substr(now, 1, 10))
now <- paste(today, " at ", substr(now, 12, 19), sep = "")
 
library(knitr)

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

## delete cached stuff
folders <- c("cache/html", "files/figure-html")
nam <- "survival2sample_"
for (i in 1:length(folders)){
  all <- list.files(path = paste(path, nam, folders[i], "/", sep = ""))
  if (length(all) > 0){file.remove(paste(path, nam, folders[i], "/", all, sep = ""))}
}


## ---- include = FALSE---------------------------------------------------------
# make the code available for download
purl("survival2sample.Rmd", quiet = TRUE)

# We have manually generated the file survival2sample_simulations.r and ran this file 
# and saved the corresponding outputs. These are input here and the computer-intensive 
# chunks below are set to eval = FALSE. This, to avoid having to run simulations 
# each time this file is compiled.

# load object "power"
load(file = "power_sim.RData")
load(file = "scen3_samplesize.RData")


## ---- include=TRUE, echo=FALSE------------------------------------------------
source("code/pval2sampleT2E.r")
source("code/power2sampleT2E.r")
source("code/milestoneTest.r")
source("code/testMedianTwoGroupsT2E.r")


## ---- include=TRUE, echo=TRUE, message = FALSE--------------------------------
## load packages
packs.html <- c("survival", "rpact", "bpcp", "survRM2", "nph")
for (i in 1:length(packs.html)){library(packs.html[i], character.only = TRUE)}


## ---- include=TRUE, echo=TRUE-------------------------------------------------
# drop-out (2.5% annually)
dropout.year <- 0.025
dropout.month <- -log(1 - dropout.year) / 12
dropout <- rep(dropout.month, 2)

# number of simulation runs
M <- 10 ^ 4

# median in control
med <- 10

# milestone timepoint
t0 <- 20

# ------------------------------------
# scenarios we consider
# ------------------------------------

# both arms exponential, no effect
scen1 <- list(1, "both arms exponential, no effect",          c(0, Inf), rep(med, 2),        
              rep(med, 2),       NA, t0, 50, 800, 380, 
              dropout, M)

# both arms exponential, effect
scen2 <- list(2, "both arms exponential, effect",             c(0, Inf), rep(med / 0.75, 2), 
              rep(med, 2),       NA, t0, 50, 800, 380, 
              dropout, M)

# NPH: piecewise exponential, delayed effect 
scen3 <- list(3, "piecewise exponential, delayed effect", c(0, 6),   c(med, med / 0.6),  
              rep(med, 2),       NA, t0, 50, 800, 440, 
              dropout, M)  

# NPH: piecewise exponential, cure
scen4 <- list(4, "piecewise exponential, cure",               c(0, 18),  c(med, Inf), 
              c(med / 0.7, Inf), NA, t0, 50, 800, 300, 
              dropout, M)  

scenarios <- list(scen1, scen2, scen3, scen4)

# add names
nams.scen <- c("scenario", "description", "piecewiseSurvivalTime", "med1", "med2", "hr", "t0", 
               "permonth", "n", "d", "dropout", "M")
for (i in 1:length(scenarios)){
  scen.i <- scenarios[[i]]
  names(scen.i) <- nams.scen
  
  # fill in derived numbers
  scen.i[["hr"]] <- scen.i[["med2"]] / scen.i[["med1"]]

  # replace
  scenarios[[i]] <- scen.i
}


## ---- echo = TRUE, results = 'asis', message = FALSE, fig.cap = "", fig.align = "center", fig.width = 7, fig.height = 10----
run.scenarios <- 1:4
par(mfrow = c(2, 2), las = 1, mar = c(4.5, 4.5, 3, 1), xaxs = "i", yaxs = "i")

for (i in 1:length(run.scenarios)){
  
  scen.i <- scenarios[[run.scenarios[i]]]
  med1 <- scen.i[["med1"]]
  med2 <- scen.i[["med2"]]

  plot(0, 0, type = "n", xlim = c(0, 70), ylim = c(0, 1), xlab = "time", 
       ylab = "event-free probability", main = paste("Scenario ", scen.i[["scenario"]], ":\n", 
                                                     scen.i[["description"]], sep = ""))
  xs <- seq(0, 500, by = 0.1)
      
  pexp1 <- getPiecewiseExponentialDistribution(xs, piecewiseSurvivalTime = 
                   scen.i[["piecewiseSurvivalTime"]], piecewiseLambda = log(2) / med1)
  pexp2 <- getPiecewiseExponentialDistribution(xs, piecewiseSurvivalTime = 
                   scen.i[["piecewiseSurvivalTime"]], piecewiseLambda = log(2) / med2)
  
  lines(xs, 1 - pexp1, col = 3, lwd = 3)
  lines(xs, 1 - pexp2, col = 2, lwd = 3)
  
  # add milestone and median
  segments(t0, 0, t0, 0.5, lty = 2)
  segments(0, 0.5, t0, 0.5, lty = 2)
  
  if (i == 1){
    text(5, 0.9, 
         paste("horizontal dotted line at median\nvertical dotted line at milestone t = ", 
               t0, sep = ""), adj = 0)
    }
}


## ---- include=TRUE, echo=TRUE, eval = TRUE------------------------------------
# significance level we will use
alpha <- 0.05



## ---- include=TRUE, echo=TRUE, eval = FALSE-----------------------------------
## # names of tests
## nams <- c("logrank", "wlogrank", "maxcombo", "milestone", "median", "RMST")
## 
## # dataframe to collect power
## power <- data.frame(matrix(NA, nrow = length(run.scenarios), ncol = length(nams)))
## colnames(power) <- nams
## 
## set.seed(2022)
## 
## for (i in 1:length(run.scenarios)){
## 
##   scen.i <- scenarios[[run.scenarios[i]]]
##   M <- scen.i[["M"]]
##   med1 <- scen.i[["med1"]]
##   med2 <- scen.i[["med2"]]
##   d <- scen.i[["d"]]
##   dropout <- scen.i[["dropout"]]
##   t0 <- scen.i[["t0"]]
##   permonth <- scen.i[["permonth"]]
##   n <- scen.i[["n"]]
## 
##   # set up matrices to collect results
##   pvals <- matrix(NA, nrow = M, ncol = length(nams))
##   colnames(pvals) <- nams
## 
##   # ----------------------------------
##   # simulate M clinical trials using rpact
##   # ----------------------------------
##   trial1 <- getSimulationSurvival(lambda1 = log(2) / med1, lambda2 = log(2) / med2,
##                                   piecewiseSurvivalTime = scen.i[["piecewiseSurvivalTime"]],
##                                   dropoutRate1 = dropout[1], dropoutRate2 = dropout[2],
##                                   dropoutTime = 12, accrualTime = c(0, n / permonth),
##                                   accrualIntensity = permonth, plannedEvents = d,
##                                   maxNumberOfIterations = M,
##                                   maxNumberOfRawDatasetsPerStage = M)
##   trial2 <- getRawData(trial1)
## 
##   for (j in 1:M){
## 
##     time <- trial2[trial2$iterationNumber == j, "timeUnderObservation"]
##     event <- as.numeric(trial2[trial2$iterationNumber == j, "event"])
##     arm <- as.numeric(trial2[trial2$iterationNumber == j, "treatmentGroup"] == 1)
##     # 0 = control, 1 = treatment
## 
##     # ----------------------------------
##     # remove patients that arrived after cutoff
##     # ----------------------------------
##     rem <- (time >= 0)
##     time <- time[rem]
##     event <- event[rem]
##     arm <- arm[rem]
## 
##     # ----------------------------------
##     # compute all tests
##     # ----------------------------------
##     pvals[j, ] <- pval2sampleT2E(time = time, event = event, arm = arm, t0 = t0,
##                                  wrho = 0, wgamma = 1, rho = c(0, 0, 1), gamma = c(0, 1, 0))
##     }
##   rownames(power)[i] <- scen.i[["description"]]
##   power[i, ] <- apply(pvals <= alpha, 2, mean)
## }


## ---- include=TRUE, echo=FALSE------------------------------------------------
knitr::kable(power)


## ---- include=TRUE, echo=TRUE-------------------------------------------------
# compute power Scenario 3 and a different number of events
s <- scenarios[[3]]

# use larger number of events
d3 <- 457


## ---- include=TRUE, echo=TRUE, eval = FALSE-----------------------------------
## # simulate power
## scen3_samplesize <- power2sampleT2E(piecewiseSurvivalTime = s$piecewiseSurvivalTime,
##                                     med1 = s$med1, med2 = s$med2, t0 = s$t0,
##                                     accrualTime = c(0, s$n / s$permonth),
##                                     accrualIntensity = s$permonth,
##                                     dropout = s$dropout, d = d3, alpha = 0.05, M = s$M,
##                                     print.progress = s$M + 1, seed = 2022, wrho = 0, wgamma = 1,
##                                     rho = c(0, 0, 1), gamma = c(0, 1, 0))
## scen3_samplesize


## ---- include=TRUE, echo=TRUE-------------------------------------------------
library(km.ci)
library(survival)

# milestone
t0 <- 24

# load dataset for illustration
data(rectum.dat)

# compute Kaplan-Meier estimator
fit <- survfit(Surv(time, status) ~ 1, data = rectum.dat)

# log-log transformed confidence intervals
fit2 <- km.ci(fit, method = "loglog")

# log-log confidence interval at t0 = 24
summary(fit2, times = t0)

# bpcp confidence interval
library(bpcp)
fit3 <- bpcp(time = rectum.dat$time, status = rectum.dat$status)
StCI(fit3, tstar = t0)


## ---- echo = TRUE, fig.width=6, fig.height=5 ,results = "markup"--------------
library(ComparisonSurv)
data(leuk2)
group <- factor(leuk2$treatment, levels = c("placebo", "6-MP"))
t0 <- 12

fixtdiff(time = leuk2$time, status = leuk2$status, group = group, testtime = t0, trans = c("cloglog"), 
         varpooled = TRUE, correct = FALSE, doall = TRUE)

Fixpoint.test(time = leuk2$time, status = leuk2$status, 
              group = ifelse(leuk2$treatment == "placebo", 0, 1), t0 = t0)


## ---- include=TRUE, echo=TRUE-------------------------------------------------
bpcp2samp(time = leuk2$time, status = leuk2$status, group = group, 
          testtime = t0, parmtype = "difference")


## ---- include = FALSE---------------------------------------------------------
## this chunk is here to delete large cached objects (likely some simulation stuff that is stored) 
## that are not needed for the html file to render properly, but which cause github 
## to reject the commit
path      <- ""

## delete cached stuff
folders <- c("cache/html")
nam <- "survival2sample_"
for (i in 1:length(folders)){
  
  all <- list.files(path = paste(path, nam, folders[i], "/", sep = ""))
  if (length(all) > 0){file.remove(paste(path, nam, folders[i], "/", all, sep = ""))}
}

