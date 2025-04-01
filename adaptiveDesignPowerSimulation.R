library(simstudy)
library(data.table)
library(tidyverse)


s_define <- function() {
  
  #--- data definition code ---#
  
  # define the outcome, should be updated to match study design 
  # ..delta[1] and ..delta[2] are placeholders for treatment effects of Low (T=2) and High (T=3) dose group
  def1 <- defDataAdd(varname = "Y", dist = "normal",
                     formula = "5 + 0*period + ..delta[1]*(Trt==2) + ..delta[2]*(Trt==3) + 0*period*Trt + e")  
  
  # define the correlated errors
  mu <- rep(0, 3)
  sigma <- rep(sqrt(3), 3)
  
  # define the treatment effect
  delta <- c(0, 0)
#  n <- 3*100
#  rho <- 0.5
  
  return(list(def1 = def1, mu = mu, sigma = sigma, delta = delta))
}

s_generate <- function(list_of_defs, argsvec) {
  
  list2env(list_of_defs, envir = environment())
  list2env(as.list(argsvec), envir = environment())
  
  #--- data generation code ---#
  
  # generate correlated data for each id and assign treatment
  
  # parameters that can be defined through argsvec for simulating different scenarios
  dtCor <- genCorData(n = n, mu = mu, sigma = sigma, rho = rho, corstr = "cs") 
    
  dtCor <- trtAssign(dtCor, nTrt = 3, balanced = TRUE, grpName = "Trt")
  
  # create longitudinal data set and generate outcome based on definition
  
  longData <- addPeriods(dtCor, nPeriods = 3, idvars = "id", timevars = c("V1","V2", "V3"), timevarName = "e")
  longData <- addColumns(def1, longData)
  
  # look at the data via plot
  
#  ggplot(data = longData, aes(x = factor(period), y = Y)) +
#    geom_line(aes(color = T, group = id)) +
#    scale_color_manual(values = c("#e38e17", "#8e17e3", "#e40e20")) +
#    xlab("Time")

  ds <- longData |> select(-e, -timeID) |> 
    pivot_wider(names_from = period, values_from = Y, names_prefix = "Y") |> 
    mutate(CFB1 = Y1-Y0, CFB2 = Y2-Y0, Trt=as.factor(Trt))
  
  return(ds)
}

s_model <- function(generated_data) {
  
  #--- model code ---#
  
  # treatment selection at interim analysis
  lmfit <- lm(CFB1 ~ Trt, data = generated_data |> filter(id <= n*frac))
  est <- summary(lmfit)$coef[-1, "Estimate"]
  select <- as.numeric(gsub("Trt", "", names(which.max(est)))) # return the selected active group id
  #pval <- summary(lmefit)$coef[2, "Pr(>|t|)"]
  
  
  # primary analysis at the end of trial
  lmfit2 <- lm(CFB2 ~ Trt, data = generated_data |> filter(Trt == 1 | Trt == select) )
  est2 <- summary(lmfit2)$coef[-1, "Estimate"]
  pval <- summary(lmfit2)$coef[2, "Pr(>|t|)"]
  
  return(data.table(est, select, est2, pval)) # model_results is a data.table
}

s_single_rep <- function(list_of_defs, argsvec) {
  
  generated_data <- s_generate(list_of_defs, argsvec)
  model_results <- s_model(generated_data)
  
  return(model_results)
}

s_replicate <- function(argsvec, nsim) {
  
  list_of_defs <- s_define()
  
  model_results <- rbindlist(
    parallel::mclapply(
      X = 1 : nsim,
      FUN = function(x) s_single_rep(list_of_defs, argsvec),
      mc.cores = 4)
  )
  
  #--- summary statistics ---#
  
  power <- model_results[, mean(pval <= 0.05)]
  prob_select <- model_results[, mean(select == 2)]
  summary_stats <- data.table(t(argsvec), prob_select, power) |> 
    mutate(ntot = n - n*frac/3) # ntot is the actual sample size at the end of trial
  
  return(summary_stats) # summary_stats is a data.table
}


# scenarios to simulate 
scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}

n <- 3*c(50, 100)
frac <- c(0.3, 0.5)
rho <- c(0.7, 0.5)

# n is the pseudo sample size when there is no treatment selection
# frac is the fraction of sample size at interim analysis

scenarios <- scenario_list(n = n, frac = frac, rho = rho)

#--- run locally ---#

summary_stats <- rbindlist(lapply(scenarios, function(a) s_replicate(a, nsim = 10000)))
summary_stats
