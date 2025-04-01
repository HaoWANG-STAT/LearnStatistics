library(simstudy)
library(data.table)
library(tidyverse)


s_define <- function() {
  
  #--- data definition code ---#
  
  # 1. define the outcome, should be updated to match study design
  # mean baseline is 4 and mean change from baseline for control and active groups are defined by a vector
  def1 <- defDataAdd(varname = "Y", dist = "normal",
                     formula = "4 + log(1 - ..reductionCon)*(TRT==0&period==1) + log(1 - ..reductionAct)*(TRT==1&period==1) + e"
                     #formula = "4 + log(1-0.1)*(TRT==0&period==1) + log(1-0.4)*(TRT==1&period==1) + e"
                     )
  
  # 2. define the stratification factors
  # probability of EOS being High is 0.5
  def2 <- defDataAdd(varname = "EOS", dist = "binary", formula = 0.5)
  # probability of NUT being High is 0.4; and prob of NUT being High given EOS High is 0.25
  # beta=log(OR), alpha=log(Odds|EOS=Low)
  def2 <- defDataAdd(def2, varname = "NUT", dist = "binary",
                 formula = "log(1) + EOS * log(1/3)", link = "logit")
  
  # 3. define the correlated errors
  mu <- rep(0, 2)
  sigma <- rep(0.96, 2)
  #reduction <- c(0.1, 0.6)
  
  return(list(def1 = def1, def2 = def2
              , mu = mu
              , sigma = sigma
              #, reduction = reduction
              )
         )
}

s_generate <- function(list_of_defs, argsvec) {
  
  
  list2env(list_of_defs, envir = environment())
  list2env(as.list(argsvec), envir = environment())
  
  #--- data generation code ---#
  
  # generate correlated data for each id and assign treatment
  
  # parameters that can be defined through argsvec for simulating different scenarios
  dtCor <- genCorData(n = n, mu = mu, sigma = sigma, rho = rho, corstr = "cs") 
  dtCor <- addColumns(def2, dtCor)
  # dtCor |> group_by(EOS, NUT) |> summarise(n=n())
  
  dtCor <- trtAssign(dtCor, nTrt = 2, balanced = TRUE, strata = c("EOS", "NUT"), grpName = "TRT")
  # dtCor |> group_by(EOS, NUT, TRT) |> summarise(n=n())
  
  # create longitudinal data set and generate outcome based on definition
  
  longData <- addPeriods(dtCor, nPeriods = 2, idvars = "id", timevars = c("V1","V2"), timevarName = "e") 
  
  longData <- addColumns(def1, longData)
  
  # look at the data via plot
  
  #  ggplot(data = longData, aes(x = factor(period), y = Y)) +
  #    geom_line(aes(color = T, group = id)) +
  #    scale_color_manual(values = c("#e38e17", "#8e17e3", "#e40e20")) +
  #    xlab("Time")
  
  ds <- longData |> select(-e, -timeID) |> 
    pivot_wider(names_from = period, values_from = Y, names_prefix = "Y") |> 
    mutate(CFB = Y1-Y0, TRT=as.factor(TRT))
  
  return(ds)
}

s_model <- function(generated_data = ds) {
  
  #--- model code ---#
  results <- 
  generated_data |> 
      bind_rows(generated_data |> mutate(EOS = -9, NUT = -9)) |>
    group_by(EOS, TRT) |> summarise(MeanBL = mean(Y0),
                               MeanEP = mean(Y1),
                               #MeanCFB = mean(CFB),
                               GMR = exp(MeanEP - MeanBL),
                               PctCFB = (1 - exp(MeanEP - MeanBL)) * 100 
    ) |> ungroup() |> 
      pivot_wider(names_from = TRT, values_from = c(MeanBL, MeanEP, GMR, PctCFB)) |> 
      mutate(Effect = PctCFB_1 - PctCFB_0,
             Outcome = ifelse(Effect > 50, "Go", ifelse(Effect < 20, "NoGo", "Eval")),
             Population = ifelse(EOS == -9, "Allcomer", "Subgroup")
             ) |> 
      select(-EOS) |>
    
   # repeat for second subgroup
   bind_rows(
     generated_data |> 
       #bind_rows(ds |> mutate(EOS = 999, NUT = 999)) |>
       group_by(NUT, TRT) |> summarise(MeanBL = mean(Y0),
                                       MeanEP = mean(Y1),
                                       #MeanCFB = mean(CFB),
                                       GMR = exp(MeanEP - MeanBL),
                                       PctCFB = (1 - exp(MeanEP - MeanBL)) * 100 
       ) |> ungroup() |> 
       pivot_wider(names_from = TRT, values_from = c(MeanBL, MeanEP, GMR, PctCFB)) |> 
       mutate(Effect = PctCFB_1 - PctCFB_0,
              Outcome = ifelse(Effect > 50, "Go", ifelse(Effect < 20, "NoGo", "Eval")),
              Population = "Subgroup"
       ) |> 
       select(-NUT)
   ) 
  
 outcome <- (results |> filter(Population=="Allcomer") |> select(Outcome) |> rename(Outcome_AC=Outcome)) |>
    bind_cols(results |> filter(Population=="Subgroup") |>  
                summarise(Outcome_Sub = sum(Outcome=="Go")) |> 
                mutate(Outcome_Sub = ifelse(Outcome_Sub > 0, "Go", "NoGo"))
    ) |> mutate(Outcome_Both = ifelse(Outcome_AC == "Go" | (Outcome_AC == "Eval" & Outcome_Sub =="Go"),  "Go", "NoGo"))
    
  return(outcome) # model_results is a data.table
}

s_single_rep <- function(list_of_defs, argsvec) {
  
  list_of_defs <- s_define()
  
  generated_data <- s_generate(list_of_defs, argsvec)
  model_results <- s_model(generated_data)
  
  return(model_results)
}

s_replicate <- function(argsvec, nsim) {
  
  list_of_defs <- s_define()
  
  coresAvailable <- parallel::detectCores()
  
  model_results <- rbindlist(
    parallel::mclapply(
      X = 1 : nsim,
      FUN = function(x) s_single_rep(list_of_defs, argsvec),
      mc.cores = coresAvailable)
  )
  
  #--- summary statistics ---#
  
  Prob1 <- model_results[, mean(Outcome_AC == "Go")] # Go if based on Allcomer data only
  Prob2 <- model_results[, mean(Outcome_Both == "Go")] # Go if based on Allcomer and Subgroup data
  #prob_select <- model_results[, mean(select == 2)]
  summary_stats <- data.table(t(argsvec), Prob1, Prob2) 
  #  mutate(ntot = n - n*frac/3) # ntot is the actual sample size at the end of trial
  
  return(summary_stats) # summary_stats is a data.table
}


# scenarios to simulate 
scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}

n <- 2*c(30)
rho <- c(0.7)
reductionCon <- c(0, 0.1)
reductionAct <- c(0.3, 0.4, 0.5)

scenarios <- scenario_list(n = n, rho = rho, reductionCon=reductionCon, reductionAct=reductionAct)

#--- run locally ---#
summary_stats <- rbindlist(lapply(scenarios, function(a) s_replicate(argsvec = a, nsim = 1000)))
summary_stats



#--- test functions ---#
list_of_defs <- s_define()
argsvec <- scenarios[[1]]
generated_data <- s_generate(list_of_defs, argsvec)
s_single_rep(list_of_defs, argsvec)
s_replicate(argsvec, nsim = 10)



