library(data.table)
library(tidyverse)
library(simstudy)



s_define <- function() {
  
  #--- data definition code ---#
  
  # 1. define the outcome
  # mean baseline is 4 and mean reduction from baseline for control and test groups are reduCtrl and reduTest
  def1 <- defDataAdd(varname = "Y", dist = "normal",
                     formula = "4 + log(1 - ..reduCtrl)*(TRT==0&period==1) + log(1 - ..reduTest)*(TRT==1&period==1) + e"
                     #formula = "4 + log(1-0.1)*(TRT==0&period==1) + log(1-0.4)*(TRT==1&period==1) + e"
                     )
  
  # 2. define the stratification factors
  # probability of EOS being High is 0.5
  def2 <- defDataAdd(varname = "EOS", dist = "binary", formula = 0.5)
  # probability of NEU being High is 0.4; and prob of NEU being High given EOS High is 0.25
  # beta=log(OR)=log(1/3), alpha=log(Odds|EOS=Low)=log(1)
  def2 <- defDataAdd(def2, varname = "NEU", dist = "binary",
                 formula = "log(1) + EOS * log(1/3)", link = "logit")
  
  return(list(def1 = def1, 
              def2 = def2)
         )
}

s_generate <- function(list_of_defs, argsvec) {
  
  
  list2env(list_of_defs, envir = environment())
  list2env(as.list(argsvec), envir = environment())
  
  #--- data generation code ---#
  
  # generate correlated data for each id and assign treatment
  
  # parameters that can be defined through argsvec for simulating different scenarios
  dtCor <- genCorData(n = n, mu = rep(0, 2), sigma = rep(sigma, 2), rho = rho, corstr = "cs") 
  dtCor <- addColumns(def2, dtCor)
  # dtCor |> group_by(EOS, NEU) |> summarise(n=n())
  
  dtCor <- trtAssign(dtCor, nTrt = 2, balanced = TRUE, strata = c("EOS", "NEU"), grpName = "TRT")
  # dtCor |> group_by(EOS, NEU, TRT) |> summarise(n=n())
  
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

s_model <- function(generated_data, thresholdSub, effType) {
  
  #--- model code ---#
  results <- 
  generated_data |> 
    bind_rows(generated_data |> mutate(EOS = -9, NEU = -9)) |>
    group_by(EOS, TRT) |> 
    summarise(MeanBL = mean(Y0),
              MeanEP = mean(Y1),
              MeanCFB = mean(CFB),
              GMR = exp(MeanEP - MeanBL),
              PctCFB = (1 - exp(MeanEP - MeanBL)) * 100,
              N = n(),
              .groups = 'keep'
              ) |>
    ungroup() |>
    pivot_wider(names_from = TRT, values_from = c(MeanBL, MeanEP, MeanCFB, GMR, PctCFB, N)) |> 
    mutate(Effect = PctCFB_1 - PctCFB_0, # effect defined as difference in change from baseline (%) between Selno and Pbo group
           Effect1 = (1 - exp(MeanCFB_1 - MeanCFB_0)) *100, # effect1 defined as 1 - GMR of change from baseline between Selno and Pbo group
           Population = case_when(EOS == -9 ~ "Allcomer",
                                  EOS == 1 ~ "EOSH",
                                  TRUE ~ "EOSL")) |> 
    select(-EOS) |>
    
   # for second subgroup
   bind_rows(
     generated_data |> 
       group_by(NEU, TRT) |> 
       summarise(MeanBL = mean(Y0),
                 MeanEP = mean(Y1),
                 MeanCFB = mean(CFB),
                 GMR = exp(MeanEP - MeanBL),
                 PctCFB = (1 - exp(MeanEP - MeanBL)) * 100,
                 N = n(),
                 .groups = 'keep'
       ) |> ungroup() |>
       pivot_wider(names_from = TRT, values_from = c(MeanBL, MeanEP, MeanCFB, GMR, PctCFB, N)) |> 
       mutate(Effect = PctCFB_1 - PctCFB_0, # effect defined as difference in change from baseline (%) between Selno and Pbo group
              Effect1 = (1 - exp(MeanCFB_1 - MeanCFB_0)) *100, # effect1 defined as GMR of change from baseline between Selno and Pbo group
              Population = case_when(NEU == -9 ~ "Allcomer",
                                     NEU == 1 ~ "NEUH",
                                     TRUE ~ "NEUL")) |>  
       select(-NEU)
   )

  # select which type of effect to use  
  if (effType == "relative") {
    results <- results |> mutate(Effect = Effect1)
  } 
    
#-------------------------------------------------------------------------------------------------------------#  
# decision <- results |> filter(Population == "Allcomer") |> select(Outcome) |> rename(Outcome_AC=Outcome) |>
#   bind_cols(
#     results |> filter(Population != "Allcomer") |> 
#       summarise(Outcome_Sub = sum(Outcome=="Go")) |> 
#       mutate(Outcome_Sub = ifelse(Outcome_Sub > 0, "Go", "NoGo"))
#    ) |> 
#   mutate(Outcome_Both = ifelse(Outcome_AC == "Go" | (Outcome_AC == "Eval" & Outcome_Sub =="Go"),  "Go", "NoGo"))
 
 
# decision <- results |> select(Population, Outcome) |> 
#   pivot_wider(names_from = Population, values_from = Outcome)
#    
#  return(list(results = results,
#              decision = decision)
#  )
 #----------------------------------------------------------------------------------------------------------------# 
 decision <- purrr::map_dfr(seq(length(thresholdSub)), ~ results) |> 
   mutate(Cutoff = rep(thresholdSub, each=nrow(results))) |> 
   
   mutate(Outcome = case_when(
     (Population == "Allcomer" & Effect > 50) | (Population != "Subgroup" & Effect > Cutoff) ~ "Go",
     Population == "Allcomer" & Effect <= 50 & Effect >= 20 ~ "Eval",
     TRUE ~ "NoGo")) |> 
    select(Cutoff, Population, Outcome) |> 
    pivot_wider(names_from = Population, values_from = Outcome)
 
  return(list(results = results, decision = decision))
 
}

# added a parameter thresholdSub for evaluating different threshold cutoff values
s_single_rep <- function(list_of_defs, argsvec, thresholdSub, effType) {
  
  generated_data <- s_generate(list_of_defs, argsvec)
  model_results <- s_model(generated_data, thresholdSub = thresholdSub, effType = effType)[["decision"]] # extract the decision part
  
  return(model_results)
}

s_replicate <- function(argsvec, thresholdSub, effType, nsim, seed=12345){
  
  set.seed(seed) # set main seed

  seeds <- sample.int(.Machine$integer.max, size = nsim)  # generate child seeds
  
  list_of_defs <- s_define()
  
  coresAvailable <- parallel::detectCores()
  
  model_results <- rbindlist(
    parallel::mclapply(
      X = 1 : nsim,
      FUN = function(x) {
        set.seed(seeds[x])
        s_single_rep(list_of_defs, argsvec, thresholdSub, effType)
        },
      mc.cores = coresAvailable)
  )
  
  #--- summary statistics ---#
  
  # currently the decision is only based on Allcomer and EOS low and high subgroup
  dProb <- model_results |> group_by(Cutoff) |> 
    summarise(
      ProbEval = mean(Allcomer == "Eval")*100,
      ProbGoAc = mean(Allcomer == "Go")*100,
      
      ProbGoAcEosl = mean(Allcomer == "Eval" & EOSL == "Go")*100,
      ProbGoAcEosh = mean(Allcomer == "Eval" & EOSH == "Go")*100,
      ProbGoAcEos = mean(Allcomer == "Go" | 
                         (Allcomer == "Eval" & (EOSL == "Go" | EOSH == "Go"))) * 100,
      
      ProbGoAcNeul = mean(Allcomer == "Eval" & NEUL == "Go")*100,
      ProbGoAcNeuh = mean(Allcomer == "Eval" & NEUH == "Go")*100,
      ProbGoAcNeu = mean(Allcomer == "Go" | 
                         (Allcomer == "Eval" & (NEUL == "Go" | NEUH == "Go"))) * 100,

      ProbGoAcEosNeul= mean(Allcomer == "Go" | 
                            (Allcomer == "Eval" & (EOSL == "Go" | EOSH == "Go")) |
                            (Allcomer == "Eval" & NEUL == "Go")) * 100,      
      ProbGoAcEosNeu= mean(Allcomer == "Go" | 
                           (Allcomer == "Eval" & (EOSL == "Go" | EOSH == "Go")) |
                           (Allcomer == "Eval" & (NEUL == "Go" | NEUH == "Go"))) * 100
    )
  
  summary_stats <- data.table(t(argsvec), dProb) 
  return(summary_stats)
}


#--- scenarios to simulate ---#

# n <- c(60)
# sigma <- 0.96
# rho <- c(0.7)
# reduCtrl <- c(0)
# reduTest <- seq(0.2, 0.5, 0.05)
# 
# scenarios <- expand.grid(n = n, sigma=sigma, rho = rho, reduCtrl=reduCtrl, reduTest=reduTest) |> 
#   asplit(MARGIN = 1)
# 
# if (file.exists(dataFile)) {
#   summary_stats <- readRDS(here::here("data", "FalseGoSubgroupSimulation.rds"))
#   message(paste0("Data read from ", here::here("data", "FalseGoSubgroupSimulation.rds"), "."))
# } else {
# summary_stats <- rbindlist(lapply(scenarios, function(a) s_replicate(argsvec = a, thresholdSub = thresholdSub, nsim = 10000, seed = 12345)))
# saveRDS(summary_stats, here::here("data", "FalseGoSubgroupSimulation.rds"))
# }



#--- test functions ---#
# argsvec <- scenarios[[1]]
# list_of_defs <- s_define()
# generated_data <- s_generate(list_of_defs, argsvec)
# s_model(generated_data)
# s_single_rep(list_of_defs, argsvec)
# s_replicate(argsvec, thresholdSub= c(50, 60), nsim = 10)
