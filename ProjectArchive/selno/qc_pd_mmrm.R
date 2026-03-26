install.packages(
  "mmrm",
  repos = c("https://openpharma.r-universe.dev", "https://cloud.r-project.org")
)
library(mmrm)
library(emmeans)
library(tidyverse)

# https://openpharma.github.io/mmrm/main/reference/emmeans_support.html

# check crp raw data
crp <-  arrow::read_parquet("/ocean/harbour/CDT71057/BP44551/CSRFinal_FA/ref/source/noncrf/LAB4_NC.parquet") |> 
  filter(LABTSTL == "C_REACTIVE_PROTEIN") |> 
  select(PATNUM, SUBJECT_STATUS, VISIT, VISIT_NAME, VISIT_NUMBER, VISIT_DATE, LABTPT, LABD, LABRESN)

# right join with adsl
adsl <- arrow::read_parquet("/ocean/harbour/CDT71057/BP44551/CSRFinal_Open/dev/data/adam/adsl.parquet") |> 
   mutate(PATNUM=as.numeric(substr(USUBJID, 9, 15)))
crp <- crp |> right_join(adsl |> select(PATNUM, FASFL), by="PATNUM") |> filter(FASFL == "Y")

crp |> filter(VISIT %in% c("D1", "D15")) |> group_by(PATNUM) |> summarise(n=n()) |> print(n = 60)

# check if there is missing visits
crp |> group_by(VISIT_NAME, VISIT_NUMBER) |> distinct(PATNUM) |> summarise(n=n()) |> arrange(VISIT_NUMBER) |> print()

# create dummy dataset for impute missing visits
df_visits <- expand.grid(
  PATNUM = unique(crp$PATNUM),
  VISIT = factor(unique(crp$VISIT), levels = c("SCRN", "PRN", "D1", "D8", "D15", "D22", "D29", "D43FU", "UNSCH"), ordered = T)
) |> arrange(PATNUM, VISIT)

# merge with crp
crp_full <- df_visits |> 
  left_join(crp, by = c("PATNUM" = "PATNUM", "VISIT" = "VISIT")) |> 
  group_by(PATNUM) |> 
  mutate(LABRESN_pre = lag(LABRESN))

# impute missing baseline/day1
crp_impute_missing_d1 <- crp_full |> 
    mutate(LABRESN_impute = case_when(
    VISIT == "D1" & is.na(LABRESN) ~ lag(LABRESN),
    TRUE ~ LABRESN)
  ) |> select(-LABRESN_pre) |> 
  filter(VISIT %in% c("D1", "D8", "D15", "D22", "D29", "D43FU")) |> 
  mutate(AVISITN=case_when(
    VISIT == "D1" ~ 0,
    VISIT == "D8" ~ 8,
    VISIT == "D15" ~ 15,
    VISIT == "D22" ~ 22,
    VISIT == "D29" ~ 29,
    VISIT == "D43FU" ~ 43
  ))


#----- cross check with adlb data, all matched ----
adlb <- arrow::read_parquet("/ocean/harbour/CDT71057/BP44551/CSRFinal_Open/dev/data/adam/adlb.parquet")
#adlb |> distinct(LBTEST, LBTESTCD, PARAMCD, PARAM) |> print(n=70)
adlb_crp <- adlb |>
  select(USUBJID, TRT01P, ITTFL, ANL01FL, AVISIT, AVISITN, LBTESTCD, LBORRES, AVAL) |> 
  filter(
    LBTESTCD == "CRP" & ITTFL == "Y" & ANL01FL == "Y"
  ) |> mutate(PATNUM=as.numeric(substr(USUBJID, 9, 15)))
  
crp_qc <- crp_impute_missing_d1 |> 
  left_join(adlb_crp |> select(PATNUM, AVISITN, AVAL), 
            by = c("PATNUM" = "PATNUM", "AVISITN" = "AVISITN")) |> 
  mutate(diff = LABRESN_impute - AVAL) |> 
  filter(!is.na(LABRESN_impute) & diff !=0)



# adpd is used for PD analysis
adpd <- arrow::read_parquet("/ocean/harbour/CDT71057/BP44551/CSRFinal_Open/dev/data/adam/adpd.parquet")
adpd |> distinct(PARAMCD, PARAM) |> print(n=70)
adpd_sub <- adpd |> 
  select(USUBJID, TRT01P, ITTFL, ANL01FL, AVISIT, AVISITN, PARAM, PARAMCD, AVAL, BASE, CHG, STRAT1R, STRAT2R, ABLFL) |> 
  filter(PARAMCD == "LOGBCRP" & ITTFL == "Y" & ANL01FL == "Y" & AVISIT != "BASELINE") |>
  #filter(PARAMCD == "LOGBCRP" & is.na(ABLFL) & !is.na(CHG)) %>%
  
  mutate(CFB = AVAL-BASE,
         AVISIT = factor(AVISIT, levels = c("Day 8", "Day 15", "Day 22", "Day 29")),
         #TRT01P = factor(TRT01P, levels = c("Placebo", "Low Dose", "High Dose"))
         STRAT1R = factor(STRAT1R, levels = c("LESS THAN 65%", "GREATER THAN OR EQUAL TO 65%")),
         STRAT2R = factor(STRAT2R, levels = c("LESS THAN 3%", "GREATER THAN OR EQUAL TO 3%"))                 
  ) |> 
  mutate(AVISIT = droplevels(AVISIT))

# fit MMRM model  
fit <- mmrm(
  formula = CHG ~ STRAT1R +STRAT2R + TRT01P * AVISIT + BASE * AVISIT + us(AVISIT | USUBJID),
  data = adpd_sub
)

df <- tidy(fit)
summary(fit)

df1 <- summary(emmeans(fit, ~ TRT01P | AVISIT))
df1a <- as.data.frame(df1) %>%
  mutate(
    emmean_1 = ((10^emmean) - 1) * 100,
    lower_cl = ((10^lower.CL) - 1) * 100,
   upper_cl = ((10^upper.CL) - 1) * 100
  )
df1a

df2 <- summary(pairs(emmeans(fit, ~ TRT01P | AVISIT), reverse = TRUE))
df2a <- as.data.frame(df2) %>%
  mutate(
    estimate_1 = ((10^estimate) - 1) * 100
   #lower_cl = ((10^lower.CL) - 1) * 100,
  # upper_cl = ((10^upper.CL) - 1) * 100
  )
df2a
#---- the MMRM analysis is correct. However the Day 43FU is not included in the model. ----