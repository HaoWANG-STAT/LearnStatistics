#---- data from Astegolimab study ----
DF<- rice::rice_read("root/clinical_studies/RO7187807/CDT70155/GB39242/data_analysis/CSR/data/sdtmv/lb.sas7bdat")
dm <- rice::rice_read("root/clinical_studies/RO7187807/CDT70155/GB39242/data_analysis/CSR/data/sdtmv/dm.sas7bdat")
                     
NEUT <- DF |> filter(LBTESTCD == "NEUT" & !(VISITNUM %in% c(1,2,20))) |> 
  left_join(dm |> select(USUBJID, ARM), by = "USUBJID") |> 
  filter(ARM != "Not Assigned")
  # filter(!is.na(LBORRES)) |> 
  # mutate(
  #   VISIT = factor(VISIT, levels = c("BASELINE", "CYCLE 1 DAY 8", "CYCLE 1 DAY 15", "CYCLE 2 DAY 1", "CYCLE 2 DAY 8", "CYCLE 2 DAY 15", "CYCLE 3 DAY 1", "CYCLE 3 DAY 8", "CYCLE 3 DAY 15", "CYCLE 4 DAY 1", "CYCLE 4 DAY 8", "CYCLE 4 DAY 15")),
  #   TRTA = factor(TRTA, levels = c("A: Placebo", "B: Low dose", "C: High dose"))
  # )

# boxplot of NEUT by visit and treatment arm
library(ggplot2)
NEUT |> filter(VISITNUM==3) |> 
  ggplot(aes(x = VISIT, y = LBSTRESN, fill = ARM)) +
  geom_boxplot() +
  labs(title = "Neutrophil Counts by Visit and Treatment Arm",
       x = "Visit",
       y = "Neutrophil Count (10^9/L)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# histogram of NEUT in different facets arm
NEUT |> filter(VISITNUM==3) |> 
  ggplot(aes(x = LBSTRESN, fill = ARM)) +
  geom_histogram(binwidth = 0.5, position = "dodge", color = "black") +
  labs(title = "Distribution of Neutrophil Counts by Treatment Arm at Visit 3",
       x = "Neutrophil Count (10^9/L)",
       y = "Frequency") +
  theme_minimal() +
  facet_wrap(~ ARM)


# mean and ci plot by visit and treatment arm
library(dplyr)
library(ggplot2)
library(ggpubr)
summary_data <- NEUT |> 
  group_by(VISITNUM, ARM) |> 
  summarise(
    mean_neut = mean(LBSTRESN, na.rm = TRUE),
    sd_neut = sd(LBSTRESN, na.rm = TRUE),
    n = n(),
    se_neut = sd_neut / sqrt(n),
    ci_lower = mean_neut - qt(0.975, df = n - 1) * se_neut,
    ci_upper = mean_neut + qt(0.975, df = n - 1) * se_neut
  ) |> 
  ungroup()
# print(summary_data)
ggplot(summary_data, aes(x = VISITNUM, y = mean_neut, color = ARM, group = ARM)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(title = "Mean Neutrophil Counts by Visit and Treatment Arm",
       x = "Visit",
       y = "Mean Neutrophil Count (10^9/L)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# mean change from baseline plot by visit and treatment arm
baseline_data <- NEUT |> 
  filter(VISITNUM == 3) |> 
  select(USUBJID, LBSTRESN) |> 
  rename(baseline_neut = LBSTRESN)
NEUT_change <- NEUT |>
  left_join(baseline_data, by = "USUBJID") |> 
  mutate(change_from_baseline = LBSTRESN - baseline_neut)
summary_change_data <- NEUT_change |>
  group_by(VISITNUM, ARM) |> 
  summarise(
    mean_change = mean(change_from_baseline, na.rm = TRUE),
    sd_change = sd(change_from_baseline, na.rm = TRUE),
    n = n(),
    se_change = sd_change / sqrt(n),
    ci_lower = mean_change - qt(0.975, df = n - 1) * se_change,
    ci_upper = mean_change + qt(0.975, df = n - 1) * se_change
  ) |> 
  ungroup()
# print(summary_change_data)
ggplot(summary_change_data, aes(x = VISITNUM, y = mean_change, color = ARM, group = ARM)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(title = "Mean Change in Neutrophil Counts from Baseline by Visit and Treatment Arm",
       x = "Visit",
       y = "Mean Change in Neutrophil Count (10^9/L)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#---- check blood eos < 300 /150 pts number ----
LAB2 <- arrow::read_parquet("/ocean/harbour/CDT71057/BP44551/CSRFinal_FA/ref/source/noncrf/LAB2_NC.parquet") |> 
  filter(VISIT_CODE =="PRN" & LABTSTL=="EOSINOPHIL_ABS") 

# read dm.sas7bdat from /ocean/harbour/CDT71057/BP44551/CSRFinal_FA/dev/data/sdtmv
# dm <- haven::read_sas("/ocean/harbour/CDT71057/BP44551/CSRFinal_FA/dev/data/sdtmv/dm.sas7bdat") |> 
#   filter(ARM %in% c("SELNOFLAST 200MG", "PLACEBO"))

# get pts got randomized
disc <- arrow::read_parquet("/ocean/harbour/CDT71057/BP44551/CSRFinal_FA/ref/source/crf/DISC1.parquet") |> 
  filter(EPOCH=="Blinded Treatment" & (is.na(DCRS) | DCRS != "FAILURE TO MEET RANDOMIZATION CRITERIA") )


LAB2_1 <- LAB2 |> 
  # merge by LAB_89.PATNUM and dis.SUBJECT
  inner_join(disc |> select(SUBJECT), by = c("PATNUM_RAW" = "SUBJECT")) 
# count LABRESN<0.15 and LABRESN<0.3
LAB2_1 |> filter(LABRESN<0.15) |> distinct(PATNUM_RAW) |> nrow()
LAB2_1 |> filter(LABRESN<0.3) |> distinct(PATNUM_RAW) |> nrow()
  

#---- check FENO data ----
# subjects with multiple FENO records on single visit: 10044, 10048
feno <- arrow::read_parquet("/ocean/harbour/CDT71057/BP44551/CSRFinal_FA/ref/source/noncrf/FNO1_NC.parquet") |> 
  filter(SUBJECT %in% c("10044", "10048")) |> 
  select(PATNUM, VISIT, LABTPT, LABTM, LABRESN)
view(feno)

# ----- check sputum neutrophil data ----
LAB8 <- arrow::read_parquet("/ocean/harbour/CDT71057/BP44551/CSRFinal_FA/ref/source/noncrf/LAB8_NC.parquet")
LAB9 <- arrow::read_parquet("/ocean/harbour/CDT71057/BP44551/CSRFinal_FA/ref/source/noncrf/LAB9_NC.parquet")
LAB_89 <- bind_rows(LAB8, LAB9) |> 
    filter(LABTSTL %in% c("SPUTUM_NEUTROPHIL_TOTAL_PER_NON_SQUAMOUS_CELLS_PCT", "SPUTUM_NEUTROPHIL_PERC"))

# sputum assessment eCRF data
sputma <- arrow::read_parquet("/ocean/harbour/CDT71057/BP44551/CSRFinal_FA/ref/source/crf/SPUTMA.parquet") |> 
  select(SUBJECT, FOLDER, SPLWT, TOTCEL)
  