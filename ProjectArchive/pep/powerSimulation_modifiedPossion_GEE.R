library(geepack)   # For GEE
library(tidyverse)
library(sandwich)   # For robust SEs
library(lmtest)     # For coeftest

# Parameters
n_clusters_per_group <- 70 #69  # Households per group
cluster_size <- 2             # Subjects per household
baseline_bxm <- 0.02          # BXM transmission rate
baseline_placebo <- 0.10      # Placebo transmission rate
sd_household <- 0.02          # Household-level SD
alpha <- 0.05                 # Significance level

# Function to simulate data and test power
simulate_power <- function() {
  # Simulate households with random effects
  households <- data.frame(
    household_id = 1:(2 * n_clusters_per_group),
    treatment = factor(rep(c("BXM", "Placebo"), each = n_clusters_per_group),
                       levels = c("BXM", "Placebo")),
    baseline = rep(c(baseline_bxm, baseline_placebo), each = n_clusters_per_group)
  )
  
  # Household random effect
  households$u <- rnorm(nrow(households), mean = 0, sd = sd_household)
  households$p <- households$baseline + households$u
  households$p <- pmin(pmax(households$p, 0), 1)  # Truncate to [0, 1]
  
  # Generate individual outcomes
  individual_data <- households %>%
    slice(rep(1:n(), each = cluster_size)) %>%
    mutate(
      y = rbinom(n(), 1, p)  # Binary outcome (0/1)
    )
  
  # Aggregate to household level: "any case" indicator
  household_data <- individual_data %>%
    group_by(household_id, treatment) %>%
    summarize(
      any_case = as.integer(sum(y) > 0),  # 1 if ≥1 case, 0 otherwise
      .groups = "drop"
    )
  
  # Create contingency table for Fisher's test
  contingency_table <- table(
    household_data$treatment,
    household_data$any_case
  )
  
  # Initialize results
  gee_sig <- NA
  poisson_indiv_sig <- NA
  household_analysis_sig <- NA
  
  # 1. GEE (clustered individual-level)
  tryCatch({
    gee_model <- geeglm(
      y ~ treatment,
      family = binomial(link = "logit"),
      data = individual_data,
      id = household_id,
      corstr = "exchangeable"
    )
    gee_p <- summary(gee_model)$coefficients["treatmentPlacebo", "Pr(>|W|)"]
    gee_sig <- gee_p < alpha
  }, error = function(e) {})
  
  # 2. Modified Poisson (unclustered individual-level)
  tryCatch({
    poisson_indiv_model <- glm(
      y ~ treatment,
      family = poisson(link = "log"),
      data = individual_data
    )
    poisson_indiv_summary <- coeftest(poisson_indiv_model, vcov = vcovHC(poisson_indiv_model, type = "HC0"))
    poisson_indiv_p <- poisson_indiv_summary["treatmentPlacebo", "Pr(>|z|)"]
    poisson_indiv_sig <- poisson_indiv_p < alpha
  }, error = function(e) {})
  
  # 3. Household-level analysis: Modified Poisson (relative risk)
  tryCatch({
    # Fit modified Poisson regression for binary outcome (any_case)
    household_model <- glm(
      any_case ~ treatment,
      family = poisson(link = "log"),  # Log link for relative risk
      data = household_data
    )
    household_summary <- coeftest(household_model, vcov = vcovHC(household_model, type = "HC0"))
    household_p <- household_summary["treatmentPlacebo", "Pr(>|z|)"]
    household_analysis_sig <- household_p < alpha
  }, error = function(e) {})
  
  # 4. Household-level Fisher's exact test
  tryCatch({
    fisher_test <- fisher.test(contingency_table)
    fisher_sig <- fisher_test$p.value < alpha
  }, error = function(e) {})
  
  return(c(
    gee = gee_sig,
    poisson_indiv = poisson_indiv_sig,
    household_poisson = household_analysis_sig,
    fisher = fisher_sig
  ))
}

# Run simulations
set.seed(123)
#results <- replicate(n_simulations, simulate_power())
coresAvailable <- parallel::detectCores()
results <- bind_rows(parallel::mclapply(1:10000, 
                                        function(i) simulate_power(),
                                        mc.cores = coresAvailable)
)

# Calculate power
power_gee <- mean(results$gee, na.rm = TRUE)
power_poisson_indiv <- mean(results$poisson_indiv, na.rm = TRUE)
power_household_poisson <- mean(results$household_poisson, na.rm = TRUE)
power_fisher <- mean(results$fisher, na.rm = TRUE)

cat("GEE Power (clustered):", round(power_gee, 3), "\n")
cat("Modified Poisson Power (individual-level):", round(power_poisson_indiv, 3), "\n")
cat("Modified Poisson Power (Household-Leve):", round(power_household_poisson, 3), "\n")
cat("Fisher's Exact Test Power:", round(power_fisher, 3))