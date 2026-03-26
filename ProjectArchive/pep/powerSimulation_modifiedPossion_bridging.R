library(geepack)   # For GEE
library(tidyverse)
library(sandwich)   # For robust SEs
library(lmtest)     # For coeftest

# Parameters
n_clusters_per_group <- 70 #69  # Households per group, total number of HHC = 2 * n_clusters_per_group * cluster_size
cluster_size <- 2             # Subjects per household
baseline_bxm <- 0.02          # BXM transmission rate
baseline_placebo <- 0.10      # Placebo transmission rate
sd_household <- 0.02          # Household-level SD, sd=0 if no within household correlation
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
  
  # Calculate RR by household
  household_rr <-  contingency_table[1,2]/contingency_table[2,2]
  
  # Initialize results
  poisson_indiv_sig <- NA
  
  
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
    
    # Get point estimate of risk ratio for the assessment of the probability of bridging
    poisson_indiv_rr <- 1/exp(poisson_indiv_summary["treatmentPlacebo", "Estimate"])
    poisson_indiv_success1 <- poisson_indiv_rr < 1
    poisson_indiv_success2 <- poisson_indiv_rr < 1 - ((1-0.14) * 0.5)
    
  }, error = function(e) {})
  
  
  
  return(c(
    poisson_indiv = poisson_indiv_sig,
    poisson_indiv_bridge1 = poisson_indiv_success1,
    poisson_indiv_bridge2 = poisson_indiv_success2,
    indiv_rr = poisson_indiv_rr,
    household_rr = household_rr
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


# calculate power or probability of bridging
power_poisson_indiv <-  mean(results$poisson_indiv, na.rm = TRUE)
prob_bridge1 <- mean(results$poisson_indiv_bridge1, na.rm = TRUE)
prob_bridge2 <- mean(results$poisson_indiv_bridge2, na.rm = TRUE)


cat("Modified Poisson Power (individual-level):", round(power_poisson_indiv, 3), "\n")
cat("Probability of bridging (RR<1) based on Modified Poisson (individual-level):", round(prob_bridge1, 3), "\n")
cat("Probability of bridging (RR<0.57) based on Modified Poisson (individual-level):", round(prob_bridge2, 3))

# compare the RR by individual vs. by household
results %>%
  mutate(diff_rr = household_rr - indiv_rr) |> 
  
  summarise(mean_diff = mean(diff_rr, na.rm = TRUE),
            sd_diff = sd(diff_rr, na.rm = TRUE),
            min_diff = min(diff_rr, na.rm = TRUE),
            max_diff = max(diff_rr, na.rm = TRUE))

results %>%
  mutate(diff_rr = household_rr - indiv_rr) |> 
  ggplot(aes(x = indiv_rr, y=diff_rr)) +
  geom_point() +
  labs(title = "Household vs Individual RR",
       x = "RR by individul",
       y = "difference in RR by household vs. individual") +
  theme_minimal()

results |> ggplot(aes(x = indiv_rr, y = household_rr)) +
  geom_point(alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  #stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(title = "Household vs Individual RR",
       x = "RR by individul",
       y = "RR by household") +
  theme_minimal()
