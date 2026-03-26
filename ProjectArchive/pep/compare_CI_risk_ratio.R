compare_rr_ci <- function(x1, n1, x2, n2, confidence_level = 0.95) {
  # --- MOVER-Wilson Method Calculation ---
  z <- qnorm(1 - (1 - confidence_level) / 2)
  
  # Wilson interval for proportion 1
  p1_hat <- x1 / n1
  term1_p1 <- x1 + z^2 / 2
  term2_p1 <- z * sqrt(x1 * (1 - x1/n1) + z^2 / 4)
  denominator1 <- n1 + z^2
  L1 <- (term1_p1 - term2_p1) / denominator1
  U1 <- (term1_p1 + term2_p1) / denominator1
  
  # Wilson interval for proportion 2
  p2_hat <- x2 / n2
  term1_p2 <- x2 + z^2 / 2
  term2_p2 <- z * sqrt(x2 * (1 - x2/n2) + z^2 / 4)
  denominator2 <- n2 + z^2
  L2 <- (term1_p2 - term2_p2) / denominator2
  U2 <- (term1_p2 + term2_p2) / denominator2
  
  # Ensure bounds are within [0, 1]
  L1 <- max(0, L1)
  U1 <- min(1, U1)
  L2 <- max(0, L2)
  U2 <- min(1, U2)
  
  # Calculate MOVER-Wilson CI
  if (L1 == 0 | U2 == 0 | L2 == 0 | U1 == 0) {
    mover_L <- NA
    mover_U <- NA
  } else {
    mover_L <- exp(log(L1) - log(U2))
    mover_U <- exp(log(U1) - log(L2))
  }
  
  # --- Log-Transformed Normal Approximation Method Calculation ---
  if (x1 == 0 | x2 == 0) {
    x1_corr <- x1 + 0.5
    n1_corr <- n1 + 1
    x2_corr <- x2 + 0.5
    n2_corr <- n2 + 1
  } else {
    x1_corr <- x1
    n1_corr <- n1
    x2_corr <- x2
    n2_corr <- n2
  }
  
  p1_hat_corr <- x1_corr / n1_corr
  p2_hat_corr <- x2_corr / n2_corr
  rr_hat <- p1_hat_corr / p2_hat_corr
  
  se_log_rr <- sqrt((1 - p1_hat) / (n1 * p1_hat) + (1 - p2_hat) / (n2 * p2_hat))
  
  log_rr_ci_L <- log(rr_hat) - z * se_log_rr
  log_rr_ci_U <- log(rr_hat) + z * se_log_rr
  
  normal_L <- exp(log_rr_ci_L)
  normal_U <- exp(log_rr_ci_U)
  
  # --- Combine and Return Results ---
  results <- matrix(
    c(rr_hat, mover_L, mover_U, rr_hat, normal_L, normal_U),
    ncol = 3, byrow = TRUE,
    dimnames = list(c("MOVER-Wilson", "Normal Approx"), c("Risk Ratio", "Lower CI", "Upper CI"))
  )
  
  return(results)
}

x1 <- 2
n1 <- 80
x2 <- 8
n2 <- 80

comparison_results <- compare_rr_ci(x1, n1, x2, n2)
print(comparison_results)