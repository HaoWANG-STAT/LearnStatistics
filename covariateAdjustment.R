# https://www.rdatagen.net/post/2025-01-28-handling-covariates-in-an-analysis-of-an-rct/
library(simstudy)
library(data.table)
#library(stargazer)
library(parallel)

# Data definitions
def_c <- 
  defData(varname = "x1", formula = 0, variance = "..s_x1^2", dist = "normal") |>
  defData(varname = "x2", formula = 0, variance = "..s_x2^2", dist = "normal") |>
  defData(varname = "A", formula = "1;1", dist = "trtAssign")

def_y <- defDataAdd(
  varname = "y", 
  formula = "5 + ..delta * A + ..b1 * x1 + ..b2 * x2", 
  variance = "..s_y^2", 
  dist = "normal"
)

# Initial parammeters
s_x1 <- 8
s_x2 <- 9
s_y <- 12
b1 <- 1.50
b2 <- 0.0
delta <- 5

# Single dataset generation
RNGkind("L'Ecuyer-CMRG")
set.seed(55)

dc <- genData(250, def_c) 
dd <- addColumns(def_y, dc)

head(dd)

dd[, .(mu_x1 = mean(x1), mu_x2 = mean(x2)), keyby = A]

calc_diff <- function(dx, rx, v) {
  
  mean_diff <- dx[get(rx)==1, mean(get(v))] - dx[get(rx)==0, mean(get(v))]
  s_pooled <- sqrt(
    (dx[get(rx)==1, (.N - 1) * var(get(v))] + 
       dx[get(rx)==0, (.N - 1) * var(get(v))] ) / dx[, .N - 2])
  
  Z_x <- mean_diff / ( s_pooled * sqrt(1/dx[get(rx)==0, .N] + 1/dx[get(rx)==1, .N]) )
  d <- mean_diff / s_pooled
  
  return(list(Z_x = Z_x, d = d))
  
}

calc_diff(dd, "A", "x1")
calc_diff(dd, "A", "x2")
dd[, .(rho_x1.y = cor(x1, y), rho_x2.y = cor(x2, y))]

# Model estimation
model_1 <- lm(data = dd, formula = y ~ A)
model_2 <- lm(data = dd, formula = y ~ A + x1)
model_3 <- lm(data = dd, formula = y ~ A + x2)
model_4 <- lm(data = dd, formula = y ~ A + x1 + x2)

# Operating characteristics (based on replicated datasets)
est_ancova <- function(dx, vars) {
  
  formula <- as.formula(paste("y ~", paste(vars, collapse = " + ")))
  model <- lm(data = dx, formula = formula)
  
  coef_summary <- summary(model)$coefficients["A", ]
  t_stat <- coef_summary["t value"]
  
  p_value <- pt(t_stat, df = model$df.residual, lower.tail = FALSE)
  ests <- data.table(t(coef_summary[1:2]), p_value)
  setnames(ests, c("est", "se", "pval"))
  
  return(ests)
  
}

replicate <- function() {
  
  dc <- genData(250, def_c) 
  dd <- addColumns(def_y, dc)
  
  est_1 <- est_ancova(dd, vars = "A")
  est_2 <- est_ancova(dd, vars = c("A", "x1"))
  est_3 <- est_ancova(dd, vars = c("A", "x2"))
  est_4 <- est_ancova(dd, vars = c("A", "x1", "x2"))
  
  return(list(est_1 = est_1, est_2 = est_2, est_3 = est_3, est_4 = est_4))
  
}
res <- mclapply(1:2000, function(x) replicate())

get.field <- function(x, field) {
  data.table(t(sapply(x, function(x) x[[field]]) ))
}

ests <- rbindlist(lapply(res, function(x) get.field(x, "est")))
sapply(ests, function(x) c(bias = mean(x) - delta, var = var(x)))

pvals <- rbindlist(lapply(res, function(x) get.field(x, "pval")))
sapply(pvals, function(x) c(mean(x < 0.025))) 

# Exploring Type 1 error rates
delta <- 0

res <- mclapply(1:2000, function(x) replicate())

pvals <- rbindlist(lapply(res, function(x) get.field(x, "pval")))
sapply(pvals, function(x) c(mean(x < 0.025))) 

# To mimic the requirement that the dataset is sampled conditional on a fixed level of covariate imbalance
replicate_2 <- function() {
  
  dd <- addColumns(def_y, dc)
  
  est_1 <- est_ancova(dd, vars = "A")
  est_2 <- est_ancova(dd, vars = c("A", "x1"))
  est_3 <- est_ancova(dd, vars = c("A", "x2"))
  est_4 <- est_ancova(dd, vars = c("A", "x1", "x2"))
  
  return(list(est_1 = est_1, est_2 = est_2, est_3 = est_3, est_4 = est_4))
  
}

dc <- genData(250, def_c) 

res <- mclapply(1:2000, function(x) replicate_2())

pvals <- rbindlist(lapply(res, function(x) get.field(x, "pval")))
sapply(pvals, function(x) c(mean(x < 0.025))) 