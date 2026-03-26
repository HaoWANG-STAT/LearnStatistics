# Estimating a risk difference (and confidence intervals) using logistic regression

set.seed(287362)

library(simstudy)
library(data.table)
library(ggplot2)
library(ggthemes)
library(parallel)

s_define <- function() {
  
  def <- defDataAdd(varname = "x1", formula = "..mu_x", variance = 8, dist = "beta") ## ..mu_x is a placeholder for the argument mu_x
  def <- defDataAdd(def, varname = "x", formula = "x1 - 0.5", dist = "nonrandom")
  def <- defDataAdd(def, varname = "y", 
                    formula = "-2 + 1 * rx + 1.5 * x",
                    dist = "binary", link="logit")
  
  return(list(def = def)) # list_of_defs is a list of simstudy data definitions
}

s_generate <- function(list_of_defs, argsvec) {
  
  list2env(list_of_defs, envir = environment())
  list2env(as.list(argsvec), envir = environment())
  
  dx <- genData(n)
  dx <- trtAssign(dx, grpName = "rx")
  dx <- addColumns(def, dx)
  
  return(dx) #  generated data is a data.table
}

s_model <- function(dx) {
  
  glmfit <- glm(y ~ rx + x, data = dx, family = "binomial")
  
  newdata <- dx[, .(rx=1, x)]
  p1 <- mean(predict(glmfit, newdata, type = "response"))
  
  newdata <- dx[, .(rx=0, x)]
  p0 <- mean(predict(glmfit, newdata, type = "response"))
  
  risk_diff <- p1 - p0
  odds_ratio <- exp(coef(glmfit)["rx"])
  
  model_results <- data.table(risk_diff, odds_ratio)
  
  return(model_results) # model_results is a data.table
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
      mc.cores = 16)
  )
  
  model_results <- cbind(t(argsvec), model_results)
  
  return(model_results) # summary_stats is a data.table
}

### Scenarios

scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}

#----

n <- 500
mu_x <- c(0.2, 0.4, 0.6, 0.8)

scenarios <- scenario_list(n = n, mu_x = mu_x)

summary_stats <- rbindlist(lapply(scenarios, function(a) s_replicate(a, nsim = 1000)))

ggplot(data = summary_stats, aes(x = risk_diff, group = mu_x)) +
  geom_density(aes(fill = factor(mu_x)), alpha = .7) +
  scale_fill_canva(palette = "Simple but bold", name = "mu_x") +
  theme(panel.grid = element_blank()) +
  xlab("estimated risk difference")

ggplot(data = summary_stats, aes(x = odds_ratio, group = mu_x)) +
  geom_density(aes(fill = factor(mu_x)), alpha = .7) +
  scale_fill_canva(palette = "Simple but bold", name = "mu_x") +
  theme(panel.grid = element_blank()) +
  xlab("estimated odds ratio")
