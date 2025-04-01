# 非劣效性试验盲态样本量重估的I类错误模拟
# 目标：验证当真实效应等于非劣效界值时，盲法样本量重估会提高I类错误率
# https://onlinelibrary.wiley.com/doi/epdf/10.1002/bimj.200610373

library(ggplot2)

set.seed(123) # 固定随机种子

### 参数设置
delta <- -0.3     # 非劣效界值 (Δ)
sigma <- 1        # 真实标准差
alpha <- 0.05     # 名义I类错误率
power_target <- 0.8 # 目标功效
n1_per_group <- 20 # 初始每组样本量（内部试点）
sim_num <- 10000  # 模拟次数

### 模拟函数：单次试验流程
simulate_trial <- function(delta, sigma, n1_per_group, power_target, adjust=TRUE) {
  # 生成试点数据（真实效应为Δ，满足H0）
  data_T <- rnorm(n1_per_group, mean = delta, sd = sigma)
  data_C <- rnorm(n1_per_group, mean = 0, sd = sigma)
  
  # 盲法合并数据（忽略组别）
  pooled_data <- c(data_T, data_C)
  n_pilot <- length(pooled_data)
  mean_pilot <- mean(pooled_data)
  var_pilot <- var(pooled_data) # 合并方差估计（含组间差异的偏倚）
  
  if (adjust) {
    # 基于合并方差重新计算样本量（使用非劣效检验公式）
    # 样本量公式：N = 2*(Z_alpha + Z_beta)^2 * sigma^2 / (delta - theta)^2
    # 此处theta=delta（H0下），但实际需根据备择假设调整，这里简化为固定效应
    # 此处假设备择效应theta=delta + 0.1（仅为计算样本量，不影响H0模拟）
    theta_alt <- delta + 0.1 # 假设的备择效应（需根据实际研究设定）
    Za <- qnorm(1 - alpha)
    Zb <- qnorm(power_target)
    N_per_group <- ceiling(2 * (Za + Zb)^2 * var_pilot / (theta_alt - delta)^2)
    N_total <- 2 * max(N_per_group, n1_per_group) # 至少不低于试点样本量
  } else {
    # 固定样本量设计（不调整）
    N_total <- 2 * n1_per_group
  }
  
  # 补充数据（如果需要）
  if (N_total > 2 * n1_per_group) {
    n_add <- N_total - 2 * n1_per_group
    data_T_add <- rnorm(n_add/2, mean = delta, sd = sigma)
    data_C_add <- rnorm(n_add/2, mean = 0, sd = sigma)
    data_T <- c(data_T, data_T_add)
    data_C <- c(data_C, data_C_add)
  }
  
  # 进行单侧t检验（非劣效性检验：H0: mu_T - mu_C <= delta）
  t_test <- t.test(data_T, data_C, alternative = "greater", mu = delta)
  p_value <- t_test$p.value
  reject <- (p_value < alpha)
  
  return(reject)
}

### 运行模拟（有/无样本量重估）
reject_fixed <- replicate(sim_num, simulate_trial(delta, sigma, n1_per_group, power_target, adjust=FALSE))
reject_adjusted <- replicate(sim_num, simulate_trial(delta, sigma, n1_per_group, power_target, adjust=TRUE))

### 计算I类错误率
alpha_act_fixed <- mean(reject_fixed)
alpha_act_adjusted <- mean(reject_adjusted)

### 输出结果
cat(paste0(
  "固定样本量设计的I类错误率: ", round(alpha_act_fixed, 4), "\n",
  "盲态重估设计的I类错误率: ", round(alpha_act_adjusted, 4), "\n"
))

### 可视化结果（拒绝率分布）
df <- data.frame(
  Design = rep(c("Fixed", "Adjusted"), each = sim_num),
  Reject = c(reject_fixed, reject_adjusted)
)
ggplot(df, aes(x = Design, fill = factor(Reject))) +
  geom_bar(position = "fill") +
  geom_hline(yintercept = alpha, linetype = "dashed", color = "red") +
  labs(y = "Rejection Rate", title = "Type I error rate (10,000 replicates)") +
  scale_fill_manual(values = c("gray", "blue"), labels = c("Not reject", "Reject")) +
  theme_minimal()
