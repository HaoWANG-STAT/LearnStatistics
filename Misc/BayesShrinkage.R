library(brms)
library(ggplot2)

# 1. Create a synthetic dataset
# We have 3 groups: Young, Elderly, and a Rare Genotype (small N)
df <- data.frame(
  group = c(rep("Young", 500), rep("Elderly", 450), rep("Rare_Genotype", 10)),
  baseline_bp = rnorm(960, 150, 10),
  reduction = c(rnorm(500, 11, 5),   # Mean 11
                rnorm(450, 9, 5),    # Mean 9
                rnorm(10, 22, 5))    # Mean 22 (Small sample, likely an outlier)
)

# 2. Fit the Bayesian Hierarchical Model
# Formula: reduction ~ 1 + (1 | group)
# This says: "The reduction is influenced by a global intercept + group-specific offsets"

fit_shrinkage <- brm(
  formula = reduction ~ baseline_bp + (1 | group),
  data = df,
  prior = c(set_prior("normal(10, 5)", class = "Intercept")),
  control = list(adapt_delta = 0.99), # Increased from default 0.8
  chains = 4, iter = 2000, warmup = 1000
)
pairs(fit_shrinkage) # Check for convergence and parameter relationships

# 3. Extract the "Shrunken" Estimates (Posterior Means)
draws <- as_draws_df(fit_shrinkage)
# Calculate the group-specific means
group_results <- coef(fit_shrinkage)$group[, , "Intercept"]
print(group_results)
# 4. Visualize the Results
# Compare Raw vs Shrunken
raw_means <- aggregate(reduction ~ group, data = df, mean)
shrunken_means <- as.data.frame(group_results)
shrunken_means$group <- rownames(shrunken_means)

comparison <- merge(raw_means, shrunken_means, by = "group")
colnames(comparison) <- c("Group", "Raw_Mean", "Shrunken_Mean", "Lower_CI", "Upper_CI", "Est_Error")

print(comparison[, c("Group", "Raw_Mean", "Shrunken_Mean")])



# 计算“成功概率” (Probability of Success)
library(dplyr)

# 提取 Rare_Genotype 的总效应（全局截距 + 该组的偏离值）
# 在 brms 中，组效应通常表示为 r_group[Rare_Genotype,Intercept]
rare_effect <- draws$b_Intercept + draws$`r_group[Rare_Genotype,Intercept]`

# 计算概率：效果 > 10 的比例
prob_greater_10 <- mean(rare_effect > 10)

print(paste0("Rare Genotype 组减压效果优于 10 mmHg 的概率是: ", round(prob_greater_10 * 100, 2), "%"))


# 可视化后验分布 (Posterior Distributions)
library(tidybayes) # 推荐安装这个包，它是贝叶斯可视化的神器

df %>%
  modelr::data_grid(group, baseline_bp = 150) %>% # 统一设定基线血压为 150 进行公平比较
  add_epred_draws(fit_shrinkage) %>% # 获得预测分布
  ggplot(aes(x = .epred, y = group, fill = group)) +
  stat_halfeye() + # 绘制“半眼图”：包含密度曲线和置信区间
  geom_vline(xintercept = 10, linetype = "dashed") +
  labs(title = "Posterior Distributions of Group Effects",
       subtitle = "Adjusted for Baseline BP = 150 mmHg",
       x = "Expected Blood Pressure Reduction (mmHg)",
       y = "Subgroup") +
  theme_minimal()
