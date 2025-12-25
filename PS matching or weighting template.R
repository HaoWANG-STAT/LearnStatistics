# https://mp.weixin.qq.com/s/4MbY3aDikFbxnUFvhKGmQA

# 测试数据构建 ----
#Generating data similar to Austin (2009) for demonstrating treatment effect estimation
gen_X <- function(n) {
  X <- matrix(rnorm(9 * n), nrow = n, ncol = 9)
  X[,5] <- as.numeric(X[,5] < .5)
  X
}

#~20% treated
gen_A <- function(X) {
  LP_A <- - 1.2 + log(2)*X[,1] - log(1.5)*X[,2] + log(2)*X[,4] - log(2.4)*X[,5] + log(2)*X[,7] - log(1.5)*X[,8]
  P_A <- plogis(LP_A)
  rbinom(nrow(X), 1, P_A)
}

# Continuous outcome
gen_Y_C <- function(A, X) {
  2*A + 2*X[,1] + 2*X[,2] + 2*X[,3] + 1*X[,4] + 2*X[,5] + 1*X[,6] + rnorm(length(A), 0, 5)
}
#Conditional:
#  MD: 2
#Marginal:
#  MD: 2

# Binary outcome
gen_Y_B <- function(A, X) {
  LP_B <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
  P_B <- plogis(LP_B)
  rbinom(length(A), 1, P_B)
}
#Conditional:
#  OR:   2.4
#  logOR: .875
#Marginal:
#  RD:    .144
#  RR:   1.54
#  logRR: .433
#  OR:   1.92
#  logOR  .655

# Survival outcome
gen_Y_S <- function(A, X) {
  LP_S <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
  sqrt(-log(runif(length(A)))*2e4*exp(-LP_S))
}
#Conditional:
#  HR:   2.4
#  logHR: .875
#Marginal:
#  HR:   1.57
#  logHR: .452

set.seed(19599)

n <- 2000
X <- gen_X(n)
A <- gen_A(X)

Y_C <- gen_Y_C(A, X)
Y_B <- gen_Y_B(A, X)
Y_S <- gen_Y_S(A, X)

d <- data.frame(A, X, Y_C, Y_B, Y_S)


# 目标人群ATT/效应值-连续变量 ----
# 加载所需的R包
library("MatchIt")        # 用于进行倾向性评分匹配分析
library("WeightIt")       # 用于进行倾向性评分加权分析
library("marginaleffects")# 用于计算边际效应

# 使用最近邻匹配法进行倾向性评分匹配
# A为处理变量，X1-X9为协变量
# 使用可变比例匹配，每个处理组样本最多匹配2-4个对照组样本
mV <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                X6 + X7 + X8 + X9,
              data = d,        # 使用数据集d
              ratio = 2,       # 基础匹配比例为1:2
              max.controls = 4)# 最多匹配4个对照
mV  # 输出匹配结果

# 提取匹配后的数据集
md <- match_data(mV)

# 查看匹配后数据集的前几行
head(md)

# 构建线性回归模型
# Y_C为因变量，A与所有协变量(X1-X9)交互作为自变量
fit1 <- lm(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + 
                        X6 + X7 + X8 + X9),
           data = md,         # 使用匹配后的数据
           weights = weights) # 使用匹配权重

# 计算处理效应
# 计算A的平均边际效应，使用匹配子类作为方差估计单位
# 仅对处理组(A==1)计算平均处理效应
avg_comparisons(fit1,
                variables = "A",     # 指定处理变量
                vcov = ~subclass,    # 使用匹配子类计算标准误
                newdata = subset(A == 1)) # 限定在处理组计算

# 计算每个处理水平下的预测值
# 使用匹配子类作为方差估计单位
# 仅对处理组(A==1)计算预测值
avg_predictions(fit1,
                variables = "A",
                vcov = ~subclass,
                newdata = subset(A == 1))

## 目标人群ATE/效应值-连续变量 ----
# 仅在最终边际效应删除人群选择行
avg_comparisons(fit1,
                variables = "A",     # 指定处理变量
                vcov = ~subclass) 

# 目标人群ATT/效应值-二分类变量-RR ----
# 构建逻辑回归模型
# Y_B为二分类因变量，A与所有协变量(X1-X9)交互作为自变量
fit2 <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + 
                         X6 + X7 + X8 + X9),
            data = md,         # 使用匹配后的数据集
            weights = weights, # 使用匹配权重
            family = quasibinomial()) # 使用准二项分布族，可以处理过度离散

# 计算处理效应的相对风险比(RR)及其置信区间
# variables="A": 指定处理变量为A
# vcov=~subclass: 使用匹配子类作为方差估计单位
# newdata=subset(A==1): 仅对处理组计算效应
# comparison="lnratioavg": 计算对数风险比的平均值
# transform="exp": 对结果取指数转换，得到相对风险比
avg_comparisons(fit2,
                variables = "A",
                vcov = ~subclass,
                newdata = subset(A == 1),
                comparison = "lnratioavg",
                transform = "exp")
# 如要计算OR，comparison设定为 "lnoravg"即可。

## 估计ATE/效应值-二分类变量-RR-weighting方法 ----
#PS weighting for the ATE with a logistic regression PS
W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + 
                X6 + X7 + X8 + X9, data = d,
              method = "glm", estimand = "ATE")
W
#Logistic regression model with covariates
fit <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + 
                                 X6 + X7 + X8 + X9),
                    data = d, weightit = W,
                    family = binomial)

#Compute effects; RR and confidence interval
avg_comparisons(fit,
                variables = "A",
                comparison = "lnratioavg",
                transform = "exp")

# Compute predicted values for each treatment level
avg_predictions(fit,
                variables = "A")

### 目标人群ATE/效应值-二分类变量-风险差RD ----
avg_comparisons(fit,
                variables = "A",
                comparison = "differenceavg")


# 目标人群ATE/效应值-生存结果-HR ----

# 对于生存结局，有几种衡量效应量的方法。在使用Cox比例风险模型时，关注的量是治疗组和对照组之间的风险比（HR）。
# 与优势比（OR）一样，风险比是不可压缩的，这意味着只有当模型中不包含其他协变量时，估计的风险比才是边际风险比的有效估计。
# 其他效应度量，例如平均生存时间的差异或给定时间后的生存概率，可以像之前描述的连续型和二分类结局一样处理。

library("survival")

# 使用Cox比例风险模型计算边际风险比
# 使用匹配后的数据集md
# 使用匹配权重weights
# 使用匹配子类subclass作为聚类单位进行稳健标准误估计
coxph(Surv(Y_S) ~ A,
      data = md,
      robust = TRUE, 
      weights = weights,
      cluster = subclass)

