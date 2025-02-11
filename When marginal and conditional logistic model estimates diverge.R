# When marginal and conditional logistic model estimates of risk difference diverge
# https://www.rdatagen.net/post/marginal-v-conditional/

set.seed(287362)

library(simstudy)
library(data.table)
library(ggplot2)
library(ggthemes)
library(parallel)

## First example ----
# Define data
def1 <- defData(varname = "clustEff", formula = 0, variance = 2, 
                id = "cID")
def1 <- defData(def1, varname = "nInd", formula = 1000, 
                dist = "noZeroPoisson")

def2 <- defDataAdd(varname = "Y0", formula = "-1 + clustEff", 
                   dist = "binary", link = "logit")
def2 <- defDataAdd(def2, varname = "Y1", 
                   formula = "-1 + clustEff + 2", 
                   dist = "binary", link = "logit")

# Generate cluster level data

set.seed(123)

dtC <- genData(n = 100, def1)

# Generate individual level data

dt <- genCluster(dtClust = dtC, cLevelVar = "cID", numIndsVar = "nInd", 
                 level1ID = "id")

dt <- addColumns(def2, dt)

dtLong <- addPeriods(dtName = dt, idvars = c("id","cID"), 
                     nPeriods = 2,timevars = c("Y0","Y1"), 
                     timevarName = "Y"
)
# Calculate average potential outcomes by exposure (which is called period)

dtMean <- dtLong[, .(Y = mean(Y)), keyby = .(period, cID)] # conditional mean
dtMMean <- dtLong[, .(Y = mean(Y)), keyby = .(period)] # marginal mean
dtMMean[, cID := 999]

ggplot(data = dtMean, aes(x=factor(period), y = Y, group= cID)) +
  # geom_jitter(width= .25, color = "grey75") +
  geom_line(color = "grey75", position=position_jitter(w=0.02, h=0.02)) +
  geom_point(data=dtMMean) +
  geom_line(data=dtMMean, size = 1, color = "red") +
  ylab("Estimated cluster probability") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.title.x = element_blank()) 

dtChange <- dt[, .(Y0 = mean(Y0), Y1 = mean(Y1)), keyby = cID]
dtChange[, riskdiff := round(Y1 - Y0, 2)]
dtChange[, loratio := round( log( (Y1 / (1-Y1)) / (Y0 / (1-Y0) )), 2)]
dtChange[, odds_ratio := round( exp(loratio), 2)]

dtChange[sample(1:100, 10, replace = F),
         .(Y0 = round(Y0,2), Y1 = round(Y1,2), riskdiff, loratio), 
         keyby=cID]

### investigate marginal vs. conditional effect estimates for OR ----
glmfit <- glm(Y ~ period, data = dtLong, family = "binomial")
newdata <- dtLong[, .(period=1, cID)] 
p1 <- mean(predict(glmfit, newdata, type = "response"))

newdata <- dtLong[, .(period=0)]
p0 <- mean(predict(glmfit, newdata, type = "response"))

c(p1, p0)

risk_diff <- p1 - p0
odds_ratio <- exp(coef(glmfit)["period"])

c(rd = risk_diff, or = odds_ratio)

# 
dtChangeMean <- dtChange[, .(RD = mean(riskdiff), logOR = mean(loratio))]
dtChangeMean[, OR := exp(logOR)]

## Second example ----
def1 <- defData(varname = "clustEff", formula = 0, variance = 2, id = "cID")
def1 <- defData(def1, varname = "nInd", formula = 100, dist = "noZeroPoisson")

# Each individual now has a measured age

def2 <- defDataAdd(varname = "age", formula = 0, variance = 2)
def2 <- defDataAdd(def2, varname = "Y", 
                   formula = "-4 + clustEff + 2*trt + 2*age", 
                   dist = "binary", link = "logit")

# Generate cluster level data

dtC <- genData(200, def1) 
dtC <- trtAssign(dtC, grpName = "trt") 

# Generate individual level data

dt <- genCluster(dtClust = dtC, cLevelVar = "cID", numIndsVar = "nInd", 
                 level1ID = "id")
dt <- addColumns(def2, dt)

glmerFit1 <- glmer(Y ~ trt + age + (1 | cID), data = dt, family="binomial")
glmFit1 <- glm(Y ~ trt + age, family = binomial, data = dt)

newCond <- expand.grid(cID = unique(dt$cID), age=seq(-4, 4, by =.1))
newCond0 <- data.table(trt = 0, newCond)
newCond1 <- data.table(trt = 1, newCond)

newMarg0 <- data.table(trt = 0, age = seq(-4, 4, by = .1))
newMarg1 <- data.table(trt = 1, age = seq(-4, 4, by = .1))

newCond0[, pCond0 := predict(glmerFit1, newCond0, type = "response")]
newCond1[, pCond1 := predict(glmerFit1, newCond1, type = "response")]

newMarg0[, pMarg0 := predict(glmFit1, newMarg0, type = "response")]
newMarg0[, pCAvg0 := predict(glmerFit1, newMarg0[,c(1,2)], 
                             re.form = NA, type="response")]

newMarg1[, pMarg1 := predict(glmFit1, newMarg1, type = "response")]
newMarg1[, pCAvg1 := predict(glmerFit1, newMarg1[,c(1,2)], 
                             re.form = NA, type="response")]

dtAvg <- data.table(age = newMarg1$age, 
                    avgMarg = newMarg1$pMarg1 - newMarg0$pMarg0, 
                    avgCond = newMarg1$pCAvg1 - newMarg0$pCAvg0
)

p1 <- ggplot(aes(x = age, y = pCond1), data=newCond1) + 
  geom_line(color="grey", aes(group = cID)) +
  geom_line(data=newMarg1, aes(x = age, y = pMarg1), color = "red", size = 1) +
  geom_line(data=newMarg1, aes(x = age, y = pCAvg1), color = "black", size = 1) +
  ggtitle("Treatment group") +
  xlab("Age") +
  ylab("Probability") +
  my_theme()

p0 <- ggplot(aes(x = age, y = pCond0), data=newCond0) + 
  geom_line(color="grey", aes(group = cID)) +
  geom_line(data=newMarg0, aes(x = age, y = pMarg0), color = "red", size = 1) +
  geom_line(data=newMarg0, aes(x = age, y = pCAvg0), color = "black", size = 1) +
  ggtitle("Control group") +
  xlab("Age") +
  ylab("Probability") +
  my_theme()

pdiff <- ggplot(data = dtAvg) + 
  geom_line(aes(x = age, y = avgMarg), color = "red", size = 1) +
  geom_line(aes(x = age, y = avgCond), color = "black", size = 1) +
  ggtitle("Risk difference") +
  xlab("Age") +
  ylab("Probability") +
  my_theme()

grid.arrange(p1, p0, pdiff)