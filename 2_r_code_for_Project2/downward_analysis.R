# R code for analysis of decay region.

# Load useful packages. ----
library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(performance)

# Set random seed. ----
set.seed(123)

# Load data. ----
myData <- read_csv("./data_for_downward_analysis.csv")

## Specify date columns. ----
myData <- myData |>
  mutate(date = mdy(date),
         positiveTestDate = mdy(positiveTestDate),
         symptomOnset = mdy(symptomOnset),
         dose1Date = mdy(dose1Date),
         dose2Date = mdy(dose2Date),
         dose3Date = mdy(dose3Date),
         dose4Date = mdy(dose4Date))

# Transform data. ----

## Create dsfd variable (days since first dose). ----
myData <- myData |>
  mutate(dsfd = as.numeric(difftime(date, dose1Date, units = 'days')))

## Create ysfd variable (years since first dose). ----
## This is done for numerical purposes, specifically, 
## to avoid convergence issues.
myData <- myData |>
    mutate(ysfd = dsfd / 365)

## Scale NAb flow a number between 0 and 1. ----
myData <- myData |>
  mutate(nabFlow = nabFlowAverage / 100)

## Windsorize NAb flow variable. ----
## * We need this as it 
##   (a) seems to give us constant variance, and
##   (b) nabFlow values are restricted in [0, 1].

### Specify adjustment term for winsorization. ----
adj <- .001 

### Windsorize NAb flow variable. ----
### First, filter out NA's.
tol <- sqrt(.Machine$double.eps)

myData <- myData |>
  filter(!is.na(nabFlow)) |>
  mutate(nabFlowAdj = ifelse(abs(nabFlow - 1) < tol, 1 - adj, nabFlow)) |>
  mutate(nabFlowAdj = ifelse(abs(nabFlowAdj - 0) < tol, adj, nabFlowAdj)) |>
  mutate(nabFlowT = log(nabFlowAdj / (1 - nabFlowAdj)))

## Center ysfd and include higher orders of it. ----
## Also include higher orders of ysfd.
myData <- myData |>
  mutate(ysfd2 = ysfd^2,
         ysfd3 = ysfd^3,
         ysfd4 = ysfd^4,
         ysfd5 = ysfd^5,
         ysfd6 = ysfd^6,
         ysfd7 = ysfd^7) |>
  mutate(ysfdc = as.vector(scale(ysfd, center = TRUE, scale = FALSE)),
         ysfdc2 = ysfdc^2,
         ysfdc3 = ysfdc^3,
         ysfdc4 = ysfdc^4,
         ysfdc5 = ysfdc^5,
         ysfdc6 = ysfdc^6,
         ysfdc7 = ysfdc^7)

## Create cleaner convalesced/naive variable. ----
## (For plotting purposes.)
myData <- myData |>
    mutate(isConvalescedLabel = ifelse(isConvalesced == "Yes", "Convalesced", "Naive"))

# Filter data. ----
# * Include subjects at least one vaccine dose, no breakthroughs,
#   no booster, only Moderna and Pfizer vaccine.
# * We chose to look only at day 37 and beyond as it focuses n the
#   downward trend.
# * Remove specific individuals as it was determined these 
#   subjects had breakthroughs despite meeting the above 
#   conditions.
daycutoff <- 37

filteredData <- myData |>
  # Filter for subjects with the vaccine.
  filter(hasDose1 == 'Yes') |>
  # Filter out breakthrough and unknown events.
  filter(cohort %in% c("1", "2", "3", "13")) |>
  # Filter for either Moderna or Pfizer for 1st dose.
  filter(dose1Type %in% c("Moderna", "Pfizer")) |>
  #filter(dose2Type %in% c("Moderna", "Pfizer")) |>  
  # Filter for observations before third dose.
  filter(is.na(dose3Date) | date < dose3Date) |>
  # Filter for observations 40 days after first dose date.
  # 36 days includes 7 convalesced subjects.
  # 37 days includes 6 convalesced subjects.
  filter(dsfd >= daycutoff) |>
  # Remove subject possible breakthroughs.
  filter(subject != "67R") |>
  filter(subject != "427R") |>
  filter(subject != "425R") |>
  filter(subject != "407R")

# Plot filtered data set. ----

myXsOrig <- seq(0, 460, by = 30)
ggplot(filteredData, aes(x = dsfd, y = nabFlow, 
                         group = subject, 
                         color = isConvalescedLabel)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_line() +  
  geom_vline(xintercept = 37, lty = 2) +    
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(labels = myXsOrig, breaks = myXsOrig) +    
  labs(x = "days since first dose", y = "NAb activity (% neutralization)", color = NULL) +
  theme_classic() +
  theme(legend.position = c(.7, .95),
        legend.box.background = element_rect(color = "black"))

# Time predictor selection. ----
model1o <- lmer(nabFlowT ~ 1 + isConvalesced + 
                 ysfdc + ysfdc:isConvalesced +
                (1 + ysfdc | subject),
               data = filteredData,
               REML = FALSE)

model2o <- lmer(nabFlowT ~ 1 + isConvalesced + 
                 ysfdc + ysfdc:isConvalesced +
                 ysfdc2 + ysfdc2:isConvalesced +
                (1 + ysfdc | subject),
               data = filteredData,
               REML = FALSE)

model3o <- lmer(nabFlowT ~ 1 + isConvalesced + 
                 ysfdc + ysfdc:isConvalesced +
                 ysfdc2 + ysfdc2:isConvalesced + 
                 ysfdc3 + ysfdc3:isConvalesced +
                (1 + ysfdc | subject),
               data = filteredData,
               REML = FALSE)

model4o <- lmer(nabFlowT ~ 1 + isConvalesced + 
                 ysfdc + ysfdc:isConvalesced +
                 ysfdc2 + ysfdc2:isConvalesced + 
                 ysfdc3 + ysfdc3:isConvalesced +
                 ysfdc4 + ysfdc4:isConvalesced +
               (1 + ysfdc | subject),
               data = filteredData,
               REML = FALSE)

model5o <- lmer(nabFlowT ~ 1 + isConvalesced + 
                 ysfdc + ysfdc:isConvalesced +
                 ysfdc2 + ysfdc2:isConvalesced + 
                 ysfdc3 + ysfdc3:isConvalesced +
                 ysfdc4 + ysfdc4:isConvalesced +
                 ysfdc5 + ysfdc5:isConvalesced +
               (1 + ysfdc | subject),
               data = filteredData,
               REML = FALSE)

model6o <- lmer(nabFlowT ~ 1 + isConvalesced + 
                 ysfdc + ysfdc:isConvalesced +
                 ysfdc2 + ysfdc2:isConvalesced + 
                 ysfdc3 + ysfdc3:isConvalesced +
                 ysfdc4 + ysfdc4:isConvalesced +
                 ysfdc5 + ysfdc5:isConvalesced +
                 ysfdc6 + ysfdc6:isConvalesced +
               (1 + ysfdc | subject),
               data = filteredData,
               REML = FALSE)

anova(model1o, model2o, model3o, model4o, model5o, model6o)
anova(model1o, model2o)
anova(model2o, model3o)
anova(model3o, model4o)
anova(model4o, model5o)

## Check model diagnostics. ----
check_heteroscedasticity(model3o)
check_autocorrelation(model3o)
check_normality(model3o)
check_outliers(model3o)
check_collinearity(model3o)
diagnosticPlots <- check_model(model6o, check = c("linearity", "homogeneity","reqq", "qq"))
diagnosticPlots

# Random effects selection. ----
modelVariances <- rep(0, 4)

model3o.0 <- lm(nabFlowT ~ 1 + isConvalesced + 
                 ysfdc + ysfdc:isConvalesced +
                 ysfdc2 + ysfdc2:isConvalesced + 
                 ysfdc3 + ysfdc3:isConvalesced,
               data = filteredData)
modelVariances[1] <- summary(model3o.0)$sigma^2

model3o.1 <- lmer(nabFlowT ~ 1 + isConvalesced + 
                 ysfdc + ysfdc:isConvalesced +
                 ysfdc2 + ysfdc2:isConvalesced + 
                 ysfdc3 + ysfdc3:isConvalesced +                   
                (1 | subject),
               data = filteredData,
               REML = TRUE)
modelVariances[2] <- summary(model3o.1)$sigma^2

model3o.2Uc <- lmer(nabFlowT ~ 1 + isConvalesced + 
                 ysfdc + ysfdc:isConvalesced +
                 ysfdc2 + ysfdc2:isConvalesced + 
                 ysfdc3 + ysfdc3:isConvalesced +
                (1 | subject) + (0 + ysfdc | subject),
               data = filteredData,
               REML = TRUE)
modelVariances[3] <- summary(model3o.2Uc)$sigma^2

model3o.2 <- lmer(nabFlowT ~ 1 + isConvalesced + 
                 ysfdc + ysfdc:isConvalesced +
                 ysfdc2 + ysfdc2:isConvalesced + 
                 ysfdc3 + ysfdc3:isConvalesced +
                (1 + ysfdc | subject),
               data = filteredData,
               REML = TRUE)
modelVariances[4] <- summary(model3o.2)$sigma^2

modelVariances
AIC(model3o.0, model3o.1, model3o.2Uc, model3o.2)

## Check model diagnostics. ----
check_heteroscedasticity(model3o.2)
check_autocorrelation(model3o.2)
check_normality(model3o.2)
check_outliers(model3o.2)
diagnosticPlots <- check_model(model3o.2, check = c("linearity", "homogeneity","reqq", "qq"))
diagnosticPlots

# Test for the relevancy of convalesced status. ----
## Test if convalesced status better explains the data. ----
modelFull <- lmer(nabFlowT ~ poly(ysfdc, 4, raw = TRUE) * isConvalesced +
                    (1 + ysfdc | subject),
                  data = filteredData, 
                  REML = FALSE)

chosenModel <- lmer(nabFlowT ~ poly(ysfdc, 4, raw = TRUE) * isConvalesced +
                    (1 + ysfdc | subject),
                  data = filteredData, 
                  REML = TRUE)

modelReduced <- lmer(nabFlowT ~ poly(ysfdc, 4, raw = TRUE) +
                    (1 + ysfdc | subject),
                  data = filteredData, 
                  REML = FALSE)

anova(modelReduced, modelFull)

chosenModel <- lmer(nabFlowT ~ 1 + isConvalesced + 
                    ysfdc + ysfdc:isConvalesced +
                    ysfdc2 + ysfdc2:isConvalesced + 
                    ysfdc3 + ysfdc3:isConvalesced +
                    (1 + ysfdc | subject),
                  data = filteredData,
                  REML = TRUE)
summary(chosenModel)
model_performance(chosenModel)

### Check diagnostics of chosen model. ----
check_heteroscedasticity(chosenModel)
check_autocorrelation(chosenModel)
check_normality(chosenModel)
check_outliers(chosenModel)
check_collinearity(chosenModel)
diagnosticPlots <- check_model(chosenModel, check = c("linearity", "homogeneity","reqq", "qq"))
diagnosticPlots

## Test if certain interactions contribute to the model. ----
modelF <- lmer(nabFlowT ~ 1 + isConvalesced + 
                    ysfdc + 
                    ysfdc2 + ysfdc2:isConvalesced + 
                    ysfdc3 + 
                    (1 + ysfdc | subject),
                  data = filteredData,
                  REML = FALSE,
                  control = lmerControl(optimizer ="Nelder_Mead"))

modelR <- lmer(nabFlowT ~ 1 + isConvalesced + 
                    ysfdc + 
                    ysfdc2 +  ysfdc2:isConvalesced +
                    ysfdc3 +
                    (1 + ysfdc | subject),
                  data = filteredData,
                  REML = FALSE,
                  control = lmerControl(optimizer ="Nelder_Mead"))
anova(modelR, modelF)

optimalModel <- lmer(nabFlowT ~ 1 + isConvalesced + 
                    ysfdc + 
                    ysfdc2 + ysfdc2:isConvalesced + 
                    ysfdc3 + 
                    (1 + ysfdc | subject),
                  data = filteredData,
                  REML = TRUE)
summary(optimalModel)

## Test if certain interactions contribute to the model. ----
m1 <- lmer(nabFlowT ~ 1 + isConvalesced + 
                    ysfdc + ysfdc:isConvalesced +
                    ysfdc2 + ysfdc2:isConvalesced + 
                    ysfdc3 + ysfdc3:isConvalesced +
                    (1 | subject) + (0 + ysfdc | subject),
                  data = filteredData,
                  REML = FALSE,
                  control = lmerControl(optimizer ="Nelder_Mead"))
summary(m1)

m2 <- lmer(nabFlowT ~ 1 + isConvalesced + 
                    ysfdc + ysfdc:isConvalesced +
                    ysfdc2 + ysfdc2:isConvalesced + 
                    ysfdc3 +
                    (1 | subject) + (0 + ysfdc | subject),
                  data = filteredData,
                  REML = FALSE,
                  control = lmerControl(optimizer ="Nelder_Mead"))
summary(m2)

m3 <- lmer(nabFlowT ~ 1 + isConvalesced + 
                    ysfdc + ysfdc:isConvalesced +
                    ysfdc2 + 
                    ysfdc3 + ysfdc3:isConvalesced +
                    (1 | subject) + (0 + ysfdc | subject),
                  data = filteredData,
                  REML = FALSE,
                  control = lmerControl(optimizer ="Nelder_Mead"))
summary(m3)

m4 <- lmer(nabFlowT ~ 1 + isConvalesced + 
                    ysfdc +
                    ysfdc2 + ysfdc2:isConvalesced + 
                    ysfdc3 +
                    (1 | subject) + (0 + ysfdc | subject),
                  data = filteredData,
                  REML = FALSE,
                  control = lmerControl(optimizer ="Nelder_Mead"))
summary(m4)

AIC(m1, m2, m3, m4)
optimalModel <- m2
summary(optimalModel)

### Check diagnostics of optimal model. ----
check_heteroscedasticity(optimalModel)
check_autocorrelation(optimalModel)
check_normality(optimalModel)
check_outliers(optimalModel)
check_collinearity(optimalModel)
diagnosticPlots <- check_model(optimalModel, check = c("linearity", "homogeneity","reqq", "qq"))
diagnosticPlots

# Plots of model fits. ----
# Testing for confidence interval.

chosenModel <- optimalModel

xPredLimits <- filteredData |> 
  select(ysfdc) |>
  range()

xPred <- seq(xPredLimits[1], xPredLimits[2], length = 100)
ConfData <- data.frame(ysfdc = rep(xPred, 2), 
                       ysfdc2 = rep(xPred^2, 2),
                       ysfdc3 = rep(xPred^3, 2),
                       isConvalesced = rep(c("Yes", "No"), each = 100),
                       isConvalescedLabel = rep(c("Convalesced", "Naive"), each = 100))

ConfData <- ConfData |>
  mutate(Estimate = predict(chosenModel, newdata = ConfData, re.form = ~0))  

myStats <- function(model.star){
  predict(model.star, newdata = ConfData, re.form = ~0)
}

bootObj <- bootMer(chosenModel, FUN = myStats, nsim = 10000)

CI <- confint(bootObj, level = 0.95, ci = 'bca')
colnames(CI) <- c("lwr", "upr")
ConfData <- cbind(ConfData, CI)

ConfData <- ConfData |>
  mutate(EstimatePN = exp(Estimate) / (1 + exp(Estimate))) |>
  mutate(lwrPN = exp(lwr) / (1 + exp(lwr))) |>
  mutate(uprPN = exp(upr) / (1 + exp(upr)))

xTicks = seq(30, 450, by = 30)
ysfdcTicks = (xTicks/365) - mean(filteredData$ysfd)

ggplot(ConfData, aes(x = ysfdc, y = Estimate, color = isConvalescedLabel)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              fill = "gray", alpha = 0.4, linetype = 0, show.legend = FALSE) +
  geom_point(data = filteredData, aes(x = ysfdc, y = nabFlowT, color = isConvalescedLabel)) +
  scale_x_continuous(breaks = ysfdcTicks, labels = xTicks) +
  labs(x = "days since first dose", y = "logit(NAb activity)", color = NULL) +
  scale_color_brewer(palette = "Dark2") +   
  theme_bw()

ggplot(ConfData, aes(x = ysfdc, y = EstimatePN, color = isConvalescedLabel)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwrPN, ymax = uprPN), 
              fill = "gray", alpha = 0.4, linetype = 0, show.legend = FALSE) +
  geom_point(data = filteredData, aes(x = ysfdc, y = nabFlow, color = isConvalescedLabel)) +
  geom_hline(yintercept = 0) +  
  scale_x_continuous(breaks = ysfdcTicks, labels = xTicks) +
  labs(x = "days since first dose", y = "NAb activity (% neutralization)", color = NULL) +
  scale_color_brewer(palette = "Dark2") +   
  theme_bw()
