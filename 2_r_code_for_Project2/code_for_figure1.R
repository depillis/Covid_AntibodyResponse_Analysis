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
daycutoff <- 0

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
