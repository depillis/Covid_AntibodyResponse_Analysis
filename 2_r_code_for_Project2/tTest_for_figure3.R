library(tidyverse)
library(car)
library(cowplot)

datasim <- read_csv("./simulatedCovariates.txt")

leveneTest(trans_A0_simulated ~ isConvalesced,
           data = datasim)

t.test(trans_A0_simulated ~ isConvalesced,
       data = datasim, 
       var.equal = TRUE)

leveneTest(trans_r1_simulated ~ isConvalesced,
       data = datasim)

t.test(trans_r1_simulated ~ isConvalesced,
       data = datasim, 
       var.equal = TRUE)

leveneTest(trans_r2_simulated ~ isConvalesced,
       data = datasim)

t.test(trans_r2_simulated ~ isConvalesced,
       data = datasim, 
       var.equal = TRUE)

leveneTest(trans_r3_simulated ~ isConvalesced,
       data = datasim)

t.test(trans_r3_simulated ~ isConvalesced,
       data = datasim, 
       var.equal = TRUE)

leveneTest(trans_r4_simulated ~ isConvalesced,
       data = datasim)

t.test(trans_r4_simulated ~ isConvalesced,
       data = datasim, 
       var.equal = FALSE)

pr1 <- ggplot(datasim, aes(x = isConvalesced, y = trans_r1_simulated)) +
        geom_boxplot(fill = "gray") +
        xlab("") +
        ylab(bquote(logit(r[1]))) +
        scale_x_discrete(labels = c("No" = "Naive", "Yes" = "Convalesced")) +
        theme_classic()

pr2 <- ggplot(datasim, aes(x = isConvalesced, y = trans_r2_simulated)) +
        geom_boxplot(fill = "gray") +
        xlab("") +
        ylab(bquote(logit(r[2]))) +
        scale_x_discrete(labels = c("No" = "Naive", "Yes" = "Convalesced")) +
        theme_classic()

pr3 <- ggplot(datasim, aes(x = isConvalesced, y = trans_r3_simulated)) +
        geom_boxplot(fill = "gray") +
        xlab("") +
        ylab(bquote(logit(r[3]))) +
        scale_x_discrete(labels = c("No" = "Naive", "Yes" = "Convalesced")) +
        theme_classic()

pr4 <- ggplot(datasim, aes(x = isConvalesced, y = trans_r4_simulated)) +
        geom_boxplot(fill = "gray") +
        xlab("") +
        ylab(bquote(logit(r[4]))) +
        scale_x_discrete(labels = c("No" = "Naive", "Yes" = "Convalesced")) +
        theme_classic()

pA0 <- ggplot(datasim, aes(x = isConvalesced, y = trans_A0_simulated)) +
        geom_boxplot(fill = "gray") +
        xlab("") +
        ylab(bquote(logit(A[0]))) +
        scale_x_discrete(labels = c("No" = "Naive", "Yes" = "Convalesced")) +
        theme_classic()

# 8 x 3 landscape
plot_grid(pA0, pr1, pr2, pr3, pr4, labels = c("(i)", "(ii)", "(iii)", "(iv)", "(v)"), 
          ncol = 3,
          hjust = 0.03)
#--------------------------------------
data <- read_csv("./covariates.txt")

leveneTest(trans_A0_mean ~ isConvalesced,
           data = data)

t.test(trans_A0_mean ~ isConvalesced,
       data = data, 
       var.equal = FALSE)

leveneTest(trans_r1_mean ~ isConvalesced,
           data = data)

t.test(trans_r1_mean ~ isConvalesced,
       data = data, 
       var.equal = TRUE)

leveneTest(trans_r2_mean ~ isConvalesced,
           data = data)

t.test(trans_r2_mean ~ isConvalesced,
       data = data, 
       var.equal = TRUE)

leveneTest(trans_r3_mean ~ isConvalesced,
           data = data)

t.test(trans_r3_mean ~ isConvalesced,
       data = data, 
       var.equal = TRUE)

leveneTest(trans_r4_mean ~ isConvalesced,
           data = data)

t.test(trans_r4_mean ~ isConvalesced,
       data = data, 
       var.equal = FALSE)

ggplot(data, aes(x = isConvalesced, y = trans_r4_mode)) +
        geom_boxplot()
#---------------------------------------
#data <- read_csv("modeCovariates.txt")

leveneTest(trans_A0_mode ~ isConvalesced,
       data = data)

t.test(trans_A0_mode ~ isConvalesced,
       data = data, 
       var.equal = FALSE)

leveneTest(trans_r1_mode ~ isConvalesced,
           data = data)

t.test(trans_r1_mode ~ isConvalesced,
       data = data, 
       var.equal = TRUE)


leveneTest(trans_r2_mode ~ isConvalesced,
       data = data)

t.test(trans_r2_mode ~ isConvalesced,
       data = data, 
       var.equal = TRUE)

leveneTest(trans_r3_mode ~ isConvalesced,
       data = data)

t.test(trans_r3_mode ~ isConvalesced,
       data = data, 
       var.equal = TRUE)

leveneTest(trans_r4_mode ~ isConvalesced,
       data = data)

t.test(trans_r4_mode ~ isConvalesced,
       data = data, 
       var.equal = FALSE)





