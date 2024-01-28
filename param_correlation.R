# look at correlation between parameters from iEpi results

library(ggplot2)
library(Hmisc)
library(corrplot)

# load one of the chains
load("iEpi_results3/Bayes_joint_1.RData")
#load("iEpi_results3/Bayes_joint_2.RData")
load("iEpi_results3/MLE_joint_10.RData")

plot(params[,"xi"], type = "l")

# while we are at this... check geweke z scores
library(coda)
geweke.diag(params[complete.cases(params),][-c(1:200),], 
            frac1 = 0.1, frac2 = 0.4)
# beta    gamma alpha.SS alpha.SI alpha.II omega.SS omega.SI omega.II       xi 
# -1.0915   1.2492   0.5163  -0.3329  -0.3851   0.9878  -0.9681  -1.0407   1.4441 
geweke.plot(as.mcmc(params[complete.cases(params),][-c(1:200),]), 
            frac1 = 0.3, frac2 = 0.2)

plot(params[-c(1:200),"xi"], type = "l", ylab = "xi samples")


# re-arrange columns and rename
params = params[,c(1,9,2:8)]

# corr plot
cormat = cor(params, use = "complete.obs")
corrplot(cormat, method = "shade", 
         diag = FALSE, type = "upper",
         tl.col = "gray30", tl.srt = 45)

