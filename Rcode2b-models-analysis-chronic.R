#2345678901234567890123456789012345678901234567890123456789012345678901234567890

# Clear the current workspace:
rm(list=ls())
gc()

# Set the working directory:
setwd("~/Git/Obesity/")

# set the data path:
data.path <- "~/Stats/BMI/data/"

# Load the survival library to fit cox proportional hazards models:
library(survival)

# A few utility functions:
lu <- function(x) length(unique(x))
su <- function(x) sort(unique(x))

# load the tidy data (which includes all with chronic illness)
load(file = paste0(data.path, "aarp-new.RData"))
n <- dim(data)[1]

# change death time of zero to something very small and keep them in the data.
# these seem to be people who died the same day they filled out the survey, 
# according to the data.
data[data[, "time"] == 0, "time"] <- 0.00001

# Sub-select only non-chronically ill men and women:
#data <- data[data[, "chronic"] == 0, ]  # removes 133,178 individuals
#n <- dim(data)[1]

# get unique values for each sex:
sex <- su(data[, "sex"])

# Store BMI_2 mean and sd:
load("ss.RData")
bmi.mean <- ss$bmi.mean
bmi.sd <- ss$bmi.sd

# Load the models:
load(file = paste0(data.path, "chronically-ill/M0.RData"))
load(file = paste0(data.path, "chronically-ill/M1.RData"))
load(file = paste0(data.path, "chronically-ill/M2.RData"))
load(file = paste0(data.path, "chronically-ill/M3.RData"))
load(file = paste0(data.path, "chronically-ill/M4.RData"))
load(file = paste0(data.path, "chronically-ill/M5.RData"))
load(file = paste0(data.path, "chronically-ill/M6.RData"))

# load the features from M2:
load(file = paste0(data.path, "chronically-ill/dt.RData"))

# load the loglik sequence from changing the exponent:
load(file = "chronically-ill/alpha_seq_ix.RData")
alpha.seq <- seq(0.1, 3.0, 0.1)
la <- length(alpha.seq)
ll.seq <- alpha.seq.ix[[1]]
coef.seq <- alpha.seq.ix[[2]]

# compute confidence intervals for optimal alpha:
ll.ci <- matrix(0, la, 2)
for (i in 1:2) ll.ci[, i] <- as.numeric(ll.seq[[i]] - max(ll.seq[[i]]) >= -2)
colnames(ll.ci) <- c("men", "women")
data.frame(alpha = alpha.seq, ll.ci)


# p-values comparing 2.0 to 1.0:
2*pnorm(-sqrt(2*(ll.seq[[1]][10] - ll.seq[[1]][20])))
2*pnorm(-sqrt(2*(ll.seq[[2]][10] - ll.seq[[2]][20])))

# p-values comparing 2.0 to optimal values:
2*pnorm(-sqrt(2*(max(ll.seq[[1]]) - ll.seq[[1]][20])))
2*pnorm(-sqrt(2*(max(ll.seq[[2]]) - ll.seq[[2]][20])))

# increase in log-likelihood by using optimal vs. 2.0
ll.seq[[1]][which.max(ll.seq[[1]])] - ll.seq[[1]][20] # men: 17.58
ll.seq[[2]][which.max(ll.seq[[2]])] - ll.seq[[2]][20] # women: 4.43

opt.alpha <- alpha.seq[c(which.max(ll.seq[[1]]), which.max(ll.seq[[2]]))]
# 1.0 and 1.4



