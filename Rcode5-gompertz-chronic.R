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
library(dplyr)

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

# get unique values for each sex:
sex <- su(data[, "sex"])
men <- data[, "sex"] == sex[1]
women <- data[, "sex"] == sex[2]

# Load BMI_2 mean and sd:
load("ss.RData")
bmi.mean <- ss$bmi.mean
bmi.sd <- ss$bmi.sd

# load a few objects from the model fitting script:
load(file = paste0(data.path, "chronically-ill/dt.RData"))
load(file = paste0(data.path, "chronically-ill/M3.RData"))

# pob.scaled:
load(file = "chronically-ill/pob_scaled.RData")
pob <- (pob.scaled - bmi.mean)/bmi.sd

# ci - the matrix containing confidence interval information
load(file = "chronically-ill/ci_derivative.RData")

# column names for the confidence interval matrix:
colnames(ci) <- c("lower", "upper", "lower2", "upper2")

# reproduce the variable num.edges, which characterizes the interval for each 
# respondent:
num.edges <- rep(2, n)
num.edges[ci[, 2] == 50] <- 1
num.edges[ci[, 3] != 0] <- 3
num.edges[ci[, 1] == 0] <- 0
num.edges[ci[, 4] != 0] <- 4

# Focus just on those with a lower and an upper value:
sel.edges <- num.edges >= 2
# excludes 117 men and 266 women (383 total)




### Compare relative risk decrease when
### (1) moving to nearest endpoint of [18.5, 24.9] interval
### (2) moving to nearest endpoint of POB interval
lp <- matrix(0, n, 3)  # lp for linear predictor
colnames(lp) <- c("current", "who", "pob")

intercept <- vector("list", 2)
for (i in 1:2) {
  intercept[[i]] <- attr(predict(M3[[i]], type = "terms"), "constant")
}

# Compute Xbeta (linear predictors) for M3 at current BMI:
for (i in 1:2) {
  lp[data[, "sex"] == sex[i], 1] <- predict(M3[[i]], newdata = dt[[i]], 
                                            type = "lp") + intercept[[i]]
}

# Next, the [18.5, 24.9] interval:
bmi.who <- data[, "bmi"]
# "WHO" for "World Health Organization" recommended interval
bmi.who[data[, "bmi"] > 25] <- 25
bmi.who[data[, "bmi"] < 18.5] <- 18.5

# Now, the POB interval:
bmi.pob <- data[, "bmi"]
too.high <- data[sel.edges, "bmi"] > ci[sel.edges, 2]
bmi.pob[sel.edges][too.high] <- ci[sel.edges, 2][too.high]
too.low <- data[sel.edges, "bmi"] < ci[sel.edges, 1]
bmi.pob[sel.edges][too.low] <- ci[sel.edges, 1][too.low]

# Look at a few examples:
dd <- data.frame(bmi = data[, "bmi"], pob.scaled, ci[, 1:2], bmi.who, bmi.pob)
dd[1:10, ]

# Measure relative risk for WHO guideline:
for (i in 1:2) {
  nd <- dt[[i]]
  nd[, "bmis"] <- (bmi.who[data[, "sex"] == sex[i]] - bmi.mean)/bmi.sd
  lp[data[, "sex"] == sex[i], 2] <- predict(M3[[i]], newdata = nd, type = "lp") + 
                                    intercept[[i]]
}

# Measure relative risk for POB intervals:
for (i in 1:2) {
  nd <- dt[[i]]
  nd[, "bmis"] <- (bmi.pob[data[, "sex"] == sex[i]] - bmi.mean)/bmi.sd
  lp[data[, "sex"] == sex[i], 3] <- predict(M3[[i]], newdata = nd, type = "lp") + 
                                    intercept[[i]]
}


### Let's look at the baseline hazard for fit:
base <- vector("list", 2)
for (i in 1:2) {
  base[[i]] <- basehaz(M3[[i]], centered = FALSE)
}

# check: using centered = FALSE in the basehaz() function 
# and adding in the intercept explicitly into the lp object above using
# the predict.coxph(..., type = "lp")
# results in the correct cumulative hazard for each individual.

# random man:
j <- 123
all.equal(survfit(M3[[1]], newdata = dt[[1]][j, ])$cumhaz, 
          base[[1]][, 1] * exp(lp[men, 1][j]))
# TRUE

# random woman:
j <- 234
all.equal(survfit(M3[[2]], newdata = dt[[2]][j, ])$cumhaz, 
          base[[2]][, 1] * exp(lp[women, 1][j]))
# TRUE

# Look at an age sequence
# trimmed at both ends to avoid noisy, small data bins that would result
# from using range(base[[i]][, 2])
age.seq <- 53:80
l.age <- length(age.seq)

# Save the cumulative hazard rates for men and women:
h <- vector("list", 2)  # cumulative hazard
h0 <- vector("list", 2)  # empirical hazard rate
for (j in 1:2) {
  h[[j]] <- numeric(l.age)
  for (i in 1:l.age) {
    h[[j]][i] <- base[[j]][min(which(base[[j]][, 2] > age.seq[i])), 1]
  }
  h0[[j]] <- diff(h[[j]])
}

# Alternatively, use hazards from life tables from NAtional vital stats report
# http://www.cdc.gov/nchs/data/nvsr/nvsr54/nvsr54_14.pdf

#lt1 <- read.csv("lifetable-2003-men.csv")  # table 2 from appendix
#h0[[1]] <- lt1[54:80, 2]

#lt2 <- read.csv("lifetable-2003-women.csv")  # table 3 from appendix
#h0[[2]] <- lt2[54:80, 2]

# cumulative hazard plot:
par(mfrow = c(1, 1))
plot(age.seq, h[[1]], las = 1)
points(age.seq, h[[2]], col = 2)

# empirical hazard plot:
plot(age.seq[-l.age], h0[[1]], type="l", las=1, xlab="Age", ylab="Hazard rate", 
     ylim = c(0, max(h0[[1]])))
points(age.seq[-l.age], h0[[1]])
lines(age.seq[-l.age], h0[[2]], col = 2)
points(age.seq[-l.age], h0[[2]], col = 2)
abline(h = 0, col = gray(0.6))

# regression to estimate gompertz parameters:
f <- vector("list", 2)
x.age <- age.seq[-l.age]
for (i in 1:2) {
  y.haz <- log(h0[[i]])
  f[[i]] <- lm(y.haz ~ x.age)
}

# plot the curve through the points:
plot(age.seq[-l.age], h0[[1]], type="l", las=1, xlab="Age", ylab="Hazard rate")
points(age.seq[-l.age], h0[[1]])
lines(age.seq[-l.age], exp(predict(f[[1]])), col = 4, lty = 2)
# a good fit

# women
plot(age.seq[-l.age], h0[[2]], type="l", las=1, xlab="Age", ylab="Hazard rate")
points(age.seq[-l.age], h0[[2]])
lines(age.seq[-l.age], exp(predict(f[[2]])), col = 4, lty = 2)
# a good fit


# log version:
par(mfrow = c(1, 2))
for (i in 1:2) {
  plot(age.seq[-l.age], log(h0[[i]]), type="l", las=1, xlab="Age", 
       ylab="log(Hazard rate)")
  points(age.seq[-l.age], log(h0[[i]]))
  lines(age.seq[-l.age], predict(f[[i]]), col = 4, lty = 2)
}


# Get the GSL library for the incomplete gamma function:
library(gsl)

# Function for the mean residual life calculation:
# from V. Poynor master's thesis:
# https://users.soe.ucsc.edu/~vpoynor/MastersRevised.pdf
mgom <- function(t, a, b) {
  upper.incomplete.gamma <- gamma_inc(a = 0, x = a/b*exp(b*t))
  exp(a/b*exp(b*t))*(1/b)*upper.incomplete.gamma
}

# men (baseline)
mgom(t = 51:70, a = exp(coef(f[[1]])[1]), b = coef(f[[1]])[2]) + 51:70

# women
mgom(t = 51:70, a = exp(coef(f[[2]])[1]), b = coef(f[[2]])[2]) + 51:70


# compute mean residual life for three different xbeta values:
mrl <- matrix(0, n, 3)
for (i in 1:3) {
  for (j in 1:2) {
  	sel <- data[, "sex"] == sex[j]
    mrl[sel, i] <- mgom(t = data[sel, "age"], 
                        a = exp(coef(f[[j]])[1])*exp(lp[sel, i]), 
                        b = coef(f[[j]])[2])
  }
}
colnames(mrl) <- colnames(lp)


# Look at a few data points:
xx <- data.frame(sex = data[, "sex"], 
                 age = data[, "age"], 
                 bmi = data[, "bmi"], 
                 mrl, 
                 stringsAsFactors = FALSE)


pdf(file = "fig_lifespan-density-chronic.pdf", width = 8, height = 6)
par(mfrow = c(1, 1))
plot(density(xx[men, "age"] + xx[men, "current"]), las = 1, 
     main = "Expected Lifespan (years)")
lines(density(xx[women, "age"] + xx[women, "current"]), col = 2)
men.mean <- mean(xx[men, "age"] + xx[men, "current"])
women.mean <- mean(xx[women, "age"] + xx[women, "current"])
abline(v = men.mean, lty = 2, col = 1)
abline(v = women.mean, lty = 2, col = 2)
legend("topleft", inset = 0.01, col = 1:2, lwd = 1, 
       legend = paste(c("Men", "Women"), round(c(men.mean, women.mean), 1)))
dev.off()

# still very high expected lifespans.
# could be result of:
#   (1) AARP filter (selection bias)
#   (2) 16% response rate (non-response bias) and



# gather MRL and sex into one data.frame, and filter to only the subjects
# with two or more edges in their POB intervals:
mm <- data.frame(mrl[sel.edges, ], sex = data[sel.edges, "sex"])

mean.who <- mm %>% group_by(sex) %>% summarize(a = 12*mean(who - current))
mean.pob <- mm %>% group_by(sex) %>% summarize(a = 12*mean(pob - current))

sd.who <- mm %>% group_by(sex) %>% summarize(a = 12*sd(who - current))
sd.pob <- mm %>% group_by(sex) %>% summarize(a = 12*sd(pob - current))


# Look at results:
mean.who
round(mean.who[, 2], 1)

sd.who
round(sd.who[, 2], 1)

mean.pob
round(mean.pob[, 2], 1)

sd.pob
round(sd.pob[, 2], 1)






