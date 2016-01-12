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

# load a few objects from the model fitting script:
load(file = paste0(data.path, "chronically-ill/M2.RData"))
load(file = paste0(data.path, "chronically-ill/dt.RData"))
load(file = paste0(data.path, "chronically-ill/M3.RData"))


### POB Calculations

# Compute POB(2) for each person in the data based on Model M3:
beta <- matrix(0, 13, 2)
for (i in 1:2) beta[, i] <- coef(M3[[i]])[c(3, 1, 2, 23:31, 77)]
rownames(beta) <- names(coef(M3[[i]])[c(3, 1, 2, 23:31, 77)])

# set up quadratic formula to compute point estimate of optimal BMI:
a <- numeric(2)
b <- numeric(2)
c <- numeric(n)
for (i in 1:2) {
  # a*x^2 + b*x + c is the derivative of the cubic risk function of bmi under M3:
  a[i] <- 3*beta[3, i]
  b[i] <- 2*beta[2, i]
  c[data[, "sex"] == sex[i]] <- beta[1, i] + 
                                beta[4, i]*dt[[i]][, "ages"] + 
                                beta[5, i]*dt[[i]][, "diabetes"] + 
                                beta[6, i]*dt[[i]][, "race"] + 
                                beta[7, i]*dt[[i]][, "edu"] + 
                                beta[8, i]*dt[[i]][, "smoking"] + 
                                beta[9, i]*dt[[i]][, "physical"] + 
                                beta[10, i]*dt[[i]][, "alcohol.factor"] + 
                                beta[11, i]*dt[[i]][, "health"] + 
                                beta[12, i]*dt[[i]][, "marriage"] + 
                                beta[13, i]*dt[[i]][, "heights"]
}

# Solve for the roots:
pob <- numeric(n)
for (i in 1:2) {
  pob[data[, "sex"] == sex[i]] <- (-b[i] + 
    sqrt(b[i]^2 - 4*a[i]*c[data[, "sex"] == sex[i]]))/(2*a[i])
}
# A few of these are NA because there is no local minimum in [15, 50].
# 84 of them, to be specific (79 men and 5 women)
# table(data[, "sex"], is.na(pob))

# convert to original scale of bmi:
pob.scaled <- pob*bmi.sd + bmi.mean

# Mean BMI_2 and POB_2 by sex:
bmi.mean.by.sex <- aggregate(data[, "bmi"], by=list(data[, "sex"]), mean)
bmi.sd.by.sex <- aggregate(data[, "bmi"], by = list(data[, "sex"]), sd)

round(bmi.mean.by.sex[, 2], 1)
round(bmi.sd.by.sex[, 2], 2)

pob.mean.by.sex <- aggregate(pob.scaled, by = list(data[, "sex"]), mean, 
                             na.rm = TRUE)
pob.sd.by.sex <- aggregate(pob.scaled, by=list(data[, "sex"]), sd, na.rm=TRUE)

round(pob.mean.by.sex[, 2], 1)
round(pob.sd.by.sex[, 2], 2)

# Derivative test to see where the slope of the relative risk curve might be 
# zero for each individual:
bmi.test <- seq(15, 50, 0.1) # grid of values for derivative test
lt <- length(bmi.test)

# Covariance matrix of the coefficients from Model M3:
Sigma <- as.list(rep(NA, 2))
for (s in 1:2) Sigma[[s]] <- M3[[s]]$var[c(3, 1, 2, 23:31, 77), 
                                         c(3, 1, 2, 23:31, 77)]

### Derivative test:
t.mat <- matrix(0, n, lt)
t1 <- Sys.time()
for (s in 1:2) {
  # subset by sex:
  sel <- data[, "sex"] == sex[s]

  # matrix to hold the linear contrasts:
  d <- cbind(rep(1, sum(sel)), 
             numeric(sum(sel)), 
             numeric(sum(sel)), 
             as.matrix(dt[[s]][, c("ages", "diabetes", "race", "edu", "smoking",
             "physical", "alcohol.factor", "health", "marriage", "heights")]))

  # Loop through the alternative BMI values and run the test, and get the t-stat:
  for (i in 1:lt) {
    if (i %% 5 == 0) print(paste0(i, "/", lt))

    # Insert into d, the linear contrast, the current bmi value from this
    # iteration of the loop
    z <- (rep.int(bmi.test[i], sum(sel)) - bmi.mean)/bmi.sd
    d[, 2] <- 2*z
    d[, 3] <- 3*z^2
    diff <- d %*% beta[, s]
    se <- d %*% Sigma[[s]] * d
    se <- sqrt(apply(se, 1, sum))

    # compute the t-statistic (on sum(sel) - 13 - 1 df)
    t.mat[sel, i] <- diff/se
  }
}
t2 <- Sys.time()
t2 - t1  # ~ 12 minutes on laptop

# For each person, compute the end points of the confidence interval:
ci.orig <- abs(t.mat) < 1.96

# detect where one value of the sequence is in the interval, and neighboring 
# values are not
edge <- ci.orig[, 2:lt] - ci.orig[, 1:(lt - 1)]

# compute number of endpoints in the intervals:
num.edges <- apply(abs(edge), 1, sum)

table(num.edges, data[, "sex"])

# allocate vectors for lower bound, upper bound (where it exists), and second
# lower bound:
lower <- numeric(n)
upper <- numeric(n)
lower2 <- numeric(n)
upper2 <- numeric(n) # if they have 4 edges

# everybody with at least one edge has a lower bound:
lower[num.edges != 0] <- bmi.test[max.col(edge[num.edges != 0, ], 
                                          ties.method="first") + 1]

# upper bound depends on whether there is only 1 edge:
upper[num.edges == 1] <- bmi.test[lt]
upper[num.edges > 1] <- bmi.test[max.col(-edge[num.edges > 1, ], 
                                         ties.method="first")]

# handle those with an additional lower bound:
lower2[num.edges > 2] <- bmi.test[max.col(edge[num.edges > 2, ], 
                                          ties.method="last") + 1]

# handle those with a final upper edge:
upper2[num.edges == 4] <- bmi.test[max.col(-edge[num.edges == 4, , drop = FALSE], 
                                           ties.method = "last")]

# pull them together into an object to save:
ci <- cbind(lower, upper, lower2, upper2)

# save a couple objects for next script:
save(pob.scaled, file = "chronically-ill/pob_scaled.RData")
save(ci, file = "chronically-ill/ci_derivative.RData")










