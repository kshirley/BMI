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

# Set up a survival object:
S <- Surv(time = data[, "age"], 
          time2 = data[, "age"] + data[, "time"], 
          event = data[, "status"])

# get unique values for each sex:
sex <- su(data[, "sex"])

# Load BMI_2 mean and sd:
load("ss.RData")
bmi.mean <- ss$bmi.mean
bmi.sd <- ss$bmi.sd



### Basic Models:

# In this section, we call the functon coxph() a number of times (for both men 
# and women), with each call requiring roughly 30 seconds to two minutes 
# to execute.

# M0: (null model)
M0 <- as.list(rep(NA, 2))
for (i in 1:2) {
  sel <- data[, "sex"] == sex[i]
  M0[[i]] <- coxph(S ~ 1, data = data, subset = sel)
}  

# M1: (Adams 2006 variables)
M1 <- as.list(rep(NA, 2))
adams.formula <- as.formula(S ~ bmis + I(bmis^2) + I(bmis^3) + race + edu + 
                            smoking + physical + alcohol.factor)
for (i in 1:2) {
  sel <- data[, "sex"] == sex[i]
  M1[[i]] <- coxph(adams.formula, data = data, subset = sel)
}  

# M2: (All main effects)
M2 <- as.list(rep(NA, 2))
all.variables.formula <- as.formula(S ~ bmis + I(bmis^2) + I(bmis^3) + race + 
                                    edu + smoking + physical + alcohol.factor + 
                                    health + marriage + diabetes + ages + 
                                    heights)
for (i in 1:2) {
  sel <- data[, "sex"] == sex[i]
  M2[[i]] <- coxph(all.variables.formula, data = data, subset = sel, 
                   model = TRUE)
}  

# Gather the Xbeta terms from Model M2 for each of the seven categorical
# variables. These are the so-called 'tied-together' variables:
terms <- as.list(rep(NA, 2))
for (i in 1:2) terms[[i]] <- predict(M2[[i]], type = "terms")

# Set up a data.frame with the 4 numerical predictors and the seven categorical 
# predictors, now expressed as numerical 'tied-together' variables:
dt <- as.list(rep(NA, 2))
for (i in 1:2) {
  sel <- data[, "sex"] == sex[i]
  dt[[i]] <- data.frame(bmis = data[sel, "bmis"], 
                        diabetes = as.numeric(data[sel, "diabetes"] == "2=yes"), 
                        ages = data[sel, "ages"], 
                        heights = data[sel, "heights"], 
                        terms[[i]][, 4:10])
}

### Interaction models

# M4: The model with two-way interactions, squared terms, but excluding 
# BMI*height:
M4 <- as.list(rep(NA, 2))
t1 <- Sys.time()
for (i in 1:2) {
  sel <- data[, "sex"] == sex[i]
  M4[[i]] <- coxph(S[sel] ~ I(bmis^2) + I(bmis^3) + (bmis + ages + diabetes + 
                   race + edu + smoking + physical + alcohol.factor + health + 
                   marriage)^2 + I(race^2) + I(edu^2) + I(smoking^2) + 
                   I(physical^2) + I(alcohol.factor^2) + I(health^2) + 
                   I(marriage^2) + I(ages^2) + I(heights^2) + 
                   heights*(ages + diabetes + race + edu + smoking + physical + 
                   alcohol.factor + health + marriage), data = dt[[i]])
}
t2 <- Sys.time()
t2 - t1  # 1.67 minutes on macbook air

# M3: Now put in the BMI*height variable to see if it's significant:
M3 <- as.list(rep(NA, 2))
t1 <- Sys.time()
for (i in 1:2) {
  sel <- data[, "sex"] == sex[i]
  M3[[i]] <- coxph(S[sel] ~ I(bmis^2) + I(bmis^3) + (bmis + ages + diabetes + 
                   race + edu + smoking + physical + alcohol.factor + health + 
                   marriage)^2 + I(race^2) + I(edu^2) + I(smoking^2) + 
                   I(physical^2) + I(alcohol.factor^2) + I(health^2) + 
                   I(marriage^2) + I(ages^2) + I(heights^2) + 
                   heights*(ages + diabetes + race + edu + smoking + physical + 
                   alcohol.factor + health + marriage) + bmis*heights, 
                   data = dt[[i]], model = TRUE)
}
t2 <- Sys.time()
t2 - t1  # 1.68 minutes on macbook air

# Grid search over alpha (the exponent in BMI) for M4 (the model that excludes 
# the BMI*Height interaction)
alpha.seq <- seq(0.1, 3.0, 0.1)
la <- length(alpha.seq)

ll.seq <- as.list(rep(NA, 2))
coef.seq <- as.list(rep(NA, 2))
t1 <- Sys.time()
for (s in 1:2) {
  sel <- data[, "sex"] == sex[s]
  ll.seq[[s]] <- numeric(la)
  coef.seq[[s]] <- as.list(rep(NA, la))
  for (i in 1:la) {
    print(i)
    tmp <- dt[[s]]
    bmi <- data[sel, "weight"]*0.453592/(data[sel, "height"]^alpha.seq[i])
    tmp[, "bmis"] <- (bmi - mean(bmi))/sd(bmi)

    # Fit the cox model
    f <- coxph(S[sel] ~ I(bmis^2) + I(bmis^3) + (bmis + ages + diabetes + 
               race + edu + smoking + physical + alcohol.factor + health + 
               marriage)^2 + I(race^2) + I(edu^2) + I(smoking^2) + 
               I(physical^2) + I(alcohol.factor^2) + I(health^2) + 
               I(marriage^2) + I(ages^2) + I(heights^2) + 
               heights*(ages + diabetes + race + edu + smoking + physical + 
               alcohol.factor + health + marriage), data = tmp)
    ll.seq[[s]][i] <- f$loglik[2]
    coef.seq[[s]][[i]] <- summary(f)$coefficients
  }
}
t2 <- Sys.time()
t2 - t1  # 40 minutes on macbook air
alpha.seq.ix <- list(ll.seq=ll.seq, coef.seq=coef.seq)

# Re-run the fit for the optimal alpha for men and women:
M5 <- as.list(rep(NA, 2))
opt.alpha <- c(alpha.seq[which.max(ll.seq[[1]])], 
               alpha.seq[which.max(ll.seq[[2]])])
# 1.1 and 1.3 for men and women, respectively

for (s in 1:2) {
  sel <- data[, "sex"] == sex[s]
  tmp <- dt[[s]]
  bmi <- data[sel, "weight"]*0.453592/(data[sel, "height"]^opt.alpha[s])
  tmp[, "bmis"] <- (bmi - mean(bmi))/sd(bmi)

  # Fit the cox model
  M5[[s]] <- coxph(S[sel] ~ I(bmis^2) + I(bmis^3) + (bmis + ages + diabetes + 
                   race + edu + smoking + physical + alcohol.factor + health + 
                   marriage)^2 + I(race^2) + I(edu^2) + I(smoking^2) + 
                   I(physical^2) + I(alcohol.factor^2) + I(health^2) + 
                   I(marriage^2) + I(ages^2) + I(heights^2) + 
                   heights*(ages + diabetes + race + edu + smoking + physical + 
                   alcohol.factor + health + marriage), data = tmp)
}

# Last, check to see that, in the context of the optimal alpha fits, the 
# inclusion of BMI*height is not significant:
M6 <- as.list(rep(NA, 2))
for (s in 1:2) {
  sel <- data[, "sex"] == sex[s]
  tmp <- dt[[s]]
  bmi <- data[sel, "weight"]*0.453592/(data[sel, "height"]^opt.alpha[s])
  tmp[, "bmis"] <- (bmi - mean(bmi))/sd(bmi)

  # Fit the cox model
  M6[[s]] <- coxph(S[sel] ~ I(bmis^2) + I(bmis^3) + (bmis + ages + diabetes + 
                   race + edu + smoking + physical + alcohol.factor + health + 
                   marriage)^2 + I(race^2) + I(edu^2) + I(smoking^2) + 
                   I(physical^2) + I(alcohol.factor^2) + I(health^2) + 
                   I(marriage^2) + I(ages^2) + I(heights^2) + 
                   heights*(ages + diabetes + race + edu + smoking + physical + 
                   alcohol.factor + health + marriage) + bmis*heights, 
                   data = tmp)
}


# Save all 7 models + the features from M2 + 
# the loglik sequence using different exponents
save(M0, file = paste0(data.path, "chronically-ill/M0.RData"))
save(M1, file = paste0(data.path, "chronically-ill/M1.RData"))
save(M2, file = paste0(data.path, "chronically-ill/M2.RData"))
save(M3, file = paste0(data.path, "chronically-ill/M3.RData"))
save(M4, file = paste0(data.path, "chronically-ill/M4.RData"))
save(M5, file = paste0(data.path, "chronically-ill/M5.RData"))
save(M6, file = paste0(data.path, "chronically-ill/M6.RData"))

# input for M3:
save(dt, file = paste0(data.path, "chronically-ill/dt.RData"))

# loglik sequence:
save(alpha.seq.ix, file = "chronically-ill/alpha_seq_ix.RData")




