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
data <- data[data[, "chronic"] == 0, ]  # removes 133,178 individuals
n <- dim(data)[1]

# get unique values for each sex:
sex <- su(data[, "sex"])

# Store BMI_2 mean and sd:
load("ss.RData")
bmi.mean <- ss$bmi.mean
bmi.sd <- ss$bmi.sd

# Load the models:
load(file = paste0(data.path, "M0.RData"))
load(file = paste0(data.path, "M1.RData"))
load(file = paste0(data.path, "M2.RData"))
load(file = paste0(data.path, "M3.RData"))
load(file = paste0(data.path, "M4.RData"))
load(file = paste0(data.path, "M5.RData"))
load(file = paste0(data.path, "M6.RData"))

# load the features from M2:
load(file = paste0(data.path, "dt.RData"))

# load the loglik sequence from changing the exponent:
load(file = "alpha_seq_ix.RData")
alpha.seq <- seq(0.1, 3.0, 0.1)
la <- length(alpha.seq)
ll.seq <- alpha.seq.ix[[1]]
coef.seq <- alpha.seq.ix[[2]]

# measure correlations between height and bmi:
for (i in 1:2) {
  xx <- cor(data[data[, "sex"] == sex[i], "height"], 
            data[data[, "sex"] == sex[i], "bmi"])
  print(xx)
}


# Look at BMI interaction effects from M3:
cbind(round(summary(M3[[1]])$coefficients, 3)[c(23:31, 77), 4:5], 
      round(summary(M3[[2]])$coefficients, 3)[c(23:31, 77), 4:5])

# compute confidence intervals for optimal alpha:
ll.ci <- matrix(0, la, 2)
for (i in 1:2) ll.ci[, i] <- as.numeric(ll.seq[[i]] - max(ll.seq[[i]]) >= -2)
colnames(ll.ci) <- c("men", "women")
data.frame(alpha = alpha.seq, ll.ci)


# Paper Figure 6 color version
pdf(file = "fig_optimal_alpha_M5.pdf", width=6, height=6)
plot(alpha.seq, ll.seq[[1]] - max(ll.seq[[1]]), las = 1, 
     main = "Optimal BMI Exponent", ylab="LL - max(LL)", xlab="BMI exponent", 
     type="n")
lines(alpha.seq[ll.ci[, 2] == 1], 
      ll.seq[[2]][ll.ci[, 2] == 1] - max(ll.seq[[2]]), 
      col = "tomato", lwd = 6)
lines(alpha.seq[ll.ci[, 1] == 1], 
      ll.seq[[1]][ll.ci[, 1] == 1] - max(ll.seq[[1]]), 
      col = gray(0.3), lwd = 6)
points(alpha.seq, ll.seq[[1]] - max(ll.seq[[1]]))
points(alpha.seq, ll.seq[[2]] - max(ll.seq[[2]]), col=2)
lines(alpha.seq, ll.seq[[1]] - max(ll.seq[[1]]))
lines(alpha.seq, ll.seq[[2]] - max(ll.seq[[2]]), col=2)
abline(v = alpha.seq[which.max(ll.seq[[1]])], lty=2)
abline(v = alpha.seq[which.max(ll.seq[[2]])], lty=2, col=2)
abline(v = 2, col=gray(0.6))
legend("bottomleft", inset=0.02, border=NULL, col=c(1, 2), 
       legend = c("Men", "Women"), lwd=5)
dev.off()


# Paper Figure 6 - Black and White version:
pdf(file = "fig_optimal_alpha_M5_blackwhite.pdf", width=6, height=6)
plot(alpha.seq, ll.seq[[1]] - max(ll.seq[[1]]), las=1, 
     main="Optimal BMI Exponent", ylab="LL - max(LL)", xlab="BMI exponent", 
     type="n")
lines(alpha.seq[ll.ci[, 2] == 1], 
      ll.seq[[2]][ll.ci[, 2] == 1] - max(ll.seq[[2]]), 
      col=1, lwd=6)
lines(alpha.seq[ll.ci[, 1] == 1], 
      ll.seq[[1]][ll.ci[, 1] == 1] - max(ll.seq[[1]]), 
      col=gray(0.6), lwd=6)
points(alpha.seq, ll.seq[[1]] - max(ll.seq[[1]]), col=gray(0.6))
points(alpha.seq, ll.seq[[2]] - max(ll.seq[[2]]), col=1)
lines(alpha.seq, ll.seq[[1]] - max(ll.seq[[1]]), col=gray(0.6))
lines(alpha.seq, ll.seq[[2]] - max(ll.seq[[2]]), col=1)
abline(v = alpha.seq[which.max(ll.seq[[1]])], lty=2, col=gray(0.6))
abline(v = alpha.seq[which.max(ll.seq[[2]])], lty=2, col=1)
abline(v = 2, col=gray(0.4), lty=1)
legend("bottomleft", inset=0.02, border=NULL, col=c(gray(0.6), 1), 
       legend = c("Men", "Women"), lwd=5)
dev.off()

# p-values comparing 2.0 to 1.0:
2*pnorm(-sqrt(2*(ll.seq[[1]][10] - ll.seq[[1]][20])))
2*pnorm(-sqrt(2*(ll.seq[[2]][10] - ll.seq[[2]][20])))

# p-values comparing 2.0 to optimal values:
2*pnorm(-sqrt(2*(max(ll.seq[[1]]) - ll.seq[[1]][20])))
2*pnorm(-sqrt(2*(max(ll.seq[[2]]) - ll.seq[[2]][20])))

# increase in log-likelihood by using optimal vs. 2.0
ll.seq[[1]][11] - ll.seq[[1]][20] # men: 17.58
ll.seq[[2]][13] - ll.seq[[2]][20] # women: 4.43

opt.alpha <- alpha.seq[c(which.max(ll.seq[[1]]), which.max(ll.seq[[2]]))]
# 1.1 and 1.3



# Look at BMI*height for M6:
summary(M6[[1]])$coefficients["bmis:heights", ]
summary(M6[[2]])$coefficients["bmis:heights", ]
# z is about 0.43 and -0.03, respectively.
# p is about 0.67 and  0.98, respectively.

# Summary of the sequence of models that were fit:
model.names <- c("Null", 
                 "Adams 2006 Variables", 
                 "All Variables", 
                 "Interactions(alpha = 2.0) - BMI x Height", 
                 "All Interactions(alpha = 2.0)", 
                 "Interactions(alpha = optimal) - BMI x Height", 
                 "All Interactions(alpha = optimal)")

df.vec <- matrix(0, 7, 2)
ll.vec  <- matrix(0, 7, 2)
for (i in 1:2) {
  df.vec[1, i] <- 1
  df.vec[2, i] <- summary(M1[[i]])$waldtest[2]
  df.vec[3, i] <- summary(M2[[i]])$waldtest[2]
  df.vec[4, i] <- summary(M4[[i]])$waldtest[2]
  df.vec[5, i] <- summary(M3[[i]])$waldtest[2]
  df.vec[6, i] <- summary(M5[[i]])$waldtest[2]
  df.vec[7, i] <- summary(M6[[i]])$waldtest[2]
  ll.vec[1, i] <- M0[[i]]$loglik[1]
  ll.vec[2, i] <- M1[[i]]$loglik[2]
  ll.vec[3, i] <- M2[[i]]$loglik[2]
  ll.vec[4, i] <- M4[[i]]$loglik[2]
  ll.vec[5, i] <- M3[[i]]$loglik[2]
  ll.vec[6, i] <- M5[[i]]$loglik[2]
  ll.vec[7, i] <- M6[[i]]$loglik[2]
}

# add in the degrees of freedom from the tied-together parameters 
# to M3, M4, M5, and M6:
df.vec[4:7, ] <- df.vec[4:7, ] + 
                 lu(data[, "race"]) - 1 + 
                 lu(data[, "edu"]) - 1 + 
                 lu(data[, "smoking"]) - 1 + 
                 lu(data[, "physical"]) - 1 + 
                 lu(data[, "alcohol.factor"]) - 1 + 
                 lu(data[, "health"]) - 1 + 
                 lu(data[, "marriage"]) - 1

model.info <- data.frame(name = paste0("M", c(0, 1, 2, 4, 3, 5, 6)), 
                         description = model.names, df=df.vec[, 1], 
                         loglik.men = round(ll.vec[, 1] - max(ll.vec[, 1]), 1), 
                         loglik.women = round(ll.vec[, 2] - max(ll.vec[, 2]), 1), 
                         stringsAsFactors = FALSE)

# Table 1 in the paper:
xtab <- model.info[c(1, 2, 3, 5, 4, 6, 7), -1]
rownames(xtab) <- paste0("M", c(0, 1, 2, 3, 4, 5, 6))

# Create the LaTeX table from R:
library(xtable)
xtable(xtab, digits = 1)

# Just the three main models for the paper:
model.info[c(1, 2, 3, 5, 4, 6, 7), ][c(3, 4, 7), ]

# All the models:
model.info[c(1, 2, 3, 5, 4, 6, 7), ]










