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
men <- data[, "sex"] == sex[1]
women <- data[, "sex"] == sex[2]

# Store BMI_2 mean and sd:
load("ss.RData")
bmi.mean <- ss$bmi.mean
bmi.sd <- ss$bmi.sd

# load a few objects from the model fitting script:
load(file = paste0(data.path, "dt.RData"))
load(file = paste0(data.path, "M3.RData"))
load(file = paste0(data.path, "M6.RData"))

# pob.scaled:
load(file = "pob_scaled.RData")
pob <- (pob.scaled - bmi.mean)/bmi.sd

# Mean BMI_2 and POB_2 by sex:
bmi.mean.by.sex <- aggregate(data[, "bmi"], by = list(data[, "sex"]), mean)
pob.mean.by.sex <- aggregate(pob.scaled, by=list(data[, "sex"]), mean,  
                             na.rm = TRUE)
bmi.sd.by.sex <- aggregate(data[, "bmi"], by = list(data[, "sex"]), sd)
pob.sd.by.sex <- aggregate(pob.scaled, by=list(data[, "sex"]), sd, na.rm=TRUE)

# confidence intervals for each person:
load(file = "ci_derivative.RData")

lower <- ci[, 1]
upper <- ci[, 2]
lower2 <- ci[, 3]
upper2 <- ci[, 4]

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

### Look at interval widths:

# mean width of lower interval
mean(ci[sel.edges, 2] - ci[sel.edges, 1])

# median width of lower interval
median(ci[sel.edges, 2] - ci[sel.edges, 1])

# standard deviation of width of lower interval
sd(ci[sel.edges, 2] - ci[sel.edges, 1])

# histogram of lower interval width
hist(ci[sel.edges, 2] - ci[sel.edges, 1])

# mean width of lower interval; men
mean(ci[sel.edges & men, 2] - ci[sel.edges & men, 1])  # 1.7

# mean width of lower interval; women
mean(ci[sel.edges & women, 2] - ci[sel.edges & women, 1])  # 2.4

# Compute number of men and women who are below the mean pob:
bmi.type <- numeric(n)

# too light:
bmi.type[sel.edges][data[sel.edges, "bmi"] < ci[sel.edges, "lower"]] <- 1

# too heavy:
bmi.type[sel.edges][data[sel.edges, "bmi"] > ci[sel.edges, "upper"]] <- 2

# within POB interval:
bmi.type[sel.edges][data[sel.edges, "bmi"] >= ci[sel.edges, "lower"] & 
                    data[sel.edges, "bmi"] <= ci[sel.edges, "upper"]] <- 3

# Look at distribution of 'bmi.type'
table(bmi.type[sel.edges], data[sel.edges, "sex"])

# Look to see how many intervals do and don't contain the mean of the pob:
pob.type <- numeric(n)
pob.vec <- pob.mean.by.sex[match(data[sel.edges, "sex"], sex), 2]

# pob is below interval
pob.type[sel.edges][ci[sel.edges, "lower"] > pob.vec] <- 1

# pob is above interval
pob.type[sel.edges][ci[sel.edges, "upper"] < pob.vec] <- 2

# pob is within interval
pob.type[sel.edges][ci[sel.edges, "lower"] <= pob.vec & 
                    ci[sel.edges, "upper"] >= pob.vec] <- 3

pob.tb <- table(pob.type[sel.edges], data[sel.edges, "sex"])
pob.tb

# proportion of men whose POB_2 interval does not contain the mean POB_2:
sum(pob.tb[1:2, 1])/sum(pob.tb[, 1])
round(sum(pob.tb[1:2, 1])/sum(pob.tb[, 1]), 2)

# proportion of women whose POB_2 interval does not contain the mean POB_2:
sum(pob.tb[1:2, 2])/sum(pob.tb[, 2])
round(sum(pob.tb[1:2, 2])/sum(pob.tb[, 2]), 2)


# mean absolute differences:
mean(abs(pob.scaled[men] - data[men, "bmi"]), na.rm = TRUE)
mean(abs(pob.scaled[women] - data[women, "bmi"]), na.rm = TRUE)

# Figure 2: Histograms of differences between BMI and POB by sex:
pdf(file = "fig_absolute_value_difference.pdf", height=5, width=10)
par(mfrow=c(1, 2))
hist(data[men, "bmi"] - pob.scaled[men], breaks = -20:27, prob = TRUE, 
     main="Men", xlab=expression(BMI2[2] - POB[2]), ylim = c(0, 0.13), las = 1)
text(-10, 0.12, expression( " " %<-% underweight))
text(15, 0.12, expression( overweight %->% " "))
hist(data[women, "bmi"] - pob.scaled[women], breaks = -20:27, prob = TRUE, 
     main="Women", xlab=expression(BMI2[2] - POB[2]), ylim=c(0, 0.13), las=1)
text(-10, 0.12, expression( " " %<-% underweight))
text(15, 0.12, expression( overweight %->% " "))
dev.off()







# Figure 1: Density estimates of POB and BMI by sex:
pdf(file="fig_pob_histograms.pdf", width=6, height=4)
leg1.text <- c("Men", "Women")
leg2.text <- c(expression('POB'[2]~'density'), expression('POB'[2]~'mean'), 
               expression('BMI'[2]~'mean'))
par(mfrow=c(1, 1))
plot(density(pob.scaled[data[, "sex"] == sex[1]], na.rm=TRUE), 
     xlim = c(20, 35), 
     type = "l", 
     las = 1, 
     xlab = expression("BMI"[2]), 
     main = "")
lines(density(pob.scaled[data[, "sex"] == sex[2]], na.rm=TRUE), col=2)
title(main=expression("Density estimates of POB"[2] ~ "for men and women"))
abline(v=bmi.mean.by.sex[, 2], lty=3, col=1:2)
abline(v=pob.mean.by.sex[, 2], lty=2, col=1:2)
legend("topleft", inset=0.02, lty=1, col=c(1, 2), legend=leg1.text, cex=0.7)
legend("topright", inset=0.02, lty=c(1, 2, 3), col=1, legend=leg2.text, cex=0.7)
dev.off()


# Figure 1 of the paper - black and white version:
pdf(file="fig_pob_histograms_blackwhite.pdf", width=6, height=4)
leg1.text <- c("Men", "Women")
leg2.text <- c(expression('POB'[2]~'density'), expression('POB'[2]~'mean'),  
               expression('BMI'[2]~'mean'))
par(mfrow=c(1, 1))
plot(density(pob.scaled[data[, "sex"] == sex[1]], na.rm=TRUE), 
     xlim = c(20, 35), 
     type = "l", 
     las = 1, 
     xlab = expression("BMI"[2]), 
     main = "", 
     col = gray(0.5))
lines(density(pob.scaled[data[, "sex"] == sex[2]], na.rm=TRUE), col=1)
title(main=expression("Density estimates of POB"[2] ~ "for men and women"))
abline(v=bmi.mean.by.sex[, 2], lty=3, col=c(gray(0.5), 1))
abline(v=pob.mean.by.sex[, 2], lty=2, col=c(gray(0.5), 1))
legend("topleft", inset = 0.02, lty = 1, col = c(gray(0.5), 1), 
       legend = leg1.text, cex = 0.7, text.col = c(gray(0.5), 1))
legend("topright", inset = 0.02, lty = c(1, 2, 3), col = 1, 
       legend = leg2.text, cex = 0.7)
dev.off()







### Plot relative risk curves for poeple every cell of table of
### (num.edges x is.na(pob))

# In other words, inspect all types of intervals crossed with POBs not being NA
set.seed(7474)
examples <- matrix(0, 5, 2)
for (i in 1:5) {
  for (j in 1:2) {
  	sel <- num.edges == (0:4)[i] & is.na(pob) == c(FALSE, TRUE)[j]
    if (sum(sel) > 0) {
      if (sum(sel) == 1) {
        examples[(j - 1)*5 + i] <- which(sel)
      } else {
        examples[(j - 1)*5 + i] <- sample(which(sel), 1)
      }
    }
  }
}

# Compute the relative risk curves for these 4 individuals:
bmi.seq <- seq(15, 50, 0.2)
lb <- length(bmi.seq)
rel.risk <- matrix(NA, lb, 7)
t1 <- Sys.time()
for (i in 1:2) {
  for (j in 1:lb) {
  	if (j %% 20 == 0) print(j)
  	# examples 2:8 are the ones with non-zero sample size
  	ind <- which(data[examples[2:8], "sex"] == sex[i])
  	x <- dt[[i]][match(examples[2:8][ind], which(data[, "sex"] == sex[i])), ]
  	x["bmis"] <- (bmi.seq[j] - bmi.mean)/bmi.sd
  	rel.risk[j, ind] <- predict(M3[[i]], newdata = x, type = "lp")
  }
}
t2 <- Sys.time()
t2 - t1 # about 15 seconds

# Figure to illustrate the possibilities.
# Make sure they all make sense:
pdf(file = "fig_rr_pob_examples_v4.pdf", height = 6, width = 12)
for (i in 1:5) {
  par(mfrow = c(1, 2))
  for (j in 1:2) {
  	sel <- num.edges == (0:4)[i] & is.na(pob) == c(FALSE, TRUE)[j]
  	if (examples[i, j] > 0) {
      col.vec <- rep(1, lb)
      if (ci[examples[i, j], 1] != 0) {
        w1 <- min(which(bmi.seq > ci[examples[i, j], 1]))
        w2 <- ifelse(any(bmi.seq > ci[examples[i, j], 2]), 
                     min(which(bmi.seq > ci[examples[i, j], 2])), 
                     lb)
        col.vec[w1:w2] <- 2
      }
      if (ci[examples[i, j], 3] != 0) {
        w1 <- min(which(bmi.seq > ci[examples[i, j], 3]))
        w2 <- ifelse(ci[examples[i, j], 4] != 0, 
                     min(which(bmi.seq > ci[examples[i, j], 4])), 
                     lb)
        col.vec[w1:w2] <- 3
      }
      plot(bmi.seq, exp(rel.risk[, match((j - 1)*5 + i, 2:8)]), 
           las = 1, xlab = "BMI2", ylab = "exp(Xbeta)", 
           col = col.vec)
      title(main = paste0("num edges = ", i - 1, " & POB ", 
                          c("exists", "is not defined")[j], ", ", sum(sel), 
                          " cases"))
      if (j == 1) {
        legend("topleft", inset = 0.02, pch=1, 
               legend = paste0("POB = ", round(pob.scaled[examples[i, j]], 1)), 
               bty = "n")
      }
    } else {
      plot(0, 0, type = "n", 
           main = paste0("num edges = ", i - 1, "; POB ", 
                         c("exists", "is not defined")[j], ", ", sum(sel), 
                         " cases"))
    }
  }
}
dev.off()




### Compute relative risk curves of some example people in the data:

# compute a pob range with width 0.1 at +/- 2 standard deviations from 
# average for men and women:
pob.range <- matrix(0, 4, 2)
for (i in 1:2){ 
  tmp <- mean(pob.scaled[data[, "sex"] == sex[i]], na.rm=TRUE) + 
         c(-2, 2)*sd(pob.scaled[data[, "sex"] == sex[i]], na.rm=TRUE)
  pob.range[(i-1)*2 + 1, ] <- tmp[1] + c(-0.05, 0.05)
  pob.range[(i-1)*2 + 2, ] <- tmp[2] + c(-0.05, 0.05)
}

# Now sample an individual in the pob range specified above (men and women):
set.seed(135)
examples <- numeric(4)
for (i in 1:2) {
  for (j in 1:2) {
  	sel <- data[, "sex"] == sex[i] & 
  	       pob.scaled > pob.range[(i-1)*2 + j, 1] & 
  	       pob.scaled < pob.range[(i-1)*2 + j, 2] & 
  	       !is.na(pob.scaled)
    print(paste(i, j, sum(sel)))
  	tmp <- sample(sum(sel), 1)
    examples[(i-1)*2 + j] <- which(sel)[tmp]
  }
}
cbind(data[examples, ], pob.scaled[examples])

# Compute the relative risk curves for these 4 individuals:
bmi.seq <- seq(18, 35, 0.2)
lb <- length(bmi.seq)
rel.risk <- matrix(NA, lb, 4)
t1 <- Sys.time()
for (i in 1:4) {
  print(i)
  for (j in 1:lb) {
    if (j %% 20 == 0) print(j)
  	x <- dt[[as.numeric(i > 2) + 1]][which(which(data[, "sex"] == sex[as.numeric(i > 2) + 1]) == examples[i]), ]
  	x["bmis"] <- (bmi.seq[j] - bmi.mean)/bmi.sd
  	rel.risk[j, i] <- predict(M3[[as.numeric(i > 2) + 1]], newdata = x)  	
  }
}
t2 <- Sys.time()
t2 - t1 # 15 seconds

# compute an average relative risk curve:
# coarser grid here for faster computation:
b.seq <- seq(18, 35, 0.5)
lbs <- length(b.seq)

avg.lp <- matrix(0, lbs, 2)
avg.risk <- matrix(0, lbs, 2)
t1 <- Sys.time()
for (i in 1:2) {
  print(i)
  x <- dt[[i]]
  for (j in 1:lbs) {
    if (j %% 10 == 0) print(j)
    x[, "bmis"] <- (b.seq[j] - bmi.mean)/bmi.sd
    tmp <- predict(M3[[i]], newdata=x)
    avg.lp[j, i] <- mean(tmp)
    avg.risk[j,i] <- mean(exp(tmp))
  }
}
t2 <- Sys.time()
t2 - t1 # 1.4 minutes

# Plot of relative risk vs. average:
title.p1 <- c("(a) Male, ", "(b) Male, ", "(c) Female, ", "(d) Female, ")
title.p2 <- as.character(round(pob.scaled[examples], 1))

# pdf of the plot:
leg.x <- c(24, 24, 24, 24)
leg.y <- c(1.65, 5.75, 1.6, 4.2)
leg.sex <- c("Men", "Women")

# create the plot:
pdf(file="fig_relrisk_examples.pdf", width=9, height=7)
par(mfrow=c(2, 2), mar=c(4, 4, 3, 1))
for (i in 1:4) {
  ylm <- range(c(exp(rel.risk[, i]), avg.risk[, as.numeric(i > 2) + 1]))
  plot(x = bmi.seq, 
       y = exp(rel.risk[, i]), 
       type="l", 
       las=1, 
       main=substitute(x~POB[2] == y, list(x=title.p1[i], y=title.p2[i])), 
       ylab="Relative Risk", 
       xlab=expression(BMI[2.0]), 
       ylim=ylm)
  w1 <- min(which(bmi.seq > ci[examples[i], 1]))
  w2 <- max(which(bmi.seq < ci[examples[i], 2]))
  lines(bmi.seq[w1:w2], exp(rel.risk[w1:w2, i]), lwd=6, col=gray(0.7))
  lines(bmi.seq, exp(rel.risk[, i]))
  lines(b.seq, avg.risk[, as.numeric(i > 2) + 1], lty=as.numeric(i>2) + 2)
  legend(leg.x[i], leg.y[i], lty=c(1, as.numeric(i>2) + 2, 1), 
         legend=c("Individual", 
                  paste0("Average (", leg.sex[as.numeric(i>2) + 1], ")"), 
                  "POB confidence interval"),
         lwd=c(1, 1, 6), col=c(1, 1, gray(0.7)), cex=0.8)
}
dev.off()



# Compute POB(2) for each person in the data based on Model M3:
beta <- matrix(0, 13, 2)
for (i in 1:2) beta[, i] <- coef(M3[[i]])[c(3, 1, 2, 23:31, 77)]
rownames(beta) <- names(coef(M3[[i]])[c(3, 1, 2, 23:31, 77)])

# Compute cubic, quadratic, and linear components of relative risk curves:
a.examples <- numeric(4)
b.examples <- numeric(4)
c.examples <- matrix(0, 11, 4)
for (ii in 1:4) {
  i <- as.numeric(ii > 2) + 1
  j <- which(which(data[, "sex"] == sex[i]) == examples[ii])
  print(paste(ii, i, j))
  c.examples[, ii] <- c(beta[1, i], 
                        beta[4, i]*dt[[i]][j, "ages"], 
                        beta[5, i]*dt[[i]][j, "diabetes"], 
                        beta[6, i]*dt[[i]][j, "race"],
                        beta[7, i]*dt[[i]][j, "edu"], 
                        beta[8, i]*dt[[i]][j, "smoking"], 
                        beta[9, i]*dt[[i]][j, "physical"], 
                        beta[10, i]*dt[[i]][j, "alcohol.factor"], 
                        beta[11, i]*dt[[i]][j, "health"], 
                        beta[12, i]*dt[[i]][j, "marriage"], 
                        beta[13, i]*dt[[i]][j, "heights"])
  a.examples[ii] <- beta[3, i]
  b.examples[ii] <- beta[2, i]
}
rownames(c.examples) <- rownames(beta)[c(1, 4:13)]

# Look at the four individuals, and order the effects of each demographic 
# variable on the minimum point of the rel risk curve:
c.examples[order(c.examples[, 1]), 1]
c.examples[order(c.examples[, 2]), 2]
c.examples[order(c.examples[, 3]), 3]
c.examples[order(c.examples[, 4]), 4]
# Here, negative numbers push the curve to the right


c.examples
apply(c.examples, 2, rank)



cbind(data[examples, ], pob.scaled[examples])

zz <- cbind(data[examples, ], pob.scaled[examples])
zz <- zz[, c("sex", "age", "height", "weight", "race", "edu", "smoking", 
             "physical", "alcohol.factor", "health", "marriage", "diabetes", 
             "bmi")]
rownames(zz) <- letters[1:4]

zz[, "age"] <- floor(zz[, "age"])
zz[, "height"] <- round(zz[, "height"], 2)
zz[, "bmi"] <- round(zz[, "bmi"], 1)



data.frame(t(zz))









# Now, stratify the population by POB_2, and show mortality rates as a 
# function of observed BMI_2:

# First, add tiny noise to BMI so that there aren't so many ties:
set.seed(123)
bmi.noise <- data[, "bmi"] + rnorm(n, 0, 0.01)

# set up sequence to compute quadratic fits through points:
x.seq <- seq(20, 42, 0.1)

# set up color vector:
library(RColorBrewer)
col.vec <- brewer.pal(5, "Set1")

# s = 1 for men
s <- 1
sel.sex <- data[, "sex"] == sex[s]
n.pob <- 5
n.bmi <- c(5, 10, 20, 10, 5)

q <- quantile(pob.scaled[sel.sex], c(0, 0.1, 0.3, 0.7, 0.9, 1), na.rm = TRUE)
pob.cuts <- cut(pob.scaled[sel.sex], q, include.lowest=TRUE)

# set up objects ahead of the loop:
total.mat <- list(0, length(n.pob))
for (i in 1:n.pob) total.mat[[i]] <- numeric(n.bmi[i])

died.mat <- list(0, length(n.pob))
for (i in 1:n.pob) died.mat[[i]] <- numeric(n.bmi[i])

bmi.cuts <- as.list(rep(NA, n.pob))
pob.bin.means <- numeric(n.pob)

for (i in 1:n.pob) {
  pob.sel <- pob.cuts == levels(pob.cuts)[i] & !is.na(pob.cuts)
  bmi.cuts[[i]]$cuts <- cut(bmi.noise[sel.sex][pob.sel], 
                            quantile(bmi.noise[sel.sex][pob.sel], 
                                     seq(0, 1, length=n.bmi[i] + 1)), 
                            include.lowest = TRUE)
  bmi.cuts[[i]]$means <- numeric(n.bmi[i])
  pob.bin.means[i] <- mean(pob.scaled[sel.sex][pob.sel], na.rm = TRUE)
  for (j in 1:n.bmi[i]) {
    sel <- bmi.cuts[[i]]$cuts == levels(bmi.cuts[[i]]$cuts)[j]
    total.mat[[i]][j] <- sum(sel)
    died.mat[[i]][j] <- sum(data[sel.sex, ][pob.sel, "status"][sel] == 1 & 
                            data[sel.sex, ][pob.sel, "time"][sel] < 12.9)
    bmi.cuts[[i]]$means[j] <- mean(bmi.noise[sel.sex][pob.sel][sel])
  }
}

# Look at death rates:
total.mat
died.mat

# range of 95% intervals for observed proportions:
2*sqrt(min(unlist(died.mat)/unlist(total.mat))*(1 - min(unlist(died.mat)/unlist(total.mat)))/total.mat[[1]][1])
2*sqrt(max(unlist(died.mat)/unlist(total.mat))*(1 - max(unlist(died.mat)/unlist(total.mat)))/total.mat[[1]][1])

## Picture for Men:
pdf(file="fig_death_vs_pob_men.pdf", width=6, height=6)
par(mfrow=c(1, 1), mar=c(5, 4, 4, 2))
plot(x = bmi.cuts[[1]]$means, 
     y = died.mat[[1]]/total.mat[[1]], 
     xlim = c(19, 42), 
     ylim = c(0, 0.43), 
     las = 1, 
     xlab = expression(Observed~BMI[2]), 
     ylab = "Proportion who died during study", 
     pch = 20, 
     xaxt = "n", 
     col = col.vec[1])
abline(h = sum(data[men, "status"])/sum(men), col=gray(0.6))
axis(1, at=19:42)
y.vec <- died.mat[[1]]/total.mat[[1]]
x.vec <- bmi.cuts[[1]]$means
f <- lm(y.vec ~ x.vec + I(x.vec^2))
lines(x.seq, coef(f)[1] + coef(f)[2]*x.seq + coef(f)[3]*x.seq^2, 
      col = col.vec[1], lty = 2)
points(pob.bin.means[1], coef(f)[1] + coef(f)[2]*pob.bin.means[1] + 
       coef(f)[3]*pob.bin.means[1]^2, col = col.vec[1], pch = 2)
for (j in 2:n.pob) {
  points(bmi.cuts[[j]]$means, died.mat[[j]]/total.mat[[j]], 
         col = col.vec[j], pch = 20)
  y.vec <- died.mat[[j]]/total.mat[[j]]
  x.vec <- bmi.cuts[[j]]$means
  f <- lm(y.vec ~ x.vec + I(x.vec^2))
  lines(x.seq, coef(f)[1] + coef(f)[2]*x.seq + coef(f)[3]*x.seq^2, 
        col=col.vec[j], lty=2)
  points(pob.bin.means[j], coef(f)[1] + coef(f)[2]*pob.bin.means[j] + 
         coef(f)[3]*pob.bin.means[j]^2, col=col.vec[j], pch=2)
}
abline(h=0, lty=1, col=gray(0.6))
legend("bottomright", inset=0.01, col=c(col.vec[5:1], gray(0.7), 1), 
       legend=c(paste0("POB in ", levels(pob.cuts)[5:1]), "Mean Death Rate", 
       "Mean POB of bin"), lty=c(rep(1, 6), NA), cex=0.7, lwd=c(rep(2, 6), NA), 
       pch=c(rep(NA, 6), 2), bg="white")
title(sub="95% intervals range from 0.7% to 1.4%, (n = 4710 for each point)")
title(main=expression("(a) P(death) vs. BMI"[2] ~ " for men"))
dev.off()



# s = 2 for women:
s <- 2
sel.sex <- data[, "sex"] == sex[s]
n.pob <- 5
n.bmi <- c(5, 10, 20, 10, 5)

q <- quantile(pob.scaled[sel.sex], c(0, 0.1, 0.3, 0.7, 0.9, 1), na.rm = TRUE)
pob.cuts <- cut(pob.scaled[sel.sex], q, include.lowest=TRUE)

# set up stuff before loop:
total.mat <- list(0, length(n.pob))
for (i in 1:n.pob) total.mat[[i]] <- numeric(n.bmi[i])

died.mat <- list(0, length(n.pob))
for (i in 1:n.pob) died.mat[[i]] <- numeric(n.bmi[i])

bmi.cuts <- as.list(rep(NA, n.pob))
pob.bin.means <- numeric(n.pob)
for (i in 1:n.pob) {
  pob.sel <- pob.cuts == levels(pob.cuts)[i] & !is.na(pob.cuts)
  bmi.cuts[[i]]$cuts <- cut(bmi.noise[sel.sex][pob.sel], 
                            quantile(bmi.noise[sel.sex][pob.sel], 
                            seq(0, 1, length=n.bmi[i] + 1)), 
                            include.lowest = TRUE)
  bmi.cuts[[i]]$means <- numeric(n.bmi[i])
  pob.bin.means[i] <- mean(pob.scaled[sel.sex][pob.sel], na.rm=TRUE)
  for (j in 1:n.bmi[i]) {
    sel <- bmi.cuts[[i]]$cuts == levels(bmi.cuts[[i]]$cuts)[j]
    total.mat[[i]][j] <- sum(sel)
    died.mat[[i]][j] <- sum(data[sel.sex, ][pob.sel, "status"][sel] == 1 & 
                            data[sel.sex, ][pob.sel, "time"][sel] < 12.9)
    bmi.cuts[[i]]$means[j] <- mean(bmi.noise[sel.sex][pob.sel][sel])
  }
}

# Look at death rates:
total.mat
died.mat

# range of 95% intervals for observed proportions:
2*sqrt(min(unlist(died.mat)/unlist(total.mat))*(1 - min(unlist(died.mat)/unlist(total.mat)))/total.mat[[1]][1])
2*sqrt(max(unlist(died.mat)/unlist(total.mat))*(1 - max(unlist(died.mat)/unlist(total.mat)))/total.mat[[1]][1])


## Women:
pdf(file="fig_death_vs_pob_women.pdf", width=6, height=6)
par(mfrow=c(1, 1), mar=c(5, 4, 4, 2))
plot(x = bmi.cuts[[1]]$means, 
     y = died.mat[[1]]/total.mat[[1]], 
     xlim = c(19, 42), 
     ylim = c(0, 0.43), 
     las = 1, 
     xlab = expression(Observed~BMI[2]), 
     ylab = "Proportion who died during study", 
     pch = 20, 
     xaxt = "n", 
     col = col.vec[1])
abline(h = sum(data[women, "status"])/sum(women), col=gray(0.6))
axis(1, at=19:42)
y.vec <- died.mat[[1]]/total.mat[[1]]
x.vec <- bmi.cuts[[1]]$means
f <- lm(y.vec ~ x.vec + I(x.vec^2))
lines(x.seq, coef(f)[1] + coef(f)[2]*x.seq + coef(f)[3]*x.seq^2, 
      col = col.vec[1], lty = 2)
points(pob.bin.means[1], coef(f)[1] + coef(f)[2]*pob.bin.means[1] + 
       coef(f)[3]*pob.bin.means[1]^2, col=col.vec[1], pch=2)
for (j in 2:n.pob) {
  points(bmi.cuts[[j]]$means, died.mat[[j]]/total.mat[[j]], col=col.vec[j], 
         pch=20)
  y.vec <- died.mat[[j]]/total.mat[[j]]
  x.vec <- bmi.cuts[[j]]$means
  f <- lm(y.vec ~ x.vec + I(x.vec^2))
  lines(x.seq, coef(f)[1] + coef(f)[2]*x.seq + coef(f)[3]*x.seq^2, 
        col=col.vec[j], lty=2)
  points(pob.bin.means[j], coef(f)[1] + coef(f)[2]*pob.bin.means[j] + 
         coef(f)[3]*pob.bin.means[j]^2, col=col.vec[j], pch=2)
}
abline(h=0, lty=1, col=gray(0.6))
legend("topright", inset=0.02, col=c(col.vec[5:1], gray(0.7), 1), 
       legend=c(paste0("POB in ", levels(pob.cuts)[5:1]), "Mean Death Rate", 
       "Mean POB of bin"), lty=c(rep(1, 6), NA), cex=0.7, lwd=c(rep(2, 6), NA), 
       pch=c(rep(NA, 6), 2), bg="white")
title(sub="95% intervals range from 0.8% to 1.5%, (n = 3439 for each point)")
title(main=expression("(b) P(death) vs. BMI"[2] ~ " for women"))
dev.off()










### Let's look at the relative risk curves as a function of BMI for a few 
### different heights:
i <- 2  # 1 for men, 2 for women:

h.vec <- quantile(dt[[i]][, "heights"], c(0.025, 0.25, 0.50, 0.75, 0.975))
height.sd <- ss$height.sd
height.mean <- ss$height.mean

inches.vec <- round((h.vec*height.sd + height.mean) * 39.37)
feet.vec <- floor(inches.vec/12)
height.labels <- paste0(feet.vec, "'", inches.vec - floor(inches.vec/12)*12, 
                        "\"")
lh <- length(h.vec)

# grab a sequence of BMI values from the data:
bmi.seq <- seq(quantile(dt[[i]][, "bmis"], 0.025), 
               quantile(dt[[i]][, "bmis"], 0.975), 
               length = 300)
lb <- length(bmi.seq)

# new data.frame to use in prediction for M3[[1]]:
nd <- data.frame(bmis = rep(bmi.seq, lh), 
                 heights = rep(h.vec, each=lb), 
                 ages = mean(dt[[i]][, "ages"]), 
                 diabetes = mean(dt[[i]][, "diabetes"]), 
                 race = 0,
                 edu = 0, 
                 smoking = 0, 
                 physical = 0, 
                 alcohol.factor = 0, 
                 health = 0, 
                 marriage = 0)
# using the mean age and diabetes status for women (from dt[[i = 2]])

# matrix holding the relative risks:
rr.mat <- matrix(predict(M3[[i]], type = "risk", newdata = nd), lb, lh)
adj <- min(rr.mat)
rr.mat <- rr.mat - adj + 1

# set up the axis:
bmi.axis <- sprintf("%.1f", seq(min(bmi.seq*bmi.sd + bmi.mean), 
                                max(bmi.seq*bmi.sd + bmi.mean), length = 6))

# compute the index of the minimum
min.bmi <- cbind(bmi.seq[apply(rr.mat, 2, which.min)], apply(rr.mat, 2, min))

# set up legend text:
leg.text <- as.list(rep(NA, lh))
for (h in 1:lh) {
  leg.text[[h]] <- bquote(.(height.labels[h]) ~ ": " ~ 'POB'[2.0] ~ 
                          " = " ~ .(round(min.bmi[h, 1] * bmi.sd + bmi.mean, 1)))
}

# Make the plot for Figure 2 (left):
pdf(file = paste0("fig_", c("men", "women")[i], "_bmi2_height.pdf"), 
    width = 7, height = 7)
plot(bmi.seq, rr.mat[, 3], xaxt="n", las=1, 
     xlab=expression("BMI"[2.0]), ylab="Relative Risk", 
     type="l", col=3, yaxt="n", ylim=c(1, 1.9))
abline(h = seq(1, 2.2, 0.2), col=gray(0.5), lwd=0.5)
axis(1, at=seq(min(bmi.seq), max(bmi.seq), length=6), labels=bmi.axis)
axis(2, at=seq(1, 2.2, 0.2), labels=seq(1, 2.2, 0.2), las=2)
for (j in c(1:5)) {lines(bmi.seq, rr.mat[, j], col=j)}
legend("topleft", inset=0.02, col=1:5, lty=1, 
       legend = sapply(leg.text, as.expression), title = "Height", bg = "white")
points(min.bmi[, 1], min.bmi[, 2], col=1:5, pch=4)
title(main=expression(paste("Estimated relative risk curves for women of different heights: ", alpha, " = 2.0")))
dev.off()






## Re-do the plot with alpha = 1.3 (optimal for women) from Model M6:

# Now we'll use M6, where BMI*height is not significant:
b <- data[, "weight"]*0.453592/(data[, "height"]^1.3)
tmp <- dt[[i]]
tmp[, "bmis"] <- (b[data[, "sex"] == sex[i]] - mean(b))/sd(b)

# grab a sequence of BMI values from the data:
bmi.seq <- seq(quantile(tmp[, "bmis"], 0.025),  
               quantile(tmp[, "bmis"], 0.975), length = 300)
lb <- length(bmi.seq)

# new data.frame:
nd <- data.frame(bmis = rep(bmi.seq, lh), 
                 heights = rep(h.vec, each = lb), 
                 ages = mean(dt[[i]][, "ages"]), 
                 diabetes = mean(dt[[i]][, "diabetes"]), 
                 race = 0,
                 edu = 0, 
                 smoking = 0, 
                 physical = 0, 
                 alcohol.factor = 0, 
                 health = 0, 
                 marriage = 0)

# matrix holding the relative risks:
rr.mat <- matrix(predict(M6[[i]], type = "risk", newdata = nd), lb, lh)
adj <- min(rr.mat)
rr.mat <- rr.mat - adj + 1

# set up the axis:
bmi.axis <- sprintf("%.1f", seq(min(bmi.seq*sd(b) + mean(b)), 
                                max(bmi.seq*sd(b) + mean(b)), length = 6))

# compute the index of the 
min.bmi <- cbind(bmi.seq[apply(rr.mat, 2, which.min)], apply(rr.mat, 2, min))

# set up legend text:
leg.text <- as.list(rep(NA, lh))
for (h in 1:lh) {
  leg.text[[h]] <- bquote(.(height.labels[h]) ~ ": " ~ 'POB'[1.3] ~ " = " 
                          ~ .(round(min.bmi[h, 1] * sd(b) + mean(b), 1)))
}

# Make the plot for Figure 2 (right):
pdf(file = paste0("fig_", c("men", "women")[i], "_bmi2a_height.pdf"), 
    width = 7, height = 7)
plot(bmi.seq, rr.mat[, 3], xaxt="n", las=1, 
     xlab=expression("BMI"[1.3]), ylab="Relative Risk", 
     type="l", col=3, yaxt="n", ylim=c(1, 1.9))
abline(h = seq(1, 2.2, 0.2), col=gray(0.5), lwd=0.5)
axis(1, at=seq(min(bmi.seq), max(bmi.seq), length=6), labels=bmi.axis)
axis(2, at=seq(1, 2.2, 0.2), labels=seq(1, 2.2, 0.2), las=2)
for (j in c(1:5)) {
  lines(bmi.seq, rr.mat[, j], col = j)
}
legend("topleft", inset=0.02, col=1:5, lty=1, 
       legend = sapply(leg.text, as.expression), title = "Height", bg = "white")
points(min.bmi[, 1], min.bmi[, 2], col=1:5, pch=4)
title(main=expression(paste("Estimated relative risk curves for women of different heights: ", alpha, " = 1.3")))
dev.off()





