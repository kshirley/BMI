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

### Comparison of POBs when including or excluding chronically ill:
### Resume from the workspace above where we are including chronically ill:

# exlcuding chronically ill:
load("pob_scaled.RData")
pob1 <- pob.scaled

# including chronically ill:
load("chronically-ill/pob_scaled.RData")
pob2 <- pob.scaled
pob2 <- pob2[data[, "chronic"] == 0]

# some comparisons:
cor(pob1, pob2, use = "complete.obs")  # 0.96

mean(pob1, na.rm = TRUE)  # 25.98
mean(pob2, na.rm = TRUE)  # 26.86
round(mean(pob2, na.rm = TRUE) - mean(pob1, na.rm = TRUE), 2)

# look at standard deviations:
sd(pob1, na.rm = TRUE)
sd(pob2, na.rm = TRUE)
sd(pob2 - pob1, na.rm = TRUE)










