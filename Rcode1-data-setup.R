#2345678901234567890123456789012345678901234567890123456789012345678901234567890

##########################
### Data Setup Section ###
##########################

# Clear the current workspace:
rm(list=ls())

# Set the working directory:
setwd("~/Git/Obesity")

# set the data path:
data.path <- "~/Stats/BMI/data/"

# Load the survival library to fit cox proportional hazards models:
library(survival)

# Load utility functions:
lu <- function(x) length(unique(x))
su <- function(x) sort(unique(x))

# read in the original raw data:
# 117 MB, 566,398 rows, 56 columns:
data <- read.csv(paste0(data.path, "rawdata.20May2013.csv"), as.is = TRUE)

# read in the variables that we're used in Adams 2006:
vars <- readLines("adams_2006_vars.txt")

# Add the variables that we might want to include in our models:
vars <- c(vars, "MARRIAGE", "HEALTH", "DIABETES")

# reduce the data to only the variables we will include in our models:
data <- data[, vars]

# Change class of height and weight to numeric:
data[, "HT_CUR"] <- as.numeric(data[, "HT_CUR"])
data[, "WEIGHT"] <- as.numeric(data[, "WEIGHT"])

# Set all those who died after 12/31/2009 to being censored (and alive) on that 
# day. Our data source told us that the 2010 death records are not complete, 
#  which meanswe must consider 12/31/2009 as the date of censoring:
died <- data[, "DOD"] != "."
dod <- as.Date(data[died, "DOD"], format="%m/%d/%Y")
sel.2010deaths <- dod > "2009-12-31"
data[died, "DOD"][sel.2010deaths] <- "."
data[died, "EXIT_DT_ALLCAUSEMORT"][sel.2010deaths] <- "12/31/2009"

# Discard data that are outliers or missing based on their height, weight, or
# caloric or alcohol consumption:
height.out <- data[, "HT_CUR"] < 1.4 | data[, "HT_CUR"] > 2.1 | 
              is.na(data[, "HT_CUR"])
weight.out <- data[, "WEIGHT"] < 70 | data[, "WEIGHT"] > 400 | 
              is.na(data[, "WEIGHT"])
calories.out <- data[, "CALORIES"] < 200 | data[, "CALORIES"] > 6000
alcohol.out <- data[, "ALCOHOL_ALC_DRINKS"] > 200  # measured in grams/day

# Count how many records are removed with each criterion:
sum(height.out)  # 8567
sum(weight.out)  # 13654
sum(height.out | weight.out)  # 16987
sum(calories.out)  # 3418
sum(alcohol.out)  # 5365

# Total records removed:
sum(height.out | weight.out | calories.out | alcohol.out)  # 24266

# Sequential records removed:
sum(height.out | weight.out)  # 16987
sum(calories.out & !(height.out | weight.out))  # 3180
sum(alcohol.out & !(height.out | weight.out | calories.out))  # 4099

# records remaining in the data:
dim(data)[1] - sum(height.out | weight.out) - 
               sum(calories.out & !(height.out | weight.out)) - 
               sum(alcohol.out & !(height.out | weight.out | calories.out))
# 542132

# remove these records:
data <- data[!(height.out | weight.out | calories.out | alcohol.out), ]

### Re-code some variables to make their names line up in order in R output:

# Sex:
sex <- data[, "SEX"]

sex[sex == 0] <- "1=male"
sex[sex == 1] <- "2=female"

# Race:
race <- data[, "RACEI"]
race[race == 1] <- "1=white"
race[race == 2] <- "2=black"
race[race == 3] <- "3=hispanic"
race[race == 4] <- "4=asian"
race[race == 5] <- "5=pacificislander"
race[race == 6] <- "6=nativeamerican"
race[race == 9] <- "5=unknown"

# Education:
edu <- data[, "EDUC"]
edu[edu == 1] <- "1=lessthan8years"
edu[edu == 2] <- "2=8-11years"
edu[edu == 3] <- "3=highschool"
edu[edu == 4] <- "4=vocation/techschool"
edu[edu == 5] <- "5=somecollege"
edu[edu == 6] <- "6=collegegrad"
edu[edu == 7] <- "7=postgrad"
edu[edu == 9] <- "9=unknown"

# Smoking:
dose.vec <- c("1-10", "11-20", "21-30", "31-40", "41-60", "60+")

smoking <- data[, "SMOKE_QUIT_DETAILED"]

smoking[smoking == 0] <- "01=nonsmoker"

smoking[smoking %in% 1:6] <- paste0(c("02", "03", "04", "05", "06", "07")[as.numeric(smoking[smoking %in% 1:6])], "=", "quit10+dose", dose.vec[as.numeric(smoking[smoking %in% 1:6])])

smoking[smoking %in% 7:12] <- paste0(c("08", "09", "10", "11", "12", "13")[as.numeric(smoking[smoking %in% 7:12]) - 6], "=", "quit5-9dose", dose.vec[as.numeric(smoking[smoking %in% 7:12]) - 6])

smoking[smoking %in% 13:18] <- paste0(c("14", "15", "16", "17", "18", "19")[as.numeric(smoking[smoking %in% 13:18]) - 12], "=", "quit1-4dose", dose.vec[as.numeric(smoking[smoking %in% 13:18]) - 12])

smoking[smoking %in% 19:24] <- paste0(c("20", "21", "22", "23", "24", "25")[as.numeric(smoking[smoking %in% 19:24]) - 18], "=", "quit<1dose", dose.vec[as.numeric(smoking[smoking %in% 19:24]) - 18])

smoking[smoking %in% 25:30] <- paste0(c("26", "27", "28", "29", "30", "31")[as.numeric(smoking[smoking %in% 25:30]) - 24], "=", "currentdose", dose.vec[as.numeric(smoking[smoking %in% 25:30]) - 24])

smoking[smoking == 31] <- "32=unknown/missing"

# Physical Activity:
physical <- data[, "PHYSIC"]
physical[physical == 0] <- "1=never"
physical[physical == 1] <- "2=rarely"
physical[physical == 2] <- "3=1-3permonth"
physical[physical == 3] <- "4=1-2perweek"
physical[physical == 4] <- "5=3-4perweek"
physical[physical == 5] <- "6=5+perweek"
physical[physical == 9] <- "7=unknown/missing"

# Alcohol:
alcohol <- data[, "ALCOHOL_ALC_DRINKS"]/14
# number of drinks per day (14 grams is one standard drink)
alcohol.factor <- ceiling(alcohol)
alcohol.factor[alcohol.factor > 7] <- 8

# Chronic conditions:
chronic <- as.numeric(data[, "SELF_CANCER"] == 1 | 
                      data[, "HEART"] == 1 | 
                      data[, "RENAL"] == 1 | 
                      data[, "EMPHYSEMA"] == 1 |
                      data[, "STROKE"] == 1)

# Health (self-reported):
health <- data[, "HEALTH"]
health[health == 1] <- "1=excellent"
health[health == 2] <- "2=verygood"
health[health == 3] <- "3=good"
health[health == 4] <- "4=fair"
health[health == 5] <- "5=poor"
health[health == 9] <- "6=unknown"

# Marriage:
marriage <- data[, "MARRIAGE"]
marriage[marriage == 1] <- "1=married"
marriage[marriage == 2] <- "2=widowed"
marriage[marriage == 3] <- "3=divorced"
marriage[marriage == 4] <- "4=separated"
marriage[marriage == 5] <- "5=nevermarried"
marriage[marriage == 9] <- "6=unknown"

# Diabetes:
diabetes <- data[, "DIABETES"]
diabetes[diabetes == 0] <- "1=no"
diabetes[diabetes == 1] <- "2=yes"

# Create the data.frame to do our analysis:
data <- data.frame(data[, c("ENTRY_DT", "DOD", "EXIT_DT_ALLCAUSEMORT", 
                            "PERSONYRS_ALLCAUSEMORT")], 
                   age = as.numeric(data[, "ENTRY_AGE"]), 
                   calories = data[, "CALORIES"],
                   weight = data[, "WEIGHT"],
                   height = data[, "HT_CUR"],
                   sex = as.factor(sex), 
                   race = as.factor(race),
                   edu = as.factor(edu),
                   smoking = as.factor(smoking),
                   physical = as.factor(physical), 
                   alcohol = as.numeric(alcohol),
                   alcohol.factor = as.factor(alcohol.factor), 
                   chronic = as.factor(chronic),
                   health = as.factor(health), 
                   marriage = as.factor(marriage),
                   diabetes = as.factor(diabetes))

# Compute number of data points:
n <- dim(data)[1]

# compute bmi and add to data.frame:
bmi <- data[, "weight"]*0.453592/(data[, "height"]^2)
data <- data.frame(data, bmi)

# filter data on extreme bmi values:
data <- data[data[, "bmi"] > 15 & data[, "bmi"] < 50, ]
# removing 1455 individuals
n <- dim(data)[1]

# Define event time:
status <- as.numeric(data[, "DOD"] != ".")  # 1 = dead, 0 = alive
time <- rep(NA, n)
time[status == 1] <- as.Date(data[status == 1, "DOD"], format="%m/%d/%Y") - 
                     as.Date(data[status == 1, "ENTRY_DT"], format="%m/%d/%Y")
time[status == 0] <- as.Date(data[status == 0, "EXIT_DT_ALLCAUSEMORT"], 
                             format="%m/%d/%Y") - 
                     as.Date(data[status == 0, "ENTRY_DT"], format="%m/%d/%Y")

# convert y (time to event in days) to years, accounting for leap years:
time <- time/365.25
data <- data.frame(data, time, status)

# Store BMI_2, age, and height means and sds:
bmi.mean <- mean(data[, "bmi"])
bmi.sd <- sd(data[, "bmi"])

age.mean <- mean(data[, "age"])
age.sd <- sd(data[, "age"])

height.mean <- mean(data[, "height"])
height.sd <- sd(data[, "height"])

# save summary stats:
ss <- list(bmi.mean = bmi.mean, 
           bmi.sd = bmi.sd, 
           height.mean = height.mean, 
           height.sd = height.sd, 
           age.mean = age.mean, 
           age.sd = age.sd)
     
save(ss, file = "ss.RData")

# compute scaled versions of bmi, age, and height (so that the Cox proportional 
# hazards models converge stably when they are fit)
# (the 's' at the end of the variable name indicates they are scaled)
bmis <- as.numeric(scale(data[, "bmi"]))
ages <- as.numeric(scale(data[, "age"]))
heights <- as.numeric(scale(data[, "height"]))

# append these scaled variables to the data.frame
data <- cbind(data, bmis, ages, heights)

# Save this data.frame:
save(data, file = paste0(data.path, "aarp-new.RData"))

