setwd("C:\\R wd\\BIPP_dCVnet")

##load necessary packages####
library(tidyverse)

# Inspect Data#####
load("data/zpsypred.RData")
# Numeric variables apart from IMD and Age were z-scored.
# The z-scores are based on this exact sample.

# Check missing values#
mice::md.pattern(dt)
mice::md.pattern(zdt)
#no missing values

#Data adjustments####
#code three factors (Sex, Ethnicity and Mother_Edu) to numeric

levels(dt$Sex)
# recode to 0,1 numeric (a dummy variable):
dt$Sex <- as.numeric(dt$Sex == "2")

levels(dt$Ethnicity)
table(dt$Ethnicity)
#  2  3  4  5  6 
# 18 16 96  4 19  

# Disclaimer - I don't have the variable dictionary to hand.

# Recode Ethnicity categories to a single dummy variable
#     which indicates non-majority ethnicity:
dt$Ethnicity <- as.numeric(dt$Ethnicity != "4")

# ~ Maternal education ----
levels(dt$Mother_Edu)
table(dt$Mother_Edu)
# 16 years 17-19 ye 19 years Still in  Unknown 
#        7       22      119        2        3  

# Unknown must be imputed. 
# Still in is too small to be useful and probably implies >19
# 16 is a small category

# Recode to a single dummy variable indicating 
#   19+ years (or still in education)
#   the reference category will be <19 years (or unknown)