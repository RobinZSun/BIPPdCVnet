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

