setwd("C:\\R wd\\BIPP_dCVnet")

# load necessary packages####
library(tidyverse)
library(lavaan)
#install dCVnet from Github
install.packages("remotes")
remotes::install_local("dCVnet_1.2.0.tar.gz", dependencies = TRUE, build_vignettes = TRUE)
library(dCVnet)

# Inspect Data#####
load("data/zpsypred.RData")
# Numeric variables apart from IMD and Age were z-scored.
# The z-scores are based on this exact sample.

# Check missing values#
mice::md.pattern(dt)
mice::md.pattern(zdt)
#no missing values

# Data adjustments####
#code three factors (Sex, Ethnicity and Mother_Edu) to numeric

levels(dt$Sex)
# recode to 0,1 numeric (a dummy variable):
dt$Sex <- as.numeric(dt$Sex == "2")

levels(dt$Ethnicity)
table(dt$Ethnicity)
#  2（Asian) 3(Black) 4(Caucasian)  5(Chinese)  6(Other)
#    18        16         96            4          19  

# Recode Ethnicity categories to a single dummy variable
#     which indicates non-majority (non-Caucasian) ethnicity:
dt$Ethnicity <- as.numeric(dt$Ethnicity != "4")

#  0(Caucasian)  1(Other)
#    96         57

levels(dt$Mother_Edu)
table(dt$Mother_Edu)
# 16 years 17-19 years 19 years Still in  Unknown 
#   7         22         119      2         3  

# Recode to a single dummy variable indicating 
#   19+ years (or still in education)
#   the reference category will be <19 years (or unknown)
dt$low_mat_edu <- as.numeric(dt$Mother_Edu == "19 years" | dt$Mother_Edu == "Still in")
with(dt, table(low_mat_edu, Mother_Edu))
#            Mother_Edu
# low_mat_edu 16 years 17-19 ye 19 years Still in Unknown
#           0        7       22        0        0       3
#           1        0        0      119        2       0 

## Separate outcome and predictors ####
# the three outcome variables are psychiatric assessments measured at 8: 
# SDQ (internalizing and externalizing)
# TMCQ (surgency, effortful control & negative affect)
# and SRS; All of them are raw scores 

# outcome data:
ovars <- c("srs8", "sdq8int", "sdq8ext", "TMCQ.SU", "TMCQ.EC", "TMCQ.NA")
odat <- dt %>% select(one_of(ovars))

# Control Variables: extract and center 
# Control Variable for outcomes:
ocon <- dt %>% 
  select(Age8) %>% 
  mutate_all(scale, center = TRUE, scale = FALSE)
# Control Variable for predictors:
pcon <- dt %>% 
  select(Age2, Age4, Sex, Ethnicity, low_mat_edu) %>% 
  mutate_all(scale, center = TRUE, scale = FALSE)

# predictor data:
pdat <- dt %>% 
  select(!all_of(ovars)) %>% 
  select(!all_of(colnames(ocon))) %>% 
  select(!all_of(colnames(pcon))) %>% 
  select(-Mother_Edu)

## Outcome dimension reduction ####
# Inspect bivariate relationships:
odat %>% 
  car::scatterplotMatrix()

# correlation matrix ('X' -> p>0.05 unadjusted):
odat %>% 
  cor() %>% 
  ggcorrplot::ggcorrplot(p.mat = ggcorrplot::cor_pmat(odat))

# Dimension reduction using Parallel Analysis
psych::fa.parallel(odat, fm = "pa")
# suggests two factors.

# efa in lavaan 
oefa <- lavaan::efa(data = map_dfc(odat, scale),
                    rotation = "oblimin",
                    nfactors = 1:3)

summary(oefa) # agrees that two factors is good.
summary(oefa[[2]])

# Solution identifies:
#   f1: Externalizing + Surgency + lack of Effortful Control
#   f2: Internalizing + SRS + Negative affect + lack of Surgency
#
# note: SU loads on both in opposite directions
#       EC loads on both in same direction, stronger to f1=externalizing.

# Extract scores for two factors:
odat_efa <- lavPredict(oefa[[2]]) %>% as.data.frame()

# check correlations:
ggcorrplot::ggcorrplot(cor(cbind(odat, odat_efa)))

# So for outcome we can either use raw variables (odat),
#                     or a two factor solution (odat_efa).

## Control Variables (e.g., sex, age, etc.)to be finalize by experts ####
# outcome adjustment
# iterate through outcome variables:
# raw outcomes (6 variables)
odat_adj <- map_dfc(odat,
                    ~ {
                      # fit model (outcomes by assessment at age8):
                      mod <- lm(y ~ ., data = data.frame(y = .x, ocon))
                      # extract residuals adjusting for CV
                      r <- resid(mod)
                      # add back the intercept == mean score:
                      r + coef(mod)[[1]]
                    }
)

# factor outcomes (f1 and f2)
odat_efa_adj <- map_dfc(odat_efa,
                        ~ {
                          # fit model (outcomes by assessment at age8):
                          mod <- lm(y ~ ., data = data.frame(y = .x, ocon))
                          # extract residuals adjusting for CV
                          r <- resid(mod)
                          # add back the intercept == mean score:
                          r + coef(mod)[[1]]
                        }
)

# predictor adjustment (to be finalise)
# Always adjust for Sex, Ethnicity, and low_mat_edu
# adjust some predictors for Age2 and Age4
# some predictors should not be adjusted for age.

# set adjustment for age4 as default:
preds_metadata <- data.frame(pred = colnames(pdat),type = "Age4")

# birth weight, IMD, gestational age should not be adjusted:
preds_metadata$type[colnames(pdat) %in% c("IMD", "bir_weight", "GA")] <- "none"

# Age 2 vars: Bayley, PARCA, Mchat
preds_metadata$type[starts_with("BY", vars = colnames(pdat))] <- "Age2"
preds_metadata$type[starts_with("PR", vars = colnames(pdat))] <- "Age2"
preds_metadata$type[starts_with("Mchat", vars = colnames(pdat))] <- "Age2"

#Adjust Age4 predictors: WPPSI, BRIEF, SRS, ECBQ, SDQ
pdat_adj <- map2_dfc(pdat,
                     preds_metadata$type,
                     ~ {
                       cvars <- switch(
                         .y,
                         none = pcon[, !(colnames(pcon) %in% c("Age2", "Age4"))],
                         Age2 = pcon[, !(colnames(pcon) %in% c("Age4"))],
                         Age4 = pcon[, !(colnames(pcon) %in% c("Age2"))]
                       )
                       # to debug / check:
                       cat(paste0("cvar_ncol: ", ncol(cvars), "\n"))
                       # fit model:
                       mod <- lm(y ~ ., data = data.frame(y = .x, cvars))
                       # extract residuals adjusting for CV
                       r <- resid(mod)
                       # add back the intercept == mean score:
                       r + coef(mod)[[1]]
                     }
)

# Write out data####

save(
  # the (pretty much) raw dataset:
  dt,
  # outcome variables with and without adjustment for control variables
  odat, odat_adj, 
  # Factors derived from outcome variables with and without adjustment for
  #    control variables:
  odat_efa, odat_efa_adj,
  # predictor variables with and without adjustment for control variables:
  pdat, pdat_adj,
  # control variables for reference:
  ocon, pcon,
  # Write above as compressed R objects to following file:
  file = "data/example_dCVnet_prepped.RData")


# dCVnet analysis start
load("data/example_dCVnet_prepped.RData")

## data description####
# # dt: the raw dataset
# # odat, odat_adj: 6 outcome variables with/without adjustment
# # odat_efa, odat_efa_adj: 2 representative outcome factors with/without adjustment
# # pdat, pdat_adj: predictors with/without adjustment
# # ocon, pcon: control variables for reference

## model strategy####
# # outcome:
#   - 5 separate models for variables in odat_adj
#   - 1 joint model (5 outcomes) for odat_adj
#   - 2 separate models for odat_efa_adj
#   - 1 joint model (2 outcomes) for odat_efa_adj
# # predictors:
#   - All predictors available in pdat_adj

## model parameters settings####
#   value of k for inner + outer loop
#     - defaults (10-fold) are sensible.
#   resolution of the grid for lambda and alpha
#     - more is always better for lambda 100 is standard, but we can run 1000.
#     - typically varying alpha makes little difference. Do 0.5

# model running####
# externalising only:
m_f1 <- dCVnet(
  y = odat_efa_adj$f1,
  data = pdat_adj,
  family = "gaussian",
  alphalist = 0.2,
  nlambda = 100,
  nrep_inner = 30,
  nrep_outer = 100
)
save(m_f1, file = "models/example_dCVnet_m_f1.RData")

# internalising only:
m_f2 <- dCVnet(
  y = odat_efa_adj$f2,
  data = pdat_adj,
  family = "gaussian",
  alphalist = 0.2,
  nlambda = 100,
  nrep_inner = 30,
  nrep_outer = 100
)
save(m_f2, file = "models/example_dCVnet_m_f2.RData")

# joint internalising and externalising factors:
m_fjoint <- dCVnet(
  y = as.matrix(odat_efa_adj),
  data = pdat_adj,
  family = "mgaussian",
  alphalist = 0.2,
  nlambda = 100,
  nrep_inner = 30,
  nrep_outer = 100
)
save(m_fjoint, file = "models/example_dCVnet_m_fjoint.RData")
