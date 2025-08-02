
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)

BacterialStructure <- read.csv("BacterialStructure.csv", fileEncoding = "latin1")
# Check data
head(BacterialStructure)

# 1. The number of Obversation
total_number <- nrow(BacterialStructure)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 343

# 2. The number of Study
unique_studyid_number <- length(unique(BacterialStructure$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 81


#### 3. Overall effect size
total_effect_model <- rma.mv(yi = RR, 
                             V = Vi, 
                             random = ~ 1 | StudyID,  # StudyID is radom factor
                             data = BacterialStructure, 
                             method = "REML")
# The results of Overall effect size
summary(total_effect_model)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -605.8240  1211.6480  1215.6480  1223.3176  1215.6834   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5351  0.7315     81     no  StudyID 
# Test for Heterogeneity:
# Q(df = 342) = 11524.9783, p-val < .0001
# Model Results:
# estimate      se     zval    pval   ci.lb   ci.ub      
#   0.9740  0.0833  11.6865  <.0001  0.8106  1.1373  *** 



#### 4. Bootstrap Method
# Custom functions
boot_fun <- function(data, indices) {
  library(metafor)
  data_boot <- data[indices, ]
  model <- rma.mv(yi = RR, V = Vi, random = ~ 1 | StudyID, data = data_boot, method = "REML")
  return(coef(model))
}
# Set parallel calculation
numCores <- detectCores() - 10
cl <- makeCluster(numCores)
clusterEvalQ(cl, library(metafor))
# 
set.seed(123)
boot_results1 <- boot(data = BacterialStructure, statistic = boot_fun, R = 1000, parallel = "snow", ncpus = numCores, cl = cl)
stopCluster(cl) 
# Calculation
for (i in 1:length(boot_results1$t0)) {
  boot_ci <- boot.ci(boot_results1, index = i, type = "bca")
  estimate <- boot_results1$t0[i]  
  cat("Estimate for coefficient", i, ":", estimate, "\n")
  if (!is.null(boot_ci$bca)) {
    cat("95% BCa CI for coefficient", i, ":", boot_ci$bca[4:5], "\n\n")
  } else {
    cat("95% BCa CI for coefficient", i, ": NA (insufficient bootstrap resamples)\n\n")
  }
}
# Estimate for coefficient 1 : 0.9739586 
# 95% BCa CI for coefficient 1 : 0.9422598 1.054942  


#### 5. Funnel Plot
simple_model <- rma(yi = RR, 
                    vi = Vi, 
                    data = BacterialStructure, 
                    method = "REML")
#### 
funnel(simple_model)
# Output  6 * 6
#### Egger's test
regtest(simple_model)
# Regression Test for Funnel Plot Asymmetry
# Model:     mixed-effects meta-regression model
# Predictor: standard error
# Test for Funnel Plot Asymmetry: z = -1.1522, p = 0.2492
# Limit Estimate (as sei -> 0):   b =  0.7686 (CI: 0.5675, 0.9698)

#  Rosenthal’s Fail-Safe N
# This method estimates how many missing studies with null effect 
# would be needed to make the overall effect non-significant
fsn_rosenthal <- fsn(x = simple_model, type = "Rosenthal")
# Print the FSN result
print(fsn_rosenthal)
# Fail-safe N Calculation Using the General Approach
# Average Effect Size:         0.6597 (with file drawer: 0.0094)
# Amount of Heterogeneity:     0.5006 (with file drawer: 0.5071)
# Observed Significance Level: <.0001 (with file drawer: 0.0500)
# Target Significance Level:   0.05
# Fail-safe N: 22392


#### 8. Subgroup analysis
### 8.1 LegumeNonlegume
BacterialStructure_filteredLegumeNonlegume <- subset(BacterialStructure, LegumeNonlegume %in% c("Legume to Non-legume", "Non-legume to Legume", "Non-legume to Non-legume"))
#
BacterialStructure_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(BacterialStructure_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredLegumeNonlegume %>%
  group_by(LegumeNonlegume) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   LegumeNonlegume          Observations Unique_StudyID
#   <fct>                           <int>          <int>
# 1 Legume to Non-legume              125             22
# 2 Non-legume to Legume               96             26
# 3 Non-legume to Non-legume          119             44

overall_model_BacterialStructure_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = BacterialStructure_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 340; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -592.8821  1185.7642  1193.7642  1209.0446  1193.8847   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5850  0.7648     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 337) = 10098.9580, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 148.1886, p-val < .0001
# Model Results:
#                                          estimate      se     zval    pval   ci.lb   ci.ub      
# LegumeNonlegumeLegume to Non-legume        1.1909  0.0984  12.1051  <.0001  0.9981  1.3837  *** 
# LegumeNonlegumeNon-legume to Legume        1.0439  0.0927  11.2576  <.0001  0.8621  1.2256  *** 
# LegumeNonlegumeNon-legume to Non-legume    0.8533  0.0921   9.2606  <.0001  0.6727  1.0339  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredLegumeNonlegume)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
    # LegumeNonlegumeLegume to Non-legume     LegumeNonlegumeNon-legume to Legume LegumeNonlegumeNon-legume to Non-legume 
    #                                 "a"                                     "b"                                     "c" 


### 8.2 AMnonAM
BacterialStructure_filteredAMnonAM <- subset(BacterialStructure, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
BacterialStructure_filteredAMnonAM$AMnonAM <- droplevels(factor(BacterialStructure_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredAMnonAM %>%
  group_by(AMnonAM) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   AMnonAM     Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 AM to AM             272             66
# 2 AM to nonAM           48             14
# 3 nonAM to AM           22              7

overall_model_BacterialStructure_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = BacterialStructure_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 342; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -595.4218  1190.8435  1198.8435  1214.1475  1198.9633   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5238  0.7237     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 339) = 11163.6392, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 156.4414, p-val < .0001
# Model Results:
#                     estimate      se     zval    pval   ci.lb   ci.ub      
# AMnonAMAM to AM       0.9781  0.0875  11.1838  <.0001  0.8067  1.1495  *** 
# AMnonAMAM to nonAM    0.7568  0.0999   7.5789  <.0001  0.5611  0.9526  *** 
# AMnonAMnonAM to AM    1.3567  0.2773   4.8919  <.0001  0.8131  1.9002  ***   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredAMnonAM)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
   # AMnonAMAM to AM AMnonAMAM to nonAM AMnonAMnonAM to AM 
   #             "a"                "b"                "a" 


### 8.3 C3C4
BacterialStructure_filteredC3C4 <- subset(BacterialStructure, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
BacterialStructure_filteredC3C4$C3C4 <- droplevels(factor(BacterialStructure_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredC3C4 %>%
  group_by(C3C4) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   C3C4     Observations Unique_StudyID
#   <fct>           <int>          <int>
# 1 C3 to C3          139             47
# 2 C3 to C4          137             33
# 3 C4 to C3           58             14
overall_model_BacterialStructure_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = BacterialStructure_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 334; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -567.2086  1134.4173  1142.4173  1157.6257  1142.5400   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5532  0.7438     79     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 331) = 10994.1890, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 150.4293, p-val < .0001
# Model Results:
#               estimate      se     zval    pval   ci.lb   ci.ub      
# C3C4C3 to C3    0.8903  0.0890  10.0034  <.0001  0.7158  1.0647  *** 
# C3C4C3 to C4    1.1117  0.0908  12.2430  <.0001  0.9337  1.2897  *** 
# C3C4C4 to C3    0.9990  0.0966  10.3461  <.0001  0.8097  1.1882  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredC3C4)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredC3C4)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# C3C4C3 to C3 C3C4C3 to C4 C3C4C4 to C3 
#          "a"          "b"          "a" 


### 8.4 Annual_Pere
BacterialStructure_filteredAnnual_Pere <- subset(BacterialStructure, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
BacterialStructure_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(BacterialStructure_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredAnnual_Pere %>%
  group_by(Annual_Pere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Annual_Pere            Observations Unique_StudyID
#   <fct>                         <int>          <int>
# 1 Annual to Annual                246             56
# 2 Annual to Perennial              29             12
# 3 Perennial to Annual              54             14
# 4 Perennial to Perennial           13              6

overall_model_BacterialStructure_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = BacterialStructure_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 342; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -582.2088  1164.4177  1174.4177  1193.5329  1174.5984   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4923  0.7017     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 338) = 9689.1585, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 192.2317, p-val < .0001
# Model Results:
#                                    estimate      se    zval    pval   ci.lb   ci.ub      
# Annual_PereAnnual to Annual          0.8602  0.0905  9.5023  <.0001  0.6827  1.0376  *** 
# Annual_PereAnnual to Perennial       1.0031  0.1079  9.2950  <.0001  0.7916  1.2146  *** 
# Annual_PerePerennial to Annual       1.1922  0.1617  7.3704  <.0001  0.8751  1.5092  *** 
# Annual_PerePerennial to Perennial    1.7127  0.1736  9.8643  <.0001  1.3724  2.0530  ***    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredAnnual_Pere)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
      # Annual_PereAnnual to Annual    Annual_PereAnnual to Perennial    Annual_PerePerennial to Annual Annual_PerePerennial to Perennial 
      #                         "a"                               "b"                              "ab"                               "c" 


### 8.6 PlantStageSubgroup
BacterialStructure_filteredPlantStageSubgroup <- subset(BacterialStructure, PlantStageSubgroup %in% c("Vegetative stage","Reproductive stage", "Maturity stage","Harvest"))
#
BacterialStructure_filteredPlantStageSubgroup$PlantStageSubgroup <- droplevels(factor(BacterialStructure_filteredPlantStageSubgroup$PlantStageSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredPlantStageSubgroup %>%
  group_by(PlantStageSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   PlantStageSubgroup Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Harvest                     138             39
# 2 Maturity stage               52             16
# 3 Reproductive stage           58             14
# 4 Vegetative stage             30              6

overall_model_BacterialStructure_filteredPlantStageSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + PlantStageSubgroup, random = ~ 1 | StudyID, data = BacterialStructure_filteredPlantStageSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredPlantStageSubgroup)
# Multivariate Meta-Analysis Model (k = 278; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -428.2889   856.5777   866.5777   884.6433   866.8016   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5967  0.7724     67     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 274) = 9719.1084, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 111.7400, p-val < .0001
# Model Results:
#                                       estimate      se     zval    pval   ci.lb   ci.ub      
# PlantStageSubgroupHarvest               1.0344  0.0990  10.4522  <.0001  0.8405  1.2284  *** 
# PlantStageSubgroupMaturity stage        0.9794  0.1120   8.7417  <.0001  0.7598  1.1990  *** 
# PlantStageSubgroupReproductive stage    0.8076  0.1160   6.9623  <.0001  0.5802  1.0349  *** 
# PlantStageSubgroupVegetative stage      1.0006  0.1112   9.0010  <.0001  0.7827  1.2184  ***   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredPlantStageSubgroup)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredPlantStageSubgroup)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
           # PlantStageSubgroupHarvest     PlantStageSubgroupMaturity stage PlantStageSubgroupReproductive stage   PlantStageSubgroupVegetative stage 
           #                       "a"                                 "ab"                                  "b"                                  "a" 

### 8.7 Bulk_Rhizosphere
BacterialStructure_filteredBulk_Rhizosphere <- subset(BacterialStructure, Bulk_Rhizosphere %in% c("Non-Rhizosphere", "Rhizosphere"))
#
BacterialStructure_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(BacterialStructure_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredBulk_Rhizosphere %>%
  group_by(Bulk_Rhizosphere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Bulk_Rhizosphere Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 Non-Rhizosphere           223             57
# 2 Rhizosphere               120             40

overall_model_BacterialStructure_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = BacterialStructure_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -605.3236  1210.6473  1216.6473  1228.1429  1216.7185  
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5345  0.7311     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 341) = 11503.3693, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 139.0494, p-val < .0001
# Model Results:
#                                  estimate      se     zval    pval   ci.lb   ci.ub      
# Bulk_RhizosphereNon-Rhizosphere    0.9601  0.0838  11.4576  <.0001  0.7958  1.1243  *** 
# Bulk_RhizosphereRhizosphere        0.9953  0.0845  11.7845  <.0001  0.8298  1.1609  ***   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredBulk_Rhizosphere)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# Bulk_RhizosphereNon-Rhizosphere     Bulk_RhizosphereRhizosphere 
#                             "a"                             "a" 


### 8.8 Soil_texture
BacterialStructure_filteredSoil_texture <- subset(BacterialStructure, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
BacterialStructure_filteredSoil_texture$Soil_texture <- droplevels(factor(BacterialStructure_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredSoil_texture %>%
  group_by(Soil_texture) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Soil_texture Observations Unique_StudyID
#   <fct>               <int>          <int>
# 1 Coarse                 62             13
# 2 Fine                   98             16
# 3 Medium                 96             28

overall_model_BacterialStructure_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = BacterialStructure_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 256; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -441.7000   883.4000   891.4000   905.5335   891.5612   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5261  0.7253     57     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 253) = 9445.6198, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 108.4918, p-val < .0001
# Model Results:
#                     estimate      se    zval    pval   ci.lb   ci.ub      
# Soil_textureCoarse    1.0390  0.2041  5.0906  <.0001  0.6390  1.4391  *** 
# Soil_textureFine      1.3058  0.1875  6.9634  <.0001  0.9383  1.6733  *** 
# Soil_textureMedium    0.8169  0.1399  5.8385  <.0001  0.5426  1.0911  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredSoil_texture)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# Soil_textureCoarse   Soil_textureFine Soil_textureMedium 
#               "ab"                "a"                "b"  

### 8.9 Tillage
BacterialStructure_filteredTillage <- subset(BacterialStructure, Tillage %in% c("Tillage", "No_tillage"))
#
BacterialStructure_filteredTillage$Tillage <- droplevels(factor(BacterialStructure_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredTillage %>%
  group_by(Tillage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Tillage    Observations Unique_StudyID
#   <fct>             <int>          <int>
# 1 No_tillage           27              3
# 2 Tillage              27              9
overall_model_BacterialStructure_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = BacterialStructure_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredTillage)
# Multivariate Meta-Analysis Model (k = 54; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -53.3605  106.7210  112.7210  118.5747  113.2210   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1353  0.3678     11     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 52) = 322.8601, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 15.6900, p-val = 0.0004
# Model Results:
#                    estimate      se    zval    pval   ci.lb   ci.ub      
# TillageNo_tillage    0.3401  0.1555  2.1874  0.0287  0.0354  0.6449    * 
# TillageTillage       0.4784  0.1208  3.9593  <.0001  0.2416  0.7152  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredTillage)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredTillage)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# TillageNo_tillage    TillageTillage 
#               "a"               "a"  


### 8.10 Straw_retention
BacterialStructure_filteredStraw_retention <- subset(BacterialStructure, Straw_retention %in% c("Retention", "No_retention"))
#
BacterialStructure_filteredStraw_retention$Straw_retention <- droplevels(factor(BacterialStructure_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredStraw_retention %>%
  group_by(Straw_retention) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Straw_retention Observations Unique_StudyID
#   <fct>                  <int>          <int>
# 1 No_retention              18              9
# 2 Retention                 29              9

overall_model_BacterialStructure_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = BacterialStructure_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 47; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -83.4495  166.8991  172.8991  178.3191  173.4844   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4513  0.6718     16     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 45) = 445.3796, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 161.0150, p-val < .0001
# Model Results:
#                              estimate      se    zval    pval    ci.lb   ci.ub      
# Straw_retentionNo_retention    1.1715  0.1779  6.5855  <.0001   0.8228  1.5202  *** 
# Straw_retentionRetention       0.1964  0.1775  1.1067  0.2684  -0.1514  0.5443     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredStraw_retention)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# Straw_retentionNo_retention    Straw_retentionRetention 
#                         "a"                         "b" 


### 8.12 RotationcyclesSubgroup
BacterialStructure_filteredRotationcyclesSubgroup <- subset(BacterialStructure, RotationcyclesSubgroup %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
BacterialStructure_filteredRotationcyclesSubgroup$RotationcyclesSubgroup <- droplevels(factor(BacterialStructure_filteredRotationcyclesSubgroup$RotationcyclesSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredRotationcyclesSubgroup %>%
  group_by(RotationcyclesSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   RotationcyclesSubgroup Observations Unique_StudyID
#   <fct>                         <int>          <int>
# 1 D1                               79             34
# 2 D1-3                             86             27
# 3 D10                              81             11
# 4 D3-5                             45              8
# 5 D5-10                            52             14
overall_model_BacterialStructure_filteredRotationcyclesSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + RotationcyclesSubgroup, random = ~ 1 | StudyID, data = BacterialStructure_filteredRotationcyclesSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredRotationcyclesSubgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -572.0276  1144.0552  1156.0552  1178.9935  1156.3090   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5725  0.7566     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 338) = 9869.1940, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 195.2356, p-val < .0001
# Model Results:
#                              estimate      se     zval    pval   ci.lb   ci.ub      
# RotationcyclesSubgroupD1       1.2850  0.0975  13.1819  <.0001  1.0939  1.4761  *** 
# RotationcyclesSubgroupD1-3     0.7140  0.0978   7.3024  <.0001  0.5223  0.9056  *** 
# RotationcyclesSubgroupD10      0.7970  0.1146   6.9558  <.0001  0.5724  1.0216  *** 
# RotationcyclesSubgroupD3-5     0.8176  0.1084   7.5401  <.0001  0.6051  1.0301  *** 
# RotationcyclesSubgroupD5-10    0.9020  0.1062   8.4921  <.0001  0.6938  1.1102  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredRotationcyclesSubgroup)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredRotationcyclesSubgroup)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
   # RotationcyclesSubgroupD1  RotationcyclesSubgroupD1-3   RotationcyclesSubgroupD10  RotationcyclesSubgroupD3-5 RotationcyclesSubgroupD5-10 
   #                      "a"                         "b"                        "bc"                        "bc"                         "c"  


### 8.13 DurationSubgroup
BacterialStructure_filteredDurationSubgroup <- subset(BacterialStructure, DurationSubgroup %in% c("D1", "D2", "D3", "D4", "D5", "D6-10", "D11-20", "D20-30", "D30"))
#
BacterialStructure_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(BacterialStructure_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredDurationSubgroup %>%
  group_by(DurationSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   DurationSubgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D1                          9              7
# 2 D11-20                     38             12
# 3 D2                         37             17
# 4 D20-30                     27              8
# 5 D3                         40             14
# 6 D30                       107              4
# 7 D4                         30             10
# 8 D5                          2              2
# 9 D6-10                      53             18

overall_model_BacterialStructure_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = BacterialStructure_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -570.8843  1141.7685  1161.7685  1199.8799  1162.4496   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5249  0.7245     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 334) = 8717.0673, p-val < .0001
# Test of Moderators (coefficients 1:9):
# QM(df = 9) = 201.4640, p-val < .0001
# Model Results:
#                         estimate      se    zval    pval    ci.lb   ci.ub      
# DurationSubgroupD1        0.9425  0.1714  5.4989  <.0001   0.6066  1.2785  *** 
# DurationSubgroupD11-20    1.0953  0.1120  9.7805  <.0001   0.8758  1.3147  *** 
# DurationSubgroupD2        1.0075  0.1169  8.6154  <.0001   0.7783  1.2367  *** 
# DurationSubgroupD20-30    1.0739  0.1276  8.4185  <.0001   0.8239  1.3239  *** 
# DurationSubgroupD3        0.6536  0.1122  5.8234  <.0001   0.4336  0.8735  *** 
# DurationSubgroupD30       0.6040  0.3666  1.6478  0.0994  -0.1144  1.3225    . 
# DurationSubgroupD4        1.5234  0.2334  6.5272  <.0001   1.0659  1.9808  *** 
# DurationSubgroupD5        1.3605  0.2864  4.7507  <.0001   0.7992  1.9218  *** 
# DurationSubgroupD6-10     0.7634  0.1251  6.1020  <.0001   0.5182  1.0086  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredDurationSubgroup)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
    # DurationSubgroupD1 DurationSubgroupD11-20     DurationSubgroupD2 DurationSubgroupD20-30     DurationSubgroupD3    DurationSubgroupD30 
    #              "abc"                   "ad"                   "ab"                   "ad"                    "c"                  "abc" 
    # DurationSubgroupD4     DurationSubgroupD5  DurationSubgroupD6-10 
    #                "d"                   "ad"                   "bc" 



### 8.14 SpeciesRichnessSubgroup
BacterialStructure_filteredSpeciesRichnessSubgroup <- subset(BacterialStructure, SpeciesRichnessSubgroup %in% c("R2", "R3", "R4"))
#
BacterialStructure_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup <- droplevels(factor(BacterialStructure_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredSpeciesRichnessSubgroup %>%
  group_by(SpeciesRichnessSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
  # SpeciesRichnessSubgroup Observations Unique_StudyID
#   <fct>                          <int>          <int>
# 1 R2                               278             71
# 2 R3                                52             15
# 3 R4                                13              3

overall_model_BacterialStructure_filteredSpeciesRichnessSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichnessSubgroup, random = ~ 1 | StudyID, data = BacterialStructure_filteredSpeciesRichnessSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredSpeciesRichnessSubgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -603.6019  1207.2037  1215.2037  1230.5195  1215.3231   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5329  0.7300     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 340) = 11193.8985, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 139.0924, p-val < .0001
# Model Results:
#                            estimate      se     zval    pval    ci.lb   ci.ub      
# SpeciesRichnessSubgroupR2    0.9981  0.0851  11.7243  <.0001   0.8313  1.1650  *** 
# SpeciesRichnessSubgroupR3    0.9469  0.0965   9.8086  <.0001   0.7577  1.1361  *** 
# SpeciesRichnessSubgroupR4    0.5422  0.4243   1.2779  0.2013  -0.2894  1.3737     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredSpeciesRichnessSubgroup)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredSpeciesRichnessSubgroup)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# SpeciesRichnessSubgroupR2 SpeciesRichnessSubgroupR3 SpeciesRichnessSubgroupR4 
#                       "a"                       "a"                      "a" 


### 8.15 Primer
BacterialStructure_filteredPrimer <- subset(BacterialStructure, Primer %in% c("V1-V3", "V1-V4", "V1-V5", "V3-V4", "V3-V5", "V4", "V4-V5", "V5-V7", "V6-V8", "V9", "Full length"))
#
BacterialStructure_filteredPrimer$Primer <- droplevels(factor(BacterialStructure_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredPrimer %>%
  group_by(Primer) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Primer      Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 Full length            3              1
# 2 V1-V4                  4              1
# 3 V1-V5                  6              1
# 4 V3-V4                124             37
# 5 V3-V5                  2              1
# 6 V4                   157             22
# 7 V4-V5                 38             14
# 8 V5-V7                  3              1

overall_model_BacterialStructure_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = BacterialStructure_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 337; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -587.7933  1175.5866  1193.5866  1227.7511  1194.1509   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5712  0.7558     78     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 329) = 9853.2122, p-val < .0001
# Test of Moderators (coefficients 1:8):
# QM(df = 8) = 119.8282, p-val < .0001
# Model Results:
#                    estimate      se    zval    pval    ci.lb   ci.ub      
# PrimerFull length    0.3710  0.7636  0.4858  0.6271  -1.1257  1.8676      
# PrimerV1-V4          1.4777  0.7623  1.9384  0.0526  -0.0164  2.9719    . 
# PrimerV1-V5          0.5688  0.7607  0.7478  0.4546  -0.9221  2.0598      
# PrimerV3-V4          0.9295  0.1270  7.3164  <.0001   0.6805  1.1785  *** 
# PrimerV3-V5          0.8899  0.7593  1.1719  0.2412  -0.5984  2.3782      
# PrimerV4             0.9536  0.1666  5.7238  <.0001   0.6271  1.2801  *** 
# PrimerV4-V5          1.0257  0.2055  4.9924  <.0001   0.6230  1.4284  *** 
# PrimerV5-V7          1.2428  0.7582  1.6392  0.1012  -0.2432  2.7287      
  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredPrimer)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredPrimer)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# PrimerFull length       PrimerV1-V4       PrimerV1-V5       PrimerV3-V4       PrimerV3-V5          PrimerV4       PrimerV4-V5       PrimerV5-V7 
#               "a"               "a"               "a"               "a"               "a"               "a"               "a"               "a" 



### 8.16 Latitude_Subgroup
BacterialStructure_filteredLatitude_Subgroup <- subset(BacterialStructure, Latitude_Subgroup %in% c("La20", "La20-40", "La40"))
#
BacterialStructure_filteredLatitude_Subgroup$Latitude_Subgroup <- droplevels(factor(BacterialStructure_filteredLatitude_Subgroup$Latitude_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredLatitude_Subgroup %>%
  group_by(Latitude_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Latitude_Subgroup Observations Unique_StudyID
#   <fct>                    <int>          <int>
# 1 La20                        31              6
# 2 La20-40                     95             43
# 3 La40                       217             32

overall_model_BacterialStructure_filteredLatitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Latitude_Subgroup, random = ~ 1 | StudyID, data = BacterialStructure_filteredLatitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredLatitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -599.7217  1199.4434  1207.4434  1222.7591  1207.5628   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5022  0.7087     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 340) = 10585.1240, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 152.0068, p-val < .0001
# Model Results:
#                           estimate      se     zval    pval   ci.lb   ci.ub      
# Latitude_SubgroupLa20       1.1831  0.3016   3.9224  <.0001  0.5919  1.7743  *** 
# Latitude_SubgroupLa20-40    1.1427  0.1110  10.2913  <.0001  0.9251  1.3604  *** 
# Latitude_SubgroupLa40       0.7104  0.1282   5.5418  <.0001  0.4591  0.9616  ***    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredLatitude_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredLatitude_Subgroup)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# Latitude_SubgroupLa20 Latitude_SubgroupLa20-40    Latitude_SubgroupLa40 
#                   "ab"                      "a"                      "b" 


### 8.17 Longitude_Subgroup
BacterialStructure_filteredLongitude_Subgroup <- subset(BacterialStructure, Longitude_Subgroup %in% c("Lo-180-0", "Lo-0-180"))
#
BacterialStructure_filteredLongitude_Subgroup$Longitude_Subgroup <- droplevels(factor(BacterialStructure_filteredLongitude_Subgroup$Longitude_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredLongitude_Subgroup %>%
  group_by(Longitude_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Longitude_Subgroup Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Lo-0-180                    217             74
# 2 Lo-180-0                    126              7

overall_model_BacterialStructure_filteredLongitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Subgroup, random = ~ 1 | StudyID, data = BacterialStructure_filteredLongitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredLongitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -600.8110  1201.6221  1207.6221  1219.1177  1207.6933   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4983  0.7059     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 341) = 10487.6706, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 152.5633, p-val < .0001
# Model Results:
#                             estimate      se     zval    pval    ci.lb   ci.ub      
# Longitude_SubgroupLo-0-180    1.0373  0.0843  12.2988  <.0001   0.8720  1.2026  *** 
# Longitude_SubgroupLo-180-0    0.3104  0.2721   1.1409  0.2539  -0.2229  0.8438      
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredLongitude_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredLongitude_Subgroup)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# Longitude_SubgroupLo-0-180 Longitude_SubgroupLo-180-0 
#                        "a"                        "b" 



### 8.18 MAPmean_Subgroup
BacterialStructure_filteredMAPmean_Subgroup <- subset(BacterialStructure, MAPmean_Subgroup %in% c("MAP600", "MAP600-1200", "MAP1200"))
#
BacterialStructure_filteredMAPmean_Subgroup$MAPmean_Subgroup <- droplevels(factor(BacterialStructure_filteredMAPmean_Subgroup$MAPmean_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredMAPmean_Subgroup %>%
  group_by(MAPmean_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MAPmean_Subgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 MAP1200                    76             25
# 2 MAP600                    114             32
# 3 MAP600-1200               153             25

overall_model_BacterialStructure_filteredMAPmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean_Subgroup, random = ~ 1 | StudyID, data = BacterialStructure_filteredMAPmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredMAPmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -603.2439  1206.4877  1214.4877  1229.8035  1214.6071   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5360  0.7321     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 340) = 9655.2324, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 137.6142, p-val < .0001
# Model Results:
#                              estimate      se    zval    pval   ci.lb   ci.ub      
# MAPmean_SubgroupMAP1200        1.1076  0.1516  7.3069  <.0001  0.8105  1.4047  *** 
# MAPmean_SubgroupMAP600         0.8991  0.1087  8.2681  <.0001  0.6860  1.1122  *** 
# MAPmean_SubgroupMAP600-1200    0.9378  0.1144  8.1972  <.0001  0.7135  1.1620  *** 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredMAPmean_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredMAPmean_Subgroup)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# MAPmean_SubgroupMAP1200      MAPmean_SubgroupMAP600 MAPmean_SubgroupMAP600-1200 
#                     "a"                         "a"                         "a" 



### 8.19 MATmean_Subgroup
BacterialStructure_filteredMATmean_Subgroup <- subset(BacterialStructure, MATmean_Subgroup %in% c("MAT8", "MAT8-15", "MAT15"))
#
BacterialStructure_filteredMATmean_Subgroup$MATmean_Subgroup <- droplevels(factor(BacterialStructure_filteredMATmean_Subgroup$MATmean_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredMATmean_Subgroup %>%
  group_by(MATmean_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MATmean_Subgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 MAT15                      90             34
# 2 MAT8                      169             34
# 3 MAT8-15                    84             14

overall_model_BacterialStructure_filteredMATmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Subgroup, random = ~ 1 | StudyID, data = BacterialStructure_filteredMATmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredMATmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -602.9996  1205.9992  1213.9992  1229.3149  1214.1186   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5326  0.7298     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 340) = 10607.8937, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 139.0385, p-val < .0001
# Model Results:
#                          estimate      se    zval    pval   ci.lb   ci.ub      
# MATmean_SubgroupMAT15      1.0953  0.1291  8.4864  <.0001  0.8423  1.3482  *** 
# MATmean_SubgroupMAT8       0.8700  0.1128  7.7159  <.0001  0.6490  1.0910  *** 
# MATmean_SubgroupMAT8-15    0.9311  0.1307  7.1236  <.0001  0.6749  1.1873  *** 



# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredMATmean_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredMATmean_Subgroup)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# MATmean_SubgroupMAT15    MATmean_SubgroupMAT8 MATmean_SubgroupMAT8-15 
#                   "a"                     "a"                     "a" 



### 8.20 DurationSubgroup
BacterialStructure_filteredDurationSubgroup <- subset(BacterialStructure, DurationSubgroup %in% c("D5", "D5-10", "D10-20", "D20-30", "D30"))
#
BacterialStructure_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(BacterialStructure_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialStructure_filteredDurationSubgroup %>%
  group_by(DurationSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   DurationSubgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D10-20                     38             12
# 2 D20-30                     27              8
# 3 D30                       107              4
# 4 D5                        118             44
# 5 D5-10                      53             18

overall_model_BacterialStructure_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = BacterialStructure_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialStructure_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -591.4051  1182.8103  1194.8103  1217.7486  1195.0641   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5881  0.7668     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 338) = 10424.8730, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 148.1719, p-val < .0001
# Model Results:
#                         estimate      se     zval    pval    ci.lb   ci.ub      
# DurationSubgroupD10-20    1.2710  0.1110  11.4458  <.0001   1.0533  1.4886  *** 
# DurationSubgroupD20-30    1.2151  0.1282   9.4773  <.0001   0.9638  1.4664  *** 
# DurationSubgroupD30       0.6040  0.3875   1.5587  0.1191  -0.1555  1.3635      
# DurationSubgroupD5        0.9561  0.0996   9.5974  <.0001   0.7609  1.1514  *** 
# DurationSubgroupD5-10     0.8262  0.1248   6.6230  <.0001   0.5817  1.0707  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialStructure_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_BacterialStructure_filteredDurationSubgroup)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# DurationSubgroupD10-20 DurationSubgroupD20-30    DurationSubgroupD30     DurationSubgroupD5  DurationSubgroupD5-10 
#                  "a"                    "a"                   "ab"                    "b"                    "b" 







#### 9. Linear Mixed Effect Model
# 
BacterialStructure$Wr <- 1 / BacterialStructure$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialStructure)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialStructure)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialStructure)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialStructure)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialStructure)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialStructure)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialStructure)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialStructure)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(Species_Richness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = BacterialStructure)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = BacterialStructure)

# Define
models <- list(
  Model1 = Model1,
  Model2 = Model2,
  Model3 = Model3,
  Model4 = Model4,
  Model5 = Model5,
  Model6 = Model6,
  Model7 = Model7,
  Model8 = Model8,
  Model9 = Model9,
  Model10 = Model10
)
# 创建一个数据框来存储每个模型的AIC、BIC和logLik
model_comparison <- data.frame(
  Model = names(models),
  AIC = sapply(models, AIC),
  BIC = sapply(models, BIC),
  logLik = sapply(models, logLik)
)
# 查看结果
print(model_comparison)
# 添加 Marginal R² 和 Conditional R²
model_comparison$Marginal_R2 <- sapply(models, function(m) r.squaredGLMM(m)[1])
model_comparison$Conditional_R2 <- sapply(models, function(m) r.squaredGLMM(m)[2])
# 查看更新后的模型比较结果
print(model_comparison)
#           Model      AIC      BIC    logLik Marginal_R2 Conditional_R2
# Model1   Model1 590.3140 613.3403 -289.1570 0.001153027     0.08294768
# Model2   Model2 589.9993 613.0257 -288.9997 0.001650687     0.08235484
# Model3   Model3 588.9369 611.9633 -288.4685 0.001393869     0.08376982
# Model4   Model4 590.5099 613.5363 -289.2550 0.001074445     0.08294885
# Model5   Model5 590.7574 613.7838 -289.3787 0.002020771     0.08143270
# Model6   Model6 590.9441 613.9705 -289.4721 0.001966584     0.08141620
# Model7   Model7 589.8024 612.8288 -288.9012 0.001710457     0.08238076
# Model8   Model8 589.1321 612.1585 -288.5661 0.001305849     0.08375449
# Model9   Model9 601.6885 640.0659 -290.8443 0.008255703     0.09263795
# Model10 Model10 604.1057 642.4830 -292.0528 0.005656132     0.09052335

##### Model 3 is the best model
summary(Model3)
# Number of obs: 343, groups:  StudyID, 81
anova(Model3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 13.9433 13.9433     1 336.73  2.7542 0.09793 .
# scale(Species_Richness)      8.1625  8.1625     1 311.87  1.6123 0.20511  
# scale(Duration)              2.9572  2.9572     1 292.42  0.5841 0.44531  

#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = BacterialStructure)
summary(ModelpH)
# Number of obs: 142, groups:  StudyID, 55
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 23.3684 23.3684     1 129.123  4.0684 0.04577 *
# scale(Species_Richness)      5.3393  5.3393     1 110.508  0.9296 0.33708  
# scale(Duration)              7.4187  7.4187     1 135.825  1.2916 0.25776  
# scale(pHCK)                 12.0505 12.0505     1  70.543  2.0980 0.15193  

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = BacterialStructure)
summary(ModelSOC)
# Number of obs: 148, groups:  StudyID, 54
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 24.5777 24.5777     1 116.043  3.3531 0.06964 .
# scale(Species_Richness)     11.9045 11.9045     1 141.665  1.6241 0.20460  
# scale(Duration)              4.5353  4.5353     1 142.917  0.6187 0.43282  
# scale(SOCCK)                 9.7392  9.7392     1  58.262  1.3287 0.25374  

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = BacterialStructure)
summary(ModelTN)
# Number of obs: 97, groups:  StudyID, 37
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 20.4799 20.4799     1 79.387  3.6683 0.05906 .
# scale(Species_Richness)     12.4602 12.4602     1 81.007  2.2318 0.13908  
# scale(Duration)             12.1912 12.1912     1 86.514  2.1836 0.14312  
# scale(TNCK)                  0.9345  0.9345     1 32.396  0.1674 0.68514  


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = BacterialStructure)
summary(ModelNO3)
# Number of obs: 50, groups:  StudyID, 19
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# scale(log(Rotation_cycles)) 34.357  34.357     1 44.437  9.7448 0.003158 **
# scale(Species_Richness)      0.883   0.883     1 29.962  0.2504 0.620454   
# scale(Duration)             25.608  25.608     1 30.104  7.2633 0.011403 * 
# scale(NO3CK)                 0.040   0.040     1 16.250  0.0113 0.916516   

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = BacterialStructure)
summary(ModelNH4)
# Number of obs: 46, groups:  StudyID, 17
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.065692 0.065692     1 23.010  0.0230 0.8808
# scale(Species_Richness)     0.176288 0.176288     1 28.805  0.0617 0.8055
# scale(Duration)             0.315747 0.315747     1 40.460  0.1106 0.7412
# scale(NH4CK)                0.241883 0.241883     1 12.951  0.0847 0.7756

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = BacterialStructure)
summary(ModelAP)
# Number of obs: 120, groups:  StudyID, 48
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 11.3925 11.3925     1 104.501  1.8555 0.1761
# scale(Species_Richness)      6.3269  6.3269     1 105.083  1.0305 0.3124
# scale(Duration)              6.9034  6.9034     1 109.573  1.1243 0.2913
# scale(APCK)                  0.4012  0.4012     1  92.439  0.0653 0.7988

#### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = BacterialStructure)
summary(ModelAK)
# Number of obs: 99, groups:  StudyID, 40
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 11.5743 11.5743     1 92.862  1.7136 0.19374  
# scale(Species_Richness)      2.8466  2.8466     1 89.967  0.4215 0.51787  
# scale(Duration)              3.6370  3.6370     1 51.790  0.5385 0.46638  
# scale(AKCK)                 24.3416 24.3416     1 87.717  3.6039 0.06093 .

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = BacterialStructure)
summary(ModelAN)
# Number of obs: 80, groups:  StudyID, 25
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 31.424  31.424     1 60.276  5.4004 0.02352 *
# scale(Species_Richness)      0.266   0.266     1 73.401  0.0457 0.83139  
# scale(Duration)             36.360  36.360     1 39.815  6.2486 0.01665 *
# scale(ANCK)                 13.339  13.339     1 33.339  2.2923 0.13944  

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = BacterialStructure)
summary(ModelLatitude)
# Number of obs: 343, groups:  StudyID, 81
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 17.5154 17.5154     1 335.91  3.4598 0.06375 .
# scale(Species_Richness)      6.4383  6.4383     1 310.10  1.2718 0.26031  
# scale(Duration)              5.4558  5.4558     1 289.13  1.0777 0.30009  
# scale(Latitude)             13.9248 13.9248     1 107.49  2.7506 0.10014  

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = BacterialStructure)
summary(ModelLongitude)
# Number of obs: 343, groups:  StudyID, 81
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# scale(log(Rotation_cycles)) 15.496  15.496     1 334.31  3.0978 0.079313 . 
# scale(Species_Richness)      6.719   6.719     1 306.62  1.3432 0.247364   
# scale(Duration)             16.827  16.827     1 306.07  3.3639 0.067609 . 
# scale(Longitude)            42.750  42.750     1  83.77  8.5461 0.004452 **

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = BacterialStructure)
summary(ModelMAPmean)
# Number of obs: 343, groups:  StudyID, 81
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 17.5273 17.5273     1 335.93  3.4675 0.06346 .
# scale(Species_Richness)      8.0925  8.0925     1 309.02  1.6010 0.20672  
# scale(Duration)              5.5616  5.5616     1 276.91  1.1003 0.29512  
# scale(MAPmean)               9.3927  9.3927     1  88.77  1.8582 0.17628   

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = BacterialStructure)
summary(ModelMATmean)
# Number of obs: 343, groups:  StudyID, 81
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)  
# scale(log(Rotation_cycles)) 17.6488 17.6488     1 336.16  3.4983 0.0623 .
# scale(Species_Richness)      7.1117  7.1117     1 310.42  1.4097 0.2360  
# scale(Duration)              5.6365  5.6365     1 284.40  1.1172 0.2914  
# scale(MATmean)              12.1347 12.1347     1  81.46  2.4053 0.1248   



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## Species_Richness
sum(!is.na(BacterialStructure$Species_Richness)) ## n = 343
p1 <- ggplot(BacterialStructure, aes(y=RR, x=Species_Richness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="Species_Richness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Species_Richness" , y="lnBacterialStructure343")
p1
pdf("Species_Richness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(BacterialStructure$Duration)) ## n = 343
p2 <- ggplot(BacterialStructure, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnBacterialStructure343")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(BacterialStructure$Rotation_cycles)) ## n = 343
p3 <- ggplot(BacterialStructure, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnBacterialStructure343")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(BacterialStructure$Latitude)) ## n = 343
p5 <- ggplot(BacterialStructure, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnBacterialStructure343")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(BacterialStructure$Longitude)) ## n = 343
p6 <- ggplot(BacterialStructure, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnBacterialStructure343")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(BacterialStructure$MAPmean)) ## n = 343
p7 <- ggplot(BacterialStructure, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnBacterialStructure343")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(BacterialStructure$MATmean)) ## n = 343
p8 <- ggplot(BacterialStructure, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnBacterialStructure343")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(BacterialStructure$pHCK)) ## n = 142
p9 <- ggplot(BacterialStructure, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnBacterialStructure142")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(BacterialStructure$SOCCK)) ## n = 148
p10 <- ggplot(BacterialStructure, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnBacterialStructure148")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(BacterialStructure$TNCK)) ## n = 97
p11 <- ggplot(BacterialStructure, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnBacterialStructure97")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(BacterialStructure$NO3CK)) ## n = 50
p12 <- ggplot(BacterialStructure, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnBacterialStructure50")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(BacterialStructure$NH4CK)) ## n = 46
p13<- ggplot(BacterialStructure, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnBacterialStructure46")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(BacterialStructure$APCK)) ## n = 120
p14 <- ggplot(BacterialStructure, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnBacterialStructure120")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(BacterialStructure$AKCK)) ## n = 99
p15 <- ggplot(BacterialStructure, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnBacterialStructure99")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(BacterialStructure$ANCK)) ## n = 80
p16 <- ggplot(BacterialStructure, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnBacterialStructure80")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(BacterialStructure$RRpH)) ## n = 142
p17 <- ggplot(BacterialStructure, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="RRpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRpH" , y="RR142")
p17
pdf("RRpH.pdf",width=8,height=8)
p17
dev.off() 

## RRSOC
sum(!is.na(BacterialStructure$RRSOC)) ## n = 148
p18 <- ggplot(BacterialStructure, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="RRSOC")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRSOC" , y="RR148")
p18
pdf("RRSOC.pdf",width=8,height=8)
p18
dev.off() 

## RRTN
sum(!is.na(BacterialStructure$RRTN)) ## n = 97
p19 <- ggplot(BacterialStructure, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="RRTN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRTN" , y="RR97")
p19
pdf("RRTN.pdf",width=8,height=8)
p19
dev.off() 

## RRNO3
sum(!is.na(BacterialStructure$RRNO3)) ## n = 50
p20 <- ggplot(BacterialStructure, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="RRNO3")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNO3" , y="RR50")
p20
pdf("RRNO3.pdf",width=8,height=8)
p20
dev.off() 

## RRNH4
sum(!is.na(BacterialStructure$RRNH4)) ## n = 46
p21 <- ggplot(BacterialStructure, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="RRNH4")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNH4" , y="RR46")
p21
pdf("RRNH4.pdf",width=8,height=8)
p21
dev.off() 

## RRAP
sum(!is.na(BacterialStructure$RRAP)) ## n = 120
p22 <- ggplot(BacterialStructure, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="RRAP")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAP" , y="RR120")
p22
pdf("RRAP.pdf",width=8,height=8)
p22
dev.off() 

## RRAK
sum(!is.na(BacterialStructure$RRAK)) ## n = 99
p23 <- ggplot(BacterialStructure, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="RRAK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAK" , y="RR99")
p23
pdf("RRAK.pdf",width=8,height=8)
p23
dev.off() 

## RRAN
sum(!is.na(BacterialStructure$RRAN)) ## n = 80
p24 <- ggplot(BacterialStructure, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialStructure", x="RRAN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAN" , y="RR80")
p24
pdf("RRAN.pdf",width=8,height=8)
p24
dev.off() 

## RRYield
sum(!is.na(BacterialStructure$RRYield)) ## n = 35
p25 <- ggplot(BacterialStructure, aes(x=RR, y=RRYield)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm, color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") +
  theme_bw() +
  theme(text = element_text(family = "serif", size=20)) +
  geom_vline(aes(xintercept=0), colour="black", linewidth=0.5, linetype="dashed") +
  labs(x="RR", y="RRYield35") +
  theme(panel.grid=element_blank()) + 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = x ~ y,  
    parse = TRUE,
    color="black",
    size = 5, # 公式字体大小
    label.x = 0.05,  # 公式位置
    label.y = 0.85
  ) +
  stat_cor(method = "spearman", size = 5) +
  scale_x_continuous(limits=c(-0.2, 4.5), expand=c(0, 0)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) 
p25
pdf("RRYield_swapped_limited.pdf", width=8, height=8)
p25
dev.off()


################## Random Forest
library(rfPermute)
library(randomForest)
# Run the random forest and calculate the permutation test of variable importance
rf_model_perm3 <- rfPermute(RR ~ Latitude + Longitude + MAPmean + MATmean + Rotation_cycles + Species_Richness + Duration +
                              pHCK + SOCCK, 
                            data = BacterialStructure, 
                            na.action = na.roughfix, 
                            importance = TRUE, 
                            ntree = 500)
# Check the importance of variables and p-values
importance(rf_model_perm3)
#                    %IncMSE %IncMSE.pval IncNodePurity IncNodePurity.pval
# Latitude         24.913303   0.00990099     29.481131         0.00990099
# MAPmean          22.823704   0.00990099     19.456575         0.00990099
# Longitude        22.564796   0.00990099     27.734771         0.00990099
# pHCK             19.127078   0.00990099     16.676285         0.08910891
# MATmean          18.386225   0.00990099     17.701305         0.00990099
# Rotation_cycles  13.881860   0.01980198     11.866828         0.53465347
# SOCCK            13.855739   0.02970297     14.531093         0.31683168
# Duration         13.793967   0.01980198     17.830272         0.00990099
# Species_Richness  7.544045   0.10891089      1.318626         1.00000000


################################### Trials sorted by effect size
library(ggplot2)
library(dplyr)
BacterialStructure <- read.csv("BacterialStructure.csv", fileEncoding = "latin1")
# 计算95% CI + 显著性分类
df_plot <- BacterialStructure %>%
  filter(!is.na(RR), !is.na(Vi)) %>%
  mutate(
    SE = sqrt(Vi),
    CI_lower = RR - 1.96 * SE,
    CI_upper = RR + 1.96 * SE,
    EffectClass = case_when(
      CI_lower > 0 ~ "Positive",
      CI_upper < 0 ~ "Negative",
      TRUE ~ "Neutral"
    )
  ) %>%
  arrange(RR) %>%
  mutate(Index = row_number())

# 自定义颜色
effect_colors <- c("Negative" = "#F7AF34",
                   "Neutral"  = "#dedede", 
                   "Positive" = "#448DCD") 

# 绘图
ggplot(df_plot, aes(x = Index, y = RR, color = EffectClass)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, alpha = 0.8) +
  scale_color_manual(values = effect_colors) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Trials sorted by effect size",
       y = "Response Ratio (RR)",
       color = "Effect") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


# 统计各类数量
table(df_plot$EffectClass)
# Neutral Positive 
#      163      180 
## 8*8


################################################ piecewiseSEM
sem_data <- read.csv("BacterialStructure.csv", fileEncoding = "latin1")
# Load required packages
library(piecewiseSEM)
library(nlme)
library(dplyr)

sem_data <- sem_data[, c("StudyID", "RRpH", "RRSOC", "RR", "MATmean", "Duration", "Species_Richness")]
sem_data <- data.frame(lapply(sem_data, function(x) as.numeric(as.character(x))))

m0_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m0_2 <- lme(RRpH ~ MATmean + Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m0_3 <- lme(RRSOC ~ MATmean + Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model0 <- psem(m0_1, m0_2, m0_3)
summary(sem_model0) ## 
# AIC
# -139.387
# Fisher's C = 1.421 with P-value = 0.492

m1_1 <- lme(RR ~ RRpH + RRSOC + Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_2 <- lme(RRpH ~  Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_3 <- lme(RRSOC ~ Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model1 <- psem(m1_1, m1_2, m1_3)
summary(sem_model1)  ## 
# AIC
# -163.493
# Fisher's C = 2.106 with P-value = 0.349

m2_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_2 <- lme(RRpH ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_3 <- lme(RRSOC ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model2 <- psem(m2_1, m2_2, m2_3)
summary(sem_model2)  ## 
# AIC
# -169.061
# Fisher's C = 1.176 with P-value = 0.556

m3_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_2 <- lme(RRpH ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_3 <- lme(RRSOC ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model3 <- psem(m3_1, m3_2, m3_3)
summary(sem_model3)  ## 
# AIC
# -150.718
# Fisher's C = 1.044 with P-value = 0.593

m4_1 <- lme(RR ~ RRpH + RRSOC + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_2 <- lme(RRpH ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_3 <- lme(RRSOC ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model4 <- psem(m4_1, m4_2, m4_3)
summary(sem_model4) ## 
# AIC
# -193.021
# Fisher's C = 1.825 with P-value = 0.401

# Structural Equation Model of sem_model4 
# 
# Call:
#   RR ~ RRpH + RRSOC + Species_Richness
# RRpH ~ Species_Richness
# RRSOC ~ Species_Richness
# 
# AIC
# -193.021
# 
# ---
#   Tests of directed separation:
#   
#   Independ.Claim Test.Type DF Crit.Value P.Value 
# RRSOC ~ RRpH + ...      coef 61     0.8449  0.4015 
# 
# --
#   Global goodness-of-fit:
#   
#   Chi-Squared = NA with P-value = NA and on 1 degrees of freedom
# Fisher's C = 1.825 with P-value = 0.401 and on 2 degrees of freedom
# 
# ---
# Coefficients:
# 
#   Response        Predictor Estimate Std.Error DF Crit.Value P.Value Std.Estimate  
#         RR             RRpH  -0.5068    1.3051 60    -0.3884  0.6991      -0.0518  
#         RR            RRSOC   0.1790    0.4448 60     0.4025  0.6888       0.0541  
#         RR Species_Richness  -0.0935    0.2552 60    -0.3663  0.7155      -0.0645  
#       RRpH Species_Richness  -0.0346    0.0134 74    -2.5768  0.0120      -0.2334 *
#      RRSOC Species_Richness  -0.0638    0.0502 74    -1.2706  0.2079      -0.1458  
# 
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05
# 
# ---
# Individual R-squared:
# 
#   Response method Marginal Conditional
#         RR   none     0.00        0.72
#       RRpH   none     0.04        0.86
#      RRSOC   none     0.02        0.65
     
m5_1 <- lme(RR ~ RRpH + RRSOC + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_2 <- lme(RRpH ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_3 <- lme(RRSOC ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model5 <- psem(m5_1, m5_2, m5_3)
summary(sem_model5) ## 
# AIC
# -172.992
# Fisher's C = 1.727 with P-value = 0.422

m6_1 <- lme(RR ~ RRpH + RRSOC +MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_2 <- lme(RRpH ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_3 <- lme(RRSOC ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model6 <- psem(m6_1, m6_2, m6_3)
summary(sem_model6)  ##
# AIC
# -180.611
# Fisher's C = 0.876 with P-value = 0.645
