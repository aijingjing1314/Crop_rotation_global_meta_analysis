
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
BacterialRichness <- read.csv("BacterialRichness.csv", fileEncoding = "latin1")
# Check data
head(BacterialRichness)

# 1. The number of Obversation
total_number <- nrow(BacterialRichness)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 372

# 2. The number of Study
unique_studyid_number <- length(unique(BacterialRichness$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 94


#### 3. Overall effect size
total_effect_model <- rma.mv(yi = RR, 
                             V = Vi, 
                             random = ~ 1 | StudyID,  # StudyID is radom factor
                             data = BacterialRichness, 
                             method = "REML")
# The results of Overall effect size
summary(total_effect_model)
# Multivariate Meta-Analysis Model (k = 372; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -950.6141  1901.2283  1905.2283  1913.0607  1905.2609   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0208  0.1443     94     no  StudyID 
# Test for Heterogeneity:
# Q(df = 371) = 10133.4023, p-val < .0001
# Model Results:
# estimate      se    zval    pval   ci.lb   ci.ub    
#   0.0385  0.0156  2.4773  0.0132  0.0080  0.0690  * 


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
boot_results1 <- boot(data = BacterialRichness, statistic = boot_fun, R = 1000, parallel = "snow", ncpus = numCores, cl = cl)
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
# Estimate for coefficient 1 : 0.03854843 
# 95% BCa CI for coefficient 1 : 0.02267992 0.05651243 


#### 5. Funnel Plot
simple_model <- rma(yi = RR, 
                    vi = Vi, 
                    data = BacterialRichness, 
                    method = "REML")
#### 
funnel(simple_model)
# Output  6 * 6
#### Egger's test
regtest(simple_model)
# Regression Test for Funnel Plot Asymmetry
# Model:     mixed-effects meta-regression model
# Predictor: standard error
# Test for Funnel Plot Asymmetry: z = 0.7864, p = 0.4316
# Limit Estimate (as sei -> 0):   b = 0.0232 (CI: -0.0146, 0.0609)

#  Rosenthal’s Fail-Safe N
# This method estimates how many missing studies with null effect 
# would be needed to make the overall effect non-significant
fsn_rosenthal <- fsn(x = simple_model, type = "Rosenthal")
# Print the FSN result
print(fsn_rosenthal)
# Fail-safe N Calculation Using the General Approach
# Average Effect Size:         0.0362 (with file drawer: 0.0100)
# Amount of Heterogeneity:     0.0306 (with file drawer: 0.0309)
# Observed Significance Level: 0.0002 (with file drawer: 0.0501)
# Target Significance Level:   0.05
# Fail-safe N: 860



#### 8. Subgroup analysis
### 8.1 LegumeNonlegume
BacterialRichness_filteredLegumeNonlegume <- subset(BacterialRichness, LegumeNonlegume %in% c("Legume to Non-legume", "Non-legume to Legume", "Non-legume to Non-legume"))
#
BacterialRichness_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(BacterialRichness_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredLegumeNonlegume %>%
  group_by(LegumeNonlegume) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   LegumeNonlegume          Observations Unique_StudyID
#   <fct>                           <int>          <int>
# 1 Legume to Non-legume              109             26
# 2 Non-legume to Legume              113             37
# 3 Non-legume to Non-legume          146             50

overall_model_BacterialRichness_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = BacterialRichness_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 368; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -959.3756  1918.7513  1926.7513  1942.3509  1926.8624   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0206  0.1434     93     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 365) = 7489.1024, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 5.8647, p-val = 0.1184
# Model Results:
#                                          estimate      se    zval    pval   ci.lb   ci.ub    
# LegumeNonlegumeLegume to Non-legume        0.0338  0.0163  2.0801  0.0375  0.0020  0.0657  * 
# LegumeNonlegumeNon-legume to Legume        0.0351  0.0159  2.2036  0.0276  0.0039  0.0663  * 
# LegumeNonlegumeNon-legume to Non-legume    0.0387  0.0160  2.4193  0.0156  0.0073  0.0700  *  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredLegumeNonlegume)
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
    #                                 "a"                                     "a"                                     "a" 

### 8.2 AMnonAM
BacterialRichness_filteredAMnonAM <- subset(BacterialRichness, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
BacterialRichness_filteredAMnonAM$AMnonAM <- droplevels(factor(BacterialRichness_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredAMnonAM %>%
  group_by(AMnonAM) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   AMnonAM     Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 AM to AM             299             80
# 2 AM to nonAM           43             14
# 3 nonAM to AM           27              8

overall_model_BacterialRichness_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = BacterialRichness_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 369; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -946.8848  1893.7695  1901.7695  1917.3801  1901.8803   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0207  0.1438     93     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 366) = 9581.9612, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 16.7796, p-val = 0.0008
# Model Results:
#                     estimate      se    zval    pval    ci.lb   ci.ub     
# AMnonAMAM to AM       0.0260  0.0164  1.5879  0.1123  -0.0061  0.0580     
# AMnonAMAM to nonAM    0.0541  0.0183  2.9579  0.0031   0.0183  0.0900  ** 
# AMnonAMnonAM to AM    0.1125  0.0529  2.1288  0.0333   0.0089  0.2161   * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredAMnonAM)
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
   #             "a"                "b"               "ab" 


### 8.3 C3C4
BacterialRichness_filteredC3C4 <- subset(BacterialRichness, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
BacterialRichness_filteredC3C4$C3C4 <- droplevels(factor(BacterialRichness_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredC3C4 %>%
  group_by(C3C4) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   C3C4     Observations Unique_StudyID
#   <fct>           <int>          <int>
# 1 C3 to C3          177             56
# 2 C3 to C4          129             38
# 3 C4 to C3           57             19

overall_model_BacterialRichness_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = BacterialRichness_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 363; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -908.3147  1816.6294  1824.6294  1840.1738  1824.7421   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0199  0.1411     92     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 360) = 6410.6722, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 62.2594, p-val < .0001
# Model Results:
#               estimate      se     zval    pval    ci.lb   ci.ub      
# C3C4C3 to C3    0.0533  0.0158   3.3792  0.0007   0.0224  0.0842  *** 
# C3C4C3 to C4    0.0431  0.0159   2.7080  0.0068   0.0119  0.0743   ** 
# C3C4C4 to C3   -0.0319  0.0179  -1.7826  0.0747  -0.0669  0.0032    . 
 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredC3C4)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredC3C4)
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
#          "a"          "a"          "b" 


### 8.4 Annual_Pere
BacterialRichness_filteredAnnual_Pere <- subset(BacterialRichness, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
BacterialRichness_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(BacterialRichness_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredAnnual_Pere %>%
  group_by(Annual_Pere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Annual_Pere            Observations Unique_StudyID
#   <fct>                         <int>          <int>
# 1 Annual to Annual                269             70
# 2 Annual to Perennial              29             14
# 3 Perennial to Annual              58             15
# 4 Perennial to Perennial           15              7

overall_model_BacterialRichness_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = BacterialRichness_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 371; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -883.2096  1766.4193  1776.4193  1795.9461  1776.5855   
# Variance Components:
#            estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0206  0.1435     93     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 367) = 9182.9155, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 150.1333, p-val < .0001
# Model Results:
#                                    estimate      se     zval    pval    ci.lb   ci.ub     
# Annual_PereAnnual to Annual          0.0433  0.0165   2.6229  0.0087   0.0109  0.0756  ** 
# Annual_PereAnnual to Perennial      -0.0253  0.0179  -1.4117  0.1581  -0.0605  0.0098     
# Annual_PerePerennial to Annual       0.0824  0.0266   3.0937  0.0020   0.0302  0.1345  ** 
# Annual_PerePerennial to Perennial   -0.0507  0.0276  -1.8375  0.0661  -0.1049  0.0034   . 

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredAnnual_Pere)
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
      #                        "a"                               "b"                               "a"                               "b" 



### 8.6 PlantStageSubgroup
BacterialRichness_filteredPlantStageSubgroup <- subset(BacterialRichness, PlantStageSubgroup %in% c("Vegetative stage","Reproductive stage", "Maturity stage","Harvest"))
#
BacterialRichness_filteredPlantStageSubgroup$PlantStageSubgroup <- droplevels(factor(BacterialRichness_filteredPlantStageSubgroup$PlantStageSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredPlantStageSubgroup %>%
  group_by(PlantStageSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   PlantStageSubgroup Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Harvest                     193             51
# 2 Maturity stage               34             15
# 3 Reproductive stage           51             15
# 4 Vegetative stage             32              7
overall_model_BacterialRichness_filteredPlantStageSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + PlantStageSubgroup, random = ~ 1 | StudyID, data = BacterialRichness_filteredPlantStageSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredPlantStageSubgroup)
# Multivariate Meta-Analysis Model (k = 310; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -606.6254  1213.2507  1223.2507  1241.8686  1223.4507   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0136  0.1168     79     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 306) = 6389.7934, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 29.1176, p-val < .0001
# Model Results:
#                                       estimate      se     zval    pval    ci.lb   ci.ub      
# PlantStageSubgroupHarvest               0.0483  0.0145   3.3414  0.0008   0.0200  0.0767  *** 
# PlantStageSubgroupMaturity stage       -0.0230  0.0208  -1.1091  0.2674  -0.0638  0.0177      
# PlantStageSubgroupReproductive stage    0.0645  0.0156   4.1279  <.0001   0.0339  0.0951  *** 
# PlantStageSubgroupVegetative stage      0.0434  0.0198   2.1896  0.0286   0.0046  0.0822    * 
 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredPlantStageSubgroup)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredPlantStageSubgroup)
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
           #                       "a"                                  "b"                                  "c"                                 "ac"  


### 8.7 Bulk_Rhizosphere
BacterialRichness_filteredBulk_Rhizosphere <- subset(BacterialRichness, Bulk_Rhizosphere %in% c("Non-Rhizosphere", "Rhizosphere"))
#
BacterialRichness_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(BacterialRichness_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredBulk_Rhizosphere %>%
  group_by(Bulk_Rhizosphere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Bulk_Rhizosphere Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 Non-Rhizosphere           241             63
# 2 Rhizosphere               131             43

overall_model_BacterialRichness_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = BacterialRichness_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 372; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -929.5064  1859.0128  1865.0128  1876.7533  1865.0783   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0211  0.1453     94     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 370) = 9978.2818, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 51.6568, p-val < .0001
# Model Results:
#                                  estimate      se    zval    pval    ci.lb   ci.ub      
# Bulk_RhizosphereNon-Rhizosphere    0.0184  0.0160  1.1553  0.2480  -0.0128  0.0497      
# Bulk_RhizosphereRhizosphere        0.0726  0.0165  4.4102  <.0001   0.0403  0.1048  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredBulk_Rhizosphere)
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
#                             "a"                             "b" 


### 8.8 Soil_texture
BacterialRichness_filteredSoil_texture <- subset(BacterialRichness, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
BacterialRichness_filteredSoil_texture$Soil_texture <- droplevels(factor(BacterialRichness_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredSoil_texture %>%
  group_by(Soil_texture) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Soil_texture Observations Unique_StudyID
#   <fct>               <int>          <int>
# 1 Coarse                 88             17
# 2 Fine                  101             18
# 3 Medium                124             34

overall_model_BacterialRichness_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = BacterialRichness_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 313; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -961.4570  1922.9140  1930.9140  1945.8603  1931.0452   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0257  0.1603     69     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 310) = 7414.2958, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 3.8625, p-val = 0.2767
# Model Results:
#                     estimate      se    zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse    0.0498  0.0396  1.2561  0.2091  -0.0279  0.1275    
# Soil_textureFine      0.0126  0.0393  0.3196  0.7493  -0.0645  0.0897    
# Soil_textureMedium    0.0419  0.0284  1.4773  0.1396  -0.0137  0.0976    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredSoil_texture)
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
#                "a"                "a"                "a" 


### 8.9 Tillage
BacterialRichness_filteredTillage <- subset(BacterialRichness, Tillage %in% c("Tillage", "No_tillage"))
#
BacterialRichness_filteredTillage$Tillage <- droplevels(factor(BacterialRichness_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredTillage %>%
  group_by(Tillage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Tillage    Observations Unique_StudyID
#   <fct>             <int>          <int>
# 1 No_tillage           38              6
# 2 Tillage              49             14

overall_model_BacterialRichness_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = BacterialRichness_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredTillage)
# Multivariate Meta-Analysis Model (k = 87; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -63.9398  127.8796  133.8796  141.2075  134.1759   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0394  0.1984     18     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 85) = 993.3860, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.0573, p-val = 0.9718
# Model Results:
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.0047  0.0500  -0.0943  0.9248  -0.1027  0.0933    
# TillageTillage      -0.0003  0.0483  -0.0053  0.9958  -0.0950  0.0945     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredTillage)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredTillage)
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
BacterialRichness_filteredStraw_retention <- subset(BacterialRichness, Straw_retention %in% c("Retention", "No_retention"))
#
BacterialRichness_filteredStraw_retention$Straw_retention <- droplevels(factor(BacterialRichness_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredStraw_retention %>%
  group_by(Straw_retention) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Straw_retention Observations Unique_StudyID
#   <fct>                  <int>          <int>
# 1 No_retention              19              8
# 2 Retention                 51             13

overall_model_BacterialRichness_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = BacterialRichness_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 70; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -42.0056   84.0111   90.0111   96.6696   90.3861   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0364  0.1907     20     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 68) = 1127.7016, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 2.5749, p-val = 0.2760
# Model Results:
#                             estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention   -0.0567  0.0511  -1.1100  0.2670  -0.1569  0.0434    
# Straw_retentionRetention       0.0097  0.0464   0.2100  0.8336  -0.0812  0.1007    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredStraw_retention)
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
#                         "a"                         "a" 


### 8.12 RotationcyclesSubgroup
BacterialRichness_filteredRotationcyclesSubgroup <- subset(BacterialRichness, RotationcyclesSubgroup %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
BacterialRichness_filteredRotationcyclesSubgroup$RotationcyclesSubgroup <- droplevels(factor(BacterialRichness_filteredRotationcyclesSubgroup$RotationcyclesSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredRotationcyclesSubgroup %>%
  group_by(RotationcyclesSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   RotationcyclesSubgroup Observations Unique_StudyID
#   <fct>                         <int>          <int>
# 1 D1                              135             45
# 2 D1-3                             83             27
# 3 D10                              38             11
# 4 D3-5                             50             11
# 5 D5-10                            66             17

overall_model_BacterialRichness_filteredRotationcyclesSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + RotationcyclesSubgroup, random = ~ 1 | StudyID, data = BacterialRichness_filteredRotationcyclesSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredRotationcyclesSubgroup)
# Multivariate Meta-Analysis Model (k = 372; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -925.5043  1851.0085  1863.0085  1886.4407  1863.2418   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0239  0.1544     94     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 367) = 9302.1302, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 67.4559, p-val < .0001
# Model Results:
#                              estimate      se     zval    pval    ci.lb   ci.ub      
# RotationcyclesSubgroupD1      -0.0012  0.0180  -0.0688  0.9452  -0.0365  0.0340      
# RotationcyclesSubgroupD1-3     0.0493  0.0182   2.6997  0.0069   0.0135  0.0850   ** 
# RotationcyclesSubgroupD10      0.1013  0.0273   3.7157  0.0002   0.0478  0.1547  *** 
# RotationcyclesSubgroupD3-5     0.0353  0.0230   1.5358  0.1246  -0.0098  0.0804      
# RotationcyclesSubgroupD5-10    0.1036  0.0215   4.8210  <.0001   0.0615  0.1457  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredRotationcyclesSubgroup)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredRotationcyclesSubgroup)
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
   #                      "a"                        "bc"                        "bd"                        "ac"                         "d" 


### 8.13 DurationSubgroup
BacterialRichness_filteredDurationSubgroup <- subset(BacterialRichness, DurationSubgroup %in% c("D1", "D2", "D3", "D4", "D5", "D6-10", "D11-20", "D20-30", "D30"))
#
BacterialRichness_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(BacterialRichness_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredDurationSubgroup %>%
  group_by(DurationSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   DurationSubgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D1                         24             11
# 2 D11-20                     35             13
# 3 D2                         52             23
# 4 D20-30                     35              8
# 5 D3                         62             17
# 6 D30                        67              5
# 7 D4                         31             11
# 8 D5                         11              3
# 9 D6-10                      55             17
overall_model_BacterialRichness_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = BacterialRichness_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 372; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -917.4135  1834.8270  1854.8270  1893.7710  1855.4520   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0287  0.1696     94     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 363) = 7601.1238, p-val < .0001
# Test of Moderators (coefficients 1:9):
# QM(df = 9) = 87.7112, p-val < .0001
# Model Results:
#                         estimate      se     zval    pval    ci.lb   ci.ub      
# DurationSubgroupD1       -0.0311  0.0313  -0.9948  0.3198  -0.0924  0.0302      
# DurationSubgroupD11-20    0.1180  0.0271   4.3481  <.0001   0.0648  0.1711  *** 
# DurationSubgroupD2       -0.0220  0.0245  -0.8988  0.3687  -0.0699  0.0259      
# DurationSubgroupD20-30    0.2548  0.0379   6.7223  <.0001   0.1805  0.3290  *** 
# DurationSubgroupD3        0.0266  0.0244   1.0909  0.2753  -0.0212  0.0743      
# DurationSubgroupD30      -0.1079  0.0768  -1.4038  0.1604  -0.2585  0.0427      
# DurationSubgroupD4        0.0052  0.0525   0.0987  0.9214  -0.0977  0.1081      
# DurationSubgroupD5        0.1459  0.0666   2.1907  0.0285   0.0154  0.2765    * 
# DurationSubgroupD6-10     0.0697  0.0383   1.8211  0.0686  -0.0053  0.1448    . 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredDurationSubgroup)
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
    #                "a"                    "b"                    "a"                    "c"                   "de"                   "ad" 
    # DurationSubgroupD4     DurationSubgroupD5  DurationSubgroupD6-10 
    #             "abde"                  "bce"                   "be" 


### 8.14 SpeciesRichnessSubgroup
BacterialRichness_filteredSpeciesRichnessSubgroup <- subset(BacterialRichness, SpeciesRichnessSubgroup %in% c("R2", "R3"))
#
BacterialRichness_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup <- droplevels(factor(BacterialRichness_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredSpeciesRichnessSubgroup %>%
  group_by(SpeciesRichnessSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   SpeciesRichnessSubgroup Observations Unique_StudyID
#   <fct>                          <int>          <int>
# 1 R2                               285             80
# 2 R3                                82             21
# 3 R4                                 5              4
overall_model_BacterialRichness_filteredSpeciesRichnessSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichnessSubgroup, random = ~ 1 | StudyID, data = BacterialRichness_filteredSpeciesRichnessSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredSpeciesRichnessSubgroup)
# Multivariate Meta-Analysis Model (k = 367; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -949.3321  1898.6642  1904.6642  1916.3638  1904.7306   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0210  0.1451     92     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 365) = 10057.1022, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 19.7715, p-val < .0001
# Model Results:
#                            estimate      se    zval    pval   ci.lb   ci.ub      
# SpeciesRichnessSubgroupR2    0.0331  0.0159  2.0844  0.0371  0.0020  0.0643    * 
# SpeciesRichnessSubgroupR3    0.0714  0.0180  3.9737  <.0001  0.0362  0.1065  ***   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredSpeciesRichnessSubgroup)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredSpeciesRichnessSubgroup)
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
# SpeciesRichnessSubgroupR2 SpeciesRichnessSubgroupR3
#                       "a"                       "b"                      


### 8.15 Primer
BacterialRichness_filteredPrimer <- subset(BacterialRichness, Primer %in% c("V1-V3", "V1-V4", "V1-V5", "V3-V4", "V3-V5", "V4", "V4-V5", "V5-V7", "V6-V8", "V9", "Full length"))
#
BacterialRichness_filteredPrimer$Primer <- droplevels(factor(BacterialRichness_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredPrimer %>%
  group_by(Primer) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#    Primer      Observations Unique_StudyID
#    <fct>              <int>          <int>
#  1 Full length            3              1
#  2 V1-V3                  7              3
#  3 V1-V4                  4              1
#  4 V1-V5                  6              1
#  5 V3-V4                166             47
#  6 V3-V5                  2              1
#  7 V4                   111             21
#  8 V4-V5                 48             14
#  9 V5-V7                  3              1
# 10 V6-V8                 16              1

overall_model_BacterialRichness_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = BacterialRichness_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 366; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -886.9844  1773.9689  1795.9689  1838.5931  1796.7363   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0191  0.1381     91     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 356) = 8571.6002, p-val < .0001
# Test of Moderators (coefficients 1:10):
# QM(df = 10) = 21.3952, p-val = 0.0185
# Model Results:
#                    estimate      se     zval    pval    ci.lb   ci.ub      
# PrimerFull length    0.0019  0.1410   0.0137  0.9891  -0.2745  0.2784      
# PrimerV1-V3          0.1206  0.0842   1.4320  0.1521  -0.0445  0.2856      
# PrimerV1-V4          0.0545  0.1384   0.3939  0.6937  -0.2167  0.3258      
# PrimerV1-V5          0.4927  0.1410   3.4933  0.0005   0.2163  0.7692  *** 
# PrimerV3-V4          0.0077  0.0210   0.3695  0.7118  -0.0333  0.0488      
# PrimerV3-V5          0.0071  0.1392   0.0511  0.9593  -0.2658  0.2800      
# PrimerV4             0.0289  0.0318   0.9074  0.3642  -0.0335  0.0912      
# PrimerV4-V5          0.0905  0.0401   2.2538  0.0242   0.0118  0.1692    * 
# PrimerV5-V7         -0.1044  0.1440  -0.7251  0.4684  -0.3867  0.1779      
# PrimerV6-V8          0.0900  0.1392   0.6467  0.5178  -0.1828  0.3628      
 
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredPrimer)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredPrimer)
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
# PrimerFull length       PrimerV1-V3       PrimerV1-V4       PrimerV1-V5       PrimerV3-V4       PrimerV3-V5          PrimerV4       PrimerV4-V5 
#               "a"               "a"               "a"               "b"               "a"               "a"               "a"               "a" 
#       PrimerV5-V7       PrimerV6-V8 
#               "a"               "a" 



### 8.16 Latitude_Subgroup
BacterialRichness_filteredLatitude_Subgroup <- subset(BacterialRichness, Latitude_Subgroup %in% c("La20", "La20-40", "La40"))
#
BacterialRichness_filteredLatitude_Subgroup$Latitude_Subgroup <- droplevels(factor(BacterialRichness_filteredLatitude_Subgroup$Latitude_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredLatitude_Subgroup %>%
  group_by(Latitude_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Latitude_Subgroup Observations Unique_StudyID
#   <fct>                    <int>          <int>
# 1 La20                        37              8
# 2 La20-40                    125             48
# 3 La40                       210             38

overall_model_BacterialRichness_filteredLatitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Latitude_Subgroup, random = ~ 1 | StudyID, data = BacterialRichness_filteredLatitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredLatitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 372; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -950.1609  1900.3218  1908.3218  1923.9649  1908.4316   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0209  0.1445     94     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 369) = 9171.3970, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 8.0946, p-val = 0.0441
# Model Results:
#                           estimate      se    zval    pval    ci.lb   ci.ub    
# Latitude_SubgroupLa20       0.1097  0.0533  2.0582  0.0396   0.0052  0.2142  * 
# Latitude_SubgroupLa20-40    0.0296  0.0219  1.3542  0.1757  -0.0133  0.0725    
# Latitude_SubgroupLa40       0.0348  0.0244  1.4229  0.1548  -0.0131  0.0826         

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredLatitude_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredLatitude_Subgroup)
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
#                   "a"                      "a"                      "a" 


### 8.17 Longitude_Subgroup
BacterialRichness_filteredLongitude_Subgroup <- subset(BacterialRichness, Longitude_Subgroup %in% c("Lo-180-0", "Lo-0-180"))
#
BacterialRichness_filteredLongitude_Subgroup$Longitude_Subgroup <- droplevels(factor(BacterialRichness_filteredLongitude_Subgroup$Longitude_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredLongitude_Subgroup %>%
  group_by(Longitude_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Longitude_Subgroup Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Lo-0-180                    254             80
# 2 Lo-180-0                    118             14

overall_model_BacterialRichness_filteredLongitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Subgroup, random = ~ 1 | StudyID, data = BacterialRichness_filteredLongitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredLongitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 372; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -950.5015  1901.0030  1907.0030  1918.7435  1907.0686   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0210  0.1450     94     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 370) = 10082.9963, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 6.3124, p-val = 0.0426
# Model Results:
#                             estimate      se    zval    pval    ci.lb   ci.ub    
# Longitude_SubgroupLo-0-180    0.0416  0.0169  2.4630  0.0138   0.0085  0.0748  * 
# Longitude_SubgroupLo-180-0    0.0204  0.0411  0.4960  0.6199  -0.0602  0.1010          

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredLongitude_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredLongitude_Subgroup)
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
#                        "a"                        "a" 



### 8.18 MAPmean_Subgroup
BacterialRichness_filteredMAPmean_Subgroup <- subset(BacterialRichness, MAPmean_Subgroup %in% c("MAP600", "MAP600-1200", "MAP1200"))
#
BacterialRichness_filteredMAPmean_Subgroup$MAPmean_Subgroup <- droplevels(factor(BacterialRichness_filteredMAPmean_Subgroup$MAPmean_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredMAPmean_Subgroup %>%
  group_by(MAPmean_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MAPmean_Subgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 MAP1200                    83             25
# 2 MAP600                    160             41
# 3 MAP600-1200               129             29

overall_model_BacterialRichness_filteredMAPmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean_Subgroup, random = ~ 1 | StudyID, data = BacterialRichness_filteredMAPmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredMAPmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 372; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -943.7043  1887.4086  1895.4086  1911.0518  1895.5185   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0241  0.1554     94     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 369) = 9427.9308, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 21.5920, p-val < .0001
# Model Results:
#                              estimate      se     zval    pval    ci.lb   ci.ub     
# MAPmean_SubgroupMAP1200        0.0734  0.0322   2.2776  0.0228   0.0102  0.1366   * 
# MAPmean_SubgroupMAP600         0.0704  0.0227   3.1050  0.0019   0.0260  0.1148  ** 
# MAPmean_SubgroupMAP600-1200   -0.0388  0.0259  -1.5012  0.1333  -0.0895  0.0119     


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredMAPmean_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredMAPmean_Subgroup)
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
#                     "a"                         "a"                         "b" 



### 8.19 MATmean_Subgroup
BacterialRichness_filteredMATmean_Subgroup <- subset(BacterialRichness, MATmean_Subgroup %in% c("MAT8", "MAT8-15", "MAT15"))
#
BacterialRichness_filteredMATmean_Subgroup$MATmean_Subgroup <- droplevels(factor(BacterialRichness_filteredMATmean_Subgroup$MATmean_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredMATmean_Subgroup %>%
  group_by(MATmean_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MATmean_Subgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 MAT15                      97             33
# 2 MAT8                      211             41
# 3 MAT8-15                    64             21

overall_model_BacterialRichness_filteredMATmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Subgroup, random = ~ 1 | StudyID, data = BacterialRichness_filteredMATmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredMATmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 372; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -936.0433  1872.0865  1880.0865  1895.7297  1880.1964   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0217  0.1474     94     no  StudyID
# Test for Residual Heterogeneity:
# QE(df = 369) = 9612.3912, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 37.1541, p-val < .0001
# Model Results:
#                          estimate      se     zval    pval    ci.lb    ci.ub      
# MATmean_SubgroupMAT15      0.0706  0.0270   2.6133  0.0090   0.0177   0.1236   ** 
# MATmean_SubgroupMAT8       0.0755  0.0220   3.4337  0.0006   0.0324   0.1186  *** 
# MATmean_SubgroupMAT8-15   -0.0819  0.0275  -2.9824  0.0029  -0.1357  -0.0281   ** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredMATmean_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredMATmean_Subgroup)
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
#                   "a"                     "a"                     "b" 



### 8.20 DurationSubgroup
BacterialRichness_filteredDurationSubgroup <- subset(BacterialRichness, DurationSubgroup %in% c("D5", "D5-10", "D10-20", "D20-30", "D30"))
#
BacterialRichness_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(BacterialRichness_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialRichness_filteredDurationSubgroup %>%
  group_by(DurationSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   DurationSubgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D10-20                     35             13
# 2 D20-30                     35              8
# 3 D30                        67              5
# 4 D5                        180             56
# 5 D5-10                      55             17

overall_model_BacterialRichness_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = BacterialRichness_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialRichness_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 372; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -930.5387  1861.0773  1873.0773  1896.5095  1873.3107   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0252  0.1588     94     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 367) = 9895.0915, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 51.2680, p-val < .0001
# Model Results:
#                         estimate      se     zval    pval    ci.lb   ci.ub      
# DurationSubgroupD10-20    0.1044  0.0245   4.2564  <.0001   0.0563  0.1525  *** 
# DurationSubgroupD20-30    0.2359  0.0364   6.4877  <.0001   0.1647  0.3072  *** 
# DurationSubgroupD30      -0.1076  0.0721  -1.4918  0.1358  -0.2489  0.0338      
# DurationSubgroupD5        0.0142  0.0195   0.7268  0.4673  -0.0241  0.0525      
# DurationSubgroupD5-10     0.0360  0.0335   1.0742  0.2827  -0.0297  0.1018     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialRichness_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_BacterialRichness_filteredDurationSubgroup)
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
#                 "a"                    "b"                    "c"                    "c"                   "ac" 




#### 9. Linear Mixed Effect Model
# 
BacterialRichness$Wr <- 1 / BacterialRichness$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialRichness)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialRichness)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialRichness)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialRichness)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialRichness)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialRichness)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialRichness)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialRichness)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(Species_Richness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = BacterialRichness)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = BacterialRichness)

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
#           Model       AIC        BIC   logLik  Marginal_R2 Conditional_R2
# Model1   Model1 -142.8175 -119.30416 77.40876 1.940666e-05   0.0008705730
# Model2   Model2 -142.5956 -119.08223 77.29780 7.435816e-06   0.0008847169
# Model3   Model3 -143.3845 -119.87117 77.69227 2.794748e-05   0.0009063557
# Model4   Model4 -142.7724 -119.25900 77.38618 1.915588e-05   0.0008712059
# Model5   Model5 -141.6076 -118.09426 76.80381 6.378652e-07   0.0008595603
# Model6   Model6 -141.5900 -118.07659 76.79498 7.067301e-07   0.0008610396
# Model7   Model7 -142.5980 -119.08460 77.29898 7.104220e-06   0.0008824813
# Model8   Model8 -143.3410 -119.82762 77.67049 2.785164e-05   0.0009075618
# Model9   Model9 -112.3750  -73.18602 66.18748 8.677546e-05   0.0010326213
# Model10 Model10 -117.1597  -77.97080 68.57987 1.068835e-04   0.0011433207

##### Model 3 is the best model
summary(Model3)
# Number of obs: 372, groups:  StudyID, 94
anova(Model3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 12.622  12.622     1 157.16  1.1725 0.2805
# scale(Species_Richness)      2.195   2.195     1 266.12  0.2039 0.6520
# scale(Duration)             13.483  13.483     1 119.55  1.2525 0.2653


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = BacterialRichness)
summary(ModelpH)
# Number of obs: 192, groups:  StudyID, 61
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)  
# scale(log(Rotation_cycles)) 11.1726 11.1726     1  94.608  2.2793 0.1344  
# scale(Species_Richness)      0.0016  0.0016     1 180.528  0.0003 0.9857  
# scale(Duration)              0.1194  0.1194     1  76.558  0.0244 0.8764  
# scale(pHCK)                 16.8217 16.8217     1  63.511  3.4317 0.0686 .

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = BacterialRichness)
summary(ModelSOC)
# Number of obs: 204, groups:  StudyID, 61
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles))  0.1258  0.1258     1 106.328  0.0215 0.8837
# scale(Species_Richness)     15.7239 15.7239     1 188.366  2.6869 0.1028
# scale(Duration)              4.2501  4.2501     1  74.836  0.7262 0.3968
# scale(SOCCK)                 0.0056  0.0056     1  46.204  0.0010 0.9755 

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = BacterialRichness)
summary(ModelTN)
# Number of obs: 123, groups:  StudyID, 44
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  3.8678  3.8678     1 66.735  0.5607 0.45660  
# scale(Species_Richness)     21.2163 21.2163     1 94.457  3.0758 0.08271 .
# scale(Duration)             17.7768 17.7768     1 56.476  2.5772 0.11399  
# scale(TNCK)                  1.1136  1.1136     1 36.322  0.1614 0.69019  


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = BacterialRichness)
summary(ModelNO3)
# Number of obs: 67, groups:  StudyID, 20
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 3.8167  3.8167     1 61.973  0.6588 0.4201
# scale(Species_Richness)     7.6536  7.6536     1 18.457  1.3211 0.2651
# scale(Duration)             2.6376  2.6376     1 61.908  0.4553 0.5024
# scale(NO3CK)                0.6555  0.6555     1 26.288  0.1131 0.7393

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = BacterialRichness)
summary(ModelNH4)
# Number of obs: 47, groups:  StudyID, 17
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 2.6687  2.6687     1 33.156  0.3389 0.5644
# scale(Species_Richness)     5.6521  5.6521     1 19.596  0.7177 0.4071
# scale(Duration)             2.7143  2.7143     1 41.985  0.3447 0.5603
# scale(NH4CK)                1.6412  1.6412     1 15.124  0.2084 0.6545

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = BacterialRichness)
summary(ModelAP)
# Number of obs: 139, groups:  StudyID, 50
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)  
# scale(log(Rotation_cycles))  0.921   0.921     1  86.382  0.1284 0.7210  
# scale(Species_Richness)      6.441   6.441     1 133.054  0.8977 0.3451  
# scale(Duration)             43.370  43.370     1  84.178  6.0448 0.0160 *
# scale(APCK)                  1.122   1.122     1 119.248  0.1564 0.6932  

#### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = BacterialRichness)
summary(ModelAK)
# Number of obs: 108, groups:  StudyID, 40
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  8.6905  8.6905     1  68.460  1.2938 0.25931  
# scale(Species_Richness)     19.5987 19.5987     1  96.480  2.9177 0.09083 .
# scale(Duration)              0.0086  0.0086     1  45.418  0.0013 0.97153  
# scale(AKCK)                  2.7345  2.7345     1 101.890  0.4071 0.52488 

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = BacterialRichness)
summary(ModelAN)
# Number of obs: 90, groups:  StudyID, 25
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles))  1.5367  1.5367     1 32.041  0.2013 0.6567
# scale(Species_Richness)     10.1710 10.1710     1 83.019  1.3324 0.2517
# scale(Duration)              2.4549  2.4549     1 25.514  0.3216 0.5756
# scale(ANCK)                  1.0757  1.0757     1 15.879  0.1409 0.7123

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = BacterialRichness)
summary(ModelLatitude)
# Number of obs: 372, groups:  StudyID, 94
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 12.6907 12.6907     1 155.650  1.1800 0.2790
# scale(Species_Richness)      2.3358  2.3358     1 274.646  0.2172 0.6416
# scale(Duration)             13.4744 13.4744     1 118.017  1.2529 0.2653
# scale(Latitude)              0.0566  0.0566     1  49.169  0.0053 0.9425

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = BacterialRichness)
summary(ModelLongitude)
# Number of obs: 372, groups:  StudyID, 94
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 12.5855 12.5855     1 155.047  1.1696 0.2812
# scale(Species_Richness)      2.3172  2.3172     1 273.197  0.2153 0.6430
# scale(Duration)             13.0221 13.0221     1 126.246  1.2101 0.2734
# scale(Longitude)             0.0767  0.0767     1  65.194  0.0071 0.9330

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = BacterialRichness)
summary(ModelMAPmean)
# Number of obs: 372, groups:  StudyID, 94
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles))  9.8799  9.8799     1 151.813  0.9169 0.3398
# scale(Species_Richness)      3.3112  3.3112     1 273.024  0.3073 0.5798
# scale(Duration)             10.8914 10.8914     1 112.748  1.0108 0.3169
# scale(MAPmean)              11.1824 11.1824     1  47.271  1.0378 0.3135

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = BacterialRichness)
summary(ModelMATmean)
# Number of obs: 372, groups:  StudyID, 94
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 11.5732 11.5732     1 156.071  1.0774 0.3009
# scale(Species_Richness)      3.0119  3.0119     1 277.634  0.2804 0.5969
# scale(Duration)             12.4804 12.4804     1 117.212  1.1619 0.2833
# scale(MATmean)               4.4301  4.4301     1  47.034  0.4124 0.5239



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## Species_Richness
sum(!is.na(BacterialRichness$Species_Richness)) ## n = 372
p1 <- ggplot(BacterialRichness, aes(y=RR, x=Species_Richness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="Species_Richness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Species_Richness" , y="lnBacterialRichness372")
p1
pdf("Species_Richness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(BacterialRichness$Duration)) ## n = 372
p2 <- ggplot(BacterialRichness, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnBacterialRichness372")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(BacterialRichness$Rotation_cycles)) ## n = 372
p3 <- ggplot(BacterialRichness, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnBacterialRichness372")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(BacterialRichness$Latitude)) ## n = 372
p5 <- ggplot(BacterialRichness, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnBacterialRichness372")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(BacterialRichness$Longitude)) ## n = 372
p6 <- ggplot(BacterialRichness, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnBacterialRichness372")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(BacterialRichness$MAPmean)) ## n = 372
p7 <- ggplot(BacterialRichness, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnBacterialRichness372")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(BacterialRichness$MATmean)) ## n = 372
p8 <- ggplot(BacterialRichness, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnBacterialRichness372")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(BacterialRichness$pHCK)) ## n = 192
p9 <- ggplot(BacterialRichness, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnBacterialRichness192")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(BacterialRichness$SOCCK)) ## n = 204
p10 <- ggplot(BacterialRichness, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnBacterialRichness204")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(BacterialRichness$TNCK)) ## n = 123
p11 <- ggplot(BacterialRichness, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnBacterialRichness123")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(BacterialRichness$NO3CK)) ## n = 67
p12 <- ggplot(BacterialRichness, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnBacterialRichness67")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(BacterialRichness$NH4CK)) ## n = 47
p13<- ggplot(BacterialRichness, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnBacterialRichness47")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(BacterialRichness$APCK)) ## n = 139
p14 <- ggplot(BacterialRichness, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnBacterialRichness139")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(BacterialRichness$AKCK)) ## n = 108
p15 <- ggplot(BacterialRichness, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnBacterialRichness108")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(BacterialRichness$ANCK)) ## n = 90
p16 <- ggplot(BacterialRichness, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnBacterialRichness90")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(BacterialRichness$RRpH)) ## n = 192
p17 <- ggplot(BacterialRichness, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="RRpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRpH" , y="RR192")
p17
pdf("RRpH.pdf",width=8,height=8)
p17
dev.off() 

## RRSOC
sum(!is.na(BacterialRichness$RRSOC)) ## n = 204
p18 <- ggplot(BacterialRichness, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="RRSOC")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRSOC" , y="RR204")
p18
pdf("RRSOC.pdf",width=8,height=8)
p18
dev.off() 

## RRTN
sum(!is.na(BacterialRichness$RRTN)) ## n = 123
p19 <- ggplot(BacterialRichness, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="RRTN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRTN" , y="RR123")
p19
pdf("RRTN.pdf",width=8,height=8)
p19
dev.off() 

## RRNO3
sum(!is.na(BacterialRichness$RRNO3)) ## n = 67
p20 <- ggplot(BacterialRichness, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="RRNO3")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNO3" , y="RR67")
p20
pdf("RRNO3.pdf",width=8,height=8)
p20
dev.off() 

## RRNH4
sum(!is.na(BacterialRichness$RRNH4)) ## n = 47
p21 <- ggplot(BacterialRichness, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="RRNH4")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNH4" , y="RR47")
p21
pdf("RRNH4.pdf",width=8,height=8)
p21
dev.off() 

## RRAP
sum(!is.na(BacterialRichness$RRAP)) ## n = 139
p22 <- ggplot(BacterialRichness, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="RRAP")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAP" , y="RR139")
p22
pdf("RRAP.pdf",width=8,height=8)
p22
dev.off() 

## RRAK
sum(!is.na(BacterialRichness$RRAK)) ## n = 108
p23 <- ggplot(BacterialRichness, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="RRAK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAK" , y="RR108")
p23
pdf("RRAK.pdf",width=8,height=8)
p23
dev.off() 

## RRAN
sum(!is.na(BacterialRichness$RRAN)) ## n = 90
p24 <- ggplot(BacterialRichness, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialRichness", x="RRAN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAN" , y="RR90")
p24
pdf("RRAN.pdf",width=8,height=8)
p24
dev.off() 

## RRYield
sum(!is.na(BacterialRichness$RRYield)) ## n = 46
p25 <- ggplot(BacterialRichness, aes(x=RR, y=RRYield)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm, color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") +
  theme_bw() +
  theme(text = element_text(family = "serif", size=20)) +
  geom_vline(aes(xintercept=0), colour="black", linewidth=0.5, linetype="dashed") +
  labs(x="RR", y="RRYield46") +
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
  scale_x_continuous(limits=c(-0.4, 0.5), expand=c(0, 0)) + 
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
                            data = BacterialRichness, 
                            na.action = na.roughfix, 
                            importance = TRUE, 
                            ntree = 500)
# Check the importance of variables and p-values
importance(rf_model_perm3)
#                    %IncMSE %IncMSE.pval IncNodePurity IncNodePurity.pval
# MATmean          21.405932   0.00990099     1.5644057         0.00990099
# pHCK             15.763589   0.01980198     0.9302883         0.60396040
# SOCCK            15.484919   0.00990099     0.9822298         0.72277228
# Latitude         14.043841   0.03960396     0.9661067         0.13861386
# Longitude        13.070510   0.02970297     0.8484593         0.44554455
# MAPmean          12.330773   0.01980198     0.6843948         0.91089109
# Duration         11.773417   0.05940594     1.1117340         0.00990099
# Rotation_cycles   5.987125   0.61386139     0.8589549         0.61386139
# Species_Richness  3.323212   0.44554455     0.1745820         0.77227723

######################## Trials sorted by effect size
library(ggplot2)
library(dplyr)
BacterialRichness <- read.csv("BacterialRichness.csv", fileEncoding = "latin1")
# 计算95% CI + 显著性分类
df_plot <- BacterialRichness %>%
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
# Negative  Neutral Positive 
#       52      220      100  
## 8*8


################################################  piecewiseSEM
sem_data <- read.csv("BacterialRichness.csv", fileEncoding = "latin1")
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
# -687.786
# Fisher's C = 11.767 with P-value = 0.003 

m1_1 <- lme(RR ~ RRpH + RRSOC + Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_2 <- lme(RRpH ~  Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_3 <- lme(RRSOC ~ Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model1 <- psem(m1_1, m1_2, m1_3)
summary(sem_model1)  ## 
# AIC
# -707.447
# Fisher's C = 12.837 with P-value = 0.002

m2_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_2 <- lme(RRpH ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_3 <- lme(RRSOC ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model2 <- psem(m2_1, m2_2, m2_3)
summary(sem_model2)  ## 
# AIC
# -719.684
# Fisher's C = 10.933 with P-value = 0.004

m3_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_2 <- lme(RRpH ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_3 <- lme(RRSOC ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model3 <- psem(m3_1, m3_2, m3_3)
summary(sem_model3)  ## 
# AIC
# -703.099
# Fisher's C = 14.449 with P-value = 0.001

m4_1 <- lme(RR ~ RRpH + RRSOC + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_2 <- lme(RRpH ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_3 <- lme(RRSOC ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model4 <- psem(m4_1, m4_2, m4_3)
summary(sem_model4) ## 
# AIC
# -741.761
# Fisher's C = 11.978 with P-value = 0.003 

# Structural Equation Model of sem_model4 
# 
# Call:
#   RR ~ RRpH + RRSOC + Species_Richness
# RRpH ~ Species_Richness
# RRSOC ~ Species_Richness
# 
# AIC
# -741.761
# 
# ---
#   Tests of directed separation:
#   
#   Independ.Claim Test.Type  DF Crit.Value P.Value   
# RRSOC ~ RRpH + ...      coef 106     3.0965  0.0025 **
#   
#   --
#   Global goodness-of-fit:
#   
#   Chi-Squared = NA with P-value = NA and on 1 degrees of freedom
# Fisher's C = 11.978 with P-value = 0.003 and on 2 degrees of freedom
# 
# ---
# Coefficients:
# 
#   Response        Predictor Estimate Std.Error  DF Crit.Value P.Value Std.Estimate  
#         RR             RRpH   0.3254    0.1804 105     1.8036  0.0742       0.1184  
#         RR            RRSOC  -0.0040    0.0741 105    -0.0533  0.9576      -0.0042  
#         RR Species_Richness   0.0064    0.0250 105     0.2566  0.7980       0.0162  
#       RRpH Species_Richness   0.0114    0.0113 114     1.0104  0.3144       0.0791  
#      RRSOC Species_Richness   0.0564    0.0274 124     2.0549  0.0420       0.1327 *
# 
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05
# 
# ---
# Individual R-squared:
# 
#   Response method Marginal Conditional
#         RR   none     0.02        0.62
#       RRpH   none     0.00        0.67
#      RRSOC   none     0.01        0.78





m5_1 <- lme(RR ~ RRpH + RRSOC + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_2 <- lme(RRpH ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_3 <- lme(RRSOC ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model5 <- psem(m5_1, m5_2, m5_3)
summary(sem_model5) ## 
# AIC
# -725.041
# Fisher's C = 14.921 with P-value = 0.001

m6_1 <- lme(RR ~ RRpH + RRSOC +MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_2 <- lme(RRpH ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_3 <- lme(RRSOC ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model6 <- psem(m6_1, m6_2, m6_3)
summary(sem_model6)  ##
# AIC
# -737.555
# Fisher's C = 13.091 with P-value = 0.001
