
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
BacterialShannon <- read.csv("BacterialShannon.csv", fileEncoding = "latin1")
# Check data
head(BacterialShannon)

# 1. The number of Obversation
total_number <- nrow(BacterialShannon)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 393

# 2. The number of Study
unique_studyid_number <- length(unique(BacterialShannon$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 93


#### 3. Overall effect size
total_effect_model <- rma.mv(yi = RR, 
                             V = Vi, 
                             random = ~ 1 | StudyID,  # StudyID is radom factor
                             data = BacterialShannon, 
                             method = "REML")
# The results of Overall effect size
summary(total_effect_model)
# Multivariate Meta-Analysis Model (k = 393; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  140.6505  -281.3010  -277.3010  -269.3585  -277.2701   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0015  0.0392     93     no  StudyID 
# Test for Heterogeneity:
# Q(df = 392) = 5347.0610, p-val < .0001
# Model Results:
# estimate      se    zval    pval   ci.lb   ci.ub     
#   0.0131  0.0047  2.7628  0.0057  0.0038  0.0224  ** 

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
boot_results1 <- boot(data = BacterialShannon, statistic = boot_fun, R = 1000, parallel = "snow", ncpus = numCores, cl = cl)
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
# Estimate for coefficient 1 : 0.01310549 
# 95% BCa CI for coefficient 1 : 0.008001568 0.01876952 


#### 5. Funnel Plot
simple_model <- rma(yi = RR, 
                    vi = Vi, 
                    data = BacterialShannon, 
                    method = "REML")
#### 
funnel(simple_model)
# Output  6 * 6
#### Egger's test
regtest(simple_model)
# Regression Test for Funnel Plot Asymmetry
# Model:     mixed-effects meta-regression model
# Predictor: standard error
# Test for Funnel Plot Asymmetry: z = 1.0446, p = 0.2962
# Limit Estimate (as sei -> 0):   b = 0.0062 (CI: -0.0015, 0.0139)

#  Rosenthal’s Fail-Safe N
# This method estimates how many missing studies with null effect 
# would be needed to make the overall effect non-significant
fsn_rosenthal <- fsn(x = simple_model, type = "Rosenthal")
# Print the FSN result
print(fsn_rosenthal)
# Fail-safe N Calculation Using the General Approach
# Average Effect Size:         0.0090 (with file drawer: 0.0034)
# Amount of Heterogeneity:     0.0015 (with file drawer: 0.0015)
# Observed Significance Level: 0.0013 (with file drawer: 0.0501)
# Target Significance Level:   0.05
# Fail-safe N: 344



#### 8. Subgroup analysis
### 8.1 LegumeNonlegume
BacterialShannon_filteredLegumeNonlegume <- subset(BacterialShannon, LegumeNonlegume %in% c("Legume to Non-legume", "Non-legume to Legume", "Non-legume to Non-legume"))
#
BacterialShannon_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(BacterialShannon_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredLegumeNonlegume %>%
  group_by(LegumeNonlegume) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   LegumeNonlegume          Observations Unique_StudyID
#   <fct>                           <int>          <int>
# 1 Legume to Non-legume              109             23
# 2 Non-legume to Legume              133             36
# 3 Non-legume to Non-legume          147             54

overall_model_BacterialShannon_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = BacterialShannon_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 389; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  138.8665  -277.7331  -269.7331  -253.9097  -269.6281   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0015  0.0389     92     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 386) = 4862.6619, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 32.8130, p-val < .0001
# Model Results:
#                                          estimate      se    zval    pval    ci.lb   ci.ub      
# LegumeNonlegumeLegume to Non-legume        0.0213  0.0055  3.9034  <.0001   0.0106  0.0321  *** 
# LegumeNonlegumeNon-legume to Legume        0.0098  0.0052  1.8693  0.0616  -0.0005  0.0200    . 
# LegumeNonlegumeNon-legume to Non-legume    0.0111  0.0052  2.1378  0.0325   0.0009  0.0213    *  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredLegumeNonlegume)
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
    #                                 "a"                                     "b"                                     "b" 


### 8.2 AMnonAM
BacterialShannon_filteredAMnonAM <- subset(BacterialShannon, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
BacterialShannon_filteredAMnonAM$AMnonAM <- droplevels(factor(BacterialShannon_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredAMnonAM %>%
  group_by(AMnonAM) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   AMnonAM     Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 AM to AM             281             74
# 2 AM to nonAM           53             17
# 3 nonAM to AM           56             11

overall_model_BacterialShannon_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = BacterialShannon_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 390; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  131.7933  -263.5867  -255.5867  -239.7530  -255.4820   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0016  0.0395     92     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 387) = 5283.4641, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 7.2372, p-val = 0.0647
# Model Results:
#                     estimate      se    zval    pval    ci.lb   ci.ub    
# AMnonAMAM to AM       0.0125  0.0051  2.4219  0.0154   0.0024  0.0225  * 
# AMnonAMAM to nonAM    0.0123  0.0058  2.1332  0.0329   0.0010  0.0236  * 
# AMnonAMnonAM to AM    0.0155  0.0133  1.1626  0.2450  -0.0106  0.0416    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredAMnonAM)
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
   #             "a"                "a"                "a" 


### 8.3 C3C4
BacterialShannon_filteredC3C4 <- subset(BacterialShannon, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
BacterialShannon_filteredC3C4$C3C4 <- droplevels(factor(BacterialShannon_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredC3C4 %>%
  group_by(C3C4) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   C3C4     Observations Unique_StudyID
#   <fct>           <int>          <int>
# 1 C3 to C3          200             60
# 2 C3 to C4          130             39
# 3 C4 to C3           54             16

overall_model_BacterialShannon_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = BacterialShannon_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 384; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  158.0868  -316.1735  -308.1735  -292.4024  -308.0672   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0014  0.0379     91     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 381) = 4755.5316, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 66.9251, p-val < .0001
# Model Results:
#               estimate      se     zval    pval    ci.lb   ci.ub      
# C3C4C3 to C3    0.0146  0.0048   3.0512  0.0023   0.0052  0.0240   ** 
# C3C4C3 to C4    0.0168  0.0049   3.4442  0.0006   0.0072  0.0264  *** 
# C3C4C4 to C3   -0.0074  0.0056  -1.3320  0.1829  -0.0184  0.0035      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredC3C4)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredC3C4)
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
BacterialShannon_filteredAnnual_Pere <- subset(BacterialShannon, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
BacterialShannon_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(BacterialShannon_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredAnnual_Pere %>%
  group_by(Annual_Pere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Annual_Pere            Observations Unique_StudyID
#   <fct>                         <int>          <int>
# 1 Annual to Annual                294             69
# 2 Annual to Perennial              35             14
# 3 Perennial to Annual              48             15
# 4 Perennial to Perennial           15              7

overall_model_BacterialShannon_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = BacterialShannon_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 392; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  144.7210  -289.4420  -279.4420  -259.6370  -279.2849   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0014  0.0380     92     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 388) = 5163.9992, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 32.4123, p-val < .0001
# Model Results:
#                                    estimate      se     zval    pval    ci.lb   ci.ub     
# Annual_PereAnnual to Annual          0.0112  0.0052   2.1559  0.0311   0.0010  0.0214   * 
# Annual_PereAnnual to Perennial      -0.0011  0.0063  -0.1797  0.8574  -0.0136  0.0113     
# Annual_PerePerennial to Annual       0.0322  0.0103   3.1313  0.0017   0.0121  0.0524  ** 
# Annual_PerePerennial to Perennial    0.0100  0.0107   0.9310  0.3518  -0.0110  0.0309    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredAnnual_Pere)
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
      #                        "ab"                               "c"                               "a"                              "bc" 


### 8.6 PlantStageSubgroup
BacterialShannon_filteredPlantStageSubgroup <- subset(BacterialShannon, PlantStageSubgroup %in% c("Vegetative stage","Reproductive stage", "Maturity stage","Harvest"))
#
BacterialShannon_filteredPlantStageSubgroup$PlantStageSubgroup <- droplevels(factor(BacterialShannon_filteredPlantStageSubgroup$PlantStageSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredPlantStageSubgroup %>%
  group_by(PlantStageSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   PlantStageSubgroup Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Harvest                     168             48
# 2 Maturity stage               41             17
# 3 Reproductive stage           81             14
# 4 Vegetative stage             44              8

overall_model_BacterialShannon_filteredPlantStageSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + PlantStageSubgroup, random = ~ 1 | StudyID, data = BacterialShannon_filteredPlantStageSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredPlantStageSubgroup)
# Multivariate Meta-Analysis Model (k = 334; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  332.6374  -665.2748  -655.2748  -636.2794  -655.0897 
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0026  0.0512     77     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 330) = 4594.8562, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 463.5185, p-val < .0001
# Model Results:
#                                       estimate      se     zval    pval    ci.lb    ci.ub      
# PlantStageSubgroupHarvest               0.0340  0.0066   5.1243  <.0001   0.0210   0.0470  *** 
# PlantStageSubgroupMaturity stage       -0.0604  0.0074  -8.2088  <.0001  -0.0748  -0.0459  *** 
# PlantStageSubgroupReproductive stage    0.0162  0.0075   2.1587  0.0309   0.0015   0.0310    * 
# PlantStageSubgroupVegetative stage     -0.0048  0.0099  -0.4857  0.6272  -0.0242   0.0146        

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredPlantStageSubgroup)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredPlantStageSubgroup)
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
           #                       "a"                                  "b"                                  "c"                                  "d" 


### 8.7 Bulk_Rhizosphere
BacterialShannon_filteredBulk_Rhizosphere <- subset(BacterialShannon, Bulk_Rhizosphere %in% c("Non-Rhizosphere", "Rhizosphere"))
#
BacterialShannon_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(BacterialShannon_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredBulk_Rhizosphere %>%
  group_by(Bulk_Rhizosphere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Bulk_Rhizosphere Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 Non-Rhizosphere           248             57
# 2 Rhizosphere               145             45

overall_model_BacterialShannon_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = BacterialShannon_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 393; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  250.3338  -500.6677  -494.6677  -482.7616  -494.6057   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0018  0.0428     93     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 391) = 5334.9825, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 232.8051, p-val < .0001
# Model Results:
#                                  estimate      se     zval    pval    ci.lb   ci.ub      
# Bulk_RhizosphereNon-Rhizosphere    0.0285  0.0052   5.4812  <.0001   0.0183  0.0387  *** 
# Bulk_RhizosphereRhizosphere       -0.0080  0.0053  -1.5135  0.1301  -0.0184  0.0024   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredBulk_Rhizosphere)
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
BacterialShannon_filteredSoil_texture <- subset(BacterialShannon, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
BacterialShannon_filteredSoil_texture$Soil_texture <- droplevels(factor(BacterialShannon_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredSoil_texture %>%
  group_by(Soil_texture) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Soil_texture Observations Unique_StudyID
#   <fct>               <int>          <int>
# 1 Coarse                 78             16
# 2 Fine                   98             17
# 3 Medium                155             32

overall_model_BacterialShannon_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = BacterialShannon_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 331; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
#  38.3915  -76.7830  -68.7830  -53.6110  -68.6592   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0014  0.0376     65     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 328) = 4306.1928, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.7554, p-val = 0.4309
# Model Results:
#                     estimate      se    zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse    0.0159  0.0103  1.5531  0.1204  -0.0042  0.0360    
# Soil_textureFine      0.0050  0.0103  0.4811  0.6305  -0.0153  0.0252    
# Soil_textureMedium    0.0027  0.0080  0.3344  0.7381  -0.0130  0.0184    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredSoil_texture)
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
BacterialShannon_filteredTillage <- subset(BacterialShannon, Tillage %in% c("Tillage", "No_tillage"))
#
BacterialShannon_filteredTillage$Tillage <- droplevels(factor(BacterialShannon_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredTillage %>%
  group_by(Tillage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Tillage    Observations Unique_StudyID
#   <fct>             <int>          <int>
# 1 No_tillage           30              4
# 2 Tillage              47             13

overall_model_BacterialShannon_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = BacterialShannon_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredTillage)
# Multivariate Meta-Analysis Model (k = 77; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#   93.4294  -186.8588  -180.8588  -173.9063  -180.5208   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0033  0.0572     15     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 75) = 352.4073, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.6378, p-val = 0.7269
# Model Results:
#                    estimate      se    zval    pval    ci.lb   ci.ub    
# TillageNo_tillage    0.0030  0.0172  0.1736  0.8622  -0.0307  0.0367    
# TillageTillage       0.0084  0.0161  0.5222  0.6015  -0.0231  0.0400   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredTillage)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredTillage)
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
BacterialShannon_filteredStraw_retention <- subset(BacterialShannon, Straw_retention %in% c("Retention", "No_retention"))
#
BacterialShannon_filteredStraw_retention$Straw_retention <- droplevels(factor(BacterialShannon_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredStraw_retention %>%
  group_by(Straw_retention) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Straw_retention Observations Unique_StudyID
#   <fct>                  <int>          <int>
# 1 No_retention              16              8
# 2 Retention                 44             11

overall_model_BacterialShannon_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = BacterialShannon_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 60; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
#  27.9223  -55.8447  -49.8447  -43.6633  -49.4002   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0030  0.0548     18     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 58) = 649.5153, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.1847, p-val = 0.5530
# Model Results:
#                              estimate      se    zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention    0.0057  0.0139  0.4083  0.6830  -0.0216  0.0330    
# Straw_retentionRetention       0.0095  0.0139  0.6820  0.4953  -0.0177  0.0367    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredStraw_retention)
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
BacterialShannon_filteredRotationcyclesSubgroup <- subset(BacterialShannon, RotationcyclesSubgroup %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
BacterialShannon_filteredRotationcyclesSubgroup$RotationcyclesSubgroup <- droplevels(factor(BacterialShannon_filteredRotationcyclesSubgroup$RotationcyclesSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredRotationcyclesSubgroup %>%
  group_by(RotationcyclesSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   RotationcyclesSubgroup Observations Unique_StudyID
#   <fct>                         <int>          <int>
# 1 D1                              124             40
# 2 D1-3                             83             25
# 3 D10                              34             10
# 4 D3-5                             70             13
# 5 D5-10                            82             19

overall_model_BacterialShannon_filteredRotationcyclesSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + RotationcyclesSubgroup, random = ~ 1 | StudyID, data = BacterialShannon_filteredRotationcyclesSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredRotationcyclesSubgroup)
# Multivariate Meta-Analysis Model (k = 393; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  143.4018  -286.8036  -274.8036  -251.0376  -274.5831   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0017  0.0410     93     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 388) = 5232.7042, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 27.9456, p-val < .0001
# Model Results:
#                              estimate      se     zval    pval    ci.lb   ci.ub     
# RotationcyclesSubgroupD1       0.0197  0.0061   3.2341  0.0012   0.0078  0.0317  ** 
# RotationcyclesSubgroupD1-3    -0.0025  0.0064  -0.3876  0.6983  -0.0151  0.0101     
# RotationcyclesSubgroupD10      0.0275  0.0128   2.1393  0.0324   0.0023  0.0526   * 
# RotationcyclesSubgroupD3-5     0.0072  0.0102   0.7030  0.4821  -0.0128  0.0271     
# RotationcyclesSubgroupD5-10    0.0207  0.0090   2.3056  0.0211   0.0031  0.0382   * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredRotationcyclesSubgroup)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredRotationcyclesSubgroup)
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
   #                      "a"                         "b"                         "a"                        "ab"                         "a" 


### 8.13 DurationSubgroup
BacterialShannon_filteredDurationSubgroup <- subset(BacterialShannon, DurationSubgroup %in% c("D1", "D2", "D3", "D4", "D5", "D6-10", "D11-20", "D20-30", "D30"))
#
BacterialShannon_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(BacterialShannon_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredDurationSubgroup %>%
  group_by(DurationSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   DurationSubgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D1                         21             10
# 2 D11-20                     76             15
# 3 D2                         45             21
# 4 D20-30                     37              9
# 5 D3                         55             16
# 6 D30                        63              4
# 7 D4                         31              9
# 8 D5                         10              2
# 9 D6-10                      55             16

overall_model_BacterialShannon_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = BacterialShannon_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 393; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  144.1662  -288.3324  -268.3324  -228.8259  -267.7426   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0016  0.0400     93     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 384) = 4650.4173, p-val < .0001
# Test of Moderators (coefficients 1:9):
# QM(df = 9) = 40.6573, p-val < .0001
# Model Results:
#                         estimate      se     zval    pval    ci.lb   ci.ub      
# DurationSubgroupD1       -0.0101  0.0132  -0.7674  0.4428  -0.0360  0.0157      
# DurationSubgroupD11-20    0.0329  0.0097   3.3779  0.0007   0.0138  0.0520  *** 
# DurationSubgroupD2        0.0243  0.0081   3.0097  0.0026   0.0085  0.0401   ** 
# DurationSubgroupD20-30    0.0228  0.0135   1.6865  0.0917  -0.0037  0.0494    . 
# DurationSubgroupD3       -0.0003  0.0078  -0.0352  0.9719  -0.0155  0.0150      
# DurationSubgroupD30      -0.0173  0.0218  -0.7930  0.4278  -0.0600  0.0254      
# DurationSubgroupD4       -0.0070  0.0139  -0.5029  0.6150  -0.0343  0.0203      
# DurationSubgroupD5       -0.0131  0.0441  -0.2970  0.7665  -0.0995  0.0733      
# DurationSubgroupD6-10     0.0334  0.0110   3.0267  0.0025   0.0118  0.0550   ** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredDurationSubgroup)
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
    #                "a"                    "b"                   "bc"                  "abc"                    "a"                   "ac" 
    # DurationSubgroupD4     DurationSubgroupD5  DurationSubgroupD6-10 
    #               "ac"                  "abc"                    "b" 



### 8.14 SpeciesRichnessSubgroup
BacterialShannon_filteredSpeciesRichnessSubgroup <- subset(BacterialShannon, SpeciesRichnessSubgroup %in% c("R2", "R3", "R4"))
#
BacterialShannon_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup <- droplevels(factor(BacterialShannon_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredSpeciesRichnessSubgroup %>%
  group_by(SpeciesRichnessSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   SpeciesRichnessSubgroup Observations Unique_StudyID
#   <fct>                          <int>          <int>
# 1 R2                               290             82
# 2 R3                                93             20
# 3 R4                                10              3

overall_model_BacterialShannon_filteredSpeciesRichnessSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichnessSubgroup, random = ~ 1 | StudyID, data = BacterialShannon_filteredSpeciesRichnessSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredSpeciesRichnessSubgroup)
# Multivariate Meta-Analysis Model (k = 393; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  138.8380  -277.6761  -269.6761  -253.8115  -269.5722   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0016  0.0398     93     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 390) = 5307.3556, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 12.6975, p-val = 0.0053
# Model Results:
#                            estimate      se     zval    pval    ci.lb   ci.ub      
# SpeciesRichnessSubgroupR2    0.0128  0.0049   2.6424  0.0082   0.0033  0.0224   ** 
# SpeciesRichnessSubgroupR3    0.0179  0.0052   3.4194  0.0006   0.0076  0.0281  *** 
# SpeciesRichnessSubgroupR4   -0.0094  0.0320  -0.2939  0.7689  -0.0721  0.0533      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredSpeciesRichnessSubgroup)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredSpeciesRichnessSubgroup)
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
#                       "a"                       "b"                      "ab" 


### 8.15 Primer
BacterialShannon_filteredPrimer <- subset(BacterialShannon, Primer %in% c("V1-V3", "V1-V4", "V1-V5", "V3-V4", "V3-V5", "V4", "V4-V5", "V5-V7", "V6-V8", "V9", "Full length"))
#
BacterialShannon_filteredPrimer$Primer <- droplevels(factor(BacterialShannon_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredPrimer %>%
  group_by(Primer) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Primer      Observations Unique_StudyID
#    <fct>              <int>          <int>
#  1 Full length            3              1
#  2 V1-V3                  7              3
#  3 V1-V4                  4              1
#  4 V1-V5                  6              1
#  5 V3-V4                191             46
#  6 V3-V5                  3              2
#  7 V4                    99             20
#  8 V4-V5                 58             15
#  9 V5-V7                  3              1
# 10 V6-V8                 16              1
# 11 V9                     1              1

overall_model_BacterialShannon_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = BacterialShannon_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 391; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  123.4230  -246.8460  -222.8460  -175.5639  -221.9958   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0017  0.0414     92     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 380) = 5225.9241, p-val < .0001
# 
# Test of Moderators (coefficients 1:11):
# QM(df = 11) = 9.5039, p-val = 0.5755
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerFull length    0.0045  0.0582   0.0770  0.9386  -0.1095  0.1185    
# PrimerV1-V3          0.0239  0.0309   0.7743  0.4388  -0.0366  0.0845    
# PrimerV1-V4          0.0023  0.0416   0.0544  0.9566  -0.0792  0.0837    
# PrimerV1-V5          0.0509  0.0505   1.0083  0.3133  -0.0481  0.1499    
# PrimerV3-V4          0.0095  0.0070   1.3546  0.1755  -0.0042  0.0232    
# PrimerV3-V5          0.0260  0.0426   0.6101  0.5418  -0.0576  0.1096    
# PrimerV4             0.0143  0.0108   1.3279  0.1842  -0.0068  0.0355    
# PrimerV4-V5          0.0153  0.0122   1.2600  0.2077  -0.0085  0.0392    
# PrimerV5-V7         -0.0095  0.0430  -0.2212  0.8249  -0.0938  0.0748    
# PrimerV6-V8          0.0240  0.0450   0.5337  0.5936  -0.0642  0.1123    
# PrimerV9             0.0773  0.0548   1.4095  0.1587  -0.0302  0.1847    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredPrimer)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredPrimer)
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
#               "a"               "a"               "a"               "a"               "a"               "a"               "a"               "a" 
#       PrimerV5-V7       PrimerV6-V8          PrimerV9 
#               "a"               "a"               "a" 


### 8.16 Latitude_Subgroup
BacterialShannon_filteredLatitude_Subgroup <- subset(BacterialShannon, Latitude_Subgroup %in% c("La20", "La20-40", "La40"))
#
BacterialShannon_filteredLatitude_Subgroup$Latitude_Subgroup <- droplevels(factor(BacterialShannon_filteredLatitude_Subgroup$Latitude_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredLatitude_Subgroup %>%
  group_by(Latitude_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Latitude_Subgroup Observations Unique_StudyID
#   <fct>                    <int>          <int>
# 1 La20                        15              5
# 2 La20-40                    129             51
# 3 La40                       249             37

overall_model_BacterialShannon_filteredLatitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Latitude_Subgroup, random = ~ 1 | StudyID, data = BacterialShannon_filteredLatitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredLatitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 393; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  137.9479  -275.8958  -267.8958  -252.0312  -267.7919   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0016  0.0398     93     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 390) = 5331.7468, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 8.4328, p-val = 0.0379
# Model Results:
#                           estimate      se    zval    pval    ci.lb   ci.ub    
# Latitude_SubgroupLa20       0.0193  0.0206  0.9368  0.3489  -0.0211  0.0596    
# Latitude_SubgroupLa20-40    0.0167  0.0065  2.5590  0.0105   0.0039  0.0294  * 
# Latitude_SubgroupLa40       0.0076  0.0076  1.0032  0.3157  -0.0072  0.0224       

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredLatitude_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredLatitude_Subgroup)
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
BacterialShannon_filteredLongitude_Subgroup <- subset(BacterialShannon, Longitude_Subgroup %in% c("Lo-180-0", "Lo-0-180"))
#
BacterialShannon_filteredLongitude_Subgroup$Longitude_Subgroup <- droplevels(factor(BacterialShannon_filteredLongitude_Subgroup$Longitude_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredLongitude_Subgroup %>%
  group_by(Longitude_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Longitude_Subgroup Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Lo-0-180                    250             81
# 2 Lo-180-0                    143             12

overall_model_BacterialShannon_filteredLongitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Subgroup, random = ~ 1 | StudyID, data = BacterialShannon_filteredLongitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredLongitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 393; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  139.6469  -279.2938  -273.2938  -261.3877  -273.2318   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0016  0.0396     93     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 391) = 5347.3159, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 7.8343, p-val = 0.0199
# Model Results:
#                             estimate      se    zval    pval    ci.lb   ci.ub     
# Longitude_SubgroupLo-0-180    0.0141  0.0051  2.7582  0.0058   0.0041  0.0240  ** 
# Longitude_SubgroupLo-180-0    0.0065  0.0137  0.4759  0.6341  -0.0203  0.0334          

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredLongitude_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredLongitude_Subgroup)
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
BacterialShannon_filteredMAPmean_Subgroup <- subset(BacterialShannon, MAPmean_Subgroup %in% c("MAP600", "MAP600-1200", "MAP1200"))
#
BacterialShannon_filteredMAPmean_Subgroup$MAPmean_Subgroup <- droplevels(factor(BacterialShannon_filteredMAPmean_Subgroup$MAPmean_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredMAPmean_Subgroup %>%
  group_by(MAPmean_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MAPmean_Subgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 MAP1200                    65             24
# 2 MAP600                    179             37
# 3 MAP600-1200               149             33

overall_model_BacterialShannon_filteredMAPmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean_Subgroup, random = ~ 1 | StudyID, data = BacterialShannon_filteredMAPmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredMAPmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 393; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  138.7218  -277.4436  -269.4436  -253.5790  -269.3397   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0015  0.0394     93     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 390) = 5310.8921, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 9.9342, p-val = 0.0191
# Model Results:
#                              estimate      se    zval    pval    ci.lb   ci.ub    
# MAPmean_SubgroupMAP1200        0.0201  0.0095  2.1236  0.0337   0.0015  0.0386  * 
# MAPmean_SubgroupMAP600         0.0049  0.0072  0.6764  0.4988  -0.0093  0.0190    
# MAPmean_SubgroupMAP600-1200    0.0185  0.0082  2.2469  0.0246   0.0024  0.0347  * 
       

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredMAPmean_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredMAPmean_Subgroup)
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
BacterialShannon_filteredMATmean_Subgroup <- subset(BacterialShannon, MATmean_Subgroup %in% c("MAT8", "MAT8-15", "MAT15"))
#
BacterialShannon_filteredMATmean_Subgroup$MATmean_Subgroup <- droplevels(factor(BacterialShannon_filteredMATmean_Subgroup$MATmean_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredMATmean_Subgroup %>%
  group_by(MATmean_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MATmean_Subgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 MAT15                      86             35
# 2 MAT8                      236             39
# 3 MAT8-15                    71             20

overall_model_BacterialShannon_filteredMATmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Subgroup, random = ~ 1 | StudyID, data = BacterialShannon_filteredMATmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredMATmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 393; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  138.8711  -277.7423  -269.7423  -253.8777  -269.6384   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0015  0.0391     93     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 390) = 5298.2915, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 10.4429, p-val = 0.0152
# Model Results:
#                          estimate      se    zval    pval    ci.lb   ci.ub     
# MATmean_SubgroupMAT15      0.0232  0.0081  2.8711  0.0041   0.0074  0.0391  ** 
# MATmean_SubgroupMAT8       0.0103  0.0070  1.4664  0.1425  -0.0035  0.0240     
# MATmean_SubgroupMAT8-15    0.0027  0.0100  0.2707  0.7866  -0.0169  0.0224     


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredMATmean_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredMATmean_Subgroup)
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
BacterialShannon_filteredDurationSubgroup <- subset(BacterialShannon, DurationSubgroup %in% c("D5", "D5-10", "D10-20", "D20-30", "D30"))
#
BacterialShannon_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(BacterialShannon_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialShannon_filteredDurationSubgroup %>%
  group_by(DurationSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   DurationSubgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D10-20                     76             15
# 2 D20-30                     37              9
# 3 D30                        63              4
# 4 D5                        162             53
# 5 D5-10                      55             16

overall_model_BacterialShannon_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = BacterialShannon_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialShannon_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 393; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
#  141.9029  -283.8058  -271.8058  -248.0398  -271.5854   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0016  0.0395     93     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 388) = 4952.9548, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 22.0112, p-val = 0.0005
# Model Results:
#                         estimate      se     zval    pval    ci.lb   ci.ub      
# DurationSubgroupD10-20    0.0343  0.0095   3.6126  0.0003   0.0157  0.0529  *** 
# DurationSubgroupD20-30    0.0229  0.0134   1.7020  0.0887  -0.0035  0.0492    . 
# DurationSubgroupD30      -0.0173  0.0215  -0.8043  0.4212  -0.0595  0.0249      
# DurationSubgroupD5        0.0039  0.0060   0.6548  0.5126  -0.0078  0.0156      
# DurationSubgroupD5-10     0.0335  0.0109   3.0655  0.0022   0.0121  0.0549   ** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialShannon_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_BacterialShannon_filteredDurationSubgroup)
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
#                    "a"                   "ab"                    "b"                    "b"                    "a" 






#### 9. Linear Mixed Effect Model
# 
BacterialShannon$Wr <- 1 / BacterialShannon$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialShannon)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialShannon)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialShannon)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialShannon)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialShannon)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialShannon)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialShannon)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialShannon)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(Species_Richness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = BacterialShannon)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = BacterialShannon)

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
#           Model       AIC       BIC   logLik  Marginal_R2 Conditional_R2
# Model1   Model1 -803.3363 -779.4935 407.6682 2.090231e-05   0.0002807545
# Model2   Model2 -801.5121 -777.6693 406.7561 2.260485e-06   0.0002643959
# Model3   Model3 -801.8640 -778.0211 406.9320 4.787966e-06   0.0002660542
# Model4   Model4 -803.2839 -779.4410 407.6419 2.084924e-05   0.0002806321
# Model5   Model5 -801.1119 -777.2690 406.5559 1.892810e-06   0.0002696911
# Model6   Model6 -801.0583 -777.2154 406.5291 1.882438e-06   0.0002696281
# Model7   Model7 -801.5667 -777.7239 406.7834 2.265774e-06   0.0002643773
# Model8   Model8 -801.8161 -777.9732 406.9080 4.801892e-06   0.0002659928
# Model9   Model9 -768.4815 -728.7434 394.2408 3.901833e-05   0.0002966788
# Model10 Model10 -766.8162 -727.0781 393.4081 4.062126e-05   0.0002971489

##### Model 1 is the best model
summary(Model1)
# Number of obs: 393, groups:  StudyID, 93
anova(Model1) 
#                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)  6.6506  6.6506     1 111.39  1.5995 0.2086
# scale(Species_Richness) 2.2897  2.2897     1 388.95  0.5507 0.4585
# scale(Duration)         7.7427  7.7427     1  92.48  1.8622 0.1757


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = BacterialShannon)
summary(ModelpH)
# Number of obs: 196, groups:  StudyID, 60
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                          Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)   4.2640  4.2640     1  40.144  0.6804 0.4143
# scale(Species_Richness)  0.0035  0.0035     1 178.803  0.0006 0.9812
# scale(Duration)         11.8803 11.8803     1  39.267  1.8958 0.1764
# scale(pHCK)              0.2800  0.2800     1  56.084  0.0447 0.8334

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = BacterialShannon)
summary(ModelSOC)
# Number of obs: 203, groups:  StudyID, 63
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                          Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)   6.6096  6.6096     1  45.902  1.0452 0.3120
# scale(Species_Richness)  2.9226  2.9226     1 175.891  0.4622 0.4975
# scale(Duration)         14.3919 14.3919     1  42.417  2.2759 0.1388
# scale(SOCCK)             0.6950  0.6950     1  37.419  0.1099 0.7421 

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = BacterialShannon)
summary(ModelTN)
# Number of obs: 126, groups:  StudyID, 42
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)  0.9512  0.9512     1 42.537  0.6292 0.4321
# scale(Species_Richness) 3.8239  3.8239     1 87.873  2.5293 0.1153
# scale(Duration)         3.4392  3.4392     1 50.089  2.2748 0.1378
# scale(TNCK)             0.6841  0.6841     1 30.829  0.4525 0.5062  


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = BacterialShannon)
summary(ModelNO3)
# Number of obs: 75, groups:  StudyID, 22
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles)  0.2532  0.2532     1 18.440  0.1802 0.67610  
# scale(Species_Richness) 0.7050  0.7050     1 49.650  0.5018 0.48201  
# scale(Duration)         1.1768  1.1768     1 19.769  0.8376 0.37111  
# scale(NO3CK)            5.9393  5.9393     1 42.319  4.2275 0.04598 *

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = BacterialShannon)
summary(ModelNH4)
# Number of obs: 57, groups:  StudyID, 20
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                          Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)  0.96774 0.96774     1 14.628  0.4912 0.4944
# scale(Species_Richness) 0.82436 0.82436     1 37.083  0.4184 0.5217
# scale(Duration)         1.87121 1.87121     1 15.964  0.9498 0.3443
# scale(NH4CK)            0.93097 0.93097     1 20.605  0.4726 0.4995

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = BacterialShannon)
summary(ModelAP)
# Number of obs: 142, groups:  StudyID, 49
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                          Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles)  15.7632 15.7632     1  36.772  1.9214 0.17405  
# scale(Species_Richness)  2.4405  2.4405     1 125.055  0.2975 0.58643  
# scale(Duration)         24.9192 24.9192     1  45.402  3.0375 0.08813 .
# scale(APCK)              0.0594  0.0594     1  61.659  0.0072 0.93244  

#### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = BacterialShannon)
summary(ModelAK)
# Number of obs: 117, groups:  StudyID, 42
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles)   7.742   7.742     1 73.544  0.8422 0.36176  
# scale(Species_Richness)  2.346   2.346     1 98.535  0.2552 0.61454  
# scale(Duration)          9.263   9.263     1 62.954  1.0078 0.31928  
# scale(AKCK)             40.594  40.594     1 72.484  4.4164 0.03907 *

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = BacterialShannon)
summary(ModelAN)
# Number of obs: 89, groups:  StudyID, 24
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)   3.606   3.606     1 60.537  0.3227 0.5721
# scale(Species_Richness)  2.756   2.756     1 77.536  0.2466 0.6209
# scale(Duration)          3.606   3.606     1 35.921  0.3227 0.5735
# scale(ANCK)             35.171  35.171     1 12.687  3.1474 0.1000

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = BacterialShannon)
summary(ModelLatitude)
# Number of obs: 393, groups:  StudyID, 93
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)  6.6573  6.6573     1 111.52  1.6020 0.2083
# scale(Species_Richness) 2.5134  2.5134     1 386.97  0.6048 0.4372
# scale(Duration)         7.8217  7.8217     1  93.53  1.8822 0.1734
# scale(Latitude)         0.3921  0.3921     1  62.40  0.0944 0.7597

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = BacterialShannon)
summary(ModelLongitude)
# Number of obs: 393, groups:  StudyID, 93
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)  7.4039  7.4039     1 119.56  1.7781 0.1849
# scale(Species_Richness) 2.3056  2.3056     1 387.96  0.5537 0.4573
# scale(Duration)         8.4181  8.4181     1 101.17  2.0216 0.1581
# scale(Longitude)        0.8185  0.8185     1 115.09  0.1966 0.6583

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = BacterialShannon)
summary(ModelMAPmean)
# Number of obs: 393, groups:  StudyID, 93
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)  6.1186  6.1186     1 112.51  1.4727 0.2275
# scale(Species_Richness) 2.4477  2.4477     1 387.59  0.5891 0.4432
# scale(Duration)         7.1923  7.1923     1  90.86  1.7311 0.1916
# scale(MAPmean)          0.2790  0.2790     1  59.85  0.0672 0.7964

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = BacterialShannon)
summary(ModelMATmean)
# Number of obs: 393, groups:  StudyID, 93
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)  5.9675  5.9675     1 111.01  1.4367 0.2332
# scale(Species_Richness) 2.4717  2.4717     1 387.59  0.5951 0.4409
# scale(Duration)         7.0381  7.0381     1  90.18  1.6944 0.1963
# scale(MATmean)          0.3918  0.3918     1  62.17  0.0943 0.7598 



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## Species_Richness
sum(!is.na(BacterialShannon$Species_Richness)) ## n = 393
p1 <- ggplot(BacterialShannon, aes(y=RR, x=Species_Richness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="Species_Richness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Species_Richness" , y="lnBacterialShannon393")
p1
pdf("Species_Richness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(BacterialShannon$Duration)) ## n = 393
p2 <- ggplot(BacterialShannon, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnBacterialShannon393")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(BacterialShannon$Rotation_cycles)) ## n = 393
p3 <- ggplot(BacterialShannon, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnBacterialShannon393")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(BacterialShannon$Latitude)) ## n = 393
p5 <- ggplot(BacterialShannon, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnBacterialShannon393")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(BacterialShannon$Longitude)) ## n = 393
p6 <- ggplot(BacterialShannon, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnBacterialShannon393")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(BacterialShannon$MAPmean)) ## n = 393
p7 <- ggplot(BacterialShannon, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnBacterialShannon393")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(BacterialShannon$MATmean)) ## n = 393
p8 <- ggplot(BacterialShannon, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnBacterialShannon393")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(BacterialShannon$pHCK)) ## n = 196
p9 <- ggplot(BacterialShannon, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnBacterialShannon196")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(BacterialShannon$SOCCK)) ## n = 203
p10 <- ggplot(BacterialShannon, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnBacterialShannon203")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(BacterialShannon$TNCK)) ## n = 126
p11 <- ggplot(BacterialShannon, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnBacterialShannon126")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(BacterialShannon$NO3CK)) ## n = 75
p12 <- ggplot(BacterialShannon, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnBacterialShannon75")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(BacterialShannon$NH4CK)) ## n = 57
p13<- ggplot(BacterialShannon, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnBacterialShannon57")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(BacterialShannon$APCK)) ## n = 142
p14 <- ggplot(BacterialShannon, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnBacterialShannon142")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(BacterialShannon$AKCK)) ## n = 117
p15 <- ggplot(BacterialShannon, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnBacterialShannon117")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(BacterialShannon$ANCK)) ## n = 89
p16 <- ggplot(BacterialShannon, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnBacterialShannon89")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(BacterialShannon$RRpH)) ## n = 196
p17 <- ggplot(BacterialShannon, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="RRpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRpH" , y="RR196")
p17
pdf("RRpH.pdf",width=8,height=8)
p17
dev.off() 

## RRSOC
sum(!is.na(BacterialShannon$RRSOC)) ## n = 203
p18 <- ggplot(BacterialShannon, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="RRSOC")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRSOC" , y="RR203")
p18
pdf("RRSOC.pdf",width=8,height=8)
p18
dev.off() 

## RRTN
sum(!is.na(BacterialShannon$RRTN)) ## n = 126
p19 <- ggplot(BacterialShannon, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="RRTN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRTN" , y="RR126")
p19
pdf("RRTN.pdf",width=8,height=8)
p19
dev.off() 

## RRNO3
sum(!is.na(BacterialShannon$RRNO3)) ## n = 75
p20 <- ggplot(BacterialShannon, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="RRNO3")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNO3" , y="RR75")
p20
pdf("RRNO3.pdf",width=8,height=8)
p20
dev.off() 

## RRNH4
sum(!is.na(BacterialShannon$RRNH4)) ## n = 57
p21 <- ggplot(BacterialShannon, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="RRNH4")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNH4" , y="RR57")
p21
pdf("RRNH4.pdf",width=8,height=8)
p21
dev.off() 

## RRAP
sum(!is.na(BacterialShannon$RRAP)) ## n = 142
p22 <- ggplot(BacterialShannon, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="RRAP")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAP" , y="RR142")
p22
pdf("RRAP.pdf",width=8,height=8)
p22
dev.off() 

## RRAK
sum(!is.na(BacterialShannon$RRAK)) ## n = 117
p23 <- ggplot(BacterialShannon, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="RRAK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAK" , y="RR117")
p23
pdf("RRAK.pdf",width=8,height=8)
p23
dev.off() 

## RRAN
sum(!is.na(BacterialShannon$RRAN)) ## n = 89
p24 <- ggplot(BacterialShannon, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialShannon", x="RRAN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAN" , y="RR89")
p24
pdf("RRAN.pdf",width=8,height=8)
p24
dev.off() 

## RRYield
sum(!is.na(BacterialShannon$RRYield)) ## n = 78
p25 <- ggplot(BacterialShannon, aes(x=RR, y=RRYield)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm, color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") +
  theme_bw() +
  theme(text = element_text(family = "serif", size=20)) +
  geom_vline(aes(xintercept=0), colour="black", linewidth=0.5, linetype="dashed") +
  labs(x="RR", y="RRYield78") +
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
  scale_x_continuous(limits=c(-0.2, 0.35), expand=c(0, 0)) + 
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
                            data = BacterialShannon, 
                            na.action = na.roughfix, 
                            importance = TRUE, 
                            ntree = 500)
# Check the importance of variables and p-values
importance(rf_model_perm3)
#                    %IncMSE %IncMSE.pval IncNodePurity IncNodePurity.pval
# MATmean          12.819434   0.14851485    0.11674771         0.05940594
# pHCK             12.770222   0.05940594    0.07654416         1.00000000
# MAPmean          12.215910   0.20792079    0.08673003         0.86138614
# SOCCK            12.194447   0.00990099    0.09767459         0.96039604
# Latitude         10.659123   0.44554455    0.08295573         0.89108911
# Duration          9.620362   0.45544554    0.10232333         0.08910891
# Longitude         9.159929   0.33663366    0.07975405         0.96039604
# Rotation_cycles   9.159644   0.38613861    0.14865542         0.04950495
# Species_Richness  4.577730   0.31683168    0.02678090         0.75247525

######################  Trials sorted by effect size
library(ggplot2)
library(dplyr)
BacterialShannon <- read.csv("BacterialShannon.csv", fileEncoding = "latin1")
# 95% CI + significant
df_plot <- BacterialShannon %>%
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
#       30      311       52 
## 8*8


################################################  piecewiseSEM
sem_data <- read.csv("BacterialShannon.csv", fileEncoding = "latin1")
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
# -968.547
# Fisher's C = 5.633 with P-value = 0.06

m1_1 <- lme(RR ~ RRpH + RRSOC + Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_2 <- lme(RRpH ~  Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_3 <- lme(RRSOC ~ Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model1 <- psem(m1_1, m1_2, m1_3)
summary(sem_model1)  ## 
# AIC
# -1000.111
# Fisher's C = 6.363 with P-value = 0.042

m2_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_2 <- lme(RRpH ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_3 <- lme(RRSOC ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model2 <- psem(m2_1, m2_2, m2_3)
summary(sem_model2)  ## 
# AIC
# -1006.716
# Fisher's C = 5.436 with P-value = 0.066

m3_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_2 <- lme(RRpH ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_3 <- lme(RRSOC ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model3 <- psem(m3_1, m3_2, m3_3)
summary(sem_model3)  ## 
# AIC
# -989.528
# Fisher's C = 4.723 with P-value = 0.094

m4_1 <- lme(RR ~ RRpH + RRSOC + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_2 <- lme(RRpH ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_3 <- lme(RRSOC ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model4 <- psem(m4_1, m4_2, m4_3)
summary(sem_model4) ## 
# AIC
# -1038.361
# Fisher's C = 6.16 with P-value = 0.046 

m5_1 <- lme(RR ~ RRpH + RRSOC + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_2 <- lme(RRpH ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_3 <- lme(RRSOC ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model5 <- psem(m5_1, m5_2, m5_3)
summary(sem_model5) ## 
# AIC
# -1020.890
# Fisher's C = 5.349 with P-value = 0.069

m6_1 <- lme(RR ~ RRpH + RRSOC +MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_2 <- lme(RRpH ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_3 <- lme(RRSOC ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model6 <- psem(m6_1, m6_2, m6_3)
summary(sem_model6)  ##
# AIC
# -1027.726
# Fisher's C = 4.496 with P-value = 0.106

# Structural Equation Model of sem_model6 
# Call:
#   RR ~ RRpH + RRSOC + MATmean
# RRpH ~ MATmean
# RRSOC ~ MATmean
# AIC
# -1027.726
# ---
#   Tests of directed separation:
#   
#   Independ.Claim Test.Type  DF Crit.Value P.Value 
# RRSOC ~ RRpH + ...      coef 100     1.6329  0.1056 
# --
#   Global goodness-of-fit:
#   Chi-Squared = NA with P-value = NA and on 1 degrees of freedom
# Fisher's C = 4.496 with P-value = 0.106 and on 2 degrees of freedom
# ---
# Coefficients:
#   Response Predictor Estimate Std.Error  DF Crit.Value P.Value Std.Estimate  
#         RR      RRpH   0.1723    0.0740  99     2.3272  0.0220       0.1904 *
#         RR     RRSOC   0.0042    0.0252  99     0.1655  0.8689       0.0142  
#         RR   MATmean   0.0006    0.0011  99     0.5793  0.5637       0.0641  
#       RRpH   MATmean   0.0033    0.0016 116     2.1334  0.0350       0.3101 *
#      RRSOC   MATmean   0.0058    0.0045 115     1.2961  0.1975       0.1765  
# 
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05
# ---
# Individual R-squared:
#   Response method Marginal Conditional
#         RR   none     0.06        0.56
#       RRpH   none     0.07        0.78
#      RRSOC   none     0.03        0.79