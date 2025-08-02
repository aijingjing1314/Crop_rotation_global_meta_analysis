
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)

FungalStructure <- read.csv("FungalStructure.csv", fileEncoding = "latin1")
# Check data
head(FungalStructure)

# 1. The number of Obversation
total_number <- nrow(FungalStructure)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 214

# 2. The number of Study
unique_studyid_number <- length(unique(FungalStructure$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 70


#### 3. Overall effect size
total_effect_model <- rma.mv(yi = RR, 
                             V = Vi, 
                             random = ~ 1 | StudyID,  # StudyID is radom factor
                             data = FungalStructure, 
                             method = "REML")
# The results of Overall effect size
summary(total_effect_model)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -695.6153  1391.2306  1395.2306  1401.9532  1395.2877   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4641  0.6812     70     no  StudyID 
# Test for Heterogeneity:
# Q(df = 213) = 7643.2533, p-val < .0001
# Model Results:
# estimate      se     zval    pval   ci.lb   ci.ub      
#   1.0138  0.0836  12.1261  <.0001  0.8500  1.1777  *** 

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
boot_results1 <- boot(data = FungalStructure, statistic = boot_fun, R = 1000, parallel = "snow", ncpus = numCores, cl = cl)
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
# Estimate for coefficient 1 : 1.01383 
# 95% BCa CI for coefficient 1 : 0.9317944 1.094906 


#### 5. Funnel Plot
simple_model <- rma(yi = RR, 
                    vi = Vi, 
                    data = FungalStructure, 
                    method = "REML")
#### 
funnel(simple_model)
# Output  6 * 6
#### Egger's test
regtest(simple_model)
# Regression Test for Funnel Plot Asymmetry
# Model:     mixed-effects meta-regression model
# Predictor: standard error
# Test for Funnel Plot Asymmetry: z = -0.0399, p = 0.9682
# Limit Estimate (as sei -> 0):   b =  1.0180 (CI: 0.7720, 1.2640)

#  Rosenthal’s Fail-Safe N
# This method estimates how many missing studies with null effect 
# would be needed to make the overall effect non-significant
fsn_rosenthal <- fsn(x = simple_model, type = "Rosenthal")
# Print the FSN result
print(fsn_rosenthal)
# Fail-safe N Calculation Using the General Approach
# Average Effect Size:         1.0136 (with file drawer: 0.0120)
# Amount of Heterogeneity:     0.6124 (with file drawer: 0.6250)
# Observed Significance Level: <.0001 (with file drawer: 0.0500)
# Target Significance Level:   0.05
# Fail-safe N: 16974


#### 8. Subgroup analysis
### 8.1 LegumeNonlegume
FungalStructure_filteredLegumeNonlegume <- subset(FungalStructure, LegumeNonlegume %in% c("Legume to Non-legume", "Non-legume to Legume", "Non-legume to Non-legume"))
#
FungalStructure_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(FungalStructure_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredLegumeNonlegume %>%
  group_by(LegumeNonlegume) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   LegumeNonlegume          Observations Unique_StudyID
#   <fct>                           <int>          <int>
# 1 Legume to Non-legume               28             17
# 2 Non-legume to Legume               81             25
# 3 Non-legume to Non-legume          105             39

overall_model_FungalStructure_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = FungalStructure_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -673.3142  1346.6285  1354.6285  1368.0359  1354.8227   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4373  0.6613     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 211) = 6626.2361, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 201.1388, p-val < .0001
# Model Results:
#                                          estimate      se     zval    pval   ci.lb   ci.ub      
# LegumeNonlegumeLegume to Non-legume        1.0415  0.0962  10.8244  <.0001  0.8529  1.2300  *** 
# LegumeNonlegumeNon-legume to Legume        0.8204  0.0936   8.7607  <.0001  0.6369  1.0040  *** 
# LegumeNonlegumeNon-legume to Non-legume    1.1140  0.0907  12.2861  <.0001  0.9363  1.2918  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredLegumeNonlegume)
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
    #                          "a"                                     "b"                                     "a" 


### 8.2 AMnonAM
FungalStructure_filteredAMnonAM <- subset(FungalStructure, AMnonAM %in% c("AM to AM", "AM to nonAM"))
#
FungalStructure_filteredAMnonAM$AMnonAM <- droplevels(factor(FungalStructure_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredAMnonAM %>%
  group_by(AMnonAM) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   AMnonAM     Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 AM to AM             175             64
# 2 AM to nonAM           33             10
overall_model_FungalStructure_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = FungalStructure_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 208; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -686.9877  1373.9755  1379.9755  1389.9591  1380.0943  
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4635  0.6808     67     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 206) = 7389.0578, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 151.1411, p-val < .0001
# Model Results:
#                     estimate      se     zval    pval   ci.lb   ci.ub      
# AMnonAMAM to AM       1.0463  0.0858  12.1932  <.0001  0.8781  1.2144  *** 
# AMnonAMAM to nonAM    0.8588  0.1071   8.0195  <.0001  0.6489  1.0687  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredAMnonAM)
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
   # AMnonAMAM to AM AMnonAMAM to nonAM 
   #         "a"                "b"           


### 8.3 C3C4
FungalStructure_filteredC3C4 <- subset(FungalStructure, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
FungalStructure_filteredC3C4$C3C4 <- droplevels(factor(FungalStructure_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredC3C4 %>%
  group_by(C3C4) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   C3C4     Observations Unique_StudyID
#   <fct>           <int>          <int>
# 1 C3 to C3          110             38
# 2 C3 to C4           49             30
# 3 C4 to C3           53             18

overall_model_FungalStructure_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = FungalStructure_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 212; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -685.0011  1370.0022  1378.0022  1391.3715  1378.1983   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4507  0.6714     69     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 209) = 6463.5477, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 161.9375, p-val < .0001
# Model Results:
#               estimate      se     zval    pval   ci.lb   ci.ub      
# C3C4C3 to C3    1.0895  0.0892  12.2148  <.0001  0.9147  1.2644  *** 
# C3C4C3 to C4    0.9820  0.0888  11.0584  <.0001  0.8080  1.1561  *** 
# C3C4C4 to C3    0.8761  0.0909   9.6415  <.0001  0.6980  1.0542  ***  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredC3C4)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredC3C4)
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
#     "a"          "a"          "b" 


### 8.4 Annual_Pere
FungalStructure_filteredAnnual_Pere <- subset(FungalStructure, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
FungalStructure_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(FungalStructure_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredAnnual_Pere %>%
  group_by(Annual_Pere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
  # Annual_Pere            Observations Unique_StudyID
#   <fct>                         <int>          <int>
# 1 Annual to Annual                159             52
# 2 Annual to Perennial              10              6
# 3 Perennial to Annual              35             10
# 4 Perennial to Perennial           10              7

overall_model_FungalStructure_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = FungalStructure_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -688.0776  1376.1552  1386.1552  1402.8907  1386.4493   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3658  0.6048     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 210) = 6502.0712, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 200.1832, p-val < .0001
# Model Results:
#                                    estimate      se     zval    pval   ci.lb   ci.ub      
# Annual_PereAnnual to Annual          0.9019  0.0829  10.8838  <.0001  0.7395  1.0643  *** 
# Annual_PereAnnual to Perennial       0.7728  0.1297   5.9578  <.0001  0.5186  1.0271  *** 
# Annual_PerePerennial to Annual       1.4795  0.1503   9.8431  <.0001  1.1849  1.7741  *** 
# Annual_PerePerennial to Perennial    1.4210  0.1626   8.7402  <.0001  1.1024  1.7397  ***   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredAnnual_Pere)
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
      #                  "a"                               "a"                               "b"                               "b" 


### 8.6 PlantStageSubgroup
FungalStructure_filteredPlantStageSubgroup <- subset(FungalStructure, PlantStageSubgroup %in% c("Vegetative stage","Reproductive stage", "Maturity stage","Harvest"))
#
FungalStructure_filteredPlantStageSubgroup$PlantStageSubgroup <- droplevels(factor(FungalStructure_filteredPlantStageSubgroup$PlantStageSubgroup))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredPlantStageSubgroup %>%
  group_by(PlantStageSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   PlantStageSubgroup Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Harvest                      96             31
# 2 Maturity stage               28             12
# 3 Reproductive stage           43             14
# 4 Vegetative stage             16              6
overall_model_FungalStructure_filteredPlantStageSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + PlantStageSubgroup, random = ~ 1 | StudyID, data = FungalStructure_filteredPlantStageSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredPlantStageSubgroup)
# Multivariate Meta-Analysis Model (k = 183; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -553.5612  1107.1223  1117.1223  1133.0593  1117.4692   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5482  0.7404     56     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 179) = 6564.2006, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 110.7108, p-val < .0001
# Model Results:
#                                       estimate      se    zval    pval   ci.lb   ci.ub      
# PlantStageSubgroupHarvest               0.9788  0.1098  8.9127  <.0001  0.7636  1.1941  *** 
# PlantStageSubgroupMaturity stage        1.0186  0.1375  7.4109  <.0001  0.7492  1.2880  *** 
# PlantStageSubgroupReproductive stage    1.2040  0.1313  9.1697  <.0001  0.9467  1.4614  *** 
# PlantStageSubgroupVegetative stage      1.0455  0.1298  8.0554  <.0001  0.7911  1.2999  ***   
#    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredPlantStageSubgroup)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredPlantStageSubgroup)
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
           #                 "a"                                  "a"                                  "a"                                  "a"


### 8.7 Bulk_Rhizosphere
FungalStructure_filteredBulk_Rhizosphere <- subset(FungalStructure, Bulk_Rhizosphere %in% c("Non-Rhizosphere", "Rhizosphere"))
#
FungalStructure_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(FungalStructure_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredBulk_Rhizosphere %>%
  group_by(Bulk_Rhizosphere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Bulk_Rhizosphere Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 Non-Rhizosphere           130             47
# 2 Rhizosphere                84             35

overall_model_FungalStructure_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = FungalStructure_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -694.8034  1389.6068  1395.6068  1405.6765  1395.7221   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4640  0.6812     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 212) = 7493.0555, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 149.8909, p-val < .0001
# Model Results:
#                                  estimate      se     zval    pval   ci.lb   ci.ub      
# Bulk_RhizosphereNon-Rhizosphere    1.0350  0.0845  12.2421  <.0001  0.8693  1.2007  *** 
# Bulk_RhizosphereRhizosphere        0.9836  0.0855  11.5019  <.0001  0.8160  1.1512  *** 
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredBulk_Rhizosphere)
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
FungalStructure_filteredSoil_texture <- subset(FungalStructure, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
FungalStructure_filteredSoil_texture$Soil_texture <- droplevels(factor(FungalStructure_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredSoil_texture %>%
  group_by(Soil_texture) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Soil_texture Observations Unique_StudyID
#   <fct>               <int>          <int>
# 1 Coarse                 48              9
# 2 Fine                   38             11
# 3 Medium                 75             26

overall_model_FungalStructure_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = FungalStructure_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 161; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -422.7736   845.5472   853.5472   865.7976   853.8086   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4626  0.6802     46     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 158) = 4618.8623, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 117.8693, p-val < .0001
# Model Results:
#                     estimate      se    zval    pval   ci.lb   ci.ub      
# Soil_textureCoarse    0.8937  0.2300  3.8850  0.0001  0.4428  1.3445  *** 
# Soil_textureFine      1.3668  0.2134  6.4044  <.0001  0.9485  1.7850  *** 
# Soil_textureMedium    1.0728  0.1365  7.8587  <.0001  0.8053  1.3404  *** 
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredSoil_texture)
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
#             "a"                "a"                "a" 


### 8.9 Tillage
FungalStructure_filteredTillage <- subset(FungalStructure, Tillage %in% c("Tillage", "No_tillage"))
#
FungalStructure_filteredTillage$Tillage <- droplevels(factor(FungalStructure_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredTillage %>%
  group_by(Tillage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Tillage    Observations Unique_StudyID
#   <fct>             <int>          <int>
# 1 No_tillage           40              6
# 2 Tillage              30              9

overall_model_FungalStructure_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = FungalStructure_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredTillage)
# Multivariate Meta-Analysis Model (k = 70; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -111.3783   222.7565   228.7565   235.4150   229.1315   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.9643  0.9820     12     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 68) = 2356.1406, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 15.6407, p-val = 0.0004
# Model Results:
#                    estimate      se    zval    pval   ci.lb   ci.ub      
# TillageNo_tillage    1.1175  0.2862  3.9054  <.0001  0.5567  1.6784  *** 
# TillageTillage       1.0796  0.2855  3.7816  0.0002  0.5201  1.6392  ***     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredTillage)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredTillage)
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
FungalStructure_filteredStraw_retention <- subset(FungalStructure, Straw_retention %in% c("Retention", "No_retention"))
#
FungalStructure_filteredStraw_retention$Straw_retention <- droplevels(factor(FungalStructure_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredStraw_retention %>%
  group_by(Straw_retention) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Straw_retention Observations Unique_StudyID
#   <fct>                  <int>          <int>
# 1 No_retention              21             10
# 2 Retention                 41             13

overall_model_FungalStructure_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = FungalStructure_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 62; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -190.9757   381.9513   387.9513   394.2344   388.3799   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.7189  0.8479     19     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 60) = 1610.5581, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 118.3348, p-val < .0001
# Model Results:
#                              estimate      se    zval    pval   ci.lb   ci.ub      
# Straw_retentionNo_retention    1.3394  0.2010  6.6645  <.0001  0.9455  1.7333  *** 
# Straw_retentionRetention       0.6684  0.1994  3.3516  0.0008  0.2775  1.0593  *** 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredStraw_retention)
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
FungalStructure_filteredRotationcyclesSubgroup <- subset(FungalStructure, RotationcyclesSubgroup %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
FungalStructure_filteredRotationcyclesSubgroup$RotationcyclesSubgroup <- droplevels(factor(FungalStructure_filteredRotationcyclesSubgroup$RotationcyclesSubgroup))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredRotationcyclesSubgroup %>%
  group_by(RotationcyclesSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   RotationcyclesSubgroup Observations Unique_StudyID
#   <fct>                         <int>          <int>
# 1 D1                               71             30
# 2 D1-3                             63             20
# 3 D10                              18              6
# 4 D3-5                             14              7
# 5 D5-10                            48             14

overall_model_FungalStructure_filteredRotationcyclesSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + RotationcyclesSubgroup, random = ~ 1 | StudyID, data = FungalStructure_filteredRotationcyclesSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredRotationcyclesSubgroup)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -632.7766  1265.5532  1277.5532  1297.6072  1277.9690   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.6604  0.8126     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 209) = 6360.0152, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 233.5805, p-val < .0001
# Model Results:
#                              estimate      se     zval    pval   ci.lb   ci.ub      
# RotationcyclesSubgroupD1       1.5380  0.1101  13.9723  <.0001  1.3222  1.7537  *** 
# RotationcyclesSubgroupD1-3     0.6843  0.1112   6.1551  <.0001  0.4664  0.9022  *** 
# RotationcyclesSubgroupD10      0.6664  0.1586   4.2018  <.0001  0.3556  0.9773  *** 
# RotationcyclesSubgroupD3-5     0.7766  0.1282   6.0579  <.0001  0.5253  1.0278  *** 
# RotationcyclesSubgroupD5-10    0.5949  0.1189   5.0048  <.0001  0.3619  0.8279  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredRotationcyclesSubgroup)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredRotationcyclesSubgroup)
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
   #            "a"                        "bc"                        "bc"                         "b"                         "c" 


### 8.13 DurationSubgroup
FungalStructure_filteredDurationSubgroup <- subset(FungalStructure, DurationSubgroup %in% c("D1", "D2", "D3", "D4", "D5", "D6-10", "D11-20", "D20-30", "D30"))
#
FungalStructure_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(FungalStructure_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredDurationSubgroup %>%
  group_by(DurationSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   DurationSubgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D1                          6              4
# 2 D11-20                     33             11
# 3 D2                         39             18
# 4 D20-30                     30              8
# 5 D3                         46             14
# 6 D30                         4              1
# 7 D4                          8              4
# 8 D5                          1              1
# 9 D6-10                      47             16

overall_model_FungalStructure_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = FungalStructure_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -603.1301  1206.2602  1226.2602  1259.4903  1227.3942   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4717  0.6868     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 205) = 5942.5953, p-val < .0001
# Test of Moderators (coefficients 1:9):
# QM(df = 9) = 321.4374, p-val < .0001
# Model Results:
#                         estimate      se     zval    pval    ci.lb   ci.ub      
# DurationSubgroupD1        1.6047  0.3671   4.3709  <.0001   0.8851  2.3242  *** 
# DurationSubgroupD11-20    0.5385  0.1117   4.8197  <.0001   0.3195  0.7575  *** 
# DurationSubgroupD2        1.6394  0.1159  14.1457  <.0001   1.4122  1.8665  *** 
# DurationSubgroupD20-30    0.6354  0.1268   5.0094  <.0001   0.3868  0.8839  *** 
# DurationSubgroupD3        0.7307  0.1097   6.6604  <.0001   0.5157  0.9458  *** 
# DurationSubgroupD30       1.9213  0.6979   2.7528  0.0059   0.5533  3.2892   ** 
# DurationSubgroupD4        1.3969  0.3471   4.0248  <.0001   0.7167  2.0772  *** 
# DurationSubgroupD5        0.6716  0.7043   0.9536  0.3403  -0.7088  2.0519      
# DurationSubgroupD6-10     0.7663  0.1287   5.9546  <.0001   0.5141  1.0185  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredDurationSubgroup)
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
    #                "a"                    "b"                    "a"                   "bc"                   "cd"                 "abcd" 
    # DurationSubgroupD4     DurationSubgroupD5  DurationSubgroupD6-10 
    #               "ad"                 "abcd"                  "bcd" 



### 8.14 SpeciesRichnessSubgroup
FungalStructure_filteredSpeciesRichnessSubgroup <- subset(FungalStructure, SpeciesRichnessSubgroup %in% c("R2", "R3"))
#
FungalStructure_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup <- droplevels(factor(FungalStructure_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredSpeciesRichnessSubgroup %>%
  group_by(SpeciesRichnessSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#  SpeciesRichnessSubgroup Observations Unique_StudyID
#   <fct>                          <int>          <int>
# 1 R2                               164             62
# 2 R3                                49             12

overall_model_FungalStructure_filteredSpeciesRichnessSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichnessSubgroup, random = ~ 1 | StudyID, data = FungalStructure_filteredSpeciesRichnessSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredSpeciesRichnessSubgroup)
# Multivariate Meta-Analysis Model (k = 213; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -679.5864  1359.1728  1365.1728  1375.2283  1365.2887   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4936  0.7025     69     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 211) = 6613.8081, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 166.8825, p-val < .0001
# Model Results:
#                            estimate      se     zval    pval   ci.lb   ci.ub      
# SpeciesRichnessSubgroupR2    1.0791  0.0876  12.3250  <.0001  0.9075  1.2506  *** 
# SpeciesRichnessSubgroupR3    0.6079  0.1136   5.3488  <.0001  0.3851  0.8306  ***   
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredSpeciesRichnessSubgroup)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredSpeciesRichnessSubgroup)
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
#                  "a"                       "b"                    

### 8.15 Primer
FungalStructure_filteredPrimer <- subset(FungalStructure, Primer %in% c("ITS1", "ITS1 + 5.8S + ITS2", "ITS2", "Full length", "18S"))
#
FungalStructure_filteredPrimer$Primer <- droplevels(factor(FungalStructure_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredPrimer %>%
  group_by(Primer) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Primer             Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 18S                           2              2
# 2 Full length                   3              1
# 3 ITS1                        180             56
# 4 ITS1 + 5.8S + ITS2            5              2
# 5 ITS2                         23              8

overall_model_FungalStructure_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = FungalStructure_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -690.6725  1381.3449  1393.3449  1413.3989  1393.7608   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4742  0.6886     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 209) = 5953.7148, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 146.3601, p-val < .0001
# Model Results:
#                           estimate      se     zval    pval    ci.lb   ci.ub      
# Primer18S                   0.6802  0.4938   1.3775  0.1684  -0.2876  1.6480      
# PrimerFull length           0.8074  0.6939   1.1635  0.2446  -0.5526  2.1673      
# PrimerITS1                  1.0187  0.0939  10.8511  <.0001   0.8347  1.2027  *** 
# PrimerITS1 + 5.8S + ITS2    1.6558  0.5002   3.3105  0.0009   0.6755  2.6362  *** 
# PrimerITS2                  0.9343  0.2462   3.7951  0.0001   0.4518  1.4169  *** 
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredPrimer)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredPrimer)
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
               # Primer18S        PrimerFull length               PrimerITS1 PrimerITS1 + 5.8S + ITS2               PrimerITS2 
               #   "a"                      "a"                      "a"                      "a"                      "a" 



### 8.16 Latitude_Subgroup
FungalStructure_filteredLatitude_Subgroup <- subset(FungalStructure, Latitude_Subgroup %in% c("La20", "La20-40", "La40"))
#
FungalStructure_filteredLatitude_Subgroup$Latitude_Subgroup <- droplevels(factor(FungalStructure_filteredLatitude_Subgroup$Latitude_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredLatitude_Subgroup %>%
  group_by(Latitude_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Latitude_Subgroup Observations Unique_StudyID
#   <fct>                    <int>          <int>
# 1 La20                        22              4
# 2 La20-40                     67             33
# 3 La40                       125             33

overall_model_FungalStructure_filteredLatitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Latitude_Subgroup, random = ~ 1 | StudyID, data = FungalStructure_filteredLatitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredLatitude_Subgroup)
# 
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -691.2557  1382.5113  1390.5113  1403.9187  1390.7055   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4523  0.6725     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 211) = 6695.1898, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 154.5822, p-val < .0001
# Model Results:
#                           estimate      se    zval    pval   ci.lb   ci.ub      
# Latitude_SubgroupLa20       1.5420  0.3463  4.4527  <.0001  0.8633  2.2208  *** 
# Latitude_SubgroupLa20-40    1.0862  0.1215  8.9415  <.0001  0.8481  1.3242  *** 
# Latitude_SubgroupLa40       0.8817  0.1191  7.4030  <.0001  0.6483  1.1151  ***   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredLatitude_Subgroup)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredLatitude_Subgroup)
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
FungalStructure_filteredLongitude_Subgroup <- subset(FungalStructure, Longitude_Subgroup %in% c("Lo-180-0", "Lo-0-180"))
#
FungalStructure_filteredLongitude_Subgroup$Longitude_Subgroup <- droplevels(factor(FungalStructure_filteredLongitude_Subgroup$Longitude_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredLongitude_Subgroup %>%
  group_by(Longitude_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Longitude_Subgroup Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Lo-0-180                    185             66
# 2 Lo-180-0                     29              4

overall_model_FungalStructure_filteredLongitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Subgroup, random = ~ 1 | StudyID, data = FungalStructure_filteredLongitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredLongitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -691.9323  1383.8647  1389.8647  1399.9344  1389.9801   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4398  0.6632     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 212) = 6048.5226, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 159.2515, p-val < .0001
# Model Results:
#                             estimate      se     zval    pval    ci.lb   ci.ub      
# Longitude_SubgroupLo-0-180    1.0575  0.0840  12.5816  <.0001   0.8927  1.2222  *** 
# Longitude_SubgroupLo-180-0    0.3255  0.3333   0.9769  0.3286  -0.3276  0.9787      
     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredLongitude_Subgroup)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredLongitude_Subgroup)
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
FungalStructure_filteredMAPmean_Subgroup <- subset(FungalStructure, MAPmean_Subgroup %in% c("MAP600", "MAP600-1200", "MAP1200"))
#
FungalStructure_filteredMAPmean_Subgroup$MAPmean_Subgroup <- droplevels(factor(FungalStructure_filteredMAPmean_Subgroup$MAPmean_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredMAPmean_Subgroup %>%
  group_by(MAPmean_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MAPmean_Subgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 MAP1200                    52             19
# 2 MAP600                    123             37
# 3 MAP600-1200                39             16
overall_model_FungalStructure_filteredMAPmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean_Subgroup, random = ~ 1 | StudyID, data = FungalStructure_filteredMAPmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredMAPmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -690.4533  1380.9065  1388.9065  1402.3139  1389.1007   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4350  0.6595     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 211) = 5621.6632, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 163.9980, p-val < .0001
# Model Results:
#                              estimate      se    zval    pval   ci.lb   ci.ub      
# MAPmean_SubgroupMAP1200        1.3345  0.1524  8.7568  <.0001  1.0358  1.6332  *** 
# MAPmean_SubgroupMAP600         0.9340  0.0982  9.5081  <.0001  0.7415  1.1266  *** 
# MAPmean_SubgroupMAP600-1200    0.8163  0.1200  6.8009  <.0001  0.5811  1.0516  *** 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredMAPmean_Subgroup)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredMAPmean_Subgroup)
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
#                  "a"                         "b"                         "b" 



### 8.19 MATmean_Subgroup
FungalStructure_filteredMATmean_Subgroup <- subset(FungalStructure, MATmean_Subgroup %in% c("MAT8", "MAT8-15", "MAT15"))
#
FungalStructure_filteredMATmean_Subgroup$MATmean_Subgroup <- droplevels(factor(FungalStructure_filteredMATmean_Subgroup$MATmean_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredMATmean_Subgroup %>%
  group_by(MATmean_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MATmean_Subgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 MAT15                      65             25
# 2 MAT8                      122             35
# 3 MAT8-15                    27             11

overall_model_FungalStructure_filteredMATmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Subgroup, random = ~ 1 | StudyID, data = FungalStructure_filteredMATmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredMATmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -693.8996  1387.7992  1395.7992  1409.2066  1395.9934   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4671  0.6834     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 211) = 6651.3726, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 146.7962, p-val < .0001
# Model Results:
#                          estimate      se    zval    pval   ci.lb   ci.ub      
# MATmean_SubgroupMAT15      1.1040  0.1425  7.7485  <.0001  0.8248  1.3833  *** 
# MATmean_SubgroupMAT8       0.9603  0.1069  8.9820  <.0001  0.7508  1.1699  *** 
# MATmean_SubgroupMAT8-15    0.9856  0.1362  7.2365  <.0001  0.7186  1.2525  *** 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredMATmean_Subgroup)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredMATmean_Subgroup)
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
#                 "a"                     "a"                     "a" 



### 8.20 DurationSubgroup
FungalStructure_filteredDurationSubgroup <- subset(FungalStructure, DurationSubgroup %in% c("D5", "D5-10", "D10-20", "D20-30"))
#
FungalStructure_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(FungalStructure_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- FungalStructure_filteredDurationSubgroup %>%
  group_by(DurationSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   DurationSubgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D10-20                     33             11
# 2 D20-30                     30              8
# 3 D5                        100             38
# 4 D5-10                      47             16

overall_model_FungalStructure_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = FungalStructure_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalStructure_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 210; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -684.2257  1368.4514  1378.4514  1395.0908  1378.7514   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4555  0.6749     69     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 206) = 7557.7328, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 156.2924, p-val < .0001
# Model Results:
#                         estimate      se     zval    pval   ci.lb   ci.ub      
# DurationSubgroupD10-20    0.7881  0.1057   7.4564  <.0001  0.5809  0.9952  *** 
# DurationSubgroupD20-30    0.8372  0.1229   6.8104  <.0001  0.5963  1.0781  *** 
# DurationSubgroupD5        1.1196  0.0958  11.6853  <.0001  0.9318  1.3073  *** 
# DurationSubgroupD5-10     0.9099  0.1259   7.2280  <.0001  0.6632  1.1567  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalStructure_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_FungalStructure_filteredDurationSubgroup)
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
# DurationSubgroupD10-20 DurationSubgroupD20-30     DurationSubgroupD5  DurationSubgroupD5-10 
#                    "a"                    "a"                    "b"                   "ab" 






#### 9. Linear Mixed Effect Model
# 
FungalStructure$Wr <- 1 / FungalStructure$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalStructure)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalStructure)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalStructure)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalStructure)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalStructure)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalStructure)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalStructure)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalStructure)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(Species_Richness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = FungalStructure)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = FungalStructure)

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
# Model1   Model1 483.8349 504.0307 -235.9174 0.001183001     0.04248893
# Model2   Model2 480.4356 500.6314 -234.2178 0.003320480     0.04394120
# Model3   Model3 482.1421 502.3379 -235.0710 0.001940473     0.04421435
# Model4   Model4 483.6731 503.8690 -235.8366 0.001279084     0.04256305
# Model5   Model5 480.1687 500.3645 -234.0843 0.004202636     0.04595296
# Model6   Model6 480.0935 500.2894 -234.0468 0.004250692     0.04599940
# Model7   Model7 480.5203 500.7162 -234.2602 0.003269015     0.04385552
# Model8   Model8 481.8972 502.0931 -234.9486 0.002073890     0.04433390
# Model9   Model9 494.5959 528.2556 -237.2979 0.004021275     0.04391391
# Model10 Model10 489.9807 523.6405 -234.9904 0.006621554     0.04703530

##### Model 6 is the best model
summary(Model6)
# Number of obs: 214, groups:  StudyID, 70
anova(Model6) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)  
# scale(Rotation_cycles)       12.301  12.301     1 104.19  1.3521 0.2476  
# scale(log(Species_Richness))  5.477   5.477     1 152.45  0.6021 0.4390  
# scale(log(Duration))         33.544  33.544     1 137.51  3.6872 0.0569 .
#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + scale(pHCK) + (1 | StudyID), weights = Wr, data = FungalStructure)
summary(ModelpH)
# Number of obs: 118, groups:  StudyID, 48
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)        0.7768  0.7768     1  64.466  0.0841 0.7728
# scale(log(Species_Richness)) 11.7302 11.7302     1 111.963  1.2692 0.2623
# scale(log(Duration))         24.6680 24.6680     1  59.855  2.6690 0.1076
# scale(pHCK)                   0.0805  0.0805     1  61.848  0.0087 0.9259

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = FungalStructure)
summary(ModelSOC)
# Number of obs: 132, groups:  StudyID, 49
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)        1.3694  1.3694     1  63.710  0.1277 0.7220
# scale(log(Species_Richness))  4.5275  4.5275     1 104.377  0.4222 0.5173
# scale(log(Duration))         17.7716 17.7716     1  60.326  1.6571 0.2029
# scale(SOCCK)                  0.0357  0.0357     1  64.764  0.0033 0.9542

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + scale(TNCK) + (1 | StudyID), weights = Wr, data = FungalStructure)
summary(ModelTN)
# Number of obs: 84, groups:  StudyID, 35
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)        0.4713  0.4713     1 46.728  0.0627 0.8033
# scale(log(Species_Richness))  1.3108  1.3108     1 67.264  0.1745 0.6775
# scale(log(Duration))         11.4528 11.4528     1 37.429  1.5244 0.2246
# scale(TNCK)                   5.6247  5.6247     1 37.327  0.7487 0.3924


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = FungalStructure)
summary(ModelNO3)
# Number of obs: 42, groups:  StudyID, 17
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
# Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)       0.1331  0.1331     1 25.291  0.0155 0.9018
# scale(log(Species_Richness)) 0.0237  0.0237     1 26.810  0.0028 0.9585
# scale(log(Duration))         3.5611  3.5611     1 24.403  0.4152 0.5253
# scale(NO3CK)                 1.7120  1.7120     1 16.892  0.1996 0.6607

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = FungalStructure)
summary(ModelNH4)
# Number of obs: 40, groups:  StudyID, 16
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)        0.1024  0.1024     1 22.847  0.0117 0.9148
# scale(log(Species_Richness))  0.0064  0.0064     1 25.995  0.0007 0.9787
# scale(log(Duration))          3.5645  3.5645     1 22.251  0.4077 0.5296
# scale(NH4CK)                 10.6123 10.6123     1 12.295  1.2139 0.2917

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + scale(APCK) + (1 | StudyID), weights = Wr, data = FungalStructure)
summary(ModelAP)
# Number of obs: 98, groups:  StudyID, 42
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles)       18.846  18.846     1 53.928  1.6505 0.20438  
# scale(log(Species_Richness))  0.244   0.244     1 64.039  0.0214 0.88418  
# scale(log(Duration))         70.292  70.292     1 52.298  6.1561 0.01635 *
# scale(APCK)                   0.133   0.133     1 69.811  0.0117 0.91426  

#### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + scale(AKCK) + (1 | StudyID), weights = Wr, data = FungalStructure)
summary(ModelAK)
# Number of obs: 76, groups:  StudyID, 35
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles)       34.303  34.303     1 35.865  2.4689 0.12490  
# scale(log(Species_Richness))  0.014   0.014     1 31.507  0.0010 0.97493  
# scale(log(Duration))         53.057  53.057     1 42.958  3.8188 0.05721 .
# scale(AKCK)                   6.362   6.362     1 45.255  0.4579 0.50204  

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + scale(ANCK) + (1 | StudyID), weights = Wr, data = FungalStructure)
summary(ModelAN)
# Number of obs: 64, groups:  StudyID, 22
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles)       51.156  51.156     1 20.466  5.3470 0.03127 *
# scale(log(Species_Richness)) 15.608  15.608     1 42.312  1.6314 0.20847  
# scale(log(Duration))         42.455  42.455     1 21.078  4.4375 0.04731 *
# scale(ANCK)                  18.189  18.189     1 20.917  1.9012 0.18251 

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + scale(Latitude) + (1 | StudyID), weights = Wr, data = FungalStructure)
summary(ModelLatitude)
# Number of obs: 214, groups:  StudyID, 70
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)        9.1794  9.1794     1 105.306  1.0058 0.3182
# scale(log(Species_Richness))  5.1862  5.1862     1 150.979  0.5683 0.4521
# scale(log(Duration))         24.3610 24.3610     1 134.968  2.6693 0.1046
# scale(Latitude)               1.4781  1.4781     1  71.038  0.1620 0.6886

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + scale(Longitude) + (1 | StudyID), weights = Wr, data = FungalStructure)
summary(ModelLongitude)
# Number of obs: 214, groups:  StudyID, 70
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles)       16.113  16.113     1 100.080  1.7608 0.18754  
# scale(log(Species_Richness))  1.537   1.537     1 149.778  0.1680 0.68249  
# scale(log(Duration))         34.204  34.204     1 131.317  3.7379 0.05534 .
# scale(Longitude)             39.647  39.647     1  51.613  4.3328 0.04236 *
### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = FungalStructure)
summary(ModelMAPmean)
# Number of obs: 214, groups:  StudyID, 70
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles)        9.5568  9.5568     1 111.738  1.0482 0.30814  
# scale(log(Species_Richness))  5.4614  5.4614     1 150.910  0.5990 0.44018  
# scale(log(Duration))         25.1343 25.1343     1 142.297  2.7567 0.09905 .
# scale(MAPmean)                0.1997  0.1997     1  81.852  0.0219 0.88271  

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + scale(MATmean) + (1 | StudyID), weights = Wr, data = FungalStructure)
summary(ModelMATmean)
# Number of obs: 214, groups:  StudyID, 70
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles)       10.3891 10.3891     1 105.53  1.1406 0.28796  
# scale(log(Species_Richness))  5.4269  5.4269     1 152.13  0.5958 0.44137  
# scale(log(Duration))         28.1444 28.1444     1 133.79  3.0900 0.08106 .
# scale(MATmean)                0.1006  0.1006     1  67.91  0.0110 0.91663  



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## Species_Richness
sum(!is.na(FungalStructure$Species_Richness)) ## n = 214
p1 <- ggplot(FungalStructure, aes(y=RR, x=Species_Richness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="Species_Richness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Species_Richness" , y="lnFungalStructure214")
p1
pdf("Species_Richness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(FungalStructure$Duration)) ## n = 214
p2 <- ggplot(FungalStructure, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnFungalStructure214")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(FungalStructure$Rotation_cycles)) ## n = 214
p3 <- ggplot(FungalStructure, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnFungalStructure214")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(FungalStructure$Latitude)) ## n = 214
p5 <- ggplot(FungalStructure, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnFungalStructure214")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(FungalStructure$Longitude)) ## n = 214
p6 <- ggplot(FungalStructure, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnFungalStructure214")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(FungalStructure$MAPmean)) ## n = 214
p7 <- ggplot(FungalStructure, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnFungalStructure214")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(FungalStructure$MATmean)) ## n = 214
p8 <- ggplot(FungalStructure, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnFungalStructure214")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(FungalStructure$pHCK)) ## n = 118
p9 <- ggplot(FungalStructure, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnFungalStructure118")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(FungalStructure$SOCCK)) ## n = 132
p10 <- ggplot(FungalStructure, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnFungalStructure132")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(FungalStructure$TNCK)) ## n = 84
p11 <- ggplot(FungalStructure, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnFungalStructure84")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(FungalStructure$NO3CK)) ## n = 42
p12 <- ggplot(FungalStructure, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnFungalStructure42")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(FungalStructure$NH4CK)) ## n = 40
p13<- ggplot(FungalStructure, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnFungalStructure40")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(FungalStructure$APCK)) ## n = 98
p14 <- ggplot(FungalStructure, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnFungalStructure98")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(FungalStructure$AKCK)) ## n = 76
p15 <- ggplot(FungalStructure, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnFungalStructure76")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(FungalStructure$ANCK)) ## n = 64
p16 <- ggplot(FungalStructure, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnFungalStructure64")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(FungalStructure$RRpH)) ## n = 118
p17 <- ggplot(FungalStructure, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="RRpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRpH" , y="RR118")
p17
pdf("RRpH.pdf",width=8,height=8)
p17
dev.off() 

## RRSOC
sum(!is.na(FungalStructure$RRSOC)) ## n = 132
p18 <- ggplot(FungalStructure, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="RRSOC")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRSOC" , y="RR132")
p18
pdf("RRSOC.pdf",width=8,height=8)
p18
dev.off() 

## RRTN
sum(!is.na(FungalStructure$RRTN)) ## n = 84
p19 <- ggplot(FungalStructure, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="RRTN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRTN" , y="RR84")
p19
pdf("RRTN.pdf",width=8,height=8)
p19
dev.off() 

## RRNO3
sum(!is.na(FungalStructure$RRNO3)) ## n = 42
p20 <- ggplot(FungalStructure, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="RRNO3")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNO3" , y="RR42")
p20
pdf("RRNO3.pdf",width=8,height=8)
p20
dev.off() 

## RRNH4
sum(!is.na(FungalStructure$RRNH4)) ## n = 40
p21 <- ggplot(FungalStructure, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="RRNH4")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNH4" , y="RR40")
p21
pdf("RRNH4.pdf",width=8,height=8)
p21
dev.off() 

## RRAP
sum(!is.na(FungalStructure$RRAP)) ## n = 98
p22 <- ggplot(FungalStructure, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="RRAP")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAP" , y="RR98")
p22
pdf("RRAP.pdf",width=8,height=8)
p22
dev.off() 

## RRAK
sum(!is.na(FungalStructure$RRAK)) ## n = 76
p23 <- ggplot(FungalStructure, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="RRAK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAK" , y="RR76")
p23
pdf("RRAK.pdf",width=8,height=8)
p23
dev.off() 

## RRAN
sum(!is.na(FungalStructure$RRAN)) ## n = 64
p24 <- ggplot(FungalStructure, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalStructure", x="RRAN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAN" , y="RR64")
p24
pdf("RRAN.pdf",width=8,height=8)
p24
dev.off() 

## RRYield
sum(!is.na(FungalStructure$RRYield)) ## n = 28
p25 <- ggplot(FungalStructure, aes(x=RR, y=RRYield)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm, color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") +
  theme_bw() +
  theme(text = element_text(family = "serif", size=20)) +
  geom_vline(aes(xintercept=0), colour="black", linewidth=0.5, linetype="dashed") +
  labs(x="RR", y="RRYield28") +
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
  scale_x_continuous(limits=c(-0.1, 3.3), expand=c(0, 0)) + 
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
                            data = FungalStructure, 
                            na.action = na.roughfix, 
                            importance = TRUE, 
                            ntree = 500)
# Check the importance of variables and p-values
importance(rf_model_perm3)
#                    %IncMSE %IncMSE.pval IncNodePurity IncNodePurity.pval
# Longitude        27.465344   0.00990099     19.721933         0.00990099
# MATmean          26.994259   0.00990099     20.178806         0.00990099
# Latitude         23.747210   0.00990099     13.845298         0.00990099
# Duration         23.004974   0.00990099     12.983891         0.00990099
# SOCCK            21.583007   0.00990099     13.523867         0.79207921
# MAPmean          20.792048   0.00990099     14.871108         0.00990099
# Rotation_cycles  20.660738   0.00990099      7.894426         0.21782178
# pHCK             15.112707   0.00990099     10.182979         0.93069307
# Species_Richness  4.631136   0.17821782      1.201202         1.00000000
# 

####################################### Trials sorted by effect size
library(ggplot2)
library(dplyr)
FungalStructure <- read.csv("FungalStructure.csv", fileEncoding = "latin1")
# 计算95% CI + 显著性分类
df_plot <- FungalStructure %>%
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
#       59      155 
## 8*8


################################################ piecewiseSEM
sem_data <- read.csv("FungalStructure.csv", fileEncoding = "latin1")
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
# -108.627
# Fisher's C = 1.107 with P-value = 0.575

m1_1 <- lme(RR ~ RRpH + RRSOC + Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_2 <- lme(RRpH ~  Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_3 <- lme(RRSOC ~ Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model1 <- psem(m1_1, m1_2, m1_3)
summary(sem_model1)  ## 
# AIC
# -134.562
# Fisher's C = 2.465 with P-value = 0.292

m2_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_2 <- lme(RRpH ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_3 <- lme(RRSOC ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model2 <- psem(m2_1, m2_2, m2_3)
summary(sem_model2)  ## 
# AIC
# -137.279
# Fisher's C = 0.803 with P-value = 0.669

m3_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_2 <- lme(RRpH ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_3 <- lme(RRSOC ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model3 <- psem(m3_1, m3_2, m3_3)
summary(sem_model3)  ## 
# AIC
# -114.193
# Fisher's C = 0.366 with P-value = 0.833

m4_1 <- lme(RR ~ RRpH + RRSOC + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_2 <- lme(RRpH ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_3 <- lme(RRSOC ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model4 <- psem(m4_1, m4_2, m4_3)
summary(sem_model4) ## 
# AIC
# -163.352
# Fisher's C = 2.098 with P-value = 0.35

# Structural Equation Model of sem_model4 
# Call:
#   RR ~ RRpH + RRSOC + Species_Richness
# RRpH ~ Species_Richness
# RRSOC ~ Species_Richness
# AIC
# -163.352
# ---
#   Tests of directed separation:
#   Independ.Claim Test.Type DF Crit.Value P.Value 
# RRSOC ~ RRpH + ...      coef 56     0.9418  0.3503 
# --
#   Global goodness-of-fit:
#   
#   Chi-Squared = NA with P-value = NA and on 1 degrees of freedom
# Fisher's C = 2.098 with P-value = 0.35 and on 2 degrees of freedom
# 
# ---
# Coefficients:
# 
#   Response        Predictor Estimate Std.Error DF Crit.Value P.Value Std.Estimate    
#         RR             RRpH   2.0290    1.6736 55     1.2124  0.2306       0.1791    
#         RR            RRSOC  -0.4154    0.4929 55    -0.8427  0.4030      -0.1107    
#         RR Species_Richness  -0.1591    0.2681 55    -0.5934  0.5553      -0.0866    
#       RRpH Species_Richness  -0.0527    0.0149 58    -3.5446  0.0008      -0.3248 ***
#      RRSOC Species_Richness  -0.0518    0.0458 71    -1.1304  0.2621      -0.1057    
# 
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05
# 
# ---
# Individual R-squared:
# 
#   Response method Marginal Conditional
#         RR   none     0.04        0.68
#       RRpH   none     0.10        0.84
#      RRSOC   none     0.01        0.64

m5_1 <- lme(RR ~ RRpH + RRSOC + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_2 <- lme(RRpH ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_3 <- lme(RRSOC ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model5 <- psem(m5_1, m5_2, m5_3)
summary(sem_model5) ## 
# AIC
# -138.203
# Fisher's C = 1.553 with P-value = 0.46

m6_1 <- lme(RR ~ RRpH + RRSOC +MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_2 <- lme(RRpH ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_3 <- lme(RRSOC ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model6 <- psem(m6_1, m6_2, m6_3)
summary(sem_model6)  ##
# AIC
# -143.527
# Fisher's C = 0.214 with P-value = 0.898

