
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)

FungalBeta <- read.csv("FungalBeta.csv", fileEncoding = "latin1")
# Check data
head(FungalBeta)

# 1. The number of Obversation
total_number <- nrow(FungalBeta)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 214

# 2. The number of Study
unique_studyid_number <- length(unique(FungalBeta$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 70


#### 3. Overall effect size
total_effect_model <- rma.mv(yi = RR, 
                             V = Vi, 
                             random = ~ 1 | StudyID,  # StudyID is radom factor
                             data = FungalBeta, 
                             method = "REML")
# The results of Overall effect size
summary(total_effect_model)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -535.3691  1070.7383  1074.7383  1081.4609  1074.7954   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3737  0.6113     70     no  StudyID 
# Test for Heterogeneity:
# Q(df = 213) = 1911.5936, p-val < .0001
# Model Results:
# estimate      se    zval    pval   ci.lb   ci.ub    
#   0.1652  0.0776  2.1296  0.0332  0.0132  0.3172  * 

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
boot_results1 <- boot(data = FungalBeta, statistic = boot_fun, R = 1000, parallel = "snow", ncpus = numCores, cl = cl)
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
# Estimate for coefficient 1 : 0.1651901 
# 95% BCa CI for coefficient 1 : 0.05470898 0.269176


#### 5. Funnel Plot
simple_model <- rma(yi = RR, 
                    vi = Vi, 
                    data = FungalBeta, 
                    method = "REML")
#### 
funnel(simple_model)
# Output  6 * 6
#### Egger's test
regtest(simple_model)
# Regression Test for Funnel Plot Asymmetry
# Model:     mixed-effects meta-regression model
# Predictor: standard error
# Test for Funnel Plot Asymmetry: z = 0.3631, p = 0.7165
# Limit Estimate (as sei -> 0):   b = 0.0966 (CI: -0.1756, 0.3689)

#  Rosenthal’s Fail-Safe N
# This method estimates how many missing studies with null effect 
# would be needed to make the overall effect non-significant
fsn_rosenthal <- fsn(x = simple_model, type = "Rosenthal")
# Print the FSN result
print(fsn_rosenthal)
# Fail-safe N Calculation Using the General Approach
# Average Effect Size:         0.1431 (with file drawer: 0.0776)
# Amount of Heterogeneity:     0.5224 (with file drawer: 0.5260)
# Observed Significance Level: 0.0077 (with file drawer: 0.0505)
# Target Significance Level:   0.05
# Fail-safe N: 164


#### 8. Subgroup analysis
### 8.1 LegumeNonlegume
FungalBeta_filteredLegumeNonlegume <- subset(FungalBeta, LegumeNonlegume %in% c("Legume to Non-legume", "Non-legume to Legume", "Non-legume to Non-legume"))
#
FungalBeta_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(FungalBeta_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredLegumeNonlegume %>%
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

overall_model_FungalBeta_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = FungalBeta_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -525.1848  1050.3696  1058.3696  1071.7771  1058.5638  
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4338  0.6586     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 211) = 1879.1779, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 24.3271, p-val < .0001
# Model Results:
#                                          estimate      se     zval    pval    ci.lb   ci.ub     
# LegumeNonlegumeLegume to Non-legume        0.1074  0.1101   0.9755  0.3293  -0.1084  0.3231     
# LegumeNonlegumeNon-legume to Legume       -0.0738  0.1035  -0.7131  0.4758  -0.2766  0.1290     
# LegumeNonlegumeNon-legume to Non-legume    0.3240  0.0986   3.2864  0.0010   0.1308  0.5172  **     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredLegumeNonlegume)
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
FungalBeta_filteredAMnonAM <- subset(FungalBeta, AMnonAM %in% c("AM to AM", "AM to nonAM"))
#
FungalBeta_filteredAMnonAM$AMnonAM <- droplevels(factor(FungalBeta_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredAMnonAM %>%
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

overall_model_FungalBeta_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = FungalBeta_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 208; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -516.3454  1032.6908  1038.6908  1048.6744  1038.8096   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3933  0.6271     67     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 206) = 1857.6966, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 18.4425, p-val < .0001
# Model Results:
#                     estimate      se    zval    pval    ci.lb   ci.ub      
# AMnonAMAM to AM       0.1319  0.0819  1.6100  0.1074  -0.0287  0.2924      
# AMnonAMAM to nonAM    0.5080  0.1213  4.1879  <.0001   0.2702  0.7457  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredAMnonAM)
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
   #          "a"                "b"            


### 8.3 C3C4
FungalBeta_filteredC3C4 <- subset(FungalBeta, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
FungalBeta_filteredC3C4$C3C4 <- droplevels(factor(FungalBeta_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredC3C4 %>%
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

overall_model_FungalBeta_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = FungalBeta_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 212; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -518.8250  1037.6499  1045.6499  1059.0193  1045.8460   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4288  0.6548     69     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 209) = 1838.5333, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 14.6895, p-val = 0.0021
# Model Results:
#               estimate      se     zval    pval    ci.lb   ci.ub     
# C3C4C3 to C3    0.3107  0.0961   3.2329  0.0012   0.1223  0.4991  ** 
# C3C4C3 to C4    0.0639  0.0961   0.6652  0.5059  -0.1245  0.2523     
# C3C4C4 to C3   -0.0053  0.0996  -0.0529  0.9578  -0.2005  0.1900      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredC3C4)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredC3C4)
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
#     "a"          "b"          "b" 


### 8.4 Annual_Pere
FungalBeta_filteredAnnual_Pere <- subset(FungalBeta, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
FungalBeta_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(FungalBeta_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredAnnual_Pere %>%
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

overall_model_FungalBeta_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = FungalBeta_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -527.3667  1054.7333  1064.7333  1081.4688  1065.0274   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3914  0.6256     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 210) = 1870.4379, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 17.7828, p-val = 0.0014
# Model Results:
#                                    estimate      se    zval    pval    ci.lb   ci.ub      
# Annual_PereAnnual to Annual          0.1595  0.0891  1.7899  0.0735  -0.0152  0.3341    . 
# Annual_PereAnnual to Perennial       0.6530  0.1549  4.2147  <.0001   0.3493  0.9567  *** 
# Annual_PerePerennial to Annual       0.0837  0.1758  0.4761  0.6340  -0.2608  0.4282      
# Annual_PerePerennial to Perennial    0.0318  0.2017  0.1578  0.8746  -0.3634  0.4271      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredAnnual_Pere)
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
      #                  "a"                               "b"                               "a"                               "a" 


### 8.6 PlantStageSubgroup
FungalBeta_filteredPlantStageSubgroup <- subset(FungalBeta, PlantStageSubgroup %in% c("Vegetative stage","Reproductive stage", "Maturity stage","Harvest"))
#
FungalBeta_filteredPlantStageSubgroup$PlantStageSubgroup <- droplevels(factor(FungalBeta_filteredPlantStageSubgroup$PlantStageSubgroup))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredPlantStageSubgroup %>%
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
overall_model_FungalBeta_filteredPlantStageSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + PlantStageSubgroup, random = ~ 1 | StudyID, data = FungalBeta_filteredPlantStageSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredPlantStageSubgroup)
# Multivariate Meta-Analysis Model (k = 183; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -390.5199   781.0397   791.0397   806.9766   791.3865   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3264  0.5713     56     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 179) = 1233.3059, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 16.3886, p-val = 0.0025
# Model Results:
#                                       estimate      se     zval    pval    ci.lb   ci.ub      
# PlantStageSubgroupHarvest               0.3407  0.0965   3.5302  0.0004   0.1515  0.5298  *** 
# PlantStageSubgroupMaturity stage       -0.0923  0.1406  -0.6564  0.5116  -0.3678  0.1833      
# PlantStageSubgroupReproductive stage   -0.0070  0.1255  -0.0556  0.9557  -0.2529  0.2390      
# PlantStageSubgroupVegetative stage      0.0188  0.1451   0.1296  0.8969  -0.2656  0.3032    
#    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredPlantStageSubgroup)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredPlantStageSubgroup)
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
           #                  "a"                                  "b"                                  "b"                                  "b"


### 8.7 Bulk_Rhizosphere
FungalBeta_filteredBulk_Rhizosphere <- subset(FungalBeta, Bulk_Rhizosphere %in% c("Non-Rhizosphere", "Rhizosphere"))
#
FungalBeta_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(FungalBeta_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredBulk_Rhizosphere %>%
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

overall_model_FungalBeta_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = FungalBeta_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -533.5857  1067.1715  1073.1715  1083.2412  1073.2869   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3776  0.6145     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 212) = 1907.0080, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 8.4401, p-val = 0.0147
# Model Results:
#                                  estimate      se    zval    pval    ci.lb   ci.ub     
# Bulk_RhizosphereNon-Rhizosphere    0.1275  0.0802  1.5900  0.1118  -0.0297  0.2848     
# Bulk_RhizosphereRhizosphere        0.2196  0.0826  2.6590  0.0078   0.0577  0.3815  ** 
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredBulk_Rhizosphere)
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
FungalBeta_filteredSoil_texture <- subset(FungalBeta, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
FungalBeta_filteredSoil_texture$Soil_texture <- droplevels(factor(FungalBeta_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredSoil_texture %>%
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

overall_model_FungalBeta_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = FungalBeta_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 161; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -383.5829   767.1658   775.1658   787.4162   775.4273   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3378  0.5812     46     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 158) = 1364.6141, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 5.3555, p-val = 0.1475
# Model Results:
#                     estimate      se    zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse    0.0505  0.2007  0.2516  0.8014  -0.3429  0.4439    
# Soil_textureFine      0.2949  0.1908  1.5452  0.1223  -0.0791  0.6689    
# Soil_textureMedium    0.2074  0.1217  1.7042  0.0883  -0.0311  0.4459  .   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredSoil_texture)
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
FungalBeta_filteredTillage <- subset(FungalBeta, Tillage %in% c("Tillage", "No_tillage"))
#
FungalBeta_filteredTillage$Tillage <- droplevels(factor(FungalBeta_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredTillage %>%
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

overall_model_FungalBeta_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = FungalBeta_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredTillage)
# ultivariate Meta-Analysis Model (k = 70; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -153.4467   306.8935   312.8935   319.5520   313.2685   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1333  0.3651     12     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 68) = 404.6188, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.7348, p-val = 0.6925
# Model Results:
#                    estimate      se    zval    pval    ci.lb   ci.ub    
# TillageNo_tillage    0.0793  0.1214  0.6532  0.5136  -0.1587  0.3173    
# TillageTillage       0.1001  0.1185  0.8448  0.3982  -0.1321  0.3324     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredTillage)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredTillage)
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
FungalBeta_filteredStraw_retention <- subset(FungalBeta, Straw_retention %in% c("Retention", "No_retention"))
#
FungalBeta_filteredStraw_retention$Straw_retention <- droplevels(factor(FungalBeta_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredStraw_retention %>%
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

overall_model_FungalBeta_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = FungalBeta_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 62; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -142.2088   284.4177   290.4177   296.7007   290.8463   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1804  0.4248     19     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 60) = 421.0801, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 60.4607, p-val < .0001
# Model Results:
#                              estimate      se     zval    pval    ci.lb    ci.ub      
# Straw_retentionNo_retention    0.4771  0.1212   3.9353  <.0001   0.2395   0.7147  *** 
# Straw_retentionRetention      -0.2332  0.1153  -2.0233  0.0430  -0.4591  -0.0073    * 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredStraw_retention)
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
FungalBeta_filteredRotationcyclesSubgroup <- subset(FungalBeta, RotationcyclesSubgroup %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
FungalBeta_filteredRotationcyclesSubgroup$RotationcyclesSubgroup <- droplevels(factor(FungalBeta_filteredRotationcyclesSubgroup$RotationcyclesSubgroup))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredRotationcyclesSubgroup %>%
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

overall_model_FungalBeta_filteredRotationcyclesSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + RotationcyclesSubgroup, random = ~ 1 | StudyID, data = FungalBeta_filteredRotationcyclesSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredRotationcyclesSubgroup)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -511.1223  1022.2446  1034.2446  1054.2986  1034.6604  
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4331  0.6581     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 209) = 1889.8568, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 50.0839, p-val < .0001
# Model Results:
#                              estimate      se     zval    pval    ci.lb   ci.ub      
# RotationcyclesSubgroupD1      -0.1103  0.1024  -1.0778  0.2811  -0.3110  0.0903      
# RotationcyclesSubgroupD1-3     0.1079  0.1072   1.0060  0.3144  -0.1023  0.3180      
# RotationcyclesSubgroupD10      0.3075  0.1780   1.7271  0.0841  -0.0415  0.6564    . 
# RotationcyclesSubgroupD3-5     0.3228  0.1380   2.3382  0.0194   0.0522  0.5933    * 
# RotationcyclesSubgroupD5-10    0.7550  0.1205   6.2656  <.0001   0.5188  0.9912  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredRotationcyclesSubgroup)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredRotationcyclesSubgroup)
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
   #               "a"                         "b"                         "b"                         "b"                         "c" 


### 8.13 DurationSubgroup
FungalBeta_filteredDurationSubgroup <- subset(FungalBeta, DurationSubgroup %in% c("D1", "D2", "D3", "D4", "D5", "D6-10", "D11-20", "D20-30", "D30"))
#
FungalBeta_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(FungalBeta_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredDurationSubgroup %>%
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

overall_model_FungalBeta_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = FungalBeta_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -509.3617  1018.7234  1038.7234  1071.9535  1039.8574   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3838  0.6195     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 205) = 1829.9106, p-val < .0001
# Test of Moderators (coefficients 1:9):
# QM(df = 9) = 46.4671, p-val < .0001
# Model Results:
#                         estimate      se     zval    pval    ci.lb    ci.ub      
# DurationSubgroupD1        0.6658  0.3543   1.8794  0.0602  -0.0285   1.3601    . 
# DurationSubgroupD11-20    0.2552  0.1237   2.0624  0.0392   0.0127   0.4978    * 
# DurationSubgroupD2       -0.1328  0.1207  -1.1005  0.2711  -0.3694   0.1037      
# DurationSubgroupD20-30    0.5877  0.1443   4.0720  <.0001   0.3048   0.8706  *** 
# DurationSubgroupD3       -0.3170  0.1161  -2.7311  0.0063  -0.5444  -0.0895   ** 
# DurationSubgroupD30      -0.7797  0.6478  -1.2037  0.2287  -2.0493   0.4899      
# DurationSubgroupD4        0.3869  0.3192   1.2120  0.2255  -0.2387   1.0124      
# DurationSubgroupD5        0.7296  0.6716   1.0865  0.2773  -0.5866   2.0459      
# DurationSubgroupD6-10     0.4736  0.1376   3.4426  0.0006   0.2040   0.7432  ***    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredDurationSubgroup)
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
    #               "ab"                    "a"                   "cd"                    "b"                    "c"                  "acd" 
    # DurationSubgroupD4     DurationSubgroupD5  DurationSubgroupD6-10 
    #              "abd"                 "abcd"                   "ab" 



### 8.14 SpeciesRichnessSubgroup
FungalBeta_filteredSpeciesRichnessSubgroup <- subset(FungalBeta, SpeciesRichnessSubgroup %in% c("R2", "R3"))
#
FungalBeta_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup <- droplevels(factor(FungalBeta_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredSpeciesRichnessSubgroup %>%
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

overall_model_FungalBeta_filteredSpeciesRichnessSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichnessSubgroup, random = ~ 1 | StudyID, data = FungalBeta_filteredSpeciesRichnessSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredSpeciesRichnessSubgroup)
# Multivariate Meta-Analysis Model (k = 213; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -533.3001  1066.6001  1072.6001  1082.6557  1072.7161   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3687  0.6072     69     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 211) = 1860.5484, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 4.0493, p-val = 0.1320
# Model Results:
#                            estimate      se    zval    pval    ci.lb   ci.ub    
# SpeciesRichnessSubgroupR2    0.1444  0.0791  1.8265  0.0678  -0.0106  0.2994  . 
# SpeciesRichnessSubgroupR3    0.1953  0.1185  1.6475  0.0994  -0.0370  0.4276  . 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredSpeciesRichnessSubgroup)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredSpeciesRichnessSubgroup)
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
#                       "a"                       "a"                   


### 8.15 Primer
FungalBeta_filteredPrimer <- subset(FungalBeta, Primer %in% c("ITS1", "ITS1 + 5.8S + ITS2", "ITS2", "Full length", "18S"))
#
FungalBeta_filteredPrimer$Primer <- droplevels(factor(FungalBeta_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredPrimer %>%
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

overall_model_FungalBeta_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = FungalBeta_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 213; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -528.3860  1056.7720  1068.7720  1088.7972  1069.1899   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3657  0.6047     69     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 208) = 1772.0930, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 10.6287, p-val = 0.0593
# Model Results:
#                           estimate      se     zval    pval    ci.lb   ci.ub     
# Primer18S                   0.4718  0.4481   1.0528  0.2924  -0.4065  1.3502     
# PrimerFull length           0.0691  0.6232   0.1109  0.9117  -1.1524  1.2907     
# PrimerITS1                  0.2283  0.0862   2.6500  0.0080   0.0595  0.3972  ** 
# PrimerITS1 + 5.8S + ITS2    0.2393  0.4605   0.5197  0.6033  -0.6632  1.1419     
# PrimerITS2                 -0.3294  0.2213  -1.4884  0.1367  -0.7631  0.1044   
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredPrimer)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredPrimer)
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
               #      "ab"                     "ab"                      "a"                     "ab"                      "b" 




### 8.16 Latitude_Subgroup
FungalBeta_filteredLatitude_Subgroup <- subset(FungalBeta, Latitude_Subgroup %in% c("La20", "La20-40", "La40"))
#
FungalBeta_filteredLatitude_Subgroup$Latitude_Subgroup <- droplevels(factor(FungalBeta_filteredLatitude_Subgroup$Latitude_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredLatitude_Subgroup %>%
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

overall_model_FungalBeta_filteredLatitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Latitude_Subgroup, random = ~ 1 | StudyID, data = FungalBeta_filteredLatitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredLatitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -532.5487  1065.0975  1073.0975  1086.5049  1073.2917   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3773  0.6142     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 211) = 1791.9874, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 5.5567, p-val = 0.1353
# Model Results:
#                           estimate      se     zval    pval    ci.lb   ci.ub    
# Latitude_SubgroupLa20      -0.1422  0.3254  -0.4369  0.6622  -0.7800  0.4956    
# Latitude_SubgroupLa20-40    0.1571  0.1149   1.3671  0.1716  -0.0681  0.3822    
# Latitude_SubgroupLa40       0.2096  0.1121   1.8700  0.0615  -0.0101  0.4293  .       

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredLatitude_Subgroup)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredLatitude_Subgroup)
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
FungalBeta_filteredLongitude_Subgroup <- subset(FungalBeta, Longitude_Subgroup %in% c("Lo-180-0", "Lo-0-180"))
#
FungalBeta_filteredLongitude_Subgroup$Longitude_Subgroup <- droplevels(factor(FungalBeta_filteredLongitude_Subgroup$Longitude_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredLongitude_Subgroup %>%
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

overall_model_FungalBeta_filteredLongitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Subgroup, random = ~ 1 | StudyID, data = FungalBeta_filteredLongitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredLongitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -533.9578  1067.9156  1073.9156  1083.9854  1074.0310   
# Variance Components
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3813  0.6175     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 212) = 1900.5228, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 4.5000, p-val = 0.1054
# Model Results:
#                             estimate      se    zval    pval    ci.lb   ci.ub    
# Longitude_SubgroupLo-0-180    0.1617  0.0808  1.9998  0.0455   0.0032  0.3201  * 
# Longitude_SubgroupLo-180-0    0.2212  0.3125  0.7078  0.4791  -0.3913  0.8336              

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredLongitude_Subgroup)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredLongitude_Subgroup)
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
FungalBeta_filteredMAPmean_Subgroup <- subset(FungalBeta, MAPmean_Subgroup %in% c("MAP600", "MAP600-1200", "MAP1200"))
#
FungalBeta_filteredMAPmean_Subgroup$MAPmean_Subgroup <- droplevels(factor(FungalBeta_filteredMAPmean_Subgroup$MAPmean_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredMAPmean_Subgroup %>%
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

overall_model_FungalBeta_filteredMAPmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean_Subgroup, random = ~ 1 | StudyID, data = FungalBeta_filteredMAPmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredMAPmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -527.7065  1055.4130  1063.4130  1076.8204  1063.6072   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4434  0.6659     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 211) = 1791.5730, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 16.6370, p-val = 0.0008
# Model Results:
#                              estimate      se     zval    pval    ci.lb   ci.ub      
# MAPmean_SubgroupMAP1200        0.0178  0.1592   0.1115  0.9112  -0.2943  0.3298      
# MAPmean_SubgroupMAP600         0.3472  0.1039   3.3413  0.0008   0.1435  0.5509  *** 
# MAPmean_SubgroupMAP600-1200   -0.1172  0.1350  -0.8681  0.3853  -0.3818  0.1474     


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredMAPmean_Subgroup)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredMAPmean_Subgroup)
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
#                  "ab"                         "a"                         "b" 



### 8.19 MATmean_Subgroup
FungalBeta_filteredMATmean_Subgroup <- subset(FungalBeta, MATmean_Subgroup %in% c("MAT8", "MAT8-15", "MAT15"))
#
FungalBeta_filteredMATmean_Subgroup$MATmean_Subgroup <- droplevels(factor(FungalBeta_filteredMATmean_Subgroup$MATmean_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredMATmean_Subgroup %>%
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

overall_model_FungalBeta_filteredMATmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Subgroup, random = ~ 1 | StudyID, data = FungalBeta_filteredMATmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredMATmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 214; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -522.2592  1044.5184  1052.5184  1065.9258  1052.7126   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4029  0.6347     70     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 211) = 1794.8957, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 27.5152, p-val < .0001
# Model Results:
#                          estimate      se     zval    pval    ci.lb    ci.ub     
# MATmean_SubgroupMAT15      0.1704  0.1368   1.2455  0.2130  -0.0978   0.4387     
# MATmean_SubgroupMAT8       0.3190  0.1042   3.0624  0.0022   0.1148   0.5231  ** 
# MATmean_SubgroupMAT8-15   -0.3622  0.1473  -2.4581  0.0140  -0.6509  -0.0734   *    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredMATmean_Subgroup)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredMATmean_Subgroup)
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
#                  "a"                     "a"                     "b"



### 8.20 DurationSubgroup
FungalBeta_filteredDurationSubgroup <- subset(FungalBeta, DurationSubgroup %in% c("D5", "D5-10", "D10-20", "D20-30"))
#
FungalBeta_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(FungalBeta_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- FungalBeta_filteredDurationSubgroup %>%
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

overall_model_FungalBeta_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = FungalBeta_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalBeta_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 210; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -512.4283  1024.8566  1034.8566  1051.4960  1035.1566   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3987  0.6314     69     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 206) = 1867.1733, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 31.1781, p-val < .0001
# Model Results:
#                         estimate      se     zval    pval    ci.lb   ci.ub      
# DurationSubgroupD10-20    0.3568  0.1218   2.9285  0.0034   0.1180  0.5956   ** 
# DurationSubgroupD20-30    0.6613  0.1443   4.5814  <.0001   0.3784  0.9442  *** 
# DurationSubgroupD5       -0.0924  0.0988  -0.9347  0.3499  -0.2861  0.1013      
# DurationSubgroupD5-10     0.5110  0.1391   3.6737  0.0002   0.2384  0.7835  *** 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalBeta_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_FungalBeta_filteredDurationSubgroup)
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
# DurationSubgroupD10-20 DurationSubgroupD20-30       DurationSubgroupD5  DurationSubgroupD5-10 
#        "a"                    "b"                       "c"                   "ab" 





#### 9. Linear Mixed Effect Model
# 
FungalBeta$Wr <- 1 / FungalBeta$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalBeta)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalBeta)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalBeta)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalBeta)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalBeta)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalBeta)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalBeta)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalBeta)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(Species_Richness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = FungalBeta)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = FungalBeta)

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
# Model1   Model1 516.4949 536.6907 -252.2474 0.001214313     0.01377533
# Model2   Model2 513.9535 534.1493 -250.9767 0.002179444     0.01579099
# Model3   Model3 512.3501 532.5460 -250.1751 0.002837368     0.01425804
# Model4   Model4 516.6061 536.8019 -252.3030 0.001181415     0.01382114
# Model5   Model5 517.0338 537.2296 -252.5169 0.001336328     0.01518480
# Model6   Model6 517.1371 537.3330 -252.5686 0.001284530     0.01523850
# Model7   Model7 513.8820 534.0779 -250.9410 0.002182481     0.01572962
# Model8   Model8 512.3906 532.5864 -250.1953 0.002860901     0.01429292
# Model9   Model9 531.6051 565.2649 -255.8026 0.003252041     0.01837489
# Model10 Model10 531.2794 564.9392 -255.6397 0.002910837     0.01638472

##### Model 3 is the best model
summary(Model3)
# Number of obs: 214, groups:  StudyID, 70
anova(Model3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 39.174  39.174     1  96.862  5.6385 0.01954 *
# scale(Species_Richness)     19.804  19.804     1  57.982  2.8504 0.09672 .
# scale(Duration)             22.312  22.312     1 123.412  3.2115 0.07557 .
#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = FungalBeta)
summary(ModelpH)
# Number of obs: 118, groups:  StudyID, 48
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 21.9323 21.9323     1 86.065  3.5819 0.06177 .
# scale(Species_Richness)      5.5110  5.5110     1 77.673  0.9000 0.34571  
# scale(Duration)             13.6230 13.6230     1 94.800  2.2249 0.13913  
# scale(pHCK)                  1.7192  1.7192     1 57.077  0.2808 0.59825  

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = FungalBeta)
summary(ModelSOC)
# Number of obs: 132, groups:  StudyID, 49
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 35.410  35.410     1 70.697  4.8412 0.03106 *
# scale(Species_Richness)     20.454  20.454     1 56.424  2.7964 0.10001  
# scale(Duration)             28.769  28.769     1 78.888  3.9332 0.05082 .
# scale(SOCCK)                 2.185   2.185     1 35.159  0.2987 0.58815  

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = FungalBeta)
summary(ModelTN)
# Number of obs: 84, groups:  StudyID, 35
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 10.2576 10.2576     1 58.656  1.9021 0.1731
# scale(Species_Richness)      1.7615  1.7615     1 58.251  0.3266 0.5698
# scale(Duration)              7.5945  7.5945     1 76.378  1.4083 0.2390
# scale(TNCK)                  2.8796  2.8796     1 53.685  0.5340 0.4681


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = FungalBeta)
summary(ModelNO3)
# Number of obs: 42, groups:  StudyID, 17
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 9.7386  9.7386     1 28.625  2.0139 0.1667
# scale(Species_Richness)     0.8078  0.8078     1 16.186  0.1671 0.6881
# scale(Duration)             2.2087  2.2087     1 30.188  0.4568 0.5043
# scale(NO3CK)                0.1082  0.1082     1 14.739  0.0224 0.8831

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = FungalBeta)
summary(ModelNH4)
# Number of obs: 40, groups:  StudyID, 16
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.01211 0.01211     1 12.418  0.0029 0.9582
# scale(Species_Richness)     0.00442 0.00442     1 17.999  0.0010 0.9746
# scale(Duration)             0.47091 0.47091     1 33.902  0.1112 0.7408
# scale(NH4CK)                1.84970 1.84970     1  7.956  0.4369 0.5273

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = FungalBeta)
summary(ModelAP)
# Number of obs: 98, groups:  StudyID, 42
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 36.079  36.079     1 64.176  6.3728 0.01407 *
# scale(Species_Richness)      0.259   0.259     1 52.132  0.0457 0.83160  
# scale(Duration)             23.054  23.054     1 76.518  4.0722 0.04710 *
# scale(APCK)                  8.479   8.479     1 68.178  1.4977 0.22524 

#### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = FungalBeta)
summary(ModelAK)
# Number of obs: 76, groups:  StudyID, 35
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 29.3969 29.3969     1 44.043  4.5301 0.03894 *
# scale(Species_Richness)     19.5598 19.5598     1 17.031  3.0142 0.10059  
# scale(Duration)             14.0369 14.0369     1 46.961  2.1631 0.14803  
# scale(AKCK)                  1.0657  1.0657     1 51.340  0.1642 0.68698  

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = FungalBeta)
summary(ModelAN)
# Number of obs: 64, groups:  StudyID, 22
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles))  1.1848  1.1848     1 13.1715  0.1600 0.6956
# scale(Species_Richness)     13.5432 13.5432     1  6.5218  1.8291 0.2212
# scale(Duration)              7.0585  7.0585     1 16.8066  0.9533 0.3427
# scale(ANCK)                  6.8823  6.8823     1  7.9111  0.9295 0.3635

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = FungalBeta)
summary(ModelLatitude)
# Number of obs: 214, groups:  StudyID, 70
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 36.287  36.287     1  95.172  5.2224 0.02452 *
# scale(Species_Richness)     11.779  11.779     1  58.352  1.6953 0.19803  
# scale(Duration)             23.665  23.665     1 118.145  3.4058 0.06747 .
# scale(Latitude)             24.063  24.063     1  37.489  3.4631 0.07061 .

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = FungalBeta)
summary(ModelLongitude)
# Number of obs: 214, groups:  StudyID, 70
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 37.956  37.956     1  95.693  5.4707 0.02142 *
# scale(Species_Richness)     18.402  18.402     1  68.267  2.6523 0.10801  
# scale(Duration)             20.685  20.685     1 124.521  2.9814 0.08671 .
# scale(Longitude)             0.075   0.075     1  27.494  0.0108 0.91805  

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = FungalBeta)
summary(ModelMAPmean)
# Number of obs: 214, groups:  StudyID, 70
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 42.725  42.725     1  88.940  6.1462 0.01506 *
# scale(Species_Richness)     17.453  17.453     1  53.754  2.5108 0.11894  
# scale(Duration)             26.461  26.461     1 109.939  3.8066 0.05360 .
# scale(MAPmean)              10.398  10.398     1  40.558  1.4958 0.22838  

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = FungalBeta)
summary(ModelMATmean)
# Number of obs: 214, groups:  StudyID, 70
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 39.435  39.435     1  92.255  5.6672 0.01934 *
# scale(Species_Richness)     12.345  12.345     1  54.062  1.7741 0.18847  
# scale(Duration)             23.724  23.724     1 116.625  3.4095 0.06736 .
# scale(MATmean)              28.690  28.690     1  37.466  4.1230 0.04945 *



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## Species_Richness
sum(!is.na(FungalBeta$Species_Richness)) ## n = 214
p1 <- ggplot(FungalBeta, aes(y=RR, x=Species_Richness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="Species_Richness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Species_Richness" , y="lnFungalBeta214")
p1
pdf("Species_Richness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(FungalBeta$Duration)) ## n = 214
p2 <- ggplot(FungalBeta, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnFungalBeta214")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(FungalBeta$Rotation_cycles)) ## n = 214
p3 <- ggplot(FungalBeta, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnFungalBeta214")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(FungalBeta$Latitude)) ## n = 214
p5 <- ggplot(FungalBeta, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnFungalBeta214")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(FungalBeta$Longitude)) ## n = 214
p6 <- ggplot(FungalBeta, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnFungalBeta214")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(FungalBeta$MAPmean)) ## n = 214
p7 <- ggplot(FungalBeta, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnFungalBeta214")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(FungalBeta$MATmean)) ## n = 214
p8 <- ggplot(FungalBeta, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnFungalBeta214")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(FungalBeta$pHCK)) ## n = 118
p9 <- ggplot(FungalBeta, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnFungalBeta118")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(FungalBeta$SOCCK)) ## n = 132
p10 <- ggplot(FungalBeta, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnFungalBeta132")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(FungalBeta$TNCK)) ## n = 84
p11 <- ggplot(FungalBeta, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnFungalBeta84")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(FungalBeta$NO3CK)) ## n = 42
p12 <- ggplot(FungalBeta, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnFungalBeta42")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(FungalBeta$NH4CK)) ## n = 40
p13<- ggplot(FungalBeta, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnFungalBeta40")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(FungalBeta$APCK)) ## n = 98
p14 <- ggplot(FungalBeta, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnFungalBeta98")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(FungalBeta$AKCK)) ## n = 76
p15 <- ggplot(FungalBeta, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnFungalBeta76")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(FungalBeta$ANCK)) ## n = 64
p16 <- ggplot(FungalBeta, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnFungalBeta64")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(FungalBeta$RRpH)) ## n = 118
p17 <- ggplot(FungalBeta, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="RRpH")+
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
sum(!is.na(FungalBeta$RRSOC)) ## n = 132
p18 <- ggplot(FungalBeta, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="RRSOC")+
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
sum(!is.na(FungalBeta$RRTN)) ## n = 84
p19 <- ggplot(FungalBeta, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="RRTN")+
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
sum(!is.na(FungalBeta$RRNO3)) ## n = 42
p20 <- ggplot(FungalBeta, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="RRNO3")+
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
sum(!is.na(FungalBeta$RRNH4)) ## n = 40
p21 <- ggplot(FungalBeta, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="RRNH4")+
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
sum(!is.na(FungalBeta$RRAP)) ## n = 98
p22 <- ggplot(FungalBeta, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="RRAP")+
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
sum(!is.na(FungalBeta$RRAK)) ## n = 76
p23 <- ggplot(FungalBeta, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="RRAK")+
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
sum(!is.na(FungalBeta$RRAN)) ## n = 64
p24 <- ggplot(FungalBeta, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalBeta", x="RRAN")+
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
sum(!is.na(FungalBeta$RRYield)) ## n = 28
p25 <- ggplot(FungalBeta, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-1.2, 1.6), expand=c(0, 0)) + 
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
                            data = FungalBeta, 
                            na.action = na.roughfix, 
                            importance = TRUE, 
                            ntree = 500)
# Check the importance of variables and p-values
importance(rf_model_perm3)
#                    %IncMSE %IncMSE.pval IncNodePurity IncNodePurity.pval
# MATmean          18.340257   0.01980198     13.852809         0.00990099
# Latitude         15.170249   0.03960396     12.516648         0.01980198
# SOCCK            14.945640   0.01980198     21.008526         0.01980198
# Longitude        14.894347   0.03960396     11.877953         0.07920792
# MAPmean          14.405017   0.03960396     10.845285         0.14851485
# Rotation_cycles  11.164019   0.02970297      8.047624         0.08910891
# pHCK             10.562609   0.11881188     11.847687         0.69306931
# Duration          9.557627   0.21782178      8.549692         0.11881188
# Species_Richness  6.593804   0.08910891      1.898006         0.62376238
# 

####################################### Trials sorted by effect size
library(ggplot2)
library(dplyr)
FungalBeta <- read.csv("FungalBeta.csv", fileEncoding = "latin1")
# 计算95% CI + 显著性分类
df_plot <- FungalBeta %>%
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
#       37      115       62 
## 8*8


################################################ piecewiseSEM
sem_data <- read.csv("FungalBeta.csv", fileEncoding = "latin1")
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
# -110.642
# Fisher's C = 1.107 with P-value = 0.575

m1_1 <- lme(RR ~ RRpH + RRSOC + Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_2 <- lme(RRpH ~  Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_3 <- lme(RRSOC ~ Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model1 <- psem(m1_1, m1_2, m1_3)
summary(sem_model1)  ## 
# AIC
# -137.114
# Fisher's C = 2.465 with P-value = 0.292

m2_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_2 <- lme(RRpH ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_3 <- lme(RRSOC ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model2 <- psem(m2_1, m2_2, m2_3)
summary(sem_model2)  ## 
# AIC
# -140.626
# Fisher's C = 0.803 with P-value = 0.669

m3_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_2 <- lme(RRpH ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_3 <- lme(RRSOC ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model3 <- psem(m3_1, m3_2, m3_3)
summary(sem_model3)  ## 
# AIC
# -113.756
# Fisher's C = 0.366 with P-value = 0.833

m4_1 <- lme(RR ~ RRpH + RRSOC + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_2 <- lme(RRpH ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_3 <- lme(RRSOC ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model4 <- psem(m4_1, m4_2, m4_3)
summary(sem_model4) ## 
# AIC
# -167.296
# Fisher's C = 2.098 with P-value = 0.35

# Structural Equation Model of sem_model4 
# 
# Call:
#   RR ~ RRpH + RRSOC + Species_Richness
# RRpH ~ Species_Richness
# RRSOC ~ Species_Richness
# 
# AIC
# -167.296
# 
# ---
#   Tests of directed separation:
#   
#   Independ.Claim Test.Type DF Crit.Value P.Value 
# RRSOC ~ RRpH + ...      coef 56     0.9418  0.3503 
# 
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
#         RR             RRpH   0.0062    1.4877 55     0.0042  0.9967       0.0005    
#         RR            RRSOC  -0.5695    0.4713 55    -1.2083  0.2321      -0.1508    
#         RR Species_Richness   0.4472    0.2519 55     1.7752  0.0814       0.2418    
#       RRpH Species_Richness  -0.0527    0.0149 58    -3.5446  0.0008      -0.3248 ***
#      RRSOC Species_Richness  -0.0518    0.0458 71    -1.1304  0.2621      -0.1057    
# 
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05
# 
# ---
# Individual R-squared:
# 
#   Response method Marginal Conditional
#         RR   none     0.07        0.46
#       RRpH   none     0.10        0.84
#      RRSOC   none     0.01        0.64

m5_1 <- lme(RR ~ RRpH + RRSOC + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_2 <- lme(RRpH ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_3 <- lme(RRSOC ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model5 <- psem(m5_1, m5_2, m5_3)
summary(sem_model5) ## 
# AIC
# -138.477
# Fisher's C = 1.553 with P-value = 0.46

m6_1 <- lme(RR ~ RRpH + RRSOC +MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_2 <- lme(RRpH ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_3 <- lme(RRSOC ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model6 <- psem(m6_1, m6_2, m6_3)
summary(sem_model6)  ##
# AIC
# -144.728
# Fisher's C = 0.214 with P-value = 0.898

