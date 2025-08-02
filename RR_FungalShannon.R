
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)

FungalShannon <- read.csv("FungalShannon.csv", fileEncoding = "latin1")
# Check data
head(FungalShannon)

# 1. The number of Obversation
total_number <- nrow(FungalShannon)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 260

# 2. The number of Study
unique_studyid_number <- length(unique(FungalShannon$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 76


#### 3. Overall effect size
total_effect_model <- rma.mv(yi = RR, 
                             V = Vi, 
                             random = ~ 1 | StudyID,  # StudyID is radom factor
                             data = FungalShannon, 
                             method = "REML")
# The results of Overall effect size
summary(total_effect_model)
# Multivariate Meta-Analysis Model (k = 260; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -134.6784   269.3567   273.3567   280.4704   273.4036   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0125  0.1118     76     no  StudyID 
# Test for Heterogeneity:
# Q(df = 259) = 2247.5607, p-val < .0001
# Model Results:
# estimate      se    zval    pval    ci.lb   ci.ub    
#   0.0006  0.0143  0.0427  0.9660  -0.0274  0.0286    

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
boot_results1 <- boot(data = FungalShannon, statistic = boot_fun, R = 1000, parallel = "snow", ncpus = numCores, cl = cl)
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
# Estimate for coefficient 1 : 0.0006103393 
# 95% BCa CI for coefficient 1 : -0.01830931 0.01660475 


#### 5. Funnel Plot
simple_model <- rma(yi = RR, 
                    vi = Vi, 
                    data = FungalShannon, 
                    method = "REML")
#### 
funnel(simple_model)
# Output  6 * 6
#### Egger's test
regtest(simple_model)
# Regression Test for Funnel Plot Asymmetry
# Model:     mixed-effects meta-regression model
# Predictor: standard error
# Test for Funnel Plot Asymmetry: z = -0.3473, p = 0.7284
# Limit Estimate (as sei -> 0):   b =  0.0198 (CI: -0.0166, 0.0562)

#  Rosenthal’s Fail-Safe N
# This method estimates how many missing studies with null effect 
# would be needed to make the overall effect non-significant
fsn_rosenthal <- fsn(x = simple_model, type = "Rosenthal")
# Print the FSN result
print(fsn_rosenthal)
# Fail-safe N Calculation Using the General Approach
# Average Effect Size:         0.0144
# Amount of Heterogeneity:     0.0199
# Observed Significance Level: 0.1562
# Target Significance Level:   0.05
# Fail-safe N: 0


#### 8. Subgroup analysis
### 8.1 LegumeNonlegume
FungalShannon_filteredLegumeNonlegume <- subset(FungalShannon, LegumeNonlegume %in% c("Legume to Non-legume", "Non-legume to Legume", "Non-legume to Non-legume"))
#
FungalShannon_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(FungalShannon_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredLegumeNonlegume %>%
  group_by(LegumeNonlegume) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   LegumeNonlegume          Observations Unique_StudyID
#   <fct>                           <int>          <int>
# 1 Legume to Non-legume               44             23
# 2 Non-legume to Legume               97             26
# 3 Non-legume to Non-legume          118             43

overall_model_FungalShannon_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = FungalShannon_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 259; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -94.2823  188.5645  196.5645  210.7452  196.7239   
# Variance Components:
#            estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0152  0.1233     76     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 256) = 2179.5220, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 87.1099, p-val < .0001
# Model Results:
#                                          estimate      se     zval    pval    ci.lb    ci.ub      
# LegumeNonlegumeLegume to Non-legume       -0.0804  0.0195  -4.1125  <.0001  -0.1187  -0.0421  *** 
# LegumeNonlegumeNon-legume to Legume        0.0304  0.0188   1.6164  0.1060  -0.0065   0.0673      
# LegumeNonlegumeNon-legume to Non-legume    0.0203  0.0179   1.1341  0.2567  -0.0148   0.0555     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredLegumeNonlegume)
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
FungalShannon_filteredAMnonAM <- subset(FungalShannon, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
FungalShannon_filteredAMnonAM$AMnonAM <- droplevels(factor(FungalShannon_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredAMnonAM %>%
  group_by(AMnonAM) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   AMnonAM     Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 AM to AM             181             67
# 2 AM to nonAM           38             13
# 3 nonAM to AM           40              6
overall_model_FungalShannon_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = FungalShannon_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 259; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -131.8031   263.6062   271.6062   285.7869   271.7655   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0124  0.1115     76     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 256) = 2193.9033, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 12.4102, p-val = 0.0061
# Model Results:
#                     estimate      se     zval    pval    ci.lb   ci.ub     
# AMnonAMAM to AM      -0.0067  0.0145  -0.4648  0.6421  -0.0350  0.0216     
# AMnonAMAM to nonAM   -0.0101  0.0156  -0.6469  0.5177  -0.0406  0.0204     
# AMnonAMnonAM to AM    0.1048  0.0333   3.1476  0.0016   0.0395  0.1701  ** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredAMnonAM)
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
   #             "a"                "a"                "b" 


### 8.3 C3C4
FungalShannon_filteredC3C4 <- subset(FungalShannon, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
FungalShannon_filteredC3C4$C3C4 <- droplevels(factor(FungalShannon_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredC3C4 %>%
  group_by(C3C4) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   C3C4     Observations Unique_StudyID
#   <fct>           <int>          <int>
# 1 C3 to C3          137             46
# 2 C3 to C4           60             35
# 3 C4 to C3           59             17

overall_model_FungalShannon_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = FungalShannon_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 256; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -108.4652   216.9304   224.9304   239.0639   225.0917   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0121  0.1100     75     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 253) = 2107.3566, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 0.7983, p-val = 0.8499
# Model Results:
#               estimate      se     zval    pval    ci.lb   ci.ub    
# C3C4C3 to C3   -0.0017  0.0152  -0.1146  0.9088  -0.0316  0.0281    
# C3C4C3 to C4    0.0037  0.0156   0.2359  0.8135  -0.0270  0.0344    
# C3C4C4 to C3   -0.0117  0.0208  -0.5633  0.5732  -0.0525  0.0291    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredC3C4)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredC3C4)
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
#          "a"          "a"          "a" 


### 8.4 Annual_Pere
FungalShannon_filteredAnnual_Pere <- subset(FungalShannon, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
FungalShannon_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(FungalShannon_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredAnnual_Pere %>%
  group_by(Annual_Pere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Annual_Pere            Observations Unique_StudyID
#   <fct>                         <int>          <int>
# 1 Annual to Annual                196             56
# 2 Annual to Perennial              23             12
# 3 Perennial to Annual              30             12
# 4 Perennial to Perennial           11              7

overall_model_FungalShannon_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = FungalShannon_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 260; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -137.0576   274.1152   284.1152   301.8411   284.3552   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0123  0.1109     76     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 256) = 2145.3197, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 2.5209, p-val = 0.6409
# Model Results:
#                                    estimate      se     zval    pval    ci.lb   ci.ub    
# Annual_PereAnnual to Annual         -0.0092  0.0158  -0.5806  0.5615  -0.0401  0.0218    
# Annual_PereAnnual to Perennial      -0.0034  0.0229  -0.1494  0.8813  -0.0483  0.0414    
# Annual_PerePerennial to Annual       0.0391  0.0291   1.3444  0.1788  -0.0179  0.0960    
# Annual_PerePerennial to Perennial    0.0239  0.0371   0.6432  0.5201  -0.0489  0.0967  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredAnnual_Pere)
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
      #                        "a"                               "a"                               "a"                               "a" 



### 8.6 PlantStageSubgroup
FungalShannon_filteredPlantStageSubgroup <- subset(FungalShannon, PlantStageSubgroup %in% c("Vegetative stage","Reproductive stage", "Maturity stage","Harvest"))
#
FungalShannon_filteredPlantStageSubgroup$PlantStageSubgroup <- droplevels(factor(FungalShannon_filteredPlantStageSubgroup$PlantStageSubgroup))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredPlantStageSubgroup %>%
  group_by(PlantStageSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   PlantStageSubgroup Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Harvest                     103             35
# 2 Maturity stage               31             13
# 3 Reproductive stage           65             14
# 4 Vegetative stage             28              7

overall_model_FungalShannon_filteredPlantStageSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + PlantStageSubgroup, random = ~ 1 | StudyID, data = FungalShannon_filteredPlantStageSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredPlantStageSubgroup)
# Multivariate Meta-Analysis Model (k = 227; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -104.8820   209.7640   219.7640   236.7998   220.0405   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0132  0.1150     62     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 223) = 1812.0469, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 47.1637, p-val < .0001
# Model Results:
#                                       estimate      se     zval    pval    ci.lb    ci.ub      
# PlantStageSubgroupHarvest               0.0424  0.0177   2.3939  0.0167   0.0077   0.0770    * 
# PlantStageSubgroupMaturity stage       -0.0998  0.0245  -4.0698  <.0001  -0.1479  -0.0518  *** 
# PlantStageSubgroupReproductive stage   -0.0126  0.0247  -0.5087  0.6110  -0.0610   0.0359      
# PlantStageSubgroupVegetative stage      0.0495  0.0198   2.5032  0.0123   0.0107   0.0882    * 
#    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredPlantStageSubgroup)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredPlantStageSubgroup)
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
           #                       "a"                                  "b"                                  "c"                                  "a" 


### 8.7 Bulk_Rhizosphere
FungalShannon_filteredBulk_Rhizosphere <- subset(FungalShannon, Bulk_Rhizosphere %in% c("Non-Rhizosphere", "Rhizosphere"))
#
FungalShannon_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(FungalShannon_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredBulk_Rhizosphere %>%
  group_by(Bulk_Rhizosphere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Bulk_Rhizosphere Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 Non-Rhizosphere           146             43
# 2 Rhizosphere               114             40

overall_model_FungalShannon_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = FungalShannon_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 260; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -135.1741   270.3483   276.3483   287.0071   276.4427   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0127  0.1127     76     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 258) = 2232.0667, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 2.2296, p-val = 0.3280
# Model Results:
#                                  estimate      se     zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon-Rhizosphere   -0.0065  0.0152  -0.4310  0.6665  -0.0363  0.0232    
# Bulk_RhizosphereRhizosphere        0.0084  0.0153   0.5473  0.5841  -0.0217  0.0384   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredBulk_Rhizosphere)
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
FungalShannon_filteredSoil_texture <- subset(FungalShannon, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
FungalShannon_filteredSoil_texture$Soil_texture <- droplevels(factor(FungalShannon_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredSoil_texture %>%
  group_by(Soil_texture) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Soil_texture Observations Unique_StudyID
#   <fct>               <int>          <int>
# 1 Coarse                 38              8
# 2 Fine                   42             12
# 3 Medium                114             28

overall_model_FungalShannon_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = FungalShannon_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 194; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -28.6951   57.3903   65.3903   78.3993   65.6053   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0102  0.1012     48     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 191) = 1284.9935, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.4300, p-val = 0.6985
# Model Results:
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse    0.0377  0.0382   0.9855  0.3244  -0.0372  0.1126    
# Soil_textureFine      0.0119  0.0331   0.3591  0.7195  -0.0529  0.0767    
# Soil_textureMedium   -0.0121  0.0211  -0.5743  0.5658  -0.0536  0.0293      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredSoil_texture)
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
FungalShannon_filteredTillage <- subset(FungalShannon, Tillage %in% c("Tillage", "No_tillage"))
#
FungalShannon_filteredTillage$Tillage <- droplevels(factor(FungalShannon_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredTillage %>%
  group_by(Tillage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Tillage    Observations Unique_StudyID
#   <fct>             <int>          <int>
# 1 No_tillage           31              5
# 2 Tillage              35             11

overall_model_FungalShannon_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = FungalShannon_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredTillage)
# Multivariate Meta-Analysis Model (k = 66; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
#  16.9106  -33.8213  -27.8213  -21.3446  -27.4213   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0099  0.0996     13     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 64) = 252.8101, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 9.9339, p-val = 0.0070
# Model Results:
#                    estimate      se     zval    pval    ci.lb    ci.ub     
# TillageNo_tillage   -0.1082  0.0350  -3.0908  0.0020  -0.1769  -0.0396  ** 
# TillageTillage      -0.0636  0.0312  -2.0376  0.0416  -0.1249  -0.0024   *  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredTillage)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredTillage)
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
#               "a"               "b" 


### 8.10 Straw_retention
FungalShannon_filteredStraw_retention <- subset(FungalShannon, Straw_retention %in% c("Retention", "No_retention"))
#
FungalShannon_filteredStraw_retention$Straw_retention <- droplevels(factor(FungalShannon_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredStraw_retention %>%
  group_by(Straw_retention) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#  Straw_retention Observations Unique_StudyID
#   <fct>                  <int>          <int>
# 1 No_retention              18             10
# 2 Retention                 35             11

overall_model_FungalShannon_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = FungalShannon_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 53; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -45.4257   90.8514   96.8514  102.6469   97.3620   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0087  0.0932     18     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 51) = 422.2661, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.9596, p-val = 0.3754
# Model Results:
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention   -0.0349  0.0253  -1.3755  0.1690  -0.0845  0.0148    
# Straw_retentionRetention      -0.0337  0.0252  -1.3366  0.1813  -0.0830  0.0157    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredStraw_retention)
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
FungalShannon_filteredRotationcyclesSubgroup <- subset(FungalShannon, RotationcyclesSubgroup %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
FungalShannon_filteredRotationcyclesSubgroup$RotationcyclesSubgroup <- droplevels(factor(FungalShannon_filteredRotationcyclesSubgroup$RotationcyclesSubgroup))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredRotationcyclesSubgroup %>%
  group_by(RotationcyclesSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   RotationcyclesSubgroup Observations Unique_StudyID
#   <fct>                         <int>          <int>
# 1 D1                               91             33
# 2 D1-3                             51             18
# 3 D10                              14              5
# 4 D3-5                             36             10
# 5 D5-10                            68             19

overall_model_FungalShannon_filteredRotationcyclesSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + RotationcyclesSubgroup, random = ~ 1 | StudyID, data = FungalShannon_filteredRotationcyclesSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredRotationcyclesSubgroup)
# Multivariate Meta-Analysis Model (k = 260; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -136.9255   273.8509   285.8509   307.0985   286.1896   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0146  0.1209     76     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 255) = 2181.2683, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 6.2570, p-val = 0.2820
# Model Results:
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# RotationcyclesSubgroupD1      -0.0230  0.0186  -1.2368  0.2162  -0.0594  0.0134    
# RotationcyclesSubgroupD1-3    -0.0208  0.0189  -1.1006  0.2711  -0.0579  0.0163    
# RotationcyclesSubgroupD10      0.0202  0.0502   0.4028  0.6871  -0.0781  0.1185    
# RotationcyclesSubgroupD3-5     0.0510  0.0264   1.9336  0.0532  -0.0007  0.1028  . 
# RotationcyclesSubgroupD5-10    0.0321  0.0233   1.3776  0.1683  -0.0136  0.0778 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredRotationcyclesSubgroup)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredRotationcyclesSubgroup)
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
   #                      "a"                         "a"                        "ab"                         "b"                         "b" 


### 8.13 DurationSubgroup
FungalShannon_filteredDurationSubgroup <- subset(FungalShannon, DurationSubgroup %in% c("D1", "D2", "D3", "D4", "D5", "D6-10", "D11-20", "D20-30", "D30"))
#
FungalShannon_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(FungalShannon_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredDurationSubgroup %>%
  group_by(DurationSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   DurationSubgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D1                         18              7
# 2 D11-20                     71             14
# 3 D2                         36             18
# 4 D20-30                     28              8
# 5 D3                         40             13
# 6 D30                         6              2
# 7 D4                          7              4
# 8 D5                         10              2
# 9 D6-10                      44             14

overall_model_FungalShannon_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = FungalShannon_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 260; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -130.2319   260.4638   280.4638   315.7184   281.3805   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0169  0.1299     76     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 251) = 2089.9598, p-val < .0001
# Test of Moderators (coefficients 1:9):
# QM(df = 9) = 20.7439, p-val = 0.0138
# Model Results:
#                         estimate      se     zval    pval    ci.lb    ci.ub    
# DurationSubgroupD1        0.0479  0.0544   0.8804  0.3787  -0.0588   0.1546    
# DurationSubgroupD11-20    0.0679  0.0283   2.4042  0.0162   0.0125   0.1233  * 
# DurationSubgroupD2       -0.0173  0.0353  -0.4900  0.6241  -0.0864   0.0518    
# DurationSubgroupD20-30    0.0277  0.0347   0.7966  0.4257  -0.0404   0.0957    
# DurationSubgroupD3       -0.0724  0.0334  -2.1683  0.0301  -0.1379  -0.0070  * 
# DurationSubgroupD30      -0.0725  0.0974  -0.7443  0.4567  -0.2633   0.1183    
# DurationSubgroupD4        0.0278  0.0665   0.4187  0.6754  -0.1025   0.1581    
# DurationSubgroupD5       -0.1240  0.1187  -1.0448  0.2961  -0.3567   0.1086    
# DurationSubgroupD6-10    -0.0108  0.0355  -0.3034  0.7616  -0.0803   0.0588    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredDurationSubgroup)
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
#     DurationSubgroupD1 DurationSubgroupD11-20     DurationSubgroupD2 DurationSubgroupD20-30     DurationSubgroupD3    DurationSubgroupD30 
#                   "ab"                    "a"                   "ab"                    "a"                    "b"                   "ab" 
#     DurationSubgroupD4     DurationSubgroupD5  DurationSubgroupD6-10 
#                   "ab"                   "ab"                   "ab" 



### 8.14 SpeciesRichnessSubgroup
FungalShannon_filteredSpeciesRichnessSubgroup <- subset(FungalShannon, SpeciesRichnessSubgroup %in% c("R2", "R3"))
#
FungalShannon_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup <- droplevels(factor(FungalShannon_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredSpeciesRichnessSubgroup %>%
  group_by(SpeciesRichnessSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   SpeciesRichnessSubgroup Observations Unique_StudyID
#   <fct>                          <int>          <int>
# 1 R2                               187             68
# 2 R3                                69             14

overall_model_FungalShannon_filteredSpeciesRichnessSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichnessSubgroup, random = ~ 1 | StudyID, data = FungalShannon_filteredSpeciesRichnessSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredSpeciesRichnessSubgroup)
# Multivariate Meta-Analysis Model (k = 256; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -137.4547   274.9094   280.9094   291.5214   281.0054 
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0126  0.1122     74     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 254) = 2221.8527, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 3.2509, p-val = 0.1968
# Model Results:
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# SpeciesRichnessSubgroupR2    0.0056  0.0145   0.3860  0.6995  -0.0228  0.0340    
# SpeciesRichnessSubgroupR3   -0.0062  0.0155  -0.3970  0.6914  -0.0366  0.0242    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredSpeciesRichnessSubgroup)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredSpeciesRichnessSubgroup)
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
#                      "a"                       "a"                      


### 8.15 Primer
FungalShannon_filteredPrimer <- subset(FungalShannon, Primer %in% c("ITS1", "ITS1 + 5.8S + ITS2", "ITS2", "Full length", "18S"))
#
FungalShannon_filteredPrimer$Primer <- droplevels(factor(FungalShannon_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredPrimer %>%
  group_by(Primer) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Primer             Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 18S                           6              4
# 2 Full length                   3              1
# 3 ITS1                        217             58
# 4 ITS1 + 5.8S + ITS2            5              2
# 5 ITS2                         29             11

overall_model_FungalShannon_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = FungalShannon_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 260; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -135.5122   271.0243   283.0243   304.2719   283.3630   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0126  0.1122     76     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 255) = 2190.1712, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 3.8932, p-val = 0.5649
# Model Results:
#                           estimate      se     zval    pval    ci.lb   ci.ub    
# Primer18S                  -0.0712  0.0809  -0.8801  0.3788  -0.2298  0.0874    
# PrimerFull length           0.1624  0.1194   1.3599  0.1738  -0.0717  0.3965    
# PrimerITS1                  0.0060  0.0162   0.3722  0.7097  -0.0258  0.0379    
# PrimerITS1 + 5.8S + ITS2    0.0337  0.0810   0.4167  0.6769  -0.1250  0.1925    
# PrimerITS2                 -0.0372  0.0380  -0.9782  0.3280  -0.1118  0.0374    
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredPrimer)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredPrimer)
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
               #       "a"                      "a"                      "a"                      "a"                      "a" 



### 8.16 Latitude_Subgroup
FungalShannon_filteredLatitude_Subgroup <- subset(FungalShannon, Latitude_Subgroup %in% c("La20-40", "La40"))
#
FungalShannon_filteredLatitude_Subgroup$Latitude_Subgroup <- droplevels(factor(FungalShannon_filteredLatitude_Subgroup$Latitude_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredLatitude_Subgroup %>%
  group_by(Latitude_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Latitude_Subgroup Observations Unique_StudyID
#   <fct>                    <int>          <int>
# 1 La20-40                     95             39
# 2 La40                       158             34

overall_model_FungalShannon_filteredLatitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Latitude_Subgroup, random = ~ 1 | StudyID, data = FungalShannon_filteredLatitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredLatitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 253; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -142.6896   285.3792   291.3792   301.9555   291.4763   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0125  0.1119     73     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 251) = 2239.7434, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 2.3154, p-val = 0.3142
# Model Results:
#                           estimate      se     zval    pval    ci.lb   ci.ub    
# Latitude_SubgroupLa20-40    0.0198  0.0207   0.9574  0.3384  -0.0207  0.0604    
# Latitude_SubgroupLa40      -0.0244  0.0206  -1.1827  0.2369  -0.0648  0.0160      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredLatitude_Subgroup)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredLatitude_Subgroup)
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
# Latitude_SubgroupLa20-40    Latitude_SubgroupLa40 
#              "a"                      "a" 


### 8.17 Longitude_Subgroup
FungalShannon_filteredLongitude_Subgroup <- subset(FungalShannon, Longitude_Subgroup %in% c("Lo-180-0", "Lo-0-180"))
#
FungalShannon_filteredLongitude_Subgroup$Longitude_Subgroup <- droplevels(factor(FungalShannon_filteredLongitude_Subgroup$Longitude_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredLongitude_Subgroup %>%
  group_by(Longitude_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Longitude_Subgroup Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Lo-0-180                    205             70
# 2 Lo-180-0                     55              6

overall_model_FungalShannon_filteredLongitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Subgroup, random = ~ 1 | StudyID, data = FungalShannon_filteredLongitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredLongitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 260; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -134.8584   269.7167   275.7167   286.3756   275.8112   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0127  0.1128     76     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 258) = 2241.0491, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.0053, p-val = 0.9974
# Model Results:
#                             estimate      se     zval    pval    ci.lb   ci.ub    
# Longitude_SubgroupLo-0-180    0.0008  0.0151   0.0553  0.9559  -0.0287  0.0304    
# Longitude_SubgroupLo-180-0   -0.0023  0.0487  -0.0473  0.9623  -0.0977  0.0931         

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredLongitude_Subgroup)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredLongitude_Subgroup)
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
FungalShannon_filteredMAPmean_Subgroup <- subset(FungalShannon, MAPmean_Subgroup %in% c("MAP600", "MAP600-1200", "MAP1200"))
#
FungalShannon_filteredMAPmean_Subgroup$MAPmean_Subgroup <- droplevels(factor(FungalShannon_filteredMAPmean_Subgroup$MAPmean_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredMAPmean_Subgroup %>%
  group_by(MAPmean_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MAPmean_Subgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 MAP1200                    43             19
# 2 MAP600                    163             38
# 3 MAP600-1200                54             21

overall_model_FungalShannon_filteredMAPmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean_Subgroup, random = ~ 1 | StudyID, data = FungalShannon_filteredMAPmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredMAPmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 260; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -136.1296   272.2592   280.2592   294.4555   280.4179   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0128  0.1130     76     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 257) = 2246.7735, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 0.1531, p-val = 0.9848
# Model Results:
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean_SubgroupMAP1200        0.0061  0.0297   0.2052  0.8374  -0.0522  0.0644    
# MAPmean_SubgroupMAP600        -0.0039  0.0185  -0.2108  0.8330  -0.0402  0.0324    
# MAPmean_SubgroupMAP600-1200    0.0049  0.0248   0.1978  0.8432  -0.0437  0.0535    



# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredMAPmean_Subgroup)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredMAPmean_Subgroup)
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
FungalShannon_filteredMATmean_Subgroup <- subset(FungalShannon, MATmean_Subgroup %in% c("MAT8", "MAT8-15", "MAT15"))
#
FungalShannon_filteredMATmean_Subgroup$MATmean_Subgroup <- droplevels(factor(FungalShannon_filteredMATmean_Subgroup$MATmean_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredMATmean_Subgroup %>%
  group_by(MATmean_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MATmean_Subgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 MAT15                      63             28
# 2 MAT8                      140             34
# 3 MAT8-15                    57             15

overall_model_FungalShannon_filteredMATmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Subgroup, random = ~ 1 | StudyID, data = FungalShannon_filteredMATmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredMATmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 260; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -135.9960   271.9920   279.9920   294.1883   280.1508   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0128  0.1129     76     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 257) = 2211.1415, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 0.2526, p-val = 0.9687
# Model Results:
#                          estimate      se     zval    pval    ci.lb   ci.ub    
# MATmean_SubgroupMAT15      0.0049  0.0243   0.2033  0.8389  -0.0427  0.0526    
# MATmean_SubgroupMAT8      -0.0056  0.0198  -0.2815  0.7783  -0.0443  0.0332    
# MATmean_SubgroupMAT8-15    0.0072  0.0271   0.2677  0.7889  -0.0458  0.0603      


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredMATmean_Subgroup)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredMATmean_Subgroup)
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
FungalShannon_filteredDurationSubgroup <- subset(FungalShannon, DurationSubgroup %in% c("D5", "D5-10", "D10-20", "D20-30"))
#
FungalShannon_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(FungalShannon_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- FungalShannon_filteredDurationSubgroup %>%
  group_by(DurationSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   DurationSubgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D10-20                     71             14
# 2 D20-30                     28              8
# 3 D5                        111             42
# 4 D5-10                      44             14

overall_model_FungalShannon_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = FungalShannon_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalShannon_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 254; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -133.2373   266.4746   276.4746   294.0819   276.7205 
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0156  0.1249     74     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 250) = 2128.9227, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 14.3144, p-val = 0.0064
# Model Results:
#                         estimate      se     zval    pval    ci.lb   ci.ub     
# DurationSubgroupD10-20    0.0800  0.0263   3.0449  0.0023   0.0285  0.1315  ** 
# DurationSubgroupD20-30    0.0334  0.0336   0.9927  0.3208  -0.0325  0.0993     
# DurationSubgroupD5       -0.0263  0.0203  -1.2979  0.1943  -0.0660  0.0134     
# DurationSubgroupD5-10    -0.0090  0.0345  -0.2610  0.7941  -0.0765  0.0586  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalShannon_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_FungalShannon_filteredDurationSubgroup)
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
#                    "a"                   "ab"                    "b"                    "b" 







#### 9. Linear Mixed Effect Model
# 
FungalShannon$Wr <- 1 / FungalShannon$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalShannon)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalShannon)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalShannon)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalShannon)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalShannon)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalShannon)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalShannon)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalShannon)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(Species_Richness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = FungalShannon)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = FungalShannon)

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
# Model1   Model1 -140.6275 -119.26337 76.31373 3.856106e-05    0.001599910
# Model2   Model2 -140.6317 -119.26756 76.31583 2.980365e-05    0.001622214
# Model3   Model3 -140.7049 -119.34079 76.35244 3.952136e-05    0.001594012
# Model4   Model4 -140.3164 -118.95231 76.15820 3.698732e-05    0.001591713
# Model5   Model5 -140.5031 -119.13899 76.25154 3.146524e-05    0.001633232
# Model6   Model6 -140.1676 -118.80348 76.08379 2.836993e-05    0.001627177
# Model7   Model7 -140.9814 -119.61726 76.49068 3.410398e-05    0.001629514
# Model8   Model8 -140.3807 -119.01664 76.19036 3.757738e-05    0.001585466
# Model9   Model9 -112.9880  -77.38119 66.49400 1.977753e-04    0.001841755
# Model10 Model10 -114.2364  -78.62959 67.11821 2.577722e-04    0.001853040

##### Model 7 is the best model
summary(Model7)
# Number of obs: 260, groups:  StudyID, 76
anova(Model7) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.4599  1.4599     1  91.785  0.3007 0.5848
# scale(Species_Richness)     5.4787  5.4787     1 251.602  1.1283 0.2892
# scale(log(Duration))        0.1687  0.1687     1  77.144  0.0347 0.8526


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + scale(pHCK) + (1 | StudyID), weights = Wr, data = FungalShannon)
summary(ModelpH)
# Number of obs: 135, groups:  StudyID, 49
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.98484 0.98484     1  38.780  0.2256 0.6375
# scale(Species_Richness)     1.04871 1.04871     1 129.561  0.2402 0.6249
# scale(log(Duration))        1.92858 1.92858     1  42.035  0.4417 0.5099
# scale(pHCK)                 0.37448 0.37448     1  49.693  0.0858 0.7709

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = FungalShannon)
summary(ModelSOC)
# Number of obs: 149, groups:  StudyID, 52
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.63254 0.63254     1  42.096  0.1285 0.7218
# scale(Species_Richness)     1.99847 1.99847     1 143.338  0.4060 0.5250
# scale(log(Duration))        2.48372 2.48372     1  41.697  0.5045 0.4815
# scale(SOCCK)                0.40761 0.40761     1  47.339  0.0828 0.7748

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + scale(TNCK) + (1 | StudyID), weights = Wr, data = FungalShannon)
summary(ModelTN)
# Number of obs: 91, groups:  StudyID, 35
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.0085  0.0085     1 26.356  0.0029 0.9574
# scale(Species_Richness)     3.8987  3.8987     1 79.459  1.3296 0.2523
# scale(log(Duration))        0.1004  0.1004     1 26.483  0.0342 0.8546
# scale(TNCK)                 1.8085  1.8085     1 40.890  0.6168 0.4368 


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = FungalShannon)
summary(ModelNO3)
# Number of obs: 41, groups:  StudyID, 17
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.07162 0.07162     1 13.709  0.0588 0.8120
# scale(Species_Richness)     0.14157 0.14157     1 27.625  0.1162 0.7357
# scale(log(Duration))        0.83251 0.83251     1 16.618  0.6836 0.4201
# scale(NO3CK)                0.29890 0.29890     1 11.301  0.2454 0.6298

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = FungalShannon)
summary(ModelNH4)
# Number of obs: 39, groups:  StudyID, 16
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 0.0779  0.0779     1  7.7005  0.0726 0.79472  
# scale(Species_Richness)     0.0417  0.0417     1 29.8001  0.0389 0.84506  
# scale(log(Duration))        0.6079  0.6079     1  7.8352  0.5665 0.47365  
# scale(NH4CK)                6.1061  6.1061     1  9.7633  5.6900 0.03883 *

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + scale(APCK) + (1 | StudyID), weights = Wr, data = FungalShannon)
summary(ModelAP)
# Number of obs: 109, groups:  StudyID, 41
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  0.9052  0.9052     1  36.654  0.2006 0.65688  
# scale(Species_Richness)      3.6275  3.6275     1 101.850  0.8038 0.37206  
# scale(log(Duration))         0.8608  0.8608     1  30.617  0.1907 0.66536  
# scale(APCK)                 18.0912 18.0912     1  42.077  4.0089 0.05173 .

#### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + scale(AKCK) + (1 | StudyID), weights = Wr, data = FungalShannon)
summary(ModelAK)
# Number of obs: 87, groups:  StudyID, 35
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.01473 0.01473     1 35.510  0.0028 0.9583
# scale(Species_Richness)     1.66544 1.66544     1 81.656  0.3133 0.5772
# scale(log(Duration))        0.50196 0.50196     1 28.878  0.0944 0.7608
# scale(AKCK)                 1.13055 1.13055     1 66.577  0.2127 0.6462

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + scale(ANCK) + (1 | StudyID), weights = Wr, data = FungalShannon)
summary(ModelAN)
# Number of obs: 72, groups:  StudyID, 20
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 2.1503  2.1503     1  9.813  0.5135 0.4903
# scale(Species_Richness)     0.2015  0.2015     1 32.408  0.0481 0.8278
# scale(log(Duration))        1.6226  1.6226     1  8.702  0.3875 0.5496
# scale(ANCK)                 4.0713  4.0713     1  6.325  0.9723 0.3603

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + scale(Latitude) + (1 | StudyID), weights = Wr, data = FungalShannon)
summary(ModelLatitude)
# Number of obs: 260, groups:  StudyID, 76
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.6093  1.6093     1  90.588  0.3316 0.5662
# scale(Species_Richness)     4.7126  4.7126     1 250.325  0.9710 0.3254
# scale(log(Duration))        0.6188  0.6188     1  79.760  0.1275 0.7220
# scale(Latitude)             5.4193  5.4193     1  82.556  1.1166 0.2937

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + scale(Longitude) + (1 | StudyID), weights = Wr, data = FungalShannon)
summary(ModelLongitude)
# Number of obs: 260, groups:  StudyID, 76
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.4596  1.4596     1  90.783  0.3006 0.5848
# scale(Species_Richness)     5.5171  5.5171     1 250.194  1.1363 0.2875
# scale(log(Duration))        0.1840  0.1840     1  75.736  0.0379 0.8462
# scale(Longitude)            0.0009  0.0009     1  52.730  0.0002 0.9890

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = FungalShannon)
summary(ModelMAPmean)
# Number of obs: 260, groups:  StudyID, 76
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.5209  1.5209     1  95.763  0.3133 0.5770
# scale(Species_Richness)     5.5640  5.5640     1 249.412  1.1461 0.2854
# scale(log(Duration))        0.2401  0.2401     1  81.244  0.0495 0.8246
# scale(MAPmean)              0.0602  0.0602     1  56.898  0.0124 0.9117

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + scale(MATmean) + (1 | StudyID), weights = Wr, data = FungalShannon)
summary(ModelMATmean)
# Number of obs: 260, groups:  StudyID, 76
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 2.9914  2.9914     1  95.120  0.6187 0.4335
# scale(Species_Richness)     6.1379  6.1379     1 249.665  1.2694 0.2610
# scale(log(Duration))        1.7180  1.7180     1  83.536  0.3553 0.5527
# scale(MATmean)              6.0500  6.0500     1  56.322  1.2513 0.2681 



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## Species_Richness
sum(!is.na(FungalShannon$Species_Richness)) ## n = 260
p1 <- ggplot(FungalShannon, aes(y=RR, x=Species_Richness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="Species_Richness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Species_Richness" , y="lnFungalShannon260")
p1
pdf("Species_Richness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(FungalShannon$Duration)) ## n = 260
p2 <- ggplot(FungalShannon, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnFungalShannon260")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(FungalShannon$Rotation_cycles)) ## n = 260
p3 <- ggplot(FungalShannon, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnFungalShannon260")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(FungalShannon$Latitude)) ## n = 393
p5 <- ggplot(FungalShannon, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnFungalShannon260")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(FungalShannon$Longitude)) ## n = 260
p6 <- ggplot(FungalShannon, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnFungalShannon260")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(FungalShannon$MAPmean)) ## n = 260
p7 <- ggplot(FungalShannon, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnFungalShannon260")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(FungalShannon$MATmean)) ## n = 260
p8 <- ggplot(FungalShannon, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnFungalShannon260")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(FungalShannon$pHCK)) ## n = 135
p9 <- ggplot(FungalShannon, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnFungalShannon135")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(FungalShannon$SOCCK)) ## n = 149
p10 <- ggplot(FungalShannon, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnFungalShannon149")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(FungalShannon$TNCK)) ## n = 91
p11 <- ggplot(FungalShannon, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnFungalShannon91")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(FungalShannon$NO3CK)) ## n = 41
p12 <- ggplot(FungalShannon, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnFungalShannon41")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(FungalShannon$NH4CK)) ## n = 39
p13<- ggplot(FungalShannon, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnFungalShannon39")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(FungalShannon$APCK)) ## n = 109
p14 <- ggplot(FungalShannon, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnFungalShannon109")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(FungalShannon$AKCK)) ## n = 87
p15 <- ggplot(FungalShannon, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnFungalShannon87")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(FungalShannon$ANCK)) ## n = 72
p16 <- ggplot(FungalShannon, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnFungalShannon72")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(FungalShannon$RRpH)) ## n = 135
p17 <- ggplot(FungalShannon, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="RRpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRpH" , y="RR135")
p17
pdf("RRpH.pdf",width=8,height=8)
p17
dev.off() 

## RRSOC
sum(!is.na(FungalShannon$RRSOC)) ## n = 149
p18 <- ggplot(FungalShannon, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="RRSOC")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRSOC" , y="RR149")
p18
pdf("RRSOC.pdf",width=8,height=8)
p18
dev.off() 

## RRTN
sum(!is.na(FungalShannon$RRTN)) ## n = 91
p19 <- ggplot(FungalShannon, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="RRTN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRTN" , y="RR91")
p19
pdf("RRTN.pdf",width=8,height=8)
p19
dev.off() 

## RRNO3
sum(!is.na(FungalShannon$RRNO3)) ## n = 41
p20 <- ggplot(FungalShannon, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="RRNO3")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNO3" , y="RR41")
p20
pdf("RRNO3.pdf",width=8,height=8)
p20
dev.off() 

## RRNH4
sum(!is.na(FungalShannon$RRNH4)) ## n = 39
p21 <- ggplot(FungalShannon, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="RRNH4")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNH4" , y="RR39")
p21
pdf("RRNH4.pdf",width=8,height=8)
p21
dev.off() 

## RRAP
sum(!is.na(FungalShannon$RRAP)) ## n = 109
p22 <- ggplot(FungalShannon, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="RRAP")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAP" , y="RR109")
p22
pdf("RRAP.pdf",width=8,height=8)
p22
dev.off() 

## RRAK
sum(!is.na(FungalShannon$RRAK)) ## n = 87
p23 <- ggplot(FungalShannon, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="RRAK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAK" , y="RR87")
p23
pdf("RRAK.pdf",width=8,height=8)
p23
dev.off() 

## RRAN
sum(!is.na(FungalShannon$RRAN)) ## n = 72
p24 <- ggplot(FungalShannon, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalShannon", x="RRAN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAN" , y="RR72")
p24
pdf("RRAN.pdf",width=8,height=8)
p24
dev.off() 

## RRYield
sum(!is.na(FungalShannon$RRYield)) ## n = 66
p25 <- ggplot(FungalShannon, aes(x=RR, y=RRYield)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm, color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") +
  theme_bw() +
  theme(text = element_text(family = "serif", size=20)) +
  geom_vline(aes(xintercept=0), colour="black", linewidth=0.5, linetype="dashed") +
  labs(x="RR", y="RRYield66") +
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
  scale_x_continuous(limits=c(-0.7, 0.58), expand=c(0, 0)) + 
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
                            data = FungalShannon, 
                            na.action = na.roughfix, 
                            importance = TRUE, 
                            ntree = 500)
# Check the importance of variables and p-values
importance(rf_model_perm3)
#                    %IncMSE %IncMSE.pval IncNodePurity IncNodePurity.pval
# MATmean          26.224042   0.00990099     1.0768070         0.00990099
# MAPmean          21.435218   0.00990099     0.9294968         0.01980198
# Longitude        15.923233   0.00990099     0.5938006         0.38613861
# Duration         15.895406   0.00990099     0.4064541         0.56435644
# Latitude         14.193159   0.02970297     0.5432386         0.41584158
# Rotation_cycles  14.009730   0.00990099     0.4332605         0.63366337
# pHCK             13.040117   0.01980198     0.4813620         0.78217822
# SOCCK             8.405968   0.05940594     0.7777300         0.28712871
# Species_Richness  4.401483   0.21782178     0.3027411         0.05940594

################################## Trials sorted by effect size
library(ggplot2)
library(dplyr)
FungalShannon <- read.csv("FungalShannon.csv", fileEncoding = "latin1")
# 计算95% CI + 显著性分类
df_plot <- FungalShannon %>%
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
#       40      173       47 
## 8*8



################################################ piecewiseSEM
sem_data <- read.csv("FungalShannon.csv", fileEncoding = "latin1")
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
# -405.693
# Fisher's C = 2.883 with P-value = 0.237

m1_1 <- lme(RR ~ RRpH + RRSOC + Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_2 <- lme(RRpH ~  Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_3 <- lme(RRSOC ~ Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model1 <- psem(m1_1, m1_2, m1_3)
summary(sem_model1)  ## 
# AIC
# -431.139
# Fisher's C = 4.233 with P-value = 0.12

m2_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_2 <- lme(RRpH ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_3 <- lme(RRSOC ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model2 <- psem(m2_1, m2_2, m2_3)
summary(sem_model2)  ## 
# AIC
# -443.162
# Fisher's C = 2.892 with P-value = 0.236

m3_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_2 <- lme(RRpH ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_3 <- lme(RRSOC ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model3 <- psem(m3_1, m3_2, m3_3)
summary(sem_model3)  ## 
# AIC
# -427.909
# Fisher's C = 2.226 with P-value = 0.329

m4_1 <- lme(RR ~ RRpH + RRSOC + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_2 <- lme(RRpH ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_3 <- lme(RRSOC ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model4 <- psem(m4_1, m4_2, m4_3)
summary(sem_model4) ## 
# AIC
# -468.973
# Fisher's C = 4.255 with P-value = 0.119

# Structural Equation Model of sem_model4 
# 
# Call:
#   RR ~ RRpH + RRSOC + Species_Richness
# RRpH ~ Species_Richness
# RRSOC ~ Species_Richness
# 
# AIC
# -468.973
# 
# ---
#   Tests of directed separation:
#   
#   Independ.Claim Test.Type DF Crit.Value P.Value 
# RRSOC ~ RRpH + ...      coef 72     1.5773  0.1191 
# 
# --
#   Global goodness-of-fit:
#   
#   Chi-Squared = NA with P-value = NA and on 1 degrees of freedom
# Fisher's C = 4.255 with P-value = 0.119 and on 2 degrees of freedom
# 
# ---
# Coefficients:
# 
#   Response        Predictor Estimate Std.Error DF Crit.Value P.Value Std.Estimate  
#         RR             RRpH   0.3973    0.3227 71     1.2312  0.2223       0.1411  
#         RR            RRSOC   0.1995    0.0943 71     2.1151  0.0379       0.2478 *
#         RR Species_Richness  -0.0238    0.0281 71    -0.8482  0.3992      -0.0774  
#       RRpH Species_Richness  -0.0116    0.0092 73    -1.2601  0.2117      -0.1064  
#      RRSOC Species_Richness   0.0282    0.0300 84     0.9374  0.3512       0.0737  
# 
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05
# 
# ---
# Individual R-squared:
# 
#   Response method Marginal Conditional
#         RR   none     0.08        0.28
#       RRpH   none     0.01        0.59
#      RRSOC   none     0.00        0.75

m5_1 <- lme(RR ~ RRpH + RRSOC + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_2 <- lme(RRpH ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_3 <- lme(RRSOC ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model5 <- psem(m5_1, m5_2, m5_3)
summary(sem_model5) ## 
# AIC
# -452.395
# Fisher's C = 3.434 with P-value = 0.18

m6_1 <- lme(RR ~ RRpH + RRSOC +MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_2 <- lme(RRpH ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_3 <- lme(RRSOC ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model6 <- psem(m6_1, m6_2, m6_3)
summary(sem_model6)  ##
# AIC
# -464.842
# Fisher's C = 1.986 with P-value = 0.37

