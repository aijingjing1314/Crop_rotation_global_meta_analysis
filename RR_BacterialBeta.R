
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)

BacterialBeta <- read.csv("BacterialBeta.csv", fileEncoding = "latin1")
# Check data
head(BacterialBeta)

# 1. The number of Obversation
total_number <- nrow(BacterialBeta)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 343

# 2. The number of Study
unique_studyid_number <- length(unique(BacterialBeta$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 81


#### 3. Overall effect size
total_effect_model <- rma.mv(yi = RR, 
                             V = Vi, 
                             random = ~ 1 | StudyID,  # StudyID is radom factor
                             data = BacterialBeta, 
                             method = "REML")
# The results of Overall effect size
summary(total_effect_model)
# Multivariate Meta-Analysis Model (k = 343; method: REML
#    logLik   Deviance        AIC        BIC       AICc   
# -927.3829  1854.7658  1858.7658  1866.4354  1858.8012   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5605  0.7487     81     no  StudyID 
# Test for Heterogeneity:
# Q(df = 342) = 3202.5814, p-val < .0001
# Model Results:
# estimate      se     zval    pval    ci.lb   ci.ub    
 # -0.0796  0.0864  -0.9215  0.3568  -0.2490  0.0898

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
boot_results1 <- boot(data = BacterialBeta, statistic = boot_fun, R = 1000, parallel = "snow", ncpus = numCores, cl = cl)
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
# Estimate for coefficient 1 : -0.0796398 
# 95% BCa CI for coefficient 1 : -0.1841696 0.01992125 


#### 5. Funnel Plot
simple_model <- rma(yi = RR, 
                    vi = Vi, 
                    data = BacterialBeta, 
                    method = "REML")
#### 
funnel(simple_model)
# Output  6 * 6
#### Egger's test
regtest(simple_model)
# Regression Test for Funnel Plot Asymmetry
# Model:     mixed-effects meta-regression model
# Predictor: standard error
# Test for Funnel Plot Asymmetry: z =  0.7967, p = 0.4256
# Limit Estimate (as sei -> 0):   b = -0.0866 (CI: -0.3250, 0.1518)

#  Rosenthal’s Fail-Safe N
# This method estimates how many missing studies with null effect 
# would be needed to make the overall effect non-significant
fsn_rosenthal <- fsn(x = simple_model, type = "Rosenthal")
# Print the FSN result
print(fsn_rosenthal)
# Fail-safe N Calculation Using the General Approach
# Average Effect Size:         0.0041
# Amount of Heterogeneity:     0.5429
# Observed Significance Level: 0.9244
# Target Significance Level:   0.05
# Fail-safe N: 0


#### 8. Subgroup analysis
### 8.1 LegumeNonlegume
BacterialBeta_filteredLegumeNonlegume <- subset(BacterialBeta, LegumeNonlegume %in% c("Legume to Non-legume", "Non-legume to Legume", "Non-legume to Non-legume"))
#
BacterialBeta_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(BacterialBeta_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredLegumeNonlegume %>%
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

overall_model_BacterialBeta_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = BacterialBeta_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 340; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -915.9787  1831.9573  1839.9573  1855.2376  1840.0778   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5848  0.7647     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 337) = 3161.7587, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 16.1028, p-val = 0.0011
# Model Results:
#                                          estimate      se     zval    pval    ci.lb    ci.ub     
# LegumeNonlegumeLegume to Non-legume       -0.3079  0.1076  -2.8620  0.0042  -0.5188  -0.0971  ** 
# LegumeNonlegumeNon-legume to Legume       -0.1771  0.0989  -1.7900  0.0735  -0.3709   0.0168   . 
# LegumeNonlegumeNon-legume to Non-legume    0.0699  0.0973   0.7187  0.4723  -0.1207   0.2605    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredLegumeNonlegume)
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
BacterialBeta_filteredAMnonAM <- subset(BacterialBeta, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
BacterialBeta_filteredAMnonAM$AMnonAM <- droplevels(factor(BacterialBeta_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredAMnonAM %>%
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

overall_model_BacterialBeta_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = BacterialBeta_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 342; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -924.2064  1848.4128  1856.4128  1871.7168  1856.5325   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5763  0.7591     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 339) = 3162.7693, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.7330, p-val = 0.4346
# Model Results:
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# AMnonAMAM to AM      -0.1180  0.0930  -1.2689  0.2045  -0.3002  0.0643    
# AMnonAMAM to nonAM   -0.0301  0.1158  -0.2601  0.7948  -0.2570  0.1968    
# AMnonAMnonAM to AM    0.1581  0.2984   0.5298  0.5963  -0.4268  0.7430      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredAMnonAM)
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
BacterialBeta_filteredC3C4 <- subset(BacterialBeta, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
BacterialBeta_filteredC3C4$C3C4 <- droplevels(factor(BacterialBeta_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredC3C4 %>%
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
overall_model_BacterialBeta_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = BacterialBeta_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 334; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -830.2398  1660.4797  1668.4797  1683.6881  1668.6024   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.7170  0.8468     79     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 331) = 2976.1610, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 111.8097, p-val < .0001
# Model Results:
#               estimate      se     zval    pval    ci.lb    ci.ub      
# C3C4C3 to C3    0.2402  0.1038   2.3139  0.0207   0.0367   0.4437    * 
# C3C4C3 to C4   -0.5298  0.1068  -4.9612  <.0001  -0.7391  -0.3205  *** 
# C3C4C4 to C3   -0.3526  0.1158  -3.0451  0.0023  -0.5796  -0.1257   **   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredC3C4)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredC3C4)
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
#          "a"          "b"          "c" 


### 8.4 Annual_Pere
BacterialBeta_filteredAnnual_Pere <- subset(BacterialBeta, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
BacterialBeta_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(BacterialBeta_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredAnnual_Pere %>%
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

overall_model_BacterialBeta_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = BacterialBeta_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 342; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -881.0362  1762.0723  1772.0723  1791.1875  1772.2530   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.6308  0.7943     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 338) = 3019.2707, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 88.8811, p-val < .0001
# Model Results:
#                                    estimate      se     zval    pval    ci.lb    ci.ub      
# Annual_PereAnnual to Annual         -0.0165  0.1039  -0.1585  0.8740  -0.2200   0.1871      
# Annual_PereAnnual to Perennial      -0.8154  0.1290  -6.3192  <.0001  -1.0683  -0.5625  *** 
# Annual_PerePerennial to Annual       0.2048  0.1884   1.0871  0.2770  -0.1644   0.5740      
# Annual_PerePerennial to Perennial   -0.3230  0.2193  -1.4729  0.1408  -0.7527   0.1068      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredAnnual_Pere)
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
      #                        "ab"                               "c"                               "a"                               "b" 


### 8.6 PlantStageSubgroup
BacterialBeta_filteredPlantStageSubgroup <- subset(BacterialBeta, PlantStageSubgroup %in% c("Vegetative stage","Reproductive stage", "Maturity stage","Harvest"))
#
BacterialBeta_filteredPlantStageSubgroup$PlantStageSubgroup <- droplevels(factor(BacterialBeta_filteredPlantStageSubgroup$PlantStageSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredPlantStageSubgroup %>%
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

overall_model_BacterialBeta_filteredPlantStageSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + PlantStageSubgroup, random = ~ 1 | StudyID, data = BacterialBeta_filteredPlantStageSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredPlantStageSubgroup)
# Multivariate Meta-Analysis Model (k = 278; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -669.3129  1338.6257  1348.6257  1366.6914  1348.8496   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5066  0.7117     67     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 274) = 2245.8851, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 5.2087, p-val = 0.2665
# Model Results:
#                                       estimate      se     zval    pval    ci.lb   ci.ub    
# PlantStageSubgroupHarvest              -0.1610  0.0957  -1.6833  0.0923  -0.3486  0.0265  . 
# PlantStageSubgroupMaturity stage       -0.1735  0.1159  -1.4967  0.1345  -0.4006  0.0537    
# PlantStageSubgroupReproductive stage   -0.0552  0.1274  -0.4334  0.6647  -0.3050  0.1945    
# PlantStageSubgroupVegetative stage     -0.0313  0.1180  -0.2652  0.7908  -0.2625  0.1999    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredPlantStageSubgroup)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredPlantStageSubgroup)
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
           #                       "a"                                  "a"                                  "a"                                  "a" 

### 8.7 Bulk_Rhizosphere
BacterialBeta_filteredBulk_Rhizosphere <- subset(BacterialBeta, Bulk_Rhizosphere %in% c("Non-Rhizosphere", "Rhizosphere"))
#
BacterialBeta_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(BacterialBeta_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredBulk_Rhizosphere %>%
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

overall_model_BacterialBeta_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = BacterialBeta_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -925.0505  1850.1011  1856.1011  1867.5967  1856.1723   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5637  0.7508     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 341) = 3201.7087, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 6.0419, p-val = 0.0488
# Model Results:
#                                  estimate      se     zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon-Rhizosphere   -0.1113  0.0878  -1.2682  0.2047  -0.2833  0.0607    
# Bulk_RhizosphereRhizosphere       -0.0324  0.0891  -0.3633  0.7164  -0.2070  0.1423    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredBulk_Rhizosphere)
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
BacterialBeta_filteredSoil_texture <- subset(BacterialBeta, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
BacterialBeta_filteredSoil_texture$Soil_texture <- droplevels(factor(BacterialBeta_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredSoil_texture %>%
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

overall_model_BacterialBeta_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = BacterialBeta_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 256; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -735.0456  1470.0913  1478.0913  1492.2248  1478.2526   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4256  0.6524     57     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 253) = 2317.8774, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.1965, p-val = 0.5326
# Model Results:
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse   -0.1590  0.1884  -0.8436  0.3989  -0.5283  0.2103    
# Soil_textureFine     -0.2093  0.1722  -1.2153  0.2243  -0.5467  0.1282    
# Soil_textureMedium    0.0115  0.1291   0.0891  0.9290  -0.2416  0.2646    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredSoil_texture)
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
BacterialBeta_filteredTillage <- subset(BacterialBeta, Tillage %in% c("Tillage", "No_tillage"))
#
BacterialBeta_filteredTillage$Tillage <- droplevels(factor(BacterialBeta_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredTillage %>%
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
overall_model_BacterialBeta_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = BacterialBeta_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredTillage)
# Multivariate Meta-Analysis Model (k = 54; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -95.3779  190.7557  196.7557  202.6095  197.2557   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1650  0.4062     11     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 52) = 251.3093, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 7.9055, p-val = 0.0192
# Model Results:
#                    estimate      se     zval    pval    ci.lb    ci.ub    
# TillageNo_tillage    0.1694  0.2018   0.8397  0.4011  -0.2260   0.5649    
# TillageTillage      -0.3139  0.1402  -2.2391  0.0252  -0.5886  -0.0391  * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredTillage)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredTillage)
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
BacterialBeta_filteredStraw_retention <- subset(BacterialBeta, Straw_retention %in% c("Retention", "No_retention"))
#
BacterialBeta_filteredStraw_retention$Straw_retention <- droplevels(factor(BacterialBeta_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredStraw_retention %>%
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

overall_model_BacterialBeta_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = BacterialBeta_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 47; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -75.8183  151.6366  157.6366  163.0566  158.2220   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5668  0.7528     16     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 45) = 354.3042, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 8.4513, p-val = 0.0146
# Model Results:
#                              estimate      se     zval    pval    ci.lb    ci.ub    
# Straw_retentionNo_retention   -0.4413  0.2046  -2.1573  0.0310  -0.8423  -0.0404  * 
# Straw_retentionRetention      -0.1287  0.2052  -0.6274  0.5304  -0.5309   0.2734      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredStraw_retention)
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
BacterialBeta_filteredRotationcyclesSubgroup <- subset(BacterialBeta, RotationcyclesSubgroup %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
BacterialBeta_filteredRotationcyclesSubgroup$RotationcyclesSubgroup <- droplevels(factor(BacterialBeta_filteredRotationcyclesSubgroup$RotationcyclesSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredRotationcyclesSubgroup %>%
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
overall_model_BacterialBeta_filteredRotationcyclesSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + RotationcyclesSubgroup, random = ~ 1 | StudyID, data = BacterialBeta_filteredRotationcyclesSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredRotationcyclesSubgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -906.9661  1813.9323  1825.9323  1848.8706  1826.1861   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.6160  0.7849     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 338) = 3166.3285, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 38.8446, p-val < .0001
# Model Results:
#                              estimate      se     zval    pval    ci.lb    ci.ub     
# RotationcyclesSubgroupD1       0.0407  0.1089   0.3734  0.7088  -0.1728   0.2542     
# RotationcyclesSubgroupD1-3    -0.2095  0.1100  -1.9054  0.0567  -0.4251   0.0060   . 
# RotationcyclesSubgroupD10     -0.1138  0.1330  -0.8555  0.3923  -0.3746   0.1469     
# RotationcyclesSubgroupD3-5    -0.3366  0.1251  -2.6896  0.0072  -0.5818  -0.0913  ** 
# RotationcyclesSubgroupD5-10    0.0414  0.1218   0.3396  0.7341  -0.1974   0.2802   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredRotationcyclesSubgroup)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredRotationcyclesSubgroup)
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
   #                      "a"                        "bc"                        "ab"                         "c"                         "a" 


### 8.13 DurationSubgroup
BacterialBeta_filteredDurationSubgroup <- subset(BacterialBeta, DurationSubgroup %in% c("D1", "D2", "D3", "D4", "D5", "D6-10", "D11-20", "D20-30", "D30"))
#
BacterialBeta_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(BacterialBeta_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredDurationSubgroup %>%
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

overall_model_BacterialBeta_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = BacterialBeta_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -899.5079  1799.0157  1819.0157  1857.1271  1819.6968   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.6061  0.7785     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 334) = 3051.7158, p-val < .0001
# Test of Moderators (coefficients 1:9):
# QM(df = 9) = 44.7577, p-val < .0001
# Model Results:
#                         estimate      se     zval    pval    ci.lb    ci.ub     
# DurationSubgroupD1        0.0011  0.2090   0.0053  0.9958  -0.4086   0.4108     
# DurationSubgroupD11-20   -0.3319  0.1362  -2.4363  0.0148  -0.5989  -0.0649   * 
# DurationSubgroupD2        0.1805  0.1378   1.3104  0.1901  -0.0895   0.4505     
# DurationSubgroupD20-30   -0.2851  0.1593  -1.7895  0.0735  -0.5973   0.0272   . 
# DurationSubgroupD3       -0.4252  0.1319  -3.2239  0.0013  -0.6837  -0.1667  ** 
# DurationSubgroupD30       0.0742  0.3966   0.1872  0.8515  -0.7031   0.8516     
# DurationSubgroupD4       -0.0745  0.2553  -0.2917  0.7705  -0.5749   0.4259     
# DurationSubgroupD5        0.3426  0.4214   0.8129  0.4163  -0.4834   1.1686     
# DurationSubgroupD6-10     0.0602  0.1550   0.3884  0.6977  -0.2435   0.3639   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredDurationSubgroup)
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
   #               "abc"                   "ad"                    "b"                  "acd"                    "d"                 "abcd" 
   #  DurationSubgroupD4     DurationSubgroupD5  DurationSubgroupD6-10 
   #              "abcd"                 "abcd"                   "bc" 



### 8.14 SpeciesRichnessSubgroup
BacterialBeta_filteredSpeciesRichnessSubgroup <- subset(BacterialBeta, SpeciesRichnessSubgroup %in% c("R2", "R3", "R4"))
#
BacterialBeta_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup <- droplevels(factor(BacterialBeta_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredSpeciesRichnessSubgroup %>%
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

overall_model_BacterialBeta_filteredSpeciesRichnessSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichnessSubgroup, random = ~ 1 | StudyID, data = BacterialBeta_filteredSpeciesRichnessSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredSpeciesRichnessSubgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -877.0342  1754.0684  1762.0684  1777.3842  1762.1878   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5717  0.7561     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 340) = 3085.7670, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 98.4309, p-val < .0001
# Model Results:
#                            estimate      se     zval    pval    ci.lb    ci.ub      
# SpeciesRichnessSubgroupR2    0.0260  0.0895   0.2908  0.7712  -0.1494   0.2014      
# SpeciesRichnessSubgroupR3   -0.6765  0.1083  -6.2493  <.0001  -0.8887  -0.4643  *** 
# SpeciesRichnessSubgroupR4   -0.3546  0.4432  -0.8002  0.4236  -1.2232   0.5140      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredSpeciesRichnessSubgroup)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredSpeciesRichnessSubgroup)
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
BacterialBeta_filteredPrimer <- subset(BacterialBeta, Primer %in% c("V1-V3", "V1-V4", "V1-V5", "V3-V4", "V3-V5", "V4", "V4-V5", "V5-V7", "V6-V8", "V9", "Full length"))
#
BacterialBeta_filteredPrimer$Primer <- droplevels(factor(BacterialBeta_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredPrimer %>%
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

overall_model_BacterialBeta_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = BacterialBeta_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 337; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -904.6002  1809.2004  1827.2004  1861.3649  1827.7647   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5106  0.7146     78     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 329) = 2780.3135, p-val < .0001
# Test of Moderators (coefficients 1:8):
# QM(df = 8) = 7.0530, p-val = 0.5309
# Model Results:
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerFull length    0.0582  0.7384   0.0788  0.9372  -1.3890  1.5054    
# PrimerV1-V4         -0.7578  0.7287  -1.0400  0.2984  -2.1860  0.6704    
# PrimerV1-V5         -1.3061  0.7195  -1.8152  0.0695  -2.7164  0.1042  . 
# PrimerV3-V4          0.0779  0.1227   0.6350  0.5254  -0.1626  0.3185    
# PrimerV3-V5          0.4106  0.7234   0.5676  0.5703  -1.0073  1.8286    
# PrimerV4            -0.1023  0.1594  -0.6416  0.5211  -0.4146  0.2101    
# PrimerV4-V5         -0.1950  0.1988  -0.9813  0.3264  -0.5846  0.1945    
# PrimerV5-V7         -0.5450  0.7217  -0.7551  0.4502  -1.9596  0.8696     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredPrimer)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredPrimer)
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
BacterialBeta_filteredLatitude_Subgroup <- subset(BacterialBeta, Latitude_Subgroup %in% c("La20", "La20-40", "La40"))
#
BacterialBeta_filteredLatitude_Subgroup$Latitude_Subgroup <- droplevels(factor(BacterialBeta_filteredLatitude_Subgroup$Latitude_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredLatitude_Subgroup %>%
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

overall_model_BacterialBeta_filteredLatitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Latitude_Subgroup, random = ~ 1 | StudyID, data = BacterialBeta_filteredLatitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredLatitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -918.6398  1837.2797  1845.2797  1860.5955  1845.3991   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4899  0.7000     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 340) = 2983.5141, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 13.6533, p-val = 0.0034
# Model Results:
#                           estimate      se     zval    pval    ci.lb    ci.ub      
# Latitude_SubgroupLa20      -1.0703  0.3015  -3.5506  0.0004  -1.6612  -0.4795  *** 
# Latitude_SubgroupLa20-40    0.0730  0.1120   0.6517  0.5146  -0.1465   0.2924      
# Latitude_SubgroupLa40      -0.1011  0.1282  -0.7887  0.4303  -0.3523   0.1501      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredLatitude_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredLatitude_Subgroup)
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
   #                   "a"                      "b"                      "b"  


### 8.17 Longitude_Subgroup
BacterialBeta_filteredLongitude_Subgroup <- subset(BacterialBeta, Longitude_Subgroup %in% c("Lo-180-0", "Lo-0-180"))
#
BacterialBeta_filteredLongitude_Subgroup$Longitude_Subgroup <- droplevels(factor(BacterialBeta_filteredLongitude_Subgroup$Longitude_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredLongitude_Subgroup %>%
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
overall_model_BacterialBeta_filteredLongitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Subgroup, random = ~ 1 | StudyID, data = BacterialBeta_filteredLongitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredLongitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -925.2382  1850.4763  1856.4763  1867.9720  1856.5475   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5674  0.7532     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 341) = 3155.0461, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.2880, p-val = 0.5252
# Model Results:
#                             estimate      se     zval    pval    ci.lb   ci.ub    
# Longitude_SubgroupLo-0-180   -0.0614  0.0911  -0.6747  0.4999  -0.2399  0.1170    
# Longitude_SubgroupLo-180-0   -0.2661  0.2916  -0.9126  0.3615  -0.8378  0.3055         

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredLongitude_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredLongitude_Subgroup)
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
BacterialBeta_filteredMAPmean_Subgroup <- subset(BacterialBeta, MAPmean_Subgroup %in% c("MAP600", "MAP600-1200", "MAP1200"))
#
BacterialBeta_filteredMAPmean_Subgroup$MAPmean_Subgroup <- droplevels(factor(BacterialBeta_filteredMAPmean_Subgroup$MAPmean_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredMAPmean_Subgroup %>%
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

overall_model_BacterialBeta_filteredMAPmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean_Subgroup, random = ~ 1 | StudyID, data = BacterialBeta_filteredMAPmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredMAPmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -923.2012  1846.4025  1854.4025  1869.7183  1854.5219   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5412  0.7356     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 340) = 3193.9567, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 4.8395, p-val = 0.1839
# Model Results:
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean_SubgroupMAP1200       -0.0277  0.1546  -0.1793  0.8577  -0.3307  0.2753    
# MAPmean_SubgroupMAP600         0.0050  0.1157   0.0429  0.9658  -0.2218  0.2318    
# MAPmean_SubgroupMAP600-1200   -0.2431  0.1248  -1.9470  0.0515  -0.4877  0.0016  . 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredMAPmean_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredMAPmean_Subgroup)
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
BacterialBeta_filteredMATmean_Subgroup <- subset(BacterialBeta, MATmean_Subgroup %in% c("MAT8", "MAT8-15", "MAT15"))
#
BacterialBeta_filteredMATmean_Subgroup$MATmean_Subgroup <- droplevels(factor(BacterialBeta_filteredMATmean_Subgroup$MATmean_Subgroup))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredMATmean_Subgroup %>%
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

overall_model_BacterialBeta_filteredMATmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Subgroup, random = ~ 1 | StudyID, data = BacterialBeta_filteredMATmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredMATmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -924.8198  1849.6396  1857.6396  1872.9554  1857.7590   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5715  0.7559     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 340) = 3164.5189, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.5194, p-val = 0.6778
# Model Results:
#                          estimate      se     zval    pval    ci.lb   ci.ub    
# MATmean_SubgroupMAT15     -0.0958  0.1354  -0.7076  0.4792  -0.3611  0.1695    
# MATmean_SubgroupMAT8      -0.0366  0.1204  -0.3040  0.7611  -0.2727  0.1994    
# MATmean_SubgroupMAT8-15   -0.1476  0.1501  -0.9830  0.3256  -0.4419  0.1467    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredMATmean_Subgroup)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredMATmean_Subgroup)
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
BacterialBeta_filteredDurationSubgroup <- subset(BacterialBeta, DurationSubgroup %in% c("D5", "D5-10", "D10-20", "D20-30", "D30"))
#
BacterialBeta_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(BacterialBeta_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- BacterialBeta_filteredDurationSubgroup %>%
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

overall_model_BacterialBeta_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = BacterialBeta_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_BacterialBeta_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 343; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -921.9530  1843.9059  1855.9059  1878.8442  1856.1597   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.5655  0.7520     81     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 338) = 3114.0856, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 3.8703, p-val = 0.5682
# Model Results:
#                         estimate      se     zval    pval    ci.lb   ci.ub    
# DurationSubgroupD10-20   -0.1804  0.1283  -1.4059  0.1597  -0.4320  0.0711    
# DurationSubgroupD20-30   -0.1742  0.1543  -1.1291  0.2589  -0.4766  0.1282    
# DurationSubgroupD30       0.0752  0.3836   0.1961  0.8445  -0.6765  0.8270    
# DurationSubgroupD5       -0.1198  0.1052  -1.1395  0.2545  -0.3259  0.0863    
# DurationSubgroupD5-10     0.0844  0.1488   0.5675  0.5703  -0.2071  0.3760  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_BacterialBeta_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_BacterialBeta_filteredDurationSubgroup)
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
#                   "a"                    "a"                    "a"                    "a"                    "a" 







#### 9. Linear Mixed Effect Model
# 
BacterialBeta$Wr <- 1 / BacterialBeta$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialBeta)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialBeta)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialBeta)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialBeta)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialBeta)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialBeta)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = BacterialBeta)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = BacterialBeta)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(Species_Richness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = BacterialBeta)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = BacterialBeta)

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
# Model1   Model1 817.5524 840.5788 -402.7762 0.005573389     0.02720555
# Model2   Model2 816.7180 839.7444 -402.3590 0.005700507     0.02641559
# Model3   Model3 817.7374 840.7638 -402.8687 0.005510675     0.02696590
# Model4   Model4 816.6949 839.7213 -402.3474 0.005862452     0.02750298
# Model5   Model5 817.9495 840.9759 -402.9747 0.005581130     0.02698276
# Model6   Model6 817.0861 840.1125 -402.5431 0.005864939     0.02727303
# Model7   Model7 817.6212 840.6476 -402.8106 0.005403633     0.02617819
# Model8   Model8 816.8648 839.8912 -402.4324 0.005819617     0.02724462
# Model9   Model9 823.0564 861.4337 -401.5282 0.005241529     0.02847913
# Model10 Model10 814.6902 853.0675 -397.3451 0.007511718     0.02844189

##### Model 10 is the best model
summary(Model10)
# Number of obs: 343, groups:  StudyID, 81
anova(Model10) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                                                                                Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# scale(log(Rotation_cycles))                                                     0.048   0.048     1 118.455  0.0073 0.9318353    
# scale(log(Species_Richness))                                                  141.228 141.228     1 101.313 21.4633 1.075e-05 ***
# scale(log(Duration))                                                            9.584   9.584     1  66.775  1.4566 0.2317357    
# scale(log(Rotation_cycles)):scale(log(Species_Richness))                       60.578  60.578     1  82.316  9.2064 0.0032267 ** 
# scale(log(Rotation_cycles)):scale(log(Duration))                               10.220  10.220     1  85.639  1.5531 0.2160741    
# scale(log(Species_Richness)):scale(log(Duration))                             104.334 104.334     1 154.854 15.8563 0.0001048 ***
# scale(log(Rotation_cycles)):scale(log(Species_Richness)):scale(log(Duration))  39.299  39.299     1 161.332  5.9725 0.0156082 *  


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) +                    scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + scale(pHCK) + (1 | StudyID), weights = Wr, data = BacterialBeta)
summary(ModelpH)
# Number of obs: 142, groups:  StudyID, 55
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                                                                                Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles))                                                    1.8126  1.8126     1  56.769  0.2717 0.6042
# scale(log(Species_Richness))                                                   1.8529  1.8529     1  96.004  0.2777 0.5994
# scale(log(Duration))                                                           7.6025  7.6025     1 120.907  1.1395 0.2879
# scale(pHCK)                                                                    2.2750  2.2750     1  32.391  0.3410 0.5633
# scale(log(Rotation_cycles)):scale(log(Species_Richness))                       0.1016  0.1016     1  87.705  0.0152 0.9020
# scale(log(Rotation_cycles)):scale(log(Duration))                              14.4775 14.4775     1 124.333  2.1699 0.1433
# scale(log(Species_Richness)):scale(log(Duration))                             10.0589 10.0589     1 132.672  1.5076 0.2217
# scale(log(Rotation_cycles)):scale(log(Species_Richness)):scale(log(Duration)) 10.9927 10.9927     1 132.962  1.6476 0.2015

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) +                    scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = BacterialBeta)
summary(ModelSOC)
# Number of obs: 148, groups:  StudyID, 54
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                                                                                Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))                                                    1.0098  1.0098     1  55.695  0.1632 0.68774  
# scale(log(Species_Richness))                                                   3.5842  3.5842     1  69.894  0.5794 0.44910  
# scale(log(Duration))                                                           0.4131  0.4131     1  51.805  0.0668 0.79711  
# scale(SOCCK)                                                                  28.3832 28.3832     1  26.705  4.5886 0.04147 *
# scale(log(Rotation_cycles)):scale(log(Species_Richness))                       1.6756  1.6756     1  63.139  0.2709 0.60456  
# scale(log(Rotation_cycles)):scale(log(Duration))                              15.6511 15.6511     1  80.424  2.5302 0.11561  
# scale(log(Species_Richness)):scale(log(Duration))                              0.0495  0.0495     1  53.789  0.0080 0.92906  
# scale(log(Rotation_cycles)):scale(log(Species_Richness)):scale(log(Duration))  0.1765  0.1765     1 129.835  0.0285 0.86614  

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) +                    scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + scale(TNCK) + (1 | StudyID), weights = Wr, data = BacterialBeta)
summary(ModelTN)
# Number of obs: 97, groups:  StudyID, 37
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                                                                               Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles))                                                   2.7310  2.7310     1 40.475  0.3653 0.5489
# scale(log(Species_Richness))                                                  3.2945  3.2945     1 33.327  0.4407 0.5114
# scale(log(Duration))                                                          0.0998  0.0998     1 57.356  0.0133 0.9084
# scale(TNCK)                                                                   1.9593  1.9593     1 12.070  0.2621 0.6179
# scale(log(Rotation_cycles)):scale(log(Species_Richness))                      0.2207  0.2207     1 42.960  0.0295 0.8644
# scale(log(Rotation_cycles)):scale(log(Duration))                              9.0863  9.0863     1 56.681  1.2154 0.2749
# scale(log(Species_Richness)):scale(log(Duration))                             0.4001  0.4001     1 86.868  0.0535 0.8176
# scale(log(Rotation_cycles)):scale(log(Species_Richness)):scale(log(Duration)) 0.0382  0.0382     1 84.130  0.0051 0.9432


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) +                    scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = BacterialBeta)
summary(ModelNO3)
# Number of obs: 50, groups:  StudyID, 19
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                                                                               Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles))                                                   0.0066  0.0066     1 26.227  0.0011 0.9736
# scale(log(Species_Richness))                                                  0.0300  0.0300     1 19.242  0.0051 0.9440
# scale(log(Duration))                                                          0.2440  0.2440     1 33.683  0.0412 0.8405
# scale(NO3CK)                                                                  1.5710  1.5710     1 18.533  0.2650 0.6128
# scale(log(Rotation_cycles)):scale(log(Species_Richness))                      0.0770  0.0770     1 22.882  0.0130 0.9102
# scale(log(Rotation_cycles)):scale(log(Duration))                              7.9196  7.9196     1 36.253  1.3357 0.2554
# scale(log(Species_Richness)):scale(log(Duration))                             0.2172  0.2172     1 40.331  0.0366 0.8492
# scale(log(Rotation_cycles)):scale(log(Species_Richness)):scale(log(Duration)) 1.8237  1.8237     1 40.500  0.3076 0.5822

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) +                    scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = BacterialBeta)
summary(ModelNH4)
# Number of obs: 46, groups:  StudyID, 17
anova(ModelNH4) 
# ype III Analysis of Variance Table with Satterthwaite's method
#                                                                               Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles))                                                   0.1360  0.1360     1  8.054  0.0213 0.8877
# scale(log(Species_Richness))                                                  0.3085  0.3085     1 11.471  0.0482 0.8300
# scale(log(Duration))                                                          0.1198  0.1198     1 12.661  0.0187 0.8933
# scale(NH4CK)                                                                  2.5675  2.5675     1  4.615  0.4013 0.5565
# scale(log(Rotation_cycles)):scale(log(Species_Richness))                      0.0195  0.0195     1 11.670  0.0031 0.9569
# scale(log(Rotation_cycles)):scale(log(Duration))                              4.7656  4.7656     1 24.343  0.7449 0.3965
# scale(log(Species_Richness)):scale(log(Duration))                             0.0679  0.0679     1 25.644  0.0106 0.9187
# scale(log(Rotation_cycles)):scale(log(Species_Richness)):scale(log(Duration)) 2.8193  2.8193     1 36.839  0.4407 0.5109

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) +                    scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + scale(APCK) + (1 | StudyID), weights = Wr, data = BacterialBeta)
summary(ModelAP)
# Number of obs: 120, groups:  StudyID, 48
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                                                                               Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles))                                                   0.2620  0.2620     1 56.963  0.0359 0.8503
# scale(log(Species_Richness))                                                  7.1148  7.1148     1 58.054  0.9759 0.3273
# scale(log(Duration))                                                          3.1311  3.1311     1 64.512  0.4295 0.5146
# scale(APCK)                                                                   0.0656  0.0656     1 65.590  0.0090 0.9247
# scale(log(Rotation_cycles)):scale(log(Species_Richness))                      0.6537  0.6537     1 73.622  0.0897 0.7654
# scale(log(Rotation_cycles)):scale(log(Duration))                              4.0463  4.0463     1 62.307  0.5550 0.4591
# scale(log(Species_Richness)):scale(log(Duration))                             2.3990  2.3990     1 62.926  0.3291 0.5683
# scale(log(Rotation_cycles)):scale(log(Species_Richness)):scale(log(Duration)) 2.3185  2.3185     1 63.107  0.3180 0.5748

#### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) +                    scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + scale(AKCK) + (1 | StudyID), weights = Wr, data = BacterialBeta)
summary(ModelAK)
# Number of obs: 99, groups:  StudyID, 40
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                                                                                Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles))                                                    0.0856  0.0856     1 24.738  0.0103 0.9200
# scale(log(Species_Richness))                                                  18.1841 18.1841     1 41.967  2.1860 0.1467
# scale(log(Duration))                                                           0.2701  0.2701     1 55.507  0.0325 0.8577
# scale(AKCK)                                                                    7.6471  7.6471     1 44.584  0.9193 0.3428
# scale(log(Rotation_cycles)):scale(log(Species_Richness))                       5.3723  5.3723     1 41.459  0.6458 0.4262
# scale(log(Rotation_cycles)):scale(log(Duration))                               0.1778  0.1778     1 44.590  0.0214 0.8844
# scale(log(Species_Richness)):scale(log(Duration))                              0.1492  0.1492     1 50.580  0.0179 0.8940
# scale(log(Rotation_cycles)):scale(log(Species_Richness)):scale(log(Duration))  0.9390  0.9390     1 44.305  0.1129 0.7385

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) +                    scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + scale(ANCK) + (1 | StudyID), weights = Wr, data = BacterialBeta)
summary(ModelAN)
# Number of obs: 80, groups:  StudyID, 25
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                                                                                Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)  
# scale(log(Rotation_cycles))                                                   22.4115 22.4115     1 42.276  4.1161 0.0488 *
# scale(log(Species_Richness))                                                   1.6707  1.6707     1 67.587  0.3069 0.5814  
# scale(log(Duration))                                                           9.3465  9.3465     1 21.607  1.7166 0.2039  
# scale(ANCK)                                                                    0.0929  0.0929     1 23.153  0.0171 0.8972  
# scale(log(Rotation_cycles)):scale(log(Species_Richness))                       2.8011  2.8011     1 71.138  0.5145 0.4756  
# scale(log(Rotation_cycles)):scale(log(Duration))                               4.5571  4.5571     1 36.353  0.8370 0.3663  
# scale(log(Species_Richness)):scale(log(Duration))                              5.6816  5.6816     1 68.864  1.0435 0.3106  

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) +                    scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + scale(Latitude) + (1 | StudyID), weights = Wr, data = BacterialBeta)
summary(ModelLatitude)
# Number of obs: 343, groups:  StudyID, 81
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                                                                                Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# scale(log(Rotation_cycles))                                                     0.512   0.512     1 139.153  0.0782 0.7801426    
# scale(log(Species_Richness))                                                  145.807 145.807     1 102.291 22.2946 7.458e-06 ***
# scale(log(Duration))                                                            2.768   2.768     1  84.648  0.4233 0.5170663    
# scale(Latitude)                                                                15.839  15.839     1  55.875  2.4219 0.1252970    
# scale(log(Rotation_cycles)):scale(log(Species_Richness))                       62.018  62.018     1  84.971  9.4828 0.0027934 ** 
# scale(log(Rotation_cycles)):scale(log(Duration))                                6.524   6.524     1  93.617  0.9975 0.3204862    
# scale(log(Species_Richness)):scale(log(Duration))                              89.011  89.011     1 157.275 13.6102 0.0003095 ***
# scale(log(Rotation_cycles)):scale(log(Species_Richness)):scale(log(Duration))  31.980  31.980     1 170.039  4.8899 0.0283484 *  

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) +                    scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + scale(Longitude) + (1 | StudyID), weights = Wr, data = BacterialBeta)
summary(ModelLongitude)
# Number of obs: 343, groups:  StudyID, 81
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                                                                                Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# scale(log(Rotation_cycles))                                                     0.089   0.089     1 119.239  0.0135 0.9075469    
# scale(log(Species_Richness))                                                  133.649 133.649     1  99.849 20.3505 1.762e-05 ***
# scale(log(Duration))                                                           10.241  10.241     1 103.499  1.5594 0.2145738    
# scale(Longitude)                                                                1.371   1.371     1  48.390  0.2088 0.6497292    
# scale(log(Rotation_cycles)):scale(log(Species_Richness))                       59.709  59.709     1  82.489  9.0918 0.0034117 ** 
# scale(log(Rotation_cycles)):scale(log(Duration))                               11.101  11.101     1 110.768  1.6904 0.1962470    
# scale(log(Species_Richness)):scale(log(Duration))                             104.036 104.036     1 162.378 15.8415 0.0001036 ***
# scale(log(Rotation_cycles)):scale(log(Species_Richness)):scale(log(Duration))  39.174  39.174     1 169.070  5.9649 0.0156225 *  

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) +                    scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = BacterialBeta)
summary(ModelMAPmean)
# Number of obs: 343, groups:  StudyID, 81
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                                                                                Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# scale(log(Rotation_cycles))                                                     0.032   0.032     1 123.618  0.0049 0.9443213    
# scale(log(Species_Richness))                                                  138.594 138.594     1 100.210 21.0900 1.276e-05 ***
# scale(log(Duration))                                                            8.479   8.479     1  68.404  1.2902 0.2599752    
# scale(MAPmean)                                                                  0.000   0.000     1  44.104  0.0000 0.9976088    
# scale(log(Rotation_cycles)):scale(log(Species_Richness))                       59.012  59.012     1  82.693  8.9799 0.0036025 ** 
# scale(log(Rotation_cycles)):scale(log(Duration))                                9.561   9.561     1  84.839  1.4550 0.2310871    
# scale(log(Species_Richness)):scale(log(Duration))                             103.249 103.249     1 153.385 15.7115 0.0001128 ***
# scale(log(Rotation_cycles)):scale(log(Species_Richness)):scale(log(Duration))  37.147  37.147     1 160.385  5.6527 0.0186069 *  

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) +                    scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + scale(MATmean) + (1 | StudyID), weights = Wr, data = BacterialBeta)
summary(ModelMATmean)
# Number of obs: 343, groups:  StudyID, 81
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                                                                                Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# scale(log(Rotation_cycles))                                                     0.226   0.226     1 132.442  0.0344  0.853114    
# scale(log(Species_Richness))                                                  136.281 136.281     1 103.640 20.7698 1.421e-05 ***
# scale(log(Duration))                                                            4.397   4.397     1  75.354  0.6701  0.415595    
# scale(MATmean)                                                                  7.482   7.482     1  43.501  1.1404  0.291468    
# scale(log(Rotation_cycles)):scale(log(Species_Richness))                       64.224  64.224     1  84.548  9.7880  0.002410 ** 
# scale(log(Rotation_cycles)):scale(log(Duration))                                7.246   7.246     1  88.536  1.1043  0.296194    
# scale(log(Species_Richness)):scale(log(Duration))                             100.626 100.626     1 156.634 15.3358  0.000134 ***
# scale(log(Rotation_cycles)):scale(log(Species_Richness)):scale(log(Duration))  31.436  31.436     1 171.249  4.7910  0.029962 *  



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## Species_Richness
sum(!is.na(BacterialBeta$Species_Richness)) ## n = 343
p1 <- ggplot(BacterialBeta, aes(y=RR, x=Species_Richness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="Species_Richness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Species_Richness" , y="lnBacterialBeta343")
p1
pdf("Species_Richness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(BacterialBeta$Duration)) ## n = 343
p2 <- ggplot(BacterialBeta, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnBacterialBeta343")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(BacterialBeta$Rotation_cycles)) ## n = 343
p3 <- ggplot(BacterialBeta, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnBacterialBeta343")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(BacterialBeta$Latitude)) ## n = 343
p5 <- ggplot(BacterialBeta, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnBacterialBeta343")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(BacterialBeta$Longitude)) ## n = 343
p6 <- ggplot(BacterialBeta, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnBacterialBeta343")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(BacterialBeta$MAPmean)) ## n = 343
p7 <- ggplot(BacterialBeta, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnBacterialBeta343")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(BacterialBeta$MATmean)) ## n = 343
p8 <- ggplot(BacterialBeta, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnBacterialBeta343")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(BacterialBeta$pHCK)) ## n = 142
p9 <- ggplot(BacterialBeta, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnBacterialBeta142")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(BacterialBeta$SOCCK)) ## n = 148
p10 <- ggplot(BacterialBeta, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnBacterialBeta148")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(BacterialBeta$TNCK)) ## n = 97
p11 <- ggplot(BacterialBeta, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnBacterialBeta97")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(BacterialBeta$NO3CK)) ## n = 50
p12 <- ggplot(BacterialBeta, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnBacterialBeta50")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(BacterialBeta$NH4CK)) ## n = 46
p13<- ggplot(BacterialBeta, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnBacterialBeta46")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(BacterialBeta$APCK)) ## n = 120
p14 <- ggplot(BacterialBeta, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnBacterialBeta120")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(BacterialBeta$AKCK)) ## n = 99
p15 <- ggplot(BacterialBeta, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnBacterialBeta99")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(BacterialBeta$ANCK)) ## n = 80
p16 <- ggplot(BacterialBeta, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnBacterialBeta80")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(BacterialBeta$RRpH)) ## n = 142
p17 <- ggplot(BacterialBeta, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="RRpH")+
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
sum(!is.na(BacterialBeta$RRSOC)) ## n = 148
p18 <- ggplot(BacterialBeta, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="RRSOC")+
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
sum(!is.na(BacterialBeta$RRTN)) ## n = 97
p19 <- ggplot(BacterialBeta, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="RRTN")+
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
sum(!is.na(BacterialBeta$RRNO3)) ## n = 50
p20 <- ggplot(BacterialBeta, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="RRNO3")+
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
sum(!is.na(BacterialBeta$RRNH4)) ## n = 46
p21 <- ggplot(BacterialBeta, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="RRNH4")+
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
sum(!is.na(BacterialBeta$RRAP)) ## n = 120
p22 <- ggplot(BacterialBeta, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="RRAP")+
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
sum(!is.na(BacterialBeta$RRAK)) ## n = 99
p23 <- ggplot(BacterialBeta, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="RRAK")+
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
sum(!is.na(BacterialBeta$RRAN)) ## n = 80
p24 <- ggplot(BacterialBeta, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacterialBeta", x="RRAN")+
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
sum(!is.na(BacterialBeta$RRYield)) ## n = 35
p25 <- ggplot(BacterialBeta, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-3.5, 1.2), expand=c(0, 0)) + 
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
                            data = BacterialBeta, 
                            na.action = na.roughfix, 
                            importance = TRUE, 
                            ntree = 500)
# Check the importance of variables and p-values
importance(rf_model_perm3)
#                   %IncMSE %IncMSE.pval IncNodePurity IncNodePurity.pval
# Latitude         21.39039   0.00990099     28.905800         0.00990099
# MAPmean          20.24181   0.00990099     18.223556         0.00990099
# Longitude        18.93459   0.00990099     18.733585         0.00990099
# SOCCK            17.59169   0.00990099     18.100948         0.35643564
# MATmean          16.58132   0.00990099     17.431452         0.00990099
# Species_Richness 14.19178   0.00990099      7.163687         0.02970297
# Rotation_cycles  13.40450   0.06930693     13.901720         0.65346535
# Duration         12.95025   0.07920792     11.149034         0.45544554
# pHCK             12.82863   0.00990099     15.776471         0.37623762


###################################### Trials sorted by effect size
library(ggplot2)
library(dplyr)
BacterialBeta <- read.csv("BacterialBeta.csv", fileEncoding = "latin1")
# 计算95% CI + 显著性分类
df_plot <- BacterialBeta %>%
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
#       84      185       74 
## 8*8


################################################ piecewiseSEM
sem_data <- read.csv("BacterialBeta.csv", fileEncoding = "latin1")
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
# -148.454
# Fisher's C = 1.421 with P-value = 0.492

m1_1 <- lme(RR ~ RRpH + RRSOC + Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_2 <- lme(RRpH ~  Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_3 <- lme(RRSOC ~ Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model1 <- psem(m1_1, m1_2, m1_3)
summary(sem_model1)  ## 
# AIC
# -174.407
# Fisher's C = 2.106 with P-value = 0.349

m2_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_2 <- lme(RRpH ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_3 <- lme(RRSOC ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model2 <- psem(m2_1, m2_2, m2_3)
summary(sem_model2)  ## 
# AIC
# -179.604
# Fisher's C = 1.176 with P-value = 0.556 

m3_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_2 <- lme(RRpH ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_3 <- lme(RRSOC ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model3 <- psem(m3_1, m3_2, m3_3)
summary(sem_model3)  ## 
# AIC
# -158.618
# Fisher's C = 1.044 with P-value = 0.593

m4_1 <- lme(RR ~ RRpH + RRSOC + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_2 <- lme(RRpH ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_3 <- lme(RRSOC ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model4 <- psem(m4_1, m4_2, m4_3)
summary(sem_model4) ## 
# AIC
# -205.484
# Fisher's C = 1.825 with P-value = 0.401

# Structural Equation Model of sem_model4 
# 
# Call:
#   RR ~ RRpH + RRSOC + Species_Richness
# RRpH ~ Species_Richness
# RRSOC ~ Species_Richness
# 
# AIC
# -205.484
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
#         RR             RRpH  -1.2445    1.1958 60    -1.0407  0.3022      -0.1132  
#         RR            RRSOC   0.2376    0.4182 60     0.5682  0.5720       0.0638  
#         RR Species_Richness   0.2497    0.2397 60     1.0418  0.3017       0.1533  
#       RRpH Species_Richness  -0.0346    0.0134 74    -2.5768  0.0120      -0.2334 *
#      RRSOC Species_Richness  -0.0638    0.0502 74    -1.2706  0.2079      -0.1458  
# 
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05
# 
# ---
# Individual R-squared:
# 
#   Response method Marginal Conditional
#         RR   none     0.04        0.66
#       RRpH   none     0.04        0.86
#      RRSOC   none     0.02        0.65

m5_1 <- lme(RR ~ RRpH + RRSOC + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_2 <- lme(RRpH ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_3 <- lme(RRSOC ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model5 <- psem(m5_1, m5_2, m5_3)
summary(sem_model5) ## 
# AIC
# -183.284
# Fisher's C = 1.727 with P-value = 0.422

m6_1 <- lme(RR ~ RRpH + RRSOC +MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_2 <- lme(RRpH ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_3 <- lme(RRSOC ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model6 <- psem(m6_1, m6_2, m6_3)
summary(sem_model6)  ##
# AIC
# -189.956
# Fisher's C = 0.876 with P-value = 0.645
