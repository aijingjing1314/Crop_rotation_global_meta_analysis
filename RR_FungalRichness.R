
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)

FungalRichness <- read.csv("FungalRichness.csv", fileEncoding = "latin1")
# Check data
head(FungalRichness)

# 1. The number of Obversation
total_number <- nrow(FungalRichness)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 267

# 2. The number of Study
unique_studyid_number <- length(unique(FungalRichness$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 80


#### 3. Overall effect size
total_effect_model <- rma.mv(yi = RR, 
                             V = Vi, 
                             random = ~ 1 | StudyID,  # StudyID is radom factor
                             data = FungalRichness, 
                             method = "REML")
# The results of Overall effect size
summary(total_effect_model)
# Multivariate Meta-Analysis Model (k = 267; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -819.2508  1638.5017  1642.5017  1649.6687  1642.5473   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0255  0.1595     80     no  StudyID 
# Test for Heterogeneity:
# Q(df = 266) = 3912.1529, p-val < .0001
# Model Results:
# estimate      se    zval    pval    ci.lb   ci.ub    
#   0.0113  0.0188  0.6023  0.5470  -0.0256  0.0483     

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
boot_results1 <- boot(data = FungalRichness, statistic = boot_fun, R = 1000, parallel = "snow", ncpus = numCores, cl = cl)
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
# Estimate for coefficient 1 : 0.01134564 
# 95% BCa CI for coefficient 1 : -0.01033163 0.03708599 


#### 5. Funnel Plot
simple_model <- rma(yi = RR, 
                    vi = Vi, 
                    data = FungalRichness, 
                    method = "REML")
#### 
funnel(simple_model)
# Output  6 * 6
#### Egger's test
regtest(simple_model)
# Regression Test for Funnel Plot Asymmetry
# Model:     mixed-effects meta-regression model
# Predictor: standard error
# Test for Funnel Plot Asymmetry: z = -1.1974, p = 0.2312
# Limit Estimate (as sei -> 0):   b =  0.0266 (CI: -0.0239, 0.0772)

#  Rosenthal’s Fail-Safe N
# This method estimates how many missing studies with null effect 
# would be needed to make the overall effect non-significant
fsn_rosenthal <- fsn(x = simple_model, type = "Rosenthal")
# Print the FSN result
print(fsn_rosenthal)
# Fail-safe N Calculation Using the General Approach
# Average Effect Size:         0.0002
# Amount of Heterogeneity:     0.0411
# Observed Significance Level: 0.9905
# Target Significance Level:   0.05
# Fail-safe N: 0


#### 8. Subgroup analysis
### 8.1 LegumeNonlegume
FungalRichness_filteredLegumeNonlegume <- subset(FungalRichness, LegumeNonlegume %in% c("Legume to Non-legume", "Non-legume to Legume", "Non-legume to Non-legume"))
#
FungalRichness_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(FungalRichness_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredLegumeNonlegume %>%
  group_by(LegumeNonlegume) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   LegumeNonlegume          Observations Unique_StudyID
#   <fct>                           <int>          <int>
# 1 Legume to Non-legume               45             23
# 2 Non-legume to Legume               93             30
# 3 Non-legume to Non-legume          128             44

overall_model_FungalRichness_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = FungalRichness_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -787.7606  1575.5212  1583.5212  1597.8099  1583.6763   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0253  0.1592     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 263) = 3729.9814, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 43.5865, p-val < .0001
# Model Results:
#                                          estimate      se     zval    pval    ci.lb    ci.ub    
# LegumeNonlegumeLegume to Non-legume       -0.0404  0.0204  -1.9767  0.0481  -0.0804  -0.0003  * 
# LegumeNonlegumeNon-legume to Legume        0.0246  0.0196   1.2533  0.2101  -0.0139   0.0631    
# LegumeNonlegumeNon-legume to Non-legume    0.0258  0.0193   1.3340  0.1822  -0.0121   0.0637        

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredLegumeNonlegume)
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
FungalRichness_filteredAMnonAM <- subset(FungalRichness, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
FungalRichness_filteredAMnonAM$AMnonAM <- droplevels(factor(FungalRichness_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredAMnonAM %>%
  group_by(AMnonAM) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   AMnonAM     Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 AM to AM             212             73
# 2 AM to nonAM           39             13
# 3 nonAM to AM           15              4
overall_model_FungalRichness_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = FungalRichness_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -755.8143  1511.6286  1519.6286  1533.9172  1519.7836   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0264  0.1626     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 263) = 3854.2241, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 137.1558, p-val < .0001
# Model Results:
#                     estimate      se     zval    pval    ci.lb   ci.ub      
# AMnonAMAM to AM      -0.0024  0.0193  -0.1251  0.9004  -0.0402  0.0353      
# AMnonAMAM to nonAM    0.0366  0.0225   1.6310  0.1029  -0.0074  0.0806      
# AMnonAMnonAM to AM    0.1827  0.0241   7.5673  <.0001   0.1354  0.2300  *** 
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredAMnonAM)
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
   #             "a"                "b"                "c" 


### 8.3 C3C4
FungalRichness_filteredC3C4 <- subset(FungalRichness, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
FungalRichness_filteredC3C4$C3C4 <- droplevels(factor(FungalRichness_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredC3C4 %>%
  group_by(C3C4) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#  C3C4     Observations Unique_StudyID
#   <fct>           <int>          <int>
# 1 C3 to C3          141             47
# 2 C3 to C4           64             34
# 3 C4 to C3           58             20

overall_model_FungalRichness_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = FungalRichness_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 263; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -767.0391  1534.0782  1542.0782  1556.3210  1542.2351   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0257  0.1602     79     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 260) = 3799.5458, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 49.0147, p-val < .0001
# Model Results:
#               estimate      se     zval    pval    ci.lb   ci.ub    
# C3C4C3 to C3    0.0421  0.0196   2.1467  0.0318   0.0037  0.0806  * 
# C3C4C3 to C4   -0.0155  0.0203  -0.7659  0.4437  -0.0553  0.0242    
# C3C4C4 to C3   -0.0330  0.0240  -1.3794  0.1678  -0.0800  0.0139     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredC3C4)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredC3C4)
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
#        "a"          "b"          "b"  


### 8.4 Annual_Pere
FungalRichness_filteredAnnual_Pere <- subset(FungalRichness, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
FungalRichness_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(FungalRichness_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredAnnual_Pere %>%
  group_by(Annual_Pere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Annual_Pere            Observations Unique_StudyID
#   <fct>                         <int>          <int>
# 1 Annual to Annual                184             60
# 2 Annual to Perennial              26             13
# 3 Perennial to Annual              46             12
# 4 Perennial to Perennial           11              7

overall_model_FungalRichness_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = FungalRichness_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 267; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -804.7896  1609.5792  1619.5792  1637.4399  1619.8126   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0261  0.1616     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 263) = 3821.1262, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 36.6051, p-val < .0001
# Model Results:
#                                    estimate      se     zval    pval    ci.lb    ci.ub     
# Annual_PereAnnual to Annual          0.0094  0.0202   0.4647  0.6421  -0.0302   0.0489     
# Annual_PereAnnual to Perennial      -0.0634  0.0244  -2.5992  0.0093  -0.1112  -0.0156  ** 
# Annual_PerePerennial to Annual       0.0937  0.0333   2.8141  0.0049   0.0284   0.1590  ** 
# Annual_PerePerennial to Perennial   -0.0216  0.0366  -0.5907  0.5547  -0.0933   0.0501     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredAnnual_Pere)
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
      #                  "a"                               "b"                               "c"                              "ab" 



### 8.6 PlantStageSubgroup
FungalRichness_filteredPlantStageSubgroup <- subset(FungalRichness, PlantStageSubgroup %in% c("Vegetative stage","Reproductive stage", "Maturity stage","Harvest"))
#
FungalRichness_filteredPlantStageSubgroup$PlantStageSubgroup <- droplevels(factor(FungalRichness_filteredPlantStageSubgroup$PlantStageSubgroup))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredPlantStageSubgroup %>%
  group_by(PlantStageSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
  # PlantStageSubgroup Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Harvest                     139             40
# 2 Maturity stage               29             12
# 3 Reproductive stage           49             16
# 4 Vegetative stage             16              6

overall_model_FungalRichness_filteredPlantStageSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + PlantStageSubgroup, random = ~ 1 | StudyID, data = FungalRichness_filteredPlantStageSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredPlantStageSubgroup)
# Multivariate Meta-Analysis Model (k = 233; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -208.8890   417.7779   427.7779   444.9465   428.0470   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0244  0.1562     67     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 229) = 2272.2491, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 32.7605, p-val < .0001
# Model Results:
#                                       estimate      se     zval    pval    ci.lb   ci.ub     
# PlantStageSubgroupHarvest              -0.0037  0.0209  -0.1745  0.8615  -0.0447  0.0374     
# PlantStageSubgroupMaturity stage       -0.0318  0.0264  -1.2065  0.2276  -0.0836  0.0199     
# PlantStageSubgroupReproductive stage    0.0290  0.0228   1.2720  0.2034  -0.0157  0.0736     
# PlantStageSubgroupVegetative stage      0.0794  0.0244   3.2571  0.0011   0.0316  0.1272  ** 
#    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredPlantStageSubgroup)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredPlantStageSubgroup)
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
           #                      "a"                                  "a"                                  "b"                                  "c" 


### 8.7 Bulk_Rhizosphere
FungalRichness_filteredBulk_Rhizosphere <- subset(FungalRichness, Bulk_Rhizosphere %in% c("Non-Rhizosphere", "Rhizosphere"))
#
FungalRichness_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(FungalRichness_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredBulk_Rhizosphere %>%
  group_by(Bulk_Rhizosphere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Bulk_Rhizosphere Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 Non-Rhizosphere           149             49
# 2 Rhizosphere               118             41

overall_model_FungalRichness_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = FungalRichness_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 267; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -803.2450  1606.4899  1612.4899  1623.2291  1612.5819   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0267  0.1635     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 265) = 3910.6529, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 35.4822, p-val < .0001
# Model Results:
#                                  estimate      se     zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon-Rhizosphere   -0.0175  0.0199  -0.8804  0.3786  -0.0564  0.0214    
# Bulk_RhizosphereRhizosphere        0.0454  0.0201   2.2589  0.0239   0.0060  0.0849  * 
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredBulk_Rhizosphere)
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
FungalRichness_filteredSoil_texture <- subset(FungalRichness, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
FungalRichness_filteredSoil_texture$Soil_texture <- droplevels(factor(FungalRichness_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredSoil_texture %>%
  group_by(Soil_texture) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Soil_texture Observations Unique_StudyID
#   <fct>               <int>          <int>
# 1 Coarse                 66             12
# 2 Fine                   41             12
# 3 Medium                 95             31

overall_model_FungalRichness_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = FungalRichness_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 202; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -593.0581  1186.1161  1194.1161  1207.2893  1194.3223  
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0256  0.1600     55     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 199) = 2684.1621, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 5.9430, p-val = 0.1144
# Model Results:
#                     estimate      se     zval    pval    ci.lb    ci.ub    
# Soil_textureCoarse    0.0581  0.0478   1.2159  0.2240  -0.0355   0.1517    
# Soil_textureFine     -0.1067  0.0505  -2.1128  0.0346  -0.2056  -0.0077  * 
# Soil_textureMedium   -0.0009  0.0300  -0.0291  0.9768  -0.0597   0.0579     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredSoil_texture)
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
#               "a"                "b"               "ab" 


### 8.9 Tillage
FungalRichness_filteredTillage <- subset(FungalRichness, Tillage %in% c("Tillage", "No_tillage"))
#
FungalRichness_filteredTillage$Tillage <- droplevels(factor(FungalRichness_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredTillage %>%
  group_by(Tillage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Tillage    Observations Unique_StudyID
#   <fct>             <int>          <int>
# 1 No_tillage           43              7
# 2 Tillage              34             10

overall_model_FungalRichness_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = FungalRichness_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredTillage)
# Multivariate Meta-Analysis Model (k = 77; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -502.3265  1004.6531  1010.6531  1017.6056  1010.9911   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0199  0.1410     14     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 75) = 1369.1851, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 15.8305, p-val = 0.0004
# Model Results:
#                    estimate      se     zval    pval    ci.lb    ci.ub     
# TillageNo_tillage   -0.1332  0.0416  -3.1988  0.0014  -0.2148  -0.0516  ** 
# TillageTillage      -0.0660  0.0406  -1.6252  0.1041  -0.1456   0.0136     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredTillage)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredTillage)
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
FungalRichness_filteredStraw_retention <- subset(FungalRichness, Straw_retention %in% c("Retention", "No_retention"))
#
FungalRichness_filteredStraw_retention$Straw_retention <- droplevels(factor(FungalRichness_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredStraw_retention %>%
  group_by(Straw_retention) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Straw_retention Observations Unique_StudyID
#   <fct>                  <int>          <int>
# 1 No_retention              19              9
# 2 Retention                 35             10

overall_model_FungalRichness_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = FungalRichness_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 54; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -36.0327   72.0655   78.0655   83.9192   78.5655  
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0104  0.1021     17     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 52) = 369.4170, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 6.0291, p-val = 0.0491
# Model Results:
#                              estimate      se     zval    pval    ci.lb    ci.ub    
# Straw_retentionNo_retention    0.0220  0.0347   0.6337  0.5263  -0.0460   0.0900    
# Straw_retentionRetention      -0.0707  0.0332  -2.1267  0.0334  -0.1359  -0.0055  *   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredStraw_retention)
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
FungalRichness_filteredRotationcyclesSubgroup <- subset(FungalRichness, RotationcyclesSubgroup %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
FungalRichness_filteredRotationcyclesSubgroup$RotationcyclesSubgroup <- droplevels(factor(FungalRichness_filteredRotationcyclesSubgroup$RotationcyclesSubgroup))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredRotationcyclesSubgroup %>%
  group_by(RotationcyclesSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   RotationcyclesSubgroup Observations Unique_StudyID
#   <fct>                         <int>          <int>
# 1 D1                              118             39
# 2 D1-3                             65             20
# 3 D10                              20              7
# 4 D3-5                             17              9
# 5 D5-10                            47             15

overall_model_FungalRichness_filteredRotationcyclesSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + RotationcyclesSubgroup, random = ~ 1 | StudyID, data = FungalRichness_filteredRotationcyclesSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredRotationcyclesSubgroup)
# Multivariate Meta-Analysis Model (k = 267; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -807.7600  1615.5201  1627.5201  1648.9302  1627.8495   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0318  0.1784     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 262) = 3822.3825, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 32.3922, p-val < .0001
# Model Results:
#                              estimate      se     zval    pval    ci.lb   ci.ub      
# RotationcyclesSubgroupD1      -0.0035  0.0240  -0.1464  0.8836  -0.0505  0.0434      
# RotationcyclesSubgroupD1-3    -0.0399  0.0246  -1.6221  0.1048  -0.0881  0.0083      
# RotationcyclesSubgroupD10      0.0018  0.0647   0.0279  0.9777  -0.1250  0.1286      
# RotationcyclesSubgroupD3-5     0.0075  0.0379   0.1987  0.8425  -0.0667  0.0818      
# RotationcyclesSubgroupD5-10    0.1142  0.0337   3.3900  0.0007   0.0482  0.1803  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredRotationcyclesSubgroup)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredRotationcyclesSubgroup)
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
   #           "a"                         "b"                       "abc"                        "ab"                         "c"  


### 8.13 DurationSubgroup
FungalRichness_filteredDurationSubgroup <- subset(FungalRichness, DurationSubgroup %in% c("D1", "D2", "D3", "D4", "D5", "D6-10", "D11-20", "D20-30", "D30"))
#
FungalRichness_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(FungalRichness_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredDurationSubgroup %>%
  group_by(DurationSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   DurationSubgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D1                         18              7
# 2 D11-20                     37             13
# 3 D2                         45             20
# 4 D20-30                     26              7
# 5 D3                         53             15
# 6 D30                        10              3
# 7 D4                         21              7
# 8 D5                         10              2
# 9 D6-10                      47             15

overall_model_FungalRichness_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = FungalRichness_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 267; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -800.3628  1600.7256  1620.7256  1656.2552  1621.6162   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0368  0.1919     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 258) = 3697.3620, p-val < .0001
# Test of Moderators (coefficients 1:9):
# QM(df = 9) = 48.6065, p-val < .0001
# Model Results:
#                         estimate      se     zval    pval    ci.lb    ci.ub      
# DurationSubgroupD1       -0.0402  0.0462  -0.8704  0.3841  -0.1307   0.0503      
# DurationSubgroupD11-20    0.1278  0.0392   3.2617  0.0011   0.0510   0.2046   ** 
# DurationSubgroupD2       -0.0032  0.0354  -0.0915  0.9271  -0.0726   0.0661      
# DurationSubgroupD20-30    0.1266  0.0467   2.7110  0.0067   0.0351   0.2181   ** 
# DurationSubgroupD3       -0.1387  0.0353  -3.9260  <.0001  -0.2080  -0.0695  *** 
# DurationSubgroupD30      -0.0104  0.1132  -0.0921  0.9266  -0.2324   0.2115      
# DurationSubgroupD4       -0.0430  0.0756  -0.5687  0.5695  -0.1911   0.1052      
# DurationSubgroupD5        0.0277  0.1472   0.1881  0.8508  -0.2609   0.3163      
# DurationSubgroupD6-10     0.0492  0.0500   0.9841  0.3251  -0.0488   0.1471      


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredDurationSubgroup)
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
    #                "a"                    "b"                    "a"                   "bc"                    "d"                 "abcd" 
    # DurationSubgroupD4     DurationSubgroupD5  DurationSubgroupD6-10 
    #              "acd"                 "abcd"                  "abc" 



### 8.14 SpeciesRichnessSubgroup
FungalRichness_filteredSpeciesRichnessSubgroup <- subset(FungalRichness, SpeciesRichnessSubgroup %in% c("R2", "R3"))
#
FungalRichness_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup <- droplevels(factor(FungalRichness_filteredSpeciesRichnessSubgroup$SpeciesRichnessSubgroup))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredSpeciesRichnessSubgroup %>%
  group_by(SpeciesRichnessSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   SpeciesRichnessSubgroup Observations Unique_StudyID
#   <fct>                          <int>          <int>
# 1 R2                               193             68
# 2 R3                                67             16

overall_model_FungalRichness_filteredSpeciesRichnessSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichnessSubgroup, random = ~ 1 | StudyID, data = FungalRichness_filteredSpeciesRichnessSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredSpeciesRichnessSubgroup)
# Multivariate Meta-Analysis Model (k = 260; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -817.6416  1635.2832  1641.2832  1651.9420  1641.3777   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0215  0.1465     77     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 258) = 3857.7137, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.3952, p-val = 0.4978
# Model Results:
#                            estimate      se    zval    pval    ci.lb   ci.ub    
# SpeciesRichnessSubgroupR2    0.0209  0.0178  1.1761  0.2396  -0.0139  0.0557    
# SpeciesRichnessSubgroupR3    0.0209  0.0194  1.0750  0.2824  -0.0172  0.0590       

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredSpeciesRichnessSubgroup)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredSpeciesRichnessSubgroup)
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
FungalRichness_filteredPrimer <- subset(FungalRichness, Primer %in% c("ITS1", "ITS1 + 5.8S + ITS2", "ITS2", "Full length", "18S"))
#
FungalRichness_filteredPrimer$Primer <- droplevels(factor(FungalRichness_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredPrimer %>%
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
# 3 ITS1                        224             61
# 4 ITS1 + 5.8S + ITS2            5              2
# 5 ITS2                         29             12

overall_model_FungalRichness_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = FungalRichness_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 267; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -816.2551  1632.5102  1644.5102  1665.9203  1644.8397   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0236  0.1536     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 262) = 3800.7750, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 10.3499, p-val = 0.0659
# Model Results:
#                           estimate      se     zval    pval    ci.lb   ci.ub     
# Primer18S                   0.2361  0.0863   2.7344  0.0063   0.0669  0.4053  ** 
# PrimerFull length           0.2447  0.1563   1.5657  0.1174  -0.0616  0.5510     
# PrimerITS1                  0.0019  0.0208   0.0933  0.9257  -0.0387  0.0426     
# PrimerITS1 + 5.8S + ITS2    0.0208  0.1111   0.1874  0.8514  -0.1970  0.2386     
# PrimerITS2                 -0.0292  0.0475  -0.6149  0.5386  -0.1223  0.0639    
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredPrimer)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredPrimer)
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
               #       "a"                     "ab"                      "b"                     "ab"                      "b" 



### 8.16 Latitude_Subgroup
FungalRichness_filteredLatitude_Subgroup <- subset(FungalRichness, Latitude_Subgroup %in% c("La20", "La20-40", "La40"))
#
FungalRichness_filteredLatitude_Subgroup$Latitude_Subgroup <- droplevels(factor(FungalRichness_filteredLatitude_Subgroup$Latitude_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredLatitude_Subgroup %>%
  group_by(Latitude_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Latitude_Subgroup Observations Unique_StudyID
#   <fct>                    <int>          <int>
# 1 La20                        25              5
# 2 La20-40                     97             37
# 3 La40                       145             38

overall_model_FungalRichness_filteredLatitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Latitude_Subgroup, random = ~ 1 | StudyID, data = FungalRichness_filteredLatitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredLatitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 267; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -818.4417  1636.8834  1644.8834  1659.1872  1645.0379   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0255  0.1597     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 264) = 3800.2414, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.6791, p-val = 0.4438
# Model Results:
#                           estimate      se     zval    pval    ci.lb   ci.ub    
# Latitude_SubgroupLa20       0.0739  0.0770   0.9591  0.3375  -0.0771  0.2248    
# Latitude_SubgroupLa20-40    0.0329  0.0280   1.1766  0.2394  -0.0219  0.0877    
# Latitude_SubgroupLa40      -0.0166  0.0271  -0.6124  0.5403  -0.0696  0.0365         

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredLatitude_Subgroup)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredLatitude_Subgroup)
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
FungalRichness_filteredLongitude_Subgroup <- subset(FungalRichness, Longitude_Subgroup %in% c("Lo-180-0", "Lo-0-180"))
#
FungalRichness_filteredLongitude_Subgroup$Longitude_Subgroup <- droplevels(factor(FungalRichness_filteredLongitude_Subgroup$Longitude_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredLongitude_Subgroup %>%
  group_by(Longitude_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Longitude_Subgroup Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Lo-0-180                    231             72
# 2 Lo-180-0                     36              8

overall_model_FungalRichness_filteredLongitude_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Subgroup, random = ~ 1 | StudyID, data = FungalRichness_filteredLongitude_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredLongitude_Subgroup)
# Multivariate Meta-Analysis Model (k = 267; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -817.8901  1635.7802  1641.7802  1652.5194  1641.8722   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0251  0.1584     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 265) = 3866.5051, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 3.3112, p-val = 0.1910
# Model Results:
#                             estimate      se     zval    pval    ci.lb   ci.ub    
# Longitude_SubgroupLo-0-180    0.0217  0.0197   1.1041  0.2696  -0.0168  0.0602    
# Longitude_SubgroupLo-180-0   -0.0885  0.0612  -1.4465  0.1480  -0.2083  0.0314           

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredLongitude_Subgroup)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredLongitude_Subgroup)
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
FungalRichness_filteredMAPmean_Subgroup <- subset(FungalRichness, MAPmean_Subgroup %in% c("MAP600", "MAP600-1200", "MAP1200"))
#
FungalRichness_filteredMAPmean_Subgroup$MAPmean_Subgroup <- droplevels(factor(FungalRichness_filteredMAPmean_Subgroup$MAPmean_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredMAPmean_Subgroup %>%
  group_by(MAPmean_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MAPmean_Subgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 MAP1200                    57             18
# 2 MAP600                    164             45
# 3 MAP600-1200                46             19

overall_model_FungalRichness_filteredMAPmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean_Subgroup, random = ~ 1 | StudyID, data = FungalRichness_filteredMAPmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredMAPmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 267; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -817.4869  1634.9738  1642.9738  1657.2776  1643.1283   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0264  0.1624     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 264) = 3818.5835, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 5.9551, p-val = 0.1138
# Model Results:
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean_SubgroupMAP1200        0.0909  0.0413   2.2032  0.0276   0.0100  0.1718  * 
# MAPmean_SubgroupMAP600        -0.0019  0.0231  -0.0825  0.9343  -0.0472  0.0433    
# MAPmean_SubgroupMAP600-1200   -0.0285  0.0308  -0.9260  0.3544  -0.0888  0.0318    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredMAPmean_Subgroup)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredMAPmean_Subgroup)
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
#                     "a"                         "b"                         "b" 



### 8.19 MATmean_Subgroup
FungalRichness_filteredMATmean_Subgroup <- subset(FungalRichness, MATmean_Subgroup %in% c("MAT8", "MAT8-15", "MAT15"))
#
FungalRichness_filteredMATmean_Subgroup$MATmean_Subgroup <- droplevels(factor(FungalRichness_filteredMATmean_Subgroup$MATmean_Subgroup))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredMATmean_Subgroup %>%
  group_by(MATmean_Subgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MATmean_Subgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 MAT15                      76             26
# 2 MAT8                      139             40
# 3 MAT8-15                    52             15

overall_model_FungalRichness_filteredMATmean_Subgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Subgroup, random = ~ 1 | StudyID, data = FungalRichness_filteredMATmean_Subgroup, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredMATmean_Subgroup)
# Multivariate Meta-Analysis Model (k = 267; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -817.9769  1635.9538  1643.9538  1658.2576  1644.1083   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0251  0.1585     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 264) = 3817.5881, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 4.9134, p-val = 0.1782
# Model Results:
#                          estimate      se     zval    pval    ci.lb   ci.ub    
# MATmean_SubgroupMAT15      0.0610  0.0334   1.8268  0.0677  -0.0044  0.1263  . 
# MATmean_SubgroupMAT8      -0.0018  0.0241  -0.0730  0.9418  -0.0491  0.0456    
# MATmean_SubgroupMAT8-15   -0.0387  0.0328  -1.1804  0.2378  -0.1030  0.0256    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredMATmean_Subgroup)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredMATmean_Subgroup)
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
#              "a"                    "ab"                     "b" 



### 8.20 DurationSubgroup
FungalRichness_filteredDurationSubgroup <- subset(FungalRichness, DurationSubgroup %in% c("D5", "D5-10", "D10-20", "D20-30", "D30"))
#
FungalRichness_filteredDurationSubgroup$DurationSubgroup <- droplevels(factor(FungalRichness_filteredDurationSubgroup$DurationSubgroup))
# The number of Observations and StudyID
group_summary <- FungalRichness_filteredDurationSubgroup %>%
  group_by(DurationSubgroup) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   DurationSubgroup Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D10-20                     37             13
# 2 D20-30                     26              7
# 3 D30                        10              3
# 4 D5                        147             46
# 5 D5-10                      47             15

overall_model_FungalRichness_filteredDurationSubgroup <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + DurationSubgroup, random = ~ 1 | StudyID, data = FungalRichness_filteredDurationSubgroup, method = "REML")
# QM and p value
summary(overall_model_FungalRichness_filteredDurationSubgroup)
# Multivariate Meta-Analysis Model (k = 267; method: REML)
#    logLik   Deviance        AIC        BIC       AICc   
# -807.0956  1614.1911  1626.1911  1647.6012  1626.5205   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0368  0.1917     80     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 262) = 3797.1359, p-val < .0001
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 31.5122, p-val < .0001
# Model Results:
#                         estimate      se     zval    pval    ci.lb    ci.ub      
# DurationSubgroupD10-20    0.1598  0.0381   4.1914  <.0001   0.0851   0.2345  *** 
# DurationSubgroupD20-30    0.1516  0.0462   3.2842  0.0010   0.0611   0.2420   ** 
# DurationSubgroupD30      -0.0104  0.1132  -0.0922  0.9265  -0.2322   0.2114      
# DurationSubgroupD5       -0.0603  0.0277  -2.1743  0.0297  -0.1146  -0.0059    * 
# DurationSubgroupD5-10     0.0522  0.0499   1.0456  0.2957  -0.0456   0.1500     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_FungalRichness_filteredDurationSubgroup)
vcov_rotation <- vcov(overall_model_FungalRichness_filteredDurationSubgroup)
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
#              "a"                    "a"                   "ab"                    "b"                    "a" 







#### 9. Linear Mixed Effect Model
# 
FungalRichness$Wr <- 1 / FungalRichness$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalRichness)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalRichness)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalRichness)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalRichness)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalRichness)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(Species_Richness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalRichness)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(Species_Richness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = FungalRichness)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = FungalRichness)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(Species_Richness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = FungalRichness)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(Species_Richness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(Species_Richness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = FungalRichness)

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
#           Model      AIC      BIC       logLik  Marginal_R2 Conditional_R2
# Model1   Model1 11.15768 32.68118   0.42115806 1.815917e-05   0.0003443979
# Model2   Model2 12.04248 33.56597  -0.02124025 6.023736e-06   0.0003508486
# Model3   Model3 12.16853 33.69202  -0.08426635 7.895655e-06   0.0003425366
# Model4   Model4 11.40502 32.92851   0.29749209 1.725918e-05   0.0003464546
# Model5   Model5 12.11724 33.64073  -0.05861846 1.068920e-05   0.0003478322
# Model6   Model6 12.22601 33.74950  -0.11300477 1.069513e-05   0.0003484742
# Model7   Model7 11.95838 33.48187   0.02081040 5.983261e-06   0.0003496655
# Model8   Model8 12.32976 33.85325  -0.16488012 7.759457e-06   0.0003443413
# Model9   Model9 45.79067 81.66315 -12.89533265 4.990850e-05   0.0003779858
# Model10 Model10 43.50815 79.38063 -11.75407357 2.166519e-05   0.0003901719

##### Model 1 is the best model
summary(Model1)
# Number of obs: 267, groups:  StudyID, 80
anova(Model1) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                          Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)   8.5700  8.5700     1  49.596  0.6897 0.4102
# scale(Species_Richness)  2.1367  2.1367     1 194.598  0.1720 0.6788
# scale(Duration)         14.0645 14.0645     1  60.988  1.1319 0.2916

#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = FungalRichness)
summary(ModelpH)
# Number of obs: 150, groups:  StudyID, 52
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                          Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)  13.0172 13.0172     1 21.954  0.8912 0.3554
# scale(Species_Richness)  0.0309  0.0309     1 66.406  0.0021 0.9635
# scale(Duration)          5.5525  5.5525     1 25.470  0.3801 0.5430
# scale(pHCK)              5.2032  5.2032     1 28.688  0.3562 0.5553

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = FungalRichness)
summary(ModelSOC)
# Number of obs: 163, groups:  StudyID, 53
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                          Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)  16.3286 16.3286     1 26.696  1.2114 0.2809
# scale(Species_Richness)  1.3888  1.3888     1 81.739  0.1030 0.7490
# scale(Duration)          8.1405  8.1405     1 32.463  0.6039 0.4427
# scale(SOCCK)             1.2900  1.2900     1 47.867  0.0957 0.7584

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = FungalRichness)
summary(ModelTN)
# Number of obs: 108, groups:  StudyID, 40
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)  6.9623  6.9623     1  57.112  1.6389 0.2056
# scale(Species_Richness) 2.0525  2.0525     1 102.215  0.4832 0.4886
# scale(Duration)         4.4091  4.4091     1  54.889  1.0379 0.3128
# scale(TNCK)             4.8242  4.8242     1  41.352  1.1356 0.2928


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = FungalRichness)
summary(ModelNO3)
# Number of obs: 38, groups:  StudyID, 14
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                          Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)  0.82111 0.82111     1 15.801  0.1946 0.6651
# scale(Species_Richness) 0.03859 0.03859     1 11.259  0.0091 0.9255
# scale(Duration)         0.00981 0.00981     1 32.565  0.0023 0.9618
# scale(NO3CK)            0.81783 0.81783     1 12.264  0.1938 0.6674

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = FungalRichness)
summary(ModelNH4)
# Number of obs: 36, groups:  StudyID, 13
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)  5.0843  5.0843     1 12.3340  1.4461 0.2517
# scale(Species_Richness) 1.2707  1.2707     1  9.8306  0.3614 0.5613
# scale(Duration)         1.2034  1.2034     1 30.9607  0.3423 0.5628
# scale(NH4CK)            3.3969  3.3969     1  8.0247  0.9662 0.3543

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = FungalRichness)
summary(ModelAP)
# Number of obs: 113, groups:  StudyID, 43
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)  8.2334  8.2334     1  36.804  1.0957 0.3020
# scale(Species_Richness) 3.3810  3.3810     1 100.609  0.4499 0.5039
# scale(Duration)         0.8442  0.8442     1  64.913  0.1123 0.7386
# scale(APCK)             9.0747  9.0747     1  34.764  1.2076 0.2794

#### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = FungalRichness)
summary(ModelAK)
# Number of obs: 92, groups:  StudyID, 36
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                          Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)   0.2318  0.2318     1 33.924  0.0274 0.8696
# scale(Species_Richness)  1.3388  1.3388     1 85.607  0.1580 0.6920
# scale(Duration)          0.7102  0.7102     1 35.533  0.0838 0.7739
# scale(AKCK)             10.8710 10.8710     1 46.048  1.2829 0.2632

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = FungalRichness)
summary(ModelAN)
# Number of obs: 74, groups:  StudyID, 21
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)  17.490  17.490     1  8.244  1.6950 0.2281
# scale(Species_Richness)  0.013   0.013     1 68.999  0.0013 0.9718
# scale(Duration)         33.505  33.505     1 10.342  3.2470 0.1007
# scale(ANCK)             16.815  16.815     1  5.783  1.6296 0.2506

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = FungalRichness)
summary(ModelLatitude)
# Number of obs: 267, groups:  StudyID, 80
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)  11.541  11.541     1  49.529  0.9293 0.3397
# scale(Species_Richness)  4.667   4.667     1 207.588  0.3758 0.5405
# scale(Duration)         16.591  16.591     1  60.457  1.3360 0.2523
# scale(Latitude)         11.840  11.840     1  28.842  0.9534 0.3370

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = FungalRichness)
summary(ModelLongitude)
# Number of obs: 267, groups:  StudyID, 80
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                          Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)   4.2225  4.2225     1  52.620  0.3403 0.5621
# scale(Species_Richness)  1.6467  1.6467     1 193.653  0.1327 0.7160
# scale(Duration)          5.2737  5.2737     1  67.573  0.4250 0.5166
# scale(Longitude)        15.1440 15.1440     1  49.588  1.2205 0.2746

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = FungalRichness)
summary(ModelMAPmean)
# Number of obs: 267, groups:  StudyID, 80
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                          Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)   8.4613  8.4613     1  45.464  0.6805 0.4137
# scale(Species_Richness)  2.9396  2.9396     1 201.310  0.2364 0.6273
# scale(Duration)         13.0270 13.0270     1  54.511  1.0478 0.3105
# scale(MAPmean)           5.3178  5.3178     1  22.406  0.4277 0.5198

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(Rotation_cycles) + scale(Species_Richness) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = FungalRichness)
summary(ModelMATmean)
# Number of obs: 267, groups:  StudyID, 80
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                          Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)   7.4403  7.4403     1  47.335  0.6012 0.4420
# scale(Species_Richness)  3.1306  3.1306     1 200.711  0.2529 0.6156
# scale(Duration)         11.6249 11.6249     1  57.433  0.9393 0.3365
# scale(MATmean)          10.5108 10.5108     1  22.404  0.8493 0.3666 



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## Species_Richness
sum(!is.na(FungalRichness$Species_Richness)) ## n = 267
p1 <- ggplot(FungalRichness, aes(y=RR, x=Species_Richness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="Species_Richness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Species_Richness" , y="lnFungalRichness267")
p1
pdf("Species_Richness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(FungalRichness$Duration)) ## n = 267
p2 <- ggplot(FungalRichness, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnFungalRichness267")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(FungalRichness$Rotation_cycles)) ## n = 267
p3 <- ggplot(FungalRichness, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnFungalRichness267")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(FungalRichness$Latitude)) ## n = 267
p5 <- ggplot(FungalRichness, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnFungalRichness267")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(FungalRichness$Longitude)) ## n = 267
p6 <- ggplot(FungalRichness, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnFungalRichness267")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(FungalRichness$MAPmean)) ## n = 267
p7 <- ggplot(FungalRichness, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnFungalRichness267")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(FungalRichness$MATmean)) ## n = 267
p8 <- ggplot(FungalRichness, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnFungalRichness267")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(FungalRichness$pHCK)) ## n = 150
p9 <- ggplot(FungalRichness, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnFungalRichness150")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(FungalRichness$SOCCK)) ## n = 163
p10 <- ggplot(FungalRichness, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnFungalRichness163")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(FungalRichness$TNCK)) ## n = 108
p11 <- ggplot(FungalRichness, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnFungalRichness108")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(FungalRichness$NO3CK)) ## n = 38
p12 <- ggplot(FungalRichness, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnFungalRichness38")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(FungalRichness$NH4CK)) ## n = 36
p13<- ggplot(FungalRichness, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnFungalRichness36")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(FungalRichness$APCK)) ## n = 113
p14 <- ggplot(FungalRichness, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnFungalRichness113")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(FungalRichness$AKCK)) ## n = 92
p15 <- ggplot(FungalRichness, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnFungalRichness92")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(FungalRichness$ANCK)) ## n = 74
p16 <- ggplot(FungalRichness, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnFungalRichness74")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(FungalRichness$RRpH)) ## n = 150
p17 <- ggplot(FungalRichness, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="RRpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRpH" , y="RR150")
p17
pdf("RRpH.pdf",width=8,height=8)
p17
dev.off() 

## RRSOC
sum(!is.na(FungalRichness$RRSOC)) ## n = 163
p18 <- ggplot(FungalRichness, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="RRSOC")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRSOC" , y="RR163")
p18
pdf("RRSOC.pdf",width=8,height=8)
p18
dev.off() 

## RRTN
sum(!is.na(FungalRichness$RRTN)) ## n = 108
p19 <- ggplot(FungalRichness, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="RRTN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRTN" , y="RR108")
p19
pdf("RRTN.pdf",width=8,height=8)
p19
dev.off() 

## RRNO3
sum(!is.na(FungalRichness$RRNO3)) ## n = 38
p20 <- ggplot(FungalRichness, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="RRNO3")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNO3" , y="RR38")
p20
pdf("RRNO3.pdf",width=8,height=8)
p20
dev.off() 

## RRNH4
sum(!is.na(FungalRichness$RRNH4)) ## n = 36
p21 <- ggplot(FungalRichness, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="RRNH4")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNH4" , y="RR36")
p21
pdf("RRNH4.pdf",width=8,height=8)
p21
dev.off() 

## RRAP
sum(!is.na(FungalRichness$RRAP)) ## n = 113
p22 <- ggplot(FungalRichness, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="RRAP")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAP" , y="RR113")
p22
pdf("RRAP.pdf",width=8,height=8)
p22
dev.off() 

## RRAK
sum(!is.na(FungalRichness$RRAK)) ## n = 92
p23 <- ggplot(FungalRichness, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="RRAK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAK" , y="RR92")
p23
pdf("RRAK.pdf",width=8,height=8)
p23
dev.off() 

## RRAN
sum(!is.na(FungalRichness$RRAN)) ## n = 74
p24 <- ggplot(FungalRichness, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungalRichness", x="RRAN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAN" , y="RR74")
p24
pdf("RRAN.pdf",width=8,height=8)
p24
dev.off() 

## RRYield
sum(!is.na(FungalRichness$RRYield)) ## n = 34
p25 <- ggplot(FungalRichness, aes(x=RR, y=RRYield)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm, color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") +
  theme_bw() +
  theme(text = element_text(family = "serif", size=20)) +
  geom_vline(aes(xintercept=0), colour="black", linewidth=0.5, linetype="dashed") +
  labs(x="RR", y="RRYield34") +
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
  scale_x_continuous(limits=c(-0.65, 0.26), expand=c(0, 0)) + 
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
                            data = FungalRichness, 
                            na.action = na.roughfix, 
                            importance = TRUE, 
                            ntree = 500)
# Check the importance of variables and p-values
importance(rf_model_perm3)
#                    %IncMSE %IncMSE.pval IncNodePurity IncNodePurity.pval
# Longitude        11.720415   0.12871287     1.1441629         0.14851485
# Duration         11.710828   0.04950495     0.5811545         0.71287129
# Rotation_cycles  11.031602   0.04950495     0.5418028         0.80198020
# Latitude         10.666564   0.25742574     0.8119515         0.59405941
# MATmean           7.754075   0.67326733     1.0583012         0.07920792
# pHCK              6.700432   0.36633663     1.4350207         0.07920792
# MAPmean           6.602336   0.66336634     0.9320210         0.34653465
# SOCCK             2.079937   0.66336634     1.1846493         0.60396040
# Species_Richness -4.251619   0.95049505     0.4762516         0.07920792

######################### Trials sorted by effect size
library(ggplot2)
library(dplyr)
FungalRichness <- read.csv("FungalRichness.csv", fileEncoding = "latin1")
# 计算95% CI + 显著性分类
df_plot <- FungalRichness %>%
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
#       59      145       63 
## 8*8


################################################ piecewiseSEM
sem_data <- read.csv("FungalRichness.csv", fileEncoding = "latin1")
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
# -352.401
# Fisher's C = 6.329 with P-value = 0.042

m1_1 <- lme(RR ~ RRpH + RRSOC + Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_2 <- lme(RRpH ~  Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m1_3 <- lme(RRSOC ~ Duration + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model1 <- psem(m1_1, m1_2, m1_3)
summary(sem_model1)  ## 
# AIC
# -371.251
# Fisher's C = 8.81 with P-value = 0.012

m2_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_2 <- lme(RRpH ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m2_3 <- lme(RRSOC ~ MATmean + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model2 <- psem(m2_1, m2_2, m2_3)
summary(sem_model2)  ## 
# AIC
# -384.963
# Fisher's C = 5.777 with P-value = 0.056

m3_1 <- lme(RR ~ RRpH + RRSOC + MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_2 <- lme(RRpH ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m3_3 <- lme(RRSOC ~ MATmean + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model3 <- psem(m3_1, m3_2, m3_3)
summary(sem_model3)  ## 
# AIC
# -365.282
# Fisher's C = 9.507 with P-value = 0.009

m4_1 <- lme(RR ~ RRpH + RRSOC + Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_2 <- lme(RRpH ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m4_3 <- lme(RRSOC ~ Species_Richness, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model4 <- psem(m4_1, m4_2, m4_3)
summary(sem_model4) ## 
# AIC
# -405.992
# Fisher's C = 8.132 with P-value = 0.017 

# Structural Equation Model of sem_model4 
# 
# Call:
#   RR ~ RRpH + RRSOC + Species_Richness
# RRpH ~ Species_Richness
# RRSOC ~ Species_Richness
# 
# AIC
# -405.992
# 
# ---
#   Tests of directed separation:
#   
#   Independ.Claim Test.Type DF Crit.Value P.Value  
# RRSOC ~ RRpH + ...      coef 80     2.4343  0.0171 *
#   
#   --
#   Global goodness-of-fit:
#   
#   Chi-Squared = NA with P-value = NA and on 1 degrees of freedom
# Fisher's C = 8.132 with P-value = 0.017 and on 2 degrees of freedom
# 
# ---
# Coefficients:
# 
#   Response        Predictor Estimate Std.Error DF Crit.Value P.Value Std.Estimate    
#         RR             RRpH   0.4847    0.3486 79     1.3904  0.1683       0.1520    
#         RR            RRSOC   0.4571    0.1287 79     3.5524  0.0006       0.4375 ***
#         RR Species_Richness  -0.0399    0.0328 79    -1.2152  0.2279      -0.1088    
#       RRpH Species_Richness   0.0105    0.0096 83     1.0941  0.2771       0.0911    
#      RRSOC Species_Richness   0.0512    0.0247 95     2.0760  0.0406       0.1459   *
# 
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05
# 
# ---
# Individual R-squared:
# 
#   Response method Marginal Conditional
#         RR   none     0.16        0.21
#       RRpH   none     0.01        0.60
#      RRSOC   none     0.02        0.76


m5_1 <- lme(RR ~ RRpH + RRSOC + Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_2 <- lme(RRpH ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m5_3 <- lme(RRSOC ~ Duration, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model5 <- psem(m5_1, m5_2, m5_3)
summary(sem_model5) ## 
# AIC
# -385.741
# Fisher's C = 11.366 with P-value = 0.003

m6_1 <- lme(RR ~ RRpH + RRSOC +MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_2 <- lme(RRpH ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
m6_3 <- lme(RRSOC ~ MATmean, random = ~1|StudyID, data = sem_data, na.action = na.exclude)
sem_model6 <- psem(m6_1, m6_2, m6_3)
summary(sem_model6)  ##
# AIC
# -401.573
# Fisher's C = 8.351 with P-value = 0.015
