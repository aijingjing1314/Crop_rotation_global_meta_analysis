
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)

FungalPD <- read.csv("FungalPD.csv", fileEncoding = "latin1")
# Check data
head(FungalPD)

# 1. The number of Obversation
total_number <- nrow(FungalPD)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 26

# 2. The number of Study
unique_studyid_number <- length(unique(FungalPD$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 10


#### 3. Overall effect size
total_effect_model <- rma.mv(yi = RR, 
                             V = Vi, 
                             random = ~ 1 | StudyID,  # StudyID is radom factor
                             data = FungalPD, 
                             method = "REML")
# The results of Overall effect size
summary(total_effect_model)
# Multivariate Meta-Analysis Model (k = 26; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -12.2765   24.5529   28.5529   30.9907   29.0984   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0199  0.1411     10     no  StudyID 
# Test for Heterogeneity:
# Q(df = 25) = 303.5335, p-val < .0001
# Model Results:
# estimate      se    zval    pval    ci.lb   ci.ub    
#   0.0590  0.0473  1.2476  0.2122  -0.0337  0.1517     

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
boot_results1 <- boot(data = FungalPD, statistic = boot_fun, R = 1000, parallel = "snow", ncpus = numCores, cl = cl)
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
# Estimate for coefficient 1 : 0.05899907 
# 95% BCa CI for coefficient 1 : -0.00971593 0.1058678  


#### 5. Funnel Plot
simple_model <- rma(yi = RR, 
                    vi = Vi, 
                    data = FungalPD, 
                    method = "REML")
#### 
funnel(simple_model)
# Output  6 * 6
#### Egger's test
regtest(simple_model)
# Regression Test for Funnel Plot Asymmetry
# Model:     mixed-effects meta-regression model
# Predictor: standard error
# Test for Funnel Plot Asymmetry: z = 1.6507, p = 0.0988
# Limit Estimate (as sei -> 0):   b = 0.0043 (CI: -0.1100, 0.1186)


#  Rosenthal’s Fail-Safe N
# This method estimates how many missing studies with null effect 
# would be needed to make the overall effect non-significant
fsn_rosenthal <- fsn(x = simple_model, type = "Rosenthal")
# Print the FSN result
print(fsn_rosenthal)
# Fail-safe N Calculation Using the General Approach
# Average Effect Size:         0.0850 (with file drawer: 0.0436)
# Amount of Heterogeneity:     0.0188 (with file drawer: 0.0202)
# Observed Significance Level: 0.0055 (with file drawer: 0.0527)
# Target Significance Level:   0.05
# Fail-safe N: 21


library(ggplot2)
library(dplyr)
FungalPD <- read.csv("FungalPD.csv", fileEncoding = "latin1")
# 计算95% CI + 显著性分类
df_plot <- FungalPD %>%
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
#        2       13       11  
## 8*8
