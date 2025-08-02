
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)

FungalSimpson <- read.csv("FungalSimpson.csv", fileEncoding = "latin1")
# Check data
head(FungalSimpson)

# 1. The number of Obversation
total_number <- nrow(FungalSimpson)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 47

# 2. The number of Study
unique_studyid_number <- length(unique(FungalSimpson$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 22


#### 3. Overall effect size
total_effect_model <- rma.mv(yi = RR, 
                             V = Vi, 
                             random = ~ 1 | StudyID,  # StudyID is radom factor
                             data = FungalSimpson, 
                             method = "REML")
# The results of Overall effect size
summary(total_effect_model)
# Multivariate Meta-Analysis Model (k = 47; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -94.8804  189.7608  193.7608  197.4181  194.0399   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1433  0.3786     22     no  StudyID 
# Test for Heterogeneity:
# Q(df = 46) = 716.4577, p-val < .0001
# Model Results:
# estimate      se    zval    pval    ci.lb   ci.ub    
#   0.0201  0.0843  0.2383  0.8117  -0.1451  0.1852      

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
boot_results1 <- boot(data = FungalSimpson, statistic = boot_fun, R = 1000, parallel = "snow", ncpus = numCores, cl = cl)
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
# Estimate for coefficient 1 : 0.02007927 
# 95% BCa CI for coefficient 1 : -0.07636727 0.1156135  


#### 5. Funnel Plot
simple_model <- rma(yi = RR, 
                    vi = Vi, 
                    data = FungalSimpson, 
                    method = "REML")
#### 
funnel(simple_model)
# Output  6 * 6
#### Egger's test
regtest(simple_model)
# Regression Test for Funnel Plot Asymmetry
# Model:     mixed-effects meta-regression model
# Predictor: standard error
# Test for Funnel Plot Asymmetry: z =  0.5060, p = 0.6128
# Limit Estimate (as sei -> 0):   b = -0.0020 (CI: -0.1511, 0.1471)

#  Rosenthal’s Fail-Safe N
# This method estimates how many missing studies with null effect 
# would be needed to make the overall effect non-significant
fsn_rosenthal <- fsn(x = simple_model, type = "Rosenthal")
# Print the FSN result
print(fsn_rosenthal)
# Fail-safe N Calculation Using the General Approach
# Average Effect Size:         0.0271
# Amount of Heterogeneity:     0.0983
# Observed Significance Level: 0.5838
# Target Significance Level:   0.05
# Fail-safe N: 0


library(ggplot2)
library(dplyr)
FungalSimpson <- read.csv("FungalSimpson.csv", fileEncoding = "latin1")
# 计算95% CI + 显著性分类
df_plot <- FungalSimpson %>%
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
#        3       32       12 
## 8*8
