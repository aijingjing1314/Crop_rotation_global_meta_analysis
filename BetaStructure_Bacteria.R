
library(dplyr)
library(tidyr)
library(purrr)
library(readxl)
library(writexl)  # For writing Excel files

# Read the Excel file (assumes the first sheet contains the data)
df <- read_excel("BetaStructure_Bacteria.xlsx")

# Define function to calculate pairwise Euclidean distances
calc_distances <- function(coords) {
  dist_mat <- as.matrix(dist(coords))
  dist_vec <- dist_mat[lower.tri(dist_mat)]
  return(dist_vec)
}

# Function to merge the mean and standard deviation of two independent groups
merge_mean_sd <- function(mean1, sd1, n1, mean2, sd2, n2) {
  n_total <- n1 + n2
  mean_total <- (mean1 * n1 + mean2 * n2) / n_total
  sd_total <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2 + 
                      (n1 * n2 / n_total) * (mean1 - mean2)^2) / (n_total - 1))
  return(list(mean = mean_total, sd = sd_total, n = n_total))
}

# Calculate distance-based metrics within and between groups
results <- df %>%
  group_by(Observation) %>%
  group_modify(~{
    data_obs <- .x
    data_ck <- data_obs %>% filter(Treatment == "CK") %>% select(Axis1, Axis2)
    data_tr <- data_obs %>% filter(Treatment == "TR") %>% select(Axis1, Axis2)
    
    dist_ck <- if(nrow(data_ck) > 1) calc_distances(data_ck) else numeric(0)
    dist_tr <- if(nrow(data_tr) > 1) calc_distances(data_tr) else numeric(0)
    
    dist_btwn <- if(nrow(data_ck) > 0 && nrow(data_tr) > 0) {
      dist_mat <- as.matrix(dist(rbind(data_ck, data_tr)))
      dist_mat[1:nrow(data_ck), (nrow(data_ck)+1):(nrow(data_ck)+nrow(data_tr))] %>% as.vector()
    } else numeric(0)
    
    D_C_mean <- if(length(dist_ck) > 0) mean(dist_ck) else NA
    D_C_sd <- if(length(dist_ck) > 0) sd(dist_ck) else NA
    D_C_n <- length(dist_ck)
    
    D_N_mean <- if(length(dist_tr) > 0) mean(dist_tr) else NA
    D_N_sd <- if(length(dist_tr) > 0) sd(dist_tr) else NA
    D_N_n <- length(dist_tr)
    
    D_B_mean <- if(length(dist_btwn) > 0) mean(dist_btwn) else NA
    D_B_sd <- if(length(dist_btwn) > 0) sd(dist_btwn) else NA
    D_B_n <- length(dist_btwn)
    
    if(!is.na(D_C_mean) && !is.na(D_N_mean)) {
      merged <- merge_mean_sd(D_C_mean, D_C_sd, D_C_n, D_N_mean, D_N_sd, D_N_n)
      D_CN_mean <- merged$mean
      D_CN_sd <- merged$sd
      D_CN_n <- merged$n
    } else {
      D_CN_mean <- NA
      D_CN_sd <- NA
      D_CN_n <- NA
    }
    
    tibble(
      D_C_mean = D_C_mean,     # Within-group distance for CK
      D_C_sd = D_C_sd,
      D_C_n = D_C_n,
      D_N_mean = D_N_mean,     # Within-group distance for TR
      D_N_sd = D_N_sd,
      D_N_n = D_N_n,
      D_B_mean = D_B_mean,     # Between-group distance (CK vs TR)
      D_B_sd = D_B_sd,
      D_B_n = D_B_n,
      D_C_N_mean = D_CN_mean,  # Combined within-group distance (CK + TR)
      D_C_N_sd = D_CN_sd,
      D_C_N_n = D_CN_n
    )
  })

# View results
print(results)

# Save results to Excel file
write_xlsx(results, "Bacterial_Beta_Distance_results.xlsx")


















