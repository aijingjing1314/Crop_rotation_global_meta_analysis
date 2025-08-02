library(readxl)
library(dplyr)
library(writexl)

# Read the data
df <- read_excel("Bacteria_Richness.xlsx")

# Force conversion of RR and Vi columns to numeric type
df <- df %>%
  mutate(across(c(Chao1_RR, Chao1_Vi, ACE_RR, ACE_Vi, ASV_RR, ASV_Vi), as.numeric))

# Fixed-effect model combination function (automatically skips NA)
combine_fixed <- function(rr, vi) {
  valid <- !is.na(rr) & !is.na(vi)
  rr <- rr[valid]
  vi <- vi[valid]
  
  if (length(rr) == 0) {
    return(c(NA, NA))
  } else if (length(rr) == 1) {
    return(c(rr[1], vi[1]))
  } else {
    w <- 1 / vi
    rr_combined <- sum(w * rr) / sum(w)
    vi_combined <- 1 / sum(w)
    return(c(rr_combined, vi_combined))
  }
}

# Apply the function to each row
result <- df %>%
  rowwise() %>%
  mutate(
    combined = list(combine_fixed(
      rr = c(Chao1_RR, ACE_RR, ASV_RR),
      vi = c(Chao1_Vi, ACE_Vi, ASV_Vi)
    )),
    RR_combined = combined[[1]],
    Vi_combined = combined[[2]]
  ) %>%
  ungroup() %>%
  select(Beta_ID, RR_combined, Vi_combined)

# Print the result
print(result)

# Save the result as an Excel file
write_xlsx(result, "FixedEffect_BacterialRichness_Combined.xlsx")


######################  Visualization
library(tidyr)
library(ggplot2)
library(patchwork)

# Read and convert again for plotting
df <- read_excel("Bacteria_Richness.xlsx")

df <- df %>%
  mutate(across(c(ASV_RR, ASV_Vi, Chao1_RR, Chao1_Vi, ACE_RR, ACE_Vi), as.numeric))

# Convert to long format
long_df <- df %>%
  select(Beta_ID, ASV = ASV_RR, ASV_Vi = ASV_Vi,
         Chao1 = Chao1_RR, Chao1_Vi,
         ACE = ACE_RR, ACE_Vi) %>%
  pivot_longer(cols = -Beta_ID, names_to = "Metric", values_to = "Value") %>%
  mutate(
    Type = ifelse(grepl("_Vi$", Metric), "Vi", "RR"),
    Metric = gsub("_Vi$", "", Metric)
  ) %>%
  pivot_wider(names_from = Type, values_from = Value) %>%
  mutate(
    SE = sqrt(as.numeric(Vi)),
    CI_lower = RR - 1.96 * SE,
    CI_upper = RR + 1.96 * SE
  )

# Determine if any CIs overlap within the same Beta_ID
significance_summary <- long_df %>%
  group_by(Beta_ID) %>%
  summarise(
    any_overlap = {
      combn(1:n(), 2, function(idx) {
        i <- idx[1]; j <- idx[2]
        ci_i <- c(CI_lower[i], CI_upper[i])
        ci_j <- c(CI_lower[j], CI_upper[j])
        # Check if CIs overlap
        !(ci_i[2] < ci_j[1] || ci_j[2] < ci_i[1])
      }) %>% any()
    }
  ) %>%
  mutate(Significance = ifelse(any_overlap, "Not significant", "Significant"))

# Merge the significance result back to the long dataframe
long_df <- long_df %>%
  left_join(significance_summary, by = "Beta_ID")

# Forest plot: points + error bars
p1 <- ggplot(long_df, aes(x = RR, y = reorder(Beta_ID, RR), color = Metric)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.3, alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("ASV" = "#FF8080", "Chao1" = "#00BFC4", "ACE" = "#5757F9")) +
  theme_bw(base_size = 12) +
  labs(x = "Response ratio", y = "Observations", color = "Metric")

# Display plot
p1

# Save plot as PDF
ggsave("Bacteria_richness_forest_plot.pdf", p1, width = 8, height = 12)

 