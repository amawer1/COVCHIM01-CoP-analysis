# RIDGE REGRESSION

library(tidyverse)
library(glmnet)

# Data preparation ---------------------------------------------------

#Setwd
setwd("C:/Users/amawer/GitHub/COVCHIM01-CoP-analysis")

file_path <- "data/Combined_nasal_serum_mastersheet_final.csv"
df <- read.csv(file_path)

# Select the 12 serum and nasal anti-spike IgG/IgA markers
anti_spike_markers <- c(
  "Serum_IgG_CoV_2_S_BL", "Serum_IgG_CoV_2_RBD_BL",
  "Serum_IgG_CoV_2_NTD_BL",
  "Serum_IgA_CoV_2_S_BL", "Serum_IgA_CoV_2_RBD_BL",
  "Serum_IgA_CoV_2_NTD_BL",
  "Nasal_IgG_CoV_2_S_BL", "Nasal_IgG_CoV_2_RBD_BL",
  "Nasal_IgG_CoV_2_NTD_BL",
  "Nasal_IgA_CoV_2_S_BL", "Nasal_IgA_CoV_2_RBD_BL",
  "Nasal_IgA_CoV_2_NTD_BL"
)

markers <- df[, anti_spike_markers, drop = FALSE] %>%
  mutate(across(everything(), as.numeric))

log_markers <- log10(markers)
scaled_markers <- scale(log_markers, center = TRUE, scale = TRUE)

# Make the predictor matrix and outcome
X_full <- data.frame(Dose_numeric = as.numeric(df$Dose), scaled_markers)
X <- model.matrix(~ ., data = X_full)[, -1, drop = FALSE]
Y <- as.numeric(df$Quantifiable_shedding)
print(table(Y))

# Run 200 iterations of ridge regression with dose un-penalised
# 4 folds and AUC to choose lambda in each run
# Alpha = 0 for ridge

n_iter <- 200
cv_auc <- numeric(n_iter)
lambda_min_values <- numeric(n_iter)
penalty_vec <- rep(1, ncol(X))
penalty_vec[colnames(X) == "Dose_numeric"] <- 0

for (i in seq_len(n_iter)) {
  set.seed(i)
  fit <- cv.glmnet(
    x = X,
    y = Y,
    family = "binomial",
    alpha = 0,
    nfolds = 4,
    standardize = FALSE,
    penalty.factor = penalty_vec,
    type.measure = "auc"
  )
  lambda_min_values[i] <- fit$lambda.min
  lambda_index <- which.min(abs(fit$lambda - fit$lambda.min))
  cv_auc[i] <- fit$cvm[lambda_index]
}

# Summarise AUC across the 200 runs
auc_summary <- data.frame(
  Iterations = n_iter,
  Mean_AUC = mean(cv_auc),
  SD_AUC = sd(cv_auc),
  Median_AUC = median(cv_auc),
  Q1 = unname(quantile(cv_auc, 0.25)),
  Q3 = unname(quantile(cv_auc, 0.75)),
  P2.5 = unname(quantile(cv_auc, 0.025)),
  P97.5 = unname(quantile(cv_auc, 0.975)),
  Minimum = min(cv_auc),
  Maximum = max(cv_auc)
)
print(auc_summary, digits = 4, row.names = FALSE)

auc_results <- data.frame(Iteration = seq_len(n_iter), CV_AUC = cv_auc,
                          Lambda_min = lambda_min_values)

# Fit the final ridge model using the median lambda from all runs
# Compare coefficients

consensus_lambda <- median(lambda_min_values)

ridge_final <- glmnet(
  x = X, y = Y, family = "binomial", alpha = 0,
  lambda = consensus_lambda, standardize = FALSE,
  penalty.factor = penalty_vec
)

coef_matrix <- as.matrix(coef(ridge_final, s = consensus_lambda))
ridge_results <- data.frame(
  Marker = rownames(coef_matrix),
  Coefficient = as.numeric(coef_matrix[, 1]),
  row.names = NULL
)

# Report dose separately from the standardised immune markers
dose_result <- ridge_results %>%
  filter(Marker == "Dose_numeric") %>%
  mutate(OR = exp(Coefficient))

#Remove dose from the plot, keep the immune markers
ridge_immune_results <- ridge_results %>%
  filter(!Marker %in% c("(Intercept)", "Dose_numeric")) %>%
  mutate(Penalised_OR = exp(Coefficient)) %>%
  arrange(desc(abs(Coefficient)))

print(dose_result, digits = 6, row.names = FALSE)
print(ridge_immune_results, digits = 6, row.names = FALSE)

# Plot the ridge coefficients for all 12 immune markers
ridge_plot_data <- ridge_immune_results %>%
  mutate(Display_name = Marker %>%
           str_replace("_CoV_2_", " ") %>%
           str_replace("_BL$", "") %>%
           str_replace_all("_", " ")) %>%
  arrange(Coefficient) %>%
  mutate(Display_name = factor(Display_name, levels = rev(Display_name)))

ridge_coefficient_plot <- ggplot(
  ridge_plot_data, aes(x = Coefficient, y = Display_name)
) +
  geom_segment(aes(x = 0, xend = Coefficient,
                   y = Display_name, yend = Display_name), linewidth = 0.7) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.6) +
  labs(x = "Standardised ridge coefficient", y = NULL) +
  theme_classic(base_size = 5) +
  theme(axis.title = element_text(size = 7))

print(ridge_coefficient_plot)

# Save tables, figures and model results
dir.create("Outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("Outputs/figures", recursive = TRUE, showWarnings = FALSE)


# Save the coefficient plot
ggsave("Outputs/figures/Anti_spike_ridge_coefficients.svg",
       plot = ridge_coefficient_plot, width = 2.5, height = 2, units = "in")

saveRDS(list(markers = anti_spike_markers, n_repeats = n_iter,
             n_folds = 4, cv_auc = cv_auc,
             lambda_min_values = lambda_min_values,
             consensus_lambda = consensus_lambda, final_model = ridge_final,
             coefficients = ridge_immune_results, dose = dose_result,
             auc_summary = auc_summary),
        "Outputs/tables/Anti_spike_ridge_analysis.rds")
