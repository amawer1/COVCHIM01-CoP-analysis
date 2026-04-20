library(tidyverse)
library(glmnet)
library(logistf)
library(pheatmap)

# File path
file_path <- "data/Combined_nasal_serum_mastersheet_final.csv"

df <- read.csv(file_path)

#Show number with quantifiable viral shedding vs without
table(df$Quantifiable_shedding)

# Filter and log-transform immune markers
markers <- df[, grepl("BL", names(df)) & !grepl("Alpha|Omicron|Delta|Beta|Gamma|OC43|229E|HKU|CoV_1|NL63", names(df))]
markers<- markers %>% mutate_if(is.character, as.numeric)

#log10 transform
log_markers <- log10(markers)

#Scale (centre and standardise) immune markers
scaled_markers <- scale(log_markers, center = TRUE, scale = TRUE)

# Combine with outcome and Dose
EN_df <- data.frame(
  Quantifiable_shedding = as.numeric(df$Quantifiable_shedding),
  Dose_numeric = as.numeric(df$Dose),
  scaled_markers
)


# Make the matrix
X_full <- data.frame(
  Dose_numeric = EN_df$Dose_numeric,
  scaled_markers
)


X <- model.matrix(~ ., data = X_full)[, -1]
Y <- as.numeric(EN_df$Quantifiable_shedding)

colnames(X)
ncol(X)
nrow(X)

#200 iterations of EN with dose as an unpenalised covariate, auc as the performance metric, 4 folds, alpha set at 0.8
n_iter <- 200
selected_markers_list <- vector("list", n_iter)
penalty_vec <- rep(1, ncol(X)) 
penalty_vec[1] <- 0             # leave Dose unpenalised

for (i in 1:n_iter) {
  set.seed(i)
    fit <- cv.glmnet(
    x = as.matrix(X),
    y = Y,
    family = "binomial",
    alpha = 0.8,
    nfolds = 4,
    standardize = FALSE, #Standardisation is carried out previously
    penalty.factor = penalty_vec,
        type.measure = "auc"
  )
  coefs <- coef(fit, s = "lambda.min")
  selected <- rownames(coefs)[coefs[, 1] != 0]
  selected_markers_list[[i]] <- setdiff(selected, "(Intercept)")
}

# Count number of failed (empty) runs
base_vars <- c("Dose_numeric")

immune_markers_only <- lapply(selected_markers_list, setdiff, y = base_vars)

failed_runs <- sum(sapply(immune_markers_only, length) == 0)
successful_runs <- n_iter - failed_runs


cat("\nNumber of runs where no immune markers were selected (failed runs):\n")
print(failed_runs)

cat("\nNumber of runs where at least one variable was selected (successful runs):\n")
print(successful_runs)


#Count number of runs where immune markers were selected
var_freq_pct <- immune_markers_only %>%
  unlist() %>%
  table() %>%
  `/`(n_iter) %>%
  `*`(100) %>%
  sort(decreasing = TRUE)

cat("\nPercentage of runs in which each immune marker was selected:\n")
print(var_freq_pct)


#Plot marker selection frequency
options(repr.plot.width = 10, repr.plot.height = 10)

var_freq_pct <- as.data.frame(var_freq_pct)
colnames(var_freq_pct)[1:2] <- c("Marker", "Frequency")

# Create the bar chart
ggplot(var_freq_pct, aes(x = reorder(Marker, Frequency), y = Frequency)) +
  geom_col(fill = "#4C72B0") +
  coord_flip() +  
  labs(
    title = "Marker Selection Frequency",
    x = "Marker",
    y = "Selection frequency (% of runs)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11)
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", Frequency)),
    hjust = -0.1,
    size = 3.5
  ) +
  ylim(0, max(var_freq_pct$Frequency) * 1.1)

#Plot selection frequency for selected markers in >40% of runs
options(repr.plot.width = 14, repr.plot.height = 14)

# The top selected markers 
EN_selected_markers <- c(
  "Serum_IgA_CoV_2_N_BL",
  "Nasal_IgG_CoV_2_S_BL",
  "Nasal_IgG_CoV_2_RBD_BL",
  "Nasal_IgA_CoV_2_N_BL",
  "ELISpot_ORF7_BL",
  "Serum_IgM_CoV_2_NTD_BL",
  "Serum_IgA_CoV_2_NTD_BL",
  "Serum_IgA_CoV_2_S_BL"
)

cor_data_EN_selected <- dplyr::select(EN_df, dplyr::all_of(EN_selected_markers))
cor_matrix_EN_selected <- cor(cor_data_EN_selected, method = "spearman", use = "pairwise.complete.obs")


# Colour palette
col_palette <- colorRampPalette(c(
  "#B2182B",  # strong negative
  "#EF8A62",  # moderate negative
  "#FFFFFF",  # zero
  "#67A9CF",  # moderate positive
  "#2166AC"   # strong positive
))(100)

# Fixed breaks from -1 to 1
breaks <- seq(-1, 1, length.out = 101)

pheatmap(
  mat = cor_matrix_EN_selected,
  color = col_palette,
  breaks = breaks,
  main = "Spearman correlation: EN selected markers",
  display_numbers = TRUE,
  number_color = "black",
  fontsize = 12,
  fontsize_number = 12,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA
)


#plot correlation heatmap for top markers
options(repr.plot.width = 14, repr.plot.height = 14)

# The top selected markers 
top_selected_markers <- c(
  "Serum_IgA_CoV_2_N_BL",
  "Nasal_IgG_CoV_2_S_BL",
  "Nasal_IgG_CoV_2_RBD_BL",
  "Nasal_IgA_CoV_2_N_BL"
)

cor_data_top <- dplyr::select(EN_df, dplyr::all_of(top_selected_markers))
cor_matrix_top <- cor(cor_data_top, method = "spearman", use = "pairwise.complete.obs")

# Colour palette
col_palette <- colorRampPalette(c(
  "#B2182B",  # strong negative
  "#EF8A62",  # moderate negative
  "#FFFFFF",  # zero
  "#67A9CF",  # moderate positive
  "#2166AC"   # strong positive
))(100)

# Fixed breaks from -1 to 1
breaks <- seq(-1, 1, length.out = 101)

pheatmap(
  mat = cor_matrix_top,
  color = col_palette,
  breaks = breaks,
  main = "Spearman correlation: top immune markers",
  display_numbers = TRUE,
  number_color = "black",
  fontsize = 12,
  fontsize_number = 12,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA
)

#Firth's regression on the selected markers

#Nasal_IgG_CoV_2_S_BL has been removed due to extreme colinearity with Nasal_IgG_CoV_2_RBD
top_model <- logistf(
  Quantifiable_shedding ~ Dose_numeric + Nasal_IgG_CoV_2_RBD_BL + Serum_IgA_CoV_2_N_BL + Nasal_IgA_CoV_2_N_BL,
  data = EN_df
)

summary(top_model)

options(repr.plot.width=5, repr.plot.height=5, repr.plot.res=300)

# Extract coefficients and CIs from logistf model
forest_df <- data.frame(
  term = names(coef(top_model)),
  estimate = coef(top_model),
  conf.low = top_model$ci.lower,
  conf.high = top_model$ci.upper,
  p.value = top_model$prob
)

# Remove intercept and convert log-odds to odds ratios
forest_df <- forest_df %>%
  filter(term != "(Intercept)") %>%
  mutate(
    OR = exp(estimate),
    OR_low = exp(conf.low),
    OR_high = exp(conf.high),
    term_label = recode(
      term,
      "Dose_numeric" = "Dose",
      "Nasal_IgG_CoV_2_RBD_BL" = "Nasal IgG CoV-2 RBD",
      "Serum_IgA_CoV_2_N_BL" = "Serum IgA CoV-2 N",
      "Nasal_IgA_CoV_2_N_BL" = "Nasal IgA CoV-2 N"
    ),
    term_label = factor(term_label, levels = rev(term_label))
  )

# Forest plot
ggplot(forest_df, aes(x = OR, y = term_label)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_errorbarh(
    aes(xmin = OR_low, xmax = OR_high),
    height = 0.18,
    linewidth = 0.8
  ) +
  geom_point(
    size = 3,
    colour = "#0072B2"
  ) +
  scale_x_log10(limits = c (0.1,11)) +
  labs(
    x = "Odds ratio for quantifiable shedding",
    y = NULL,
    title = "Firth logistic regression model of top immune markers and dose"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black")
  )

