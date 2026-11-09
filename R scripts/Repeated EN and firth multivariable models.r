
# REPEATED ELASTIC NET ================================

#Load libraries
library(tidyverse)
library(glmnet)
library(logistf)
library(ComplexHeatmap)
library(grid)
library(ggplot2)
library(pROC)
library(purrr)


# Data preparation --------------------------------------------

# File path
file_path <- "data/Combined_nasal_serum_mastersheet_final.csv"

df <- read.csv(file_path)

dir.create("Outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("Outputs/tables", recursive = TRUE, showWarnings = FALSE)


# Show number with quantifiable viral shedding vs without
table(df$Quantifiable_shedding)


# Filter for 35 candidate SARS-CoV-2 specific immune markers

markers <- df[,grepl("BL", names(df)) & !grepl("ORF6|ORF10|Alpha|Omicron|Delta|Beta|Gamma|OC43|229E|HKU|CoV_1|NL63",names(df))]

# Uncomment to exclude Nasal IgG RBD and S from the analysis
#markers <- df[,grepl("BL", names(df)) & !grepl("Nasal_IgG_CoV_2_S|Nasal_IgG_CoV_2_RBD|ORF6|ORF10|Alpha|Omicron|Delta|Beta|Gamma|OC43|229E|HKU|CoV_1|NL63",names(df))]

markers <- markers %>%
  mutate_if(is.character, as.numeric)


# Log10 transform
log_markers <- log10(markers)


# Scale (centre and standardise) immune markers
scaled_markers <- scale(
  log_markers,
  center = TRUE,
  scale = TRUE
)


# Combine with outcome and dose
EN_df <- data.frame(
  Quantifiable_shedding = as.numeric(df$Quantifiable_shedding),
  Dose_numeric = as.numeric(df$Dose),
  scaled_markers
)


# Make the scaled markers and dose into a matrix
X_full <- data.frame(
  Dose_numeric = EN_df$Dose_numeric,
  scaled_markers
)

X <- model.matrix(~ ., data = X_full)[, -1]

Y <- as.numeric(EN_df$Quantifiable_shedding)

colnames(X)
ncol(X)
nrow(X)


# Elastic net -------------------------------------------------
# 200 iterations of Elastic Net with dose as an unpenalised covariate, AUC as the performance metric, 4 folds, alpha = 0.8.

n_iter <- 200

selected_markers_list <- vector("list", n_iter)

# Store the cross-validated AUC from each iteration
en_cv_auc <- numeric(n_iter)


penalty_vec <- rep(1, ncol(X))

penalty_vec[1] <- 0  # Leave dose unpenalised (included in every run)


for (i in 1:n_iter) {
  
  set.seed(i)
  
  fit <- cv.glmnet(
    x = as.matrix(X),
    y = Y,
    family = "binomial",
    alpha = 0.8,
    nfolds = 4,
    standardize = FALSE, # Done previously
    penalty.factor = penalty_vec,
    type.measure = "auc"
  )
  
  
  # Extract markers selected at lambda.min
  coefs <- coef(fit, s = "lambda.min")
  
  selected <- rownames(coefs)[coefs[, 1] != 0]
  
  selected_markers_list[[i]] <- setdiff(
    selected,
    "(Intercept)"
  )
  
  
  # Extract the CV AUC at lambda.min
  lambda_index <- which.min(
    abs(fit$lambda - fit$lambda.min)
  )
  
  en_cv_auc[i] <- fit$cvm[lambda_index]
  
}


# Summary statistics across the 200 iterations

auc_summary <- data.frame(
  
  Iterations = n_iter,
  
  Mean_AUC = mean(en_cv_auc),
  
  SD_AUC = sd(en_cv_auc),
  
  Median_AUC = median(en_cv_auc),
  
  Q1 = unname(quantile(en_cv_auc, 0.25)),
  
  Q3 = unname(quantile(en_cv_auc, 0.75)),
  
  P2.5 = unname(quantile(en_cv_auc, 0.025)),
  
  P97.5 = unname(quantile(en_cv_auc, 0.975)),
  
  Minimum = min(en_cv_auc),
  
  Maximum = max(en_cv_auc)
  
)

print(
  auc_summary,
  digits = 2,
  row.names = FALSE
)


cat("\nMedian CV AUC (IQR): ")

cat(
  sprintf(
    "%.2f (%.2f–%.2f)\n",
    median(en_cv_auc),
    quantile(en_cv_auc, 0.25),
    quantile(en_cv_auc, 0.75)
  )
)


cat("\nMean elastic-net CV AUC (95% empirical interval): ")

cat(
  sprintf(
    "%.2f (%.2f–%.2f)\n",
    mean(en_cv_auc),
    quantile(en_cv_auc, 0.025),
    quantile(en_cv_auc, 0.975)
  )
)



#Marker selection frequency analysis ----------------------------

# Count number of failed (empty) runs

base_vars <- c("Dose_numeric")

immune_markers_only <- lapply(
  selected_markers_list,
  setdiff,
  y = base_vars
)


failed_runs <- sum(
  sapply(immune_markers_only, length) == 0
)

successful_runs <- n_iter - failed_runs


cat("\nNumber of runs where no immune markers were selected (failed runs)\n")

print(failed_runs)


cat("\nNumber of runs where at least one immune marker was selected (successful runs)\n")

print(successful_runs)


# Count number of runs where immune markers were selected

var_freq_pct <- immune_markers_only %>%
  unlist() %>%
  table() %>%
  `/`(n_iter) %>%
  `*`(100) %>%
  sort(decreasing = TRUE)


cat("\nPercentage of runs in which each immune marker was selected\n")

print(var_freq_pct)


#Plot selection frequency -------------------------------------

var_freq_pct <- as.data.frame(var_freq_pct)


colnames(var_freq_pct)[1:2] <- c(
  "Marker",
  "Frequency"
)


#replace variable names with antigen names
var_freq_pct <- var_freq_pct %>%
  mutate(
    Marker = Marker %>%
      stringr::str_remove_all("_BL|CoV_2") %>%
      stringr::str_replace_all("_", " ") 
  )

#Take top 10 markers
top10_markers <- var_freq_pct %>%
  slice_max(Frequency, n = 10, with_ties = FALSE)



#Set ggplot theme

panel_theme <- theme_classic(base_size = 5) +
  theme(
    plot.title = element_blank(),
    text = element_text(colour = "black"),
    
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 5),
    
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 7),
    
    plot.margin = margin(5, 5, 5, 5, unit = "mm")
  )



#Plot
EN_bar_chart <- ggplot(
  top10_markers,
  aes(x = reorder(Marker, Frequency), y = Frequency)
) +
  
  geom_col(
    fill = "#4C72B0",
    width = 0.8
  ) +
  
  # Percentage labels 
  geom_text(
    aes(label = sprintf("%.1f%%", Frequency)),
    hjust = -0.15,
    size = 5,
    size.unit ="pt"
  ) +
  
  coord_flip(clip = "off") +
  
  scale_y_continuous(
    limits = c(0, 120),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    title = NULL,
    x = NULL,
    y = "Selection frequency (%)"
  ) +
  
  panel_theme +
  
  theme(
    panel.grid.major.y = element_blank()
  )

print(EN_bar_chart)

# Export to outputs
ggsave(
  "Outputs/figures/EN_top10.svg",
  plot = EN_bar_chart,
  width =90,
  height = 72,
  units = "mm"
)




# 7. Correlation heatmap for the top markers ---------------------------------

top_selected_markers <- c(
  "Serum_IgA_CoV_2_N_BL",
  "Nasal_IgG_CoV_2_S_BL",
  "Nasal_IgG_CoV_2_RBD_BL",
  "Nasal_IgA_CoV_2_N_BL"
)

marker_labels <- c(
  "Serum IgA N",
  "Nasal IgG S",
  "Nasal IgG RBD",
  "Nasal IgA N"
)


cor_data_top <- dplyr::select(
  EN_df,
  dplyr::all_of(top_selected_markers)
)


cor_matrix_top <- cor(
  cor_data_top,
  method = "spearman",
  use = "pairwise.complete.obs"
)

# Replace matrix labels with descriptive marker names
rownames(cor_matrix_top) <- marker_labels
colnames(cor_matrix_top) <- marker_labels


# Colour palette

col_palette <- colorRampPalette(
  c(
    "#B2182B",
    "#EF8A62",
    "#FFFFFF",
    "#67A9CF",
    "#2166AC"
  )
)(100)


# Fixed breaks from -1 to 1

breaks <- seq(
  -1,
  1,
  length.out = 101
)



#plot
EN_top_heatmap <- ComplexHeatmap::pheatmap(
  mat = cor_matrix_top,
  color = col_palette,
  breaks = breaks,
  
  name = "Spearman's \u03c1",
  main = NULL,
  
  display_numbers = TRUE,
  number_color = "black",
  
  fontsize = 5,
  fontsize_row = 5,
  fontsize_col = 5,
  fontsize_number = 5,
  
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontsize = 7),
    labels_gp = grid::gpar(fontsize = 5)
  ),
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  # set dendrogram
  treeheight_row = 20,
  treeheight_col = 20,
  
  border_color = NA
)


#Save as svg
svg(
  "Outputs/figures/EN_correlations.svg",
  width = 90 / 25.4,
  height = 70 / 25.4,
  pointsize = 8
)

grid::grid.newpage()

ComplexHeatmap::draw(
  EN_top_heatmap,
  heatmap_legend_side = "right",
  padding = grid::unit(c(2, 2, 2, 2), "mm")
)

dev.off()


# Display heatmap in RStudio

dev.new(width = 5, height = 5)

ComplexHeatmap::draw(
  EN_top_heatmap,
  heatmap_legend_side = "right",
  newpage = TRUE
)


# Firth regression on selected markers --------------------------------------

# Nasal IgG S has been removed due to extreme collinearity
# with Nasal IgG RBD.

top_model <- logistf(
  Quantifiable_shedding ~
    Dose_numeric +
    Nasal_IgG_CoV_2_RBD_BL +
    Serum_IgA_CoV_2_N_BL +
    Nasal_IgA_CoV_2_N_BL,
  
  data = EN_df
)


summary(top_model)



#  Extract coefficients and confidence intervals

top_model_df <- data.frame(
  term = names(coef(top_model)),
  estimate = coef(top_model),
  conf.low = top_model$ci.lower,
  conf.high = top_model$ci.upper,
  p.value = top_model$prob
)


# Remove intercept and convert log-odds to odds ratios

top_model_df <- top_model_df %>%
  
  filter(term != "(Intercept)") %>%
  
  mutate(
    
    OR = exp(estimate),
    
    OR_low = exp(conf.low),
    
    OR_high = exp(conf.high),
    
    term_label = recode(
      term,
      
      "Dose_numeric" = "Dose",
      
      "Nasal_IgG_CoV_2_RBD_BL" =
        "Nasal IgG RBD",
      
      "Serum_IgA_CoV_2_N_BL" =
        "Serum IgA N",
      
      "Nasal_IgA_CoV_2_N_BL" =
        "Nasal IgA N"
    ),
    
    term_label = factor(
      term_label,
      levels = rev(term_label)
    )
  )

print(top_model_df)


# Firth forest plot

firth_forest <- ggplot(
  top_model_df,
  aes(
    x = OR,
    y = term_label
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey40"
  ) +
  
  geom_errorbarh(
    aes(
      xmin = OR_low,
      xmax = OR_high
    ),
    height = 0.18,
    linewidth = 0.5
  ) +
  
  geom_point(
    size = 2,
    colour = "#0072B2"
  ) +
  
  scale_x_log10(
    limits = c(0.1, 11)
  ) +
  
  labs(
    x = "Adjusted odds ratio",
    y = NULL
  ) +
  
  theme_classic(base_size = 5) +
  
  theme(
    plot.title = element_text(
      size = 7,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      colour = "black"
    ),
    
    axis.text.x = element_text(
      colour = "black"
    )
  )


print(firth_forest)


# Save svg
ggsave(
  filename = "Outputs/figures/firth_forest.svg",
  plot = firth_forest,
  width = 3,
  height = 2,
  units = "in"
)



# Cross-validated  ROC analysis -------------------------------------


# Define models

cv_models <- list(
  
  "Dose only" =
    Quantifiable_shedding ~ Dose_numeric,
  
  "Dose + Serum IgA N" =
    Quantifiable_shedding ~
    Dose_numeric +
    Serum_IgA_CoV_2_N_BL,
  
  "Dose + Nasal IgG RBD" =
    Quantifiable_shedding ~
    Dose_numeric +
    Nasal_IgG_CoV_2_RBD_BL,
  
  "Dose + Nasal IgA N" =
    Quantifiable_shedding ~
    Dose_numeric +
    Nasal_IgA_CoV_2_N_BL,
  
  "Full model" =
    Quantifiable_shedding ~
    Dose_numeric +
    Serum_IgA_CoV_2_N_BL +
    Nasal_IgG_CoV_2_RBD_BL +
    Nasal_IgA_CoV_2_N_BL
)

#Settings and storage

n_repeats <- 200
k <- 4

model_names <- names(cv_models)

y <- EN_df$Quantifiable_shedding

stopifnot(
  !anyNA(EN_df[, all.vars(cv_models[["Full model"]])]),
  all(y %in% c(0, 1)),
  all(table(y) >= k)
)

# One AUROC per model per CV repetition

firth_cv_auc <- matrix(
  NA_real_,
  nrow = n_repeats,
  ncol = length(cv_models),
  dimnames = list(NULL, model_names)
)

# Store all out-of-fold predictions

oof_predictions <- array(
  NA_real_,
  dim = c(nrow(EN_df), length(cv_models), n_repeats),
  dimnames = list(NULL, model_names, NULL)
)


# Repeated four-fold cross-validation

for (r in seq_len(n_repeats)) {
  
  set.seed(20261003 + r)
  
  # Generate folds, identical allocations are used for all 5 models.
  
  
  fold_id <- sample(
    rep(seq_len(k), length.out = nrow(EN_df))
  )
  
  
  # Fit the 5 models
  
  for (m in model_names) {
    
    predictions <- rep(NA_real_, nrow(EN_df))
    
    for (fold in seq_len(k)) {
      
      train_idx <- which(fold_id != fold)
      test_idx  <- which(fold_id == fold)
      
      # Firth regression: training data only
      
      fit <- tryCatch(
        
        withCallingHandlers(
          
          logistf(
            formula = cv_models[[m]],
            data = EN_df[train_idx, , drop = FALSE],
            pl = FALSE,
            
            control = logistf.control(
              maxit = 1000,
              maxstep = 2
            )
          ),
          
          warning = function(w) {
            
            if (grepl(
              "maximum number of iterations|did not converge|non.converg",
              conditionMessage(w),
              ignore.case = TRUE
            )) {
              stop(conditionMessage(w))
            }
          }
        ),
        
        error = function(e) {
          
          stop(
            sprintf(
              "CV failed: repetition %d, model '%s', fold %d: %s",
              r, m, fold, conditionMessage(e)
            ),
            call. = FALSE
          )
        }
      )
      
      # Predict held-out participants
      
      predictions[test_idx] <- as.numeric(
        predict(
          fit,
          newdata = EN_df[test_idx, , drop = FALSE],
          type = "response"
        )
      )
    }
    
    stopifnot(
      all(is.finite(predictions)),
      all(predictions >= 0 & predictions <= 1)
    )
    
    # Save out-of-fold predictions
    
    oof_predictions[, m, r] <- predictions
    
    # Empirical ROC and AUROC for this repetition
    
    roc_object <- pROC::roc(
      response = y,
      predictor = predictions,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )
    
    firth_cv_auc[r, m] <- as.numeric(
      pROC::auc(roc_object)
    )
  }
  
  if (r %% 25 == 0) {
    cat("Completed", r, "/", n_repeats, "\n")
  }
}

stopifnot(
  all(is.finite(firth_cv_auc)),
  all(is.finite(oof_predictions))
)


# Summarise cross-validated AUROC


cv_summary <- data.frame(
  
  Model = model_names,
  
  Mean_AUC = colMeans(firth_cv_auc),
  
  P2.5 = apply(
    firth_cv_auc, 2, quantile, probs = 0.025
  ),
  
  P97.5 = apply(
    firth_cv_auc, 2, quantile, probs = 0.975
  ),
  
  Median_AUC = apply(firth_cv_auc, 2, median),
  
  Q1 = apply(firth_cv_auc, 2, quantile, probs = 0.25),
  
  Q3 = apply(firth_cv_auc, 2, quantile, probs = 0.75),
  
  row.names = NULL
)

# Formatted reporting
cv_summary$AUC_95_interval <- sprintf(
  "%.2f (%.2f–%.2f)",
  cv_summary$Mean_AUC,
  cv_summary$P2.5,
  cv_summary$P97.5
)

print(
  cv_summary[, c("Model", "AUC_95_interval")],
  row.names = FALSE
)


# 5 model CV AUC figure ----------------------------------------------------


# All five models
plot_models <- names(cv_models)



# Smooth each empirical CV ROC and average sensitivities

fpr_grid <- seq(0, 1, length.out = 201)

cv_roc_df <- purrr::map_dfr(
  
  plot_models,
  
  function(m) {
    
    sensitivity_matrix <- sapply(
      seq_len(n_repeats),
      
      function(r) {
        
        # Empirical ROC from out-of-fold predictions
        roc_object <- pROC::roc(
          response = y,
          predictor = oof_predictions[, m, r],
          levels = c(0, 1),
          direction = "<",
          quiet = TRUE
        )
        
        # Binormal smoothing
        sm <- pROC::smooth(
          roc_object,
          method = "binormal"
        )
        
        sensitivity <- pROC::coords(
          sm,
          x = 1 - fpr_grid,
          input = "specificity",
          ret = "sensitivity",
          transpose = FALSE
        )
        
        as.numeric(unlist(sensitivity, use.names = FALSE))
      }
    )
    
    mean_sensitivity <- rowMeans(sensitivity_matrix)
    
    mean_sensitivity[1] <- 0
    mean_sensitivity[length(mean_sensitivity)] <- 1
    
    tibble(
      Model = m,
      FPR = fpr_grid,
      Sensitivity = mean_sensitivity
    )
  }
)


# Set Colours and line types

roc_colours <- c(
  "Dose only" = "grey45",
  "Dose + Serum IgA N" = "#E69F00",
  "Dose + Nasal IgG RBD" = "#0072B2",
  "Dose + Nasal IgA N" = "#CC79A7",
  "Full model" = "#D55E00" 
)

roc_linetypes <- c(
  "Dose only" = "dashed",
  "Dose + Serum IgA N" = "solid",
  "Dose + Nasal IgG RBD" = "solid",
  "Dose + Nasal IgA N" = "solid",
  "Full model" = "solid"
)



# Set legend

legend_labels <- setNames(
  
  sprintf(
    "%s: AUC = %.2f (%.2f–%.2f)",
    cv_summary$Model,
    cv_summary$Mean_AUC,
    cv_summary$P2.5,
    cv_summary$P97.5
  ),
  
  cv_summary$Model
)

# Plot

CV_ROC_plot_5 <- ggplot(
  cv_roc_df,
  aes(
    x = FPR,
    y = Sensitivity,
    colour = Model,
    linetype = Model
  )
) +
  
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dotted",
    linewidth = 0.25,
    colour = "grey80"
  ) +
  
  geom_line(linewidth = 0.5) +
  
  scale_colour_manual(
    values = roc_colours,
    breaks = plot_models,
    labels = legend_labels[plot_models]
  ) +
  
  scale_linetype_manual(
    values = roc_linetypes,
    breaks = plot_models,
    labels = legend_labels[plot_models]
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 1, 0.2)
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 1, 0.2)
  ) +
  
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    colour = NULL,
    linetype = NULL
  ) +
  
  theme_classic(base_size = 5) +
  
  theme(
    legend.position = "right",
    text = element_text(colour = "black"),
    legend.text = element_text(size = 5),
    legend.key.width = grid::unit(1.5, "cm"),
    
    axis.title = element_text(size = 7),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.2
    )
  ) +
  
  guides(
    colour = guide_legend(ncol = 1),
    linetype = guide_legend(ncol = 1)
  )

print(CV_ROC_plot_5)


# Export SVG

ggsave(
  "Outputs/figures/Firth_CV_ROC_5_models.svg",
  plot = CV_ROC_plot_5,
  width = 6,
  height = 3
)


# Save complete analytical objects

saveRDS(
  list(
    firth_cv_auc = firth_cv_auc,
    oof_predictions = oof_predictions,
    cv_roc_df = cv_roc_df,
    cv_summary = cv_summary
  ),
  "Outputs/tables/Firth_CV_results.rds"
)


# Firth penalised likelihood-ratio tests


# Fit all five models using the complete cohort

firth_models <- lapply(
  cv_models,
  function(f) {
    
    logistf(
      formula = f,
      data = EN_df,
      pl = FALSE,
      
      control = logistf.control(
        maxit = 1000,
        maxstep = 2
      )
    )
  }
)

# Extract global and dose-conditional PLRT p-values

plrt_summary <- purrr::imap_dfr(
  firth_models,
  
  function(fit, model_name) {
    
    coefficient_names <- names(coef(fit))
    
    # Global PLRT:
    # test all predictors jointly against zero,
    # retaining the intercept.
    
    global_terms <- which(
      coefficient_names != "(Intercept)"
    )
    
    global_p <- logistftest(
      fit,
      test = global_terms,
      control = logistf.control(
        maxit = 1000,
        maxstep = 2
      )
    )$prob
    
    # Dose-conditional PLRT:
    # test immune marker coefficient(s) jointly,
    # retaining dose and intercept.
    
    immune_terms <- which(
      !coefficient_names %in%
        c("(Intercept)", "Dose_numeric")
    )
    
    dose_conditional_p <- if (length(immune_terms) == 0) {
      
      NA_real_
      
    } else {
      
      logistftest(
        fit,
        test = immune_terms,
        control = logistf.control(
          maxit = 1000,
          maxstep = 2
        )
      )$prob
    }
    
    tibble(
      Model = model_name,
      PLRT_vs_intercept = as.numeric(global_p),
      PLRT_vs_dose = as.numeric(dose_conditional_p)
    )
  }
)

print(
  as.data.frame(plrt_summary),
  digits = 4,
  row.names = FALSE
)
