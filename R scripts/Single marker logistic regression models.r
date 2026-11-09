library(tidyverse)
library(ggrepel)
library(writexl)

#Setwd
file_path <- "data/Combined_nasal_serum_mastersheet_final.csv"
df <- read.csv(file_path)

# Show counts with and without quantifiable shedding
print(table(df$Quantifiable_shedding))

# Select baseline markers, excluding seasonal coronaviruses and variants
immune_markers <- df[, grepl("BL", names(df)) & !grepl("ORF6|ORF10|Alpha|Delta|Omicron|Beta|Gamma|OC43|229E|HKU|CoV_1|NL63", names(df))]
immune_markers <- immune_markers %>% mutate(across(where(is.character), as.numeric))

# Log10-transform and standardise the markers
scaled_markers <- scale(log10(immune_markers))

# Combine markers with shedding status and dose
LR_df <- data.frame(
  Quantifiable_shedding = as.numeric(df$Quantifiable_shedding),
  Dose = as.numeric(df$Dose),
  scaled_markers
)

# Fit each marker separately, with and without adjustment for dose
marker_names <- colnames(scaled_markers)
results_list <- vector("list", length(marker_names))

for (i in seq_along(marker_names)) {
  
  marker <- marker_names[i]
  
  # Unadjusted model
  crude_formula <- as.formula(
    paste("Quantifiable_shedding ~", marker)
  )
  crude_fit <- glm(crude_formula, data = LR_df, family = binomial())
  
  crude_lrt <- drop1(crude_fit, test = "LRT")
  crude_p <- crude_lrt[rownames(crude_lrt) == marker, "Pr(>Chi)"]
  
  crude_beta <- coef(crude_fit)[marker]
  crude_ci <- suppressMessages(confint(crude_fit))
  
  # Adjust for dose
  adj_formula <- as.formula(
    paste("Quantifiable_shedding ~ Dose +", marker)
  )
  adj_fit <- glm(adj_formula, data = LR_df, family = binomial())
  
  adj_lrt <- drop1(adj_fit, test = "LRT")
  adj_p <- adj_lrt[rownames(adj_lrt) == marker, "Pr(>Chi)"]
  
  adj_beta <- coef(adj_fit)[marker]
  adj_ci <- suppressMessages(confint(adj_fit))
  
  # Store results
  results_list[[i]] <- data.frame(
    Marker = marker,
    
    crude_OR = exp(crude_beta),
    crude_low = exp(crude_ci[marker, 1]),
    crude_high = exp(crude_ci[marker, 2]),
    crude_p = as.numeric(crude_p),
    
    adj_OR = exp(adj_beta),
    adj_low = exp(adj_ci[marker, 1]),
    adj_high = exp(adj_ci[marker, 2]),
    adj_p = as.numeric(adj_p)
  )
}

# Combine results
results_df <- do.call(rbind, results_list)

# Format odds ratios and confidence intervals for the table
results_df$crude_CI <- sprintf("%.2f (%.2f–%.2f)",
                               results_df$crude_OR,
                               results_df$crude_low,
                               results_df$crude_high)

results_df$adj_CI <- sprintf("%.2f (%.2f–%.2f)",
                             results_df$adj_OR,
                             results_df$adj_low,
                             results_df$adj_high)

# Apply Benjamini-Hochberg correction to the dose-adjusted p-values
results_df$cor_p_BH <- p.adjust(results_df$adj_p, method = "BH")

# Order by adjusted p-value
results_df <- results_df[order(results_df$adj_p), ]

# Prepare the results table
OR_table_display <- results_df[, c("Marker", "crude_CI", "crude_p", "adj_CI", "adj_p","cor_p_BH" )]
names(OR_table_display) <- c("Marker", "Crude OR (95% CI)", "Crude p", "Adj OR (95% CI)", "Adj p", "Corrected adjusted p")

print(OR_table_display)

dir.create("Outputs/tables", recursive = TRUE, showWarnings = FALSE)

write_xlsx(
  OR_table_display,
  path = "Outputs/tables/single_marker_logistic_regression_OR_table.xlsx"
)

# Set significance colours
significance_colours <- c(
  "BH corrected < 0.05" = "red",
  "Uncorrected p < 0.05" = "orange",
  "NS" = "grey50"
)

plot_df <- results_df %>%
  mutate(
    Marker_label = Marker %>%
      stringr::str_remove_all("_BL|CoV_2_") %>%
      stringr::str_replace_all("_", " "),
    
    sig_group = case_when(
      cor_p_BH < 0.05 ~ "BH corrected < 0.05",
      adj_p  < 0.05 ~ "Uncorrected p < 0.05",
      TRUE            ~ "NS"
    ),
    sig_group = factor(
      sig_group,
      levels = c(
        "BH corrected < 0.05",
        "Uncorrected p < 0.05",
        "NS"
      )
    ),
    neglogp = -log10(adj_p),
    adj_logOR = log(adj_OR)
  )

# Plot adjusted log odds ratios against LRT p-values
volcano_plot <- ggplot(plot_df, aes(x = adj_logOR, y = neglogp, label = Marker_label)) +
  geom_point(aes(color = sig_group), size = 2) +
  ggrepel::geom_text_repel(
    data = subset(plot_df, sig_group != "NS"),
    size = 7 / ggplot2::.pt,
    max.overlaps = Inf,
    force = 10
  ) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  scale_color_manual(values = significance_colours) +
  labs(
    x = "Adjusted log(OR)",
    y = "-log10(LRT p-value)",
    color = "Significance"
  ) +
  theme_classic(base_size = 5) +
  theme(
    axis.text = element_text(size = 5),
    text = element_text(colour = "black"),
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.background = element_rect(
      fill = "white",
      colour = "black")
  )

print(volcano_plot)

# Save the figure
dir.create("Outputs/figures", recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = "Outputs/figures/single_marker_volcano_plot.svg",
  plot = volcano_plot,
  width = 7,
  height = 5,
  units = "in"
)

# Keep the markers in adjusted p-value order for the forest plot
forest_df <- plot_df %>%
  mutate(
    Marker_label = factor(
      Marker_label,
      levels = rev(unique(Marker_label))
    )
  )

# Plot adjusted odds ratios and confidence intervals
forest_plot <- ggplot(forest_df, aes(y = Marker_label, x = adj_OR, xmin = adj_low, xmax = adj_high, color = sig_group)) +
  geom_point(size = 2) +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_log10(
    name = "Adjusted Odds Ratio (profile likelihood 95% CI)"
  ) +
  scale_color_manual(
    values = significance_colours,
    name = "Significance"
  ) +
  labs(y = NULL) +
  theme_classic(base_size = 5) +
  theme(
    text = element_text(colour = "black"),
    
    axis.text.x = element_text(size = 7, colour = "black"),
    axis.text.y = element_text(size = 7, colour = "black"),
    axis.title.x = element_text(size = 7, colour = "black"),
    
    legend.text = element_text(size = 5, colour = "black"),
    legend.title = element_text(size = 7, colour = "black"),
    
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    
    legend.background = element_rect(
      fill = "white",
      colour = "black"
    )
  )

print(forest_plot)

# Save the figure
ggsave(
  filename = "Outputs/figures/forest_plot.svg",
  plot = forest_plot,
  width = 7,
  height = 8,
  units = "in"
)
