library(tidyverse)
library(ggrepel)

# File path
file_path <-  "C:\\Users\\amawer\\GitHub\\COVCHIM01-CoP-analysis\\data\\Combined_nasal_serum_mastersheet_final.csv"

df <- read.csv(file_path)


#Show number of infected vs uninfected
table(df$Quantifiable_shedding)

# Filter immune variables, to remove seasonl CoVs and VoCs
immune_markers <- df[, grepl("BL", names(df)) & !grepl("Alpha|Delta|Omicron|Beta|Gamma|OC43|229E|HKU|CoV_1|NL63", names(df))]
immune_markers <- immune_markers %>% mutate_if(is.character, as.numeric)

#log10 transform and scale variables for comparison
scaled_markers <- scale(log10(immune_markers))

head(scaled_markers)

#Create the data frame
LR_df <- data.frame(
  Quantifiable_shedding = as.numeric(df$Quantifiable_shedding),
  Dose = as.numeric(df$Dose),
  scaled_markers
)

head(LR_df)

# Markers to scan
immune_markers <- colnames(scaled_markers)

# Empty list to store results
results_list <- vector("list", length(immune_markers))

# Loop through markers
for (i in seq_along(immune_markers)) {
  
  marker <- immune_markers[i]
  
  # Crude model
  crude_formula <- as.formula(
    paste("Quantifiable_shedding ~", marker)
  )
  crude_fit <- glm(crude_formula, data = LR_df, family = binomial())
  
  crude_lrt <- drop1(crude_fit, test = "LRT")
  crude_p <- crude_lrt[rownames(crude_lrt) == marker, "Pr(>Chi)"]
  
  crude_beta <- coef(crude_fit)[marker]
  crude_ci <- suppressMessages(confint(crude_fit))
  
  # Adjusted model
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

# Add formatted columns
results_df$crude_CI <- sprintf("%.2f (%.2f–%.2f)",
                               results_df$crude_OR,
                               results_df$crude_low,
                               results_df$crude_high)

results_df$adj_CI <- sprintf("%.2f (%.2f–%.2f)",
                             results_df$adj_OR,
                             results_df$adj_low,
                             results_df$adj_high)

# Optional BH correction on adjusted p-values
results_df$cor_p_BH <- p.adjust(results_df$adj_p, method = "BH")

# Order by adjusted p-value
results_df <- results_df[order(results_df$adj_p), ]

# Final display table
OR_table_display <- results_df[, c("Marker", "crude_CI", "crude_p", "adj_CI", "adj_p","cor_p_BH" )]
names(OR_table_display) <- c("Marker", "Crude OR (95% CI)", "Crude p", "Adj OR (95% CI)", "Adj p", "Corrected adjusted p")

OR_table_display

plot_df <- results_df %>%
  mutate(
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
  )%>%
arrange(adj_p)

# Set plot size (adjust width and height as needed)
options(repr.plot.width = 7, repr.plot.height = 8)

ggplot(plot_df, aes(x = adj_logOR, y = neglogp, label = Marker)) +
  geom_point(aes(color = sig_group), size = 2) +
  ggrepel::geom_text_repel(
    data = subset(plot_df, sig_group != "NS"),
    max.overlaps = Inf,
    force = 3
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  scale_color_manual(values = c(
    "BH corrected < 0.05" = "red",
    "Uncorrected p < 0.05" = "orange",
    "NS" = "grey50"
  )) +
  labs(
    x = "Adjusted log(OR)",
    y = "-log10(LRT p-value)",
    title = "Per-marker logistic regression adjusted for dose",
    subtitle = "(BL = baseline, prior to inoculation)",
    color = "Significance"
  ) +
theme_classic() +
  theme(
    legend.background = element_rect(fill = "white", colour = "black")
  )

# Set plot size (adjust width and height as needed)
options(repr.plot.width = 7, repr.plot.height = 8)

#To order markers by p value
forest_df <- plot_df %>%
  arrange(adj_p) %>%
  mutate(Marker = factor(Marker, levels = rev(unique(Marker))))

#plot
ggplot(forest_df, aes(y = Marker, x = adj_OR, xmin = adj_low, xmax = adj_high, color = sig_group)) +
  geom_point(size = 2) +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_log10(
    name = "Adjusted Odds Ratio (profile likelihood 95% CI)"
  ) +
  scale_color_manual(
    values = c(
      "BH corrected < 0.05" = "red",
      "Uncorrected p < 0.05" = "orange",
      "NS" = "grey50"
    ),
    name = "Significance"   
  ) +
  labs(
    y = NULL,
    title = "Per-marker logistic regression (profile likelihood CIs)",
    subtitle = "(BL = baseline, prior to inoculation)"
  ) +
  theme_classic() +
  theme(
    legend.background = element_rect(fill = "white", colour = "black")
  )
