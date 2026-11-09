library(tidyverse)
library(ComplexHeatmap)
library(grid)
library(patchwork)

# Load the data
file_path <- "data/Combined_nasal_serum_mastersheet_final.csv"

df <- read.csv(file_path)

# Select baseline immune markers
markers <- df[, grepl("BL", names(df)) & !grepl("ORF10|ORF6|Alpha|Omicron|Delta|Beta|Gamma|OC43|229E|HKU|CoV_1|NL63", names(df))]
markers <- markers %>% mutate(across(where(is.character), as.numeric))

# Log10-transform the markers
log_markers <- log10(markers)

# Calculate Spearman correlations
R <- cor(
  log_markers,
  method = "spearman",
  use = "pairwise.complete.obs"
)

# Group by response type
marker_type <- dplyr::case_when(
  grepl("ACE_2|MNA", colnames(R)) ~ "Neutralisation",
  
  grepl("IgG", colnames(R)) &
    grepl("RBD|CoV_2_S|NTD", colnames(R)) &
    !grepl("_N_", colnames(R)) ~ "Anti-spike IgG antibody",
  
  grepl("IgA", colnames(R)) &
    grepl("RBD|CoV_2_S|NTD", colnames(R)) &
    !grepl("_N_", colnames(R)) ~ "Anti-spike IgA antibody",
  
  grepl("IgG", colnames(R)) &
    grepl("_N_", colnames(R)) ~ "Anti-nucleocapsid IgG antibody",
  
  grepl("IgA", colnames(R)) &
    grepl("_N_", colnames(R)) ~ "Anti-nucleocapsid IgA antibody",
  
  grepl("IgM", colnames(R)) ~ "IgM antibody",
  
  grepl("ELISpot", colnames(R)) ~ "T cell",
  
  TRUE ~ "Other"
)

names(marker_type) <- colnames(R)

# Colours
type_cols <- c(
  "Neutralisation" = "#9ECAE1",
  "Anti-spike IgG antibody" = "#08519C",
  "Anti-spike IgA antibody" = "#4292C6",
  "Anti-nucleocapsid IgG antibody" = "#A63603",
  "Anti-nucleocapsid IgA antibody" = "#FDAE6B",
  "IgM antibody" = "#737373",
  "T cell" = "#AA4499"
)

# Correlation colour scale
col_fun <- circlize::colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  c("#B2182B", "#EF8A62", "#FFFFFF", "#67A9CF", "#2166AC")
)

# Top annotation
top_ha <- HeatmapAnnotation(
  Type = marker_type,
  col = list(Type = type_cols),
  show_annotation_name = FALSE,
  simple_anno_size = unit(4, "mm")
)

# Set the marker type legend order
legend_order <- c(
  "Anti-spike IgG antibody",
  "Anti-spike IgA antibody",
  "Neutralisation",
  "Anti-nucleocapsid IgG antibody",
  "Anti-nucleocapsid IgA antibody",
  "T cell",
  "IgM antibody"
)

left_ha <- rowAnnotation(
  Type = marker_type,
  col = list(Type = type_cols),
  show_annotation_name = FALSE,
  simple_anno_size = unit(4, "mm"),
  annotation_legend_param = list(
    Type = list(
      at = legend_order,
      title = "Marker Type",
      title_gp = gpar(fontsize = 7),
      labels_gp = gpar(fontsize = 5)
    )
  )
)

# Plot the correlation heatmap
ht <- Heatmap(
  R,
  name = "Spearman's \u03c1",
  col = col_fun,
  
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = as.dist(1 - R),
  clustering_distance_columns = as.dist(1 - R),
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  
  top_annotation = top_ha,
  left_annotation = left_ha,
  
  show_row_names = FALSE,
  show_column_names = FALSE,
  
  border = FALSE,
  rect_gp = gpar(col = NA),
  
  heatmap_legend_param = list(
    title = "Spearman's \u03c1",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1"),
    direction = "horizontal",
    title_gp = gpar(fontsize = 7),
    labels_gp = gpar(fontsize = 7)
  )
)



# Save the figure
dir.create("Outputs/figures", recursive = TRUE, showWarnings = FALSE)

svg(
  file = "Outputs/figures/heatmap.svg",
  width = 8,
  height = 5
)

# Draw the heatmap into the SVG file
ComplexHeatmap::draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE,
  newpage = TRUE
)

dev.off()

# Display the heatmap in R
ComplexHeatmap::draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE,
  newpage = TRUE
)

# Match nasal and serum markers for the scatter plots

matched_marker_pairs <- tibble(
  original_name = colnames(log_markers)
) %>%
  mutate(
    compartment = str_extract(original_name, "Nasal|Serum"),
    marker_id = original_name %>%
      str_remove("Nasal|Serum") %>%
      str_replace_all("__+", "_") %>%
      str_remove("^_") %>%
      str_remove("_$")
  ) %>%
  filter(!is.na(compartment)) %>%
  pivot_wider(
    names_from = compartment,
    values_from = original_name
  ) %>%
  filter(!is.na(Nasal), !is.na(Serum)) %>%
  select(marker_id, Nasal, Serum)%>%
  mutate(
    isotype = str_extract(marker_id, "IgA|IgG|IgM"),
    antigen_order = case_when(
      str_detect(marker_id, "_S_BL$")   ~ 1,
      str_detect(marker_id, "_RBD_BL$") ~ 2,
      str_detect(marker_id, "_N_BL$")   ~ 3,
      str_detect(marker_id, "_NTD_BL$") ~ 4,
      TRUE ~ 5
    )
  ) %>%
  arrange(isotype, antigen_order) %>%
  select(-isotype, -antigen_order)

format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))
}

spearman_label <- function(data, label) {
  cor_result <- suppressWarnings(
    cor.test(data$nasal_value, data$serum_value, method = "spearman")
  )
  
  paste0(
    label, ", \u03C1 = ", round(unname(cor_result$estimate), 2),
    ", p = ", format_p(cor_result$p.value)
  )
}

make_matched_spearman_scatter <- function(marker_id, Nasal, Serum, data) {
  
  isotype <- str_extract(marker_id, "IgA|IgG|IgM")
  
  antigen <- case_when(
    str_detect(marker_id, "_RBD_BL$") ~ "RBD",
    str_detect(marker_id, "_NTD_BL$") ~ "NTD",
    str_detect(marker_id, "_S_BL$")   ~ "S",
    str_detect(marker_id, "_N_BL$")   ~ "N",
    TRUE ~ marker_id
  )
  
  plot_title <- paste(isotype, antigen)
  
  # Set axis limits for each isotype
  
  axis_settings <- switch(
    isotype,
    
    "IgA" = list(
      xlim = c(0, 5.5),
      ylim = c(-0.5, 4.5)
    ),
    
    "IgG" = list(
      xlim = c(0, 6),
      ylim = c(-0.5, 4.5)
    ),
    
    "IgM" = list(
      xlim = c(0, 4.5),
      ylim = c(-1, 2.5)
    )
  )
  
  # Add vaccination status to the paired measurements
  plot_data <- data %>%
    select(
      Vaccinated,
      nasal_value = all_of(Nasal),
      serum_value = all_of(Serum)
    ) %>%
    drop_na() %>%
    mutate(
      Vaccination_status = factor(
        if_else(Vaccinated == 1, "Vaccinated", "Unvaccinated"),
        levels = c("Unvaccinated", "Vaccinated")
      )
    )
  
  # Show correlations for all participants and vaccinated participants
  subtitle_text <- paste(
    spearman_label(plot_data, "All volunteers"),
    spearman_label(
      filter(plot_data, Vaccination_status == "Vaccinated"),
      "Vaccinated only"
    ),
    sep = "\n"
  )
  
  ggplot(
    plot_data,
    aes(
      x = serum_value,
      y = nasal_value,
      colour = Vaccination_status,
      shape = Vaccination_status
    )
  ) +
    geom_point(size = 1.5, alpha = 0.85) +
    
    coord_cartesian(
      xlim = axis_settings$xlim,
      ylim = axis_settings$ylim
    ) +
    
    scale_colour_manual(
      values = c(
        "Unvaccinated" = "black",
        "Vaccinated" = "#2166AC"
      )
    ) +
    
    scale_shape_manual(
      values = c(
        "Unvaccinated" = 16,
        "Vaccinated" = 15
      )
    ) +
    
    labs(
      title = plot_title,
      subtitle = subtitle_text,
      x = "Serum",
      y = "Nasal",
      colour = "Participant group",
      shape = "Participant group"
    ) +
    
    theme_classic(base_size = 5) +
    
    theme(
      plot.title = element_text(face = "bold", size = 7),
      plot.subtitle = element_text(size = 5),
      axis.text = element_text(size = 5),
      axis.title = element_text(size = 5),
      legend.position = "bottom",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 5)
    )
}
# Use the log-transformed markers for the scatter plots
log_markers_with_vax <- bind_cols(
  df %>% select(Vaccinated),
  log_markers
)

matched_scatter_plots_coloured <- pmap(
  matched_marker_pairs,
  make_matched_spearman_scatter,
  data = log_markers_with_vax
)

names(matched_scatter_plots_coloured) <- matched_marker_pairs$marker_id

matched_scatter_plots_combined <- wrap_plots(matched_scatter_plots_coloured, ncol = 4) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 5)
  )

print(matched_scatter_plots_combined)

# Save the figure
ggsave(
  filename = "Outputs/figures/matched_scatter_plots_combined.svg",
  plot = matched_scatter_plots_combined,
  width = 8,
  height = 5,
  units = "in"
)
