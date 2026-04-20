library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(patchwork)

# File path
file_path <- "data/Combined_nasal_serum_mastersheet_final.csv"


df <- read.csv(file_path)

# Filter and log-transform immune variables
markers <- df %>%
  select(
    matches("BL", ignore.case = FALSE),
    -matches("Alpha|Omicron|Delta|Beta|Gamma|OC43|229E|HKU|CoV_1|NL63")
  ) %>%
  dplyr::mutate_if(is.character, as.numeric)

#log transform
log_markers <- log10(markers)

head(log_markers)

# Spearman correlation matrix and heatmap
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
    !grepl("_N_", colnames(R)) ~ "Anti-spike IgG",
  
  grepl("IgA", colnames(R)) &
    grepl("RBD|CoV_2_S|NTD", colnames(R)) &
    !grepl("_N_", colnames(R)) ~ "Anti-spike IgA",
  
  grepl("IgG", colnames(R)) &
    grepl("_N_", colnames(R)) ~ "Anti-nucleocapsid IgG",
  
  grepl("IgA", colnames(R)) &
    grepl("_N_", colnames(R)) ~ "Anti-nucleocapsid IgA",
  
  grepl("IgM", colnames(R)) ~ "IgM antibody",
  
  grepl("ELISpot|^T_", colnames(R)) ~ "T cell",
  
  TRUE ~ "Other"
)

names(marker_type) <- colnames(R)

# Colours
type_cols <- c(
  "Neutralisation" = "#1f78b4",
  "Anti-spike IgG" = "#1b9e77",
  "Anti-spike IgA" = "#66c2a5",
  "Anti-nucleocapsid IgG" = "#e31a1c",
  "Anti-nucleocapsid IgA" = "#fb9a99",
  "IgM antibody" = "#ff7f00",
  "T cell" = "#6a3d9a",
  "Other" = "grey70"
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

# Left annotation
left_ha <- rowAnnotation(
  Type = marker_type,
  col = list(Type = type_cols),
  show_annotation_name = FALSE,
  simple_anno_size = unit(4, "mm"),
  annotation_legend_param = list(
    Type = list(
      title = "Marker Type",
      title_gp = gpar(fontsize = 18),
      labels_gp = gpar(fontsize = 10)
    )
  )
)

# Heatmap object
ht <- Heatmap(
  R,
  name = "Spearman rho",
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
  
  column_title = "Spearman correlation of baseline immune markers",
  column_title_gp = gpar(fontsize = 30, fontface = "bold"),
  
  heatmap_legend_param = list(
    title = "Spearman rho",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1"),
    direction = "horizontal",
    title_gp = gpar(fontsize = 18),
    labels_gp = gpar(fontsize = 18)
  )
)

# Explicit draw for Jupyter
grid::grid.newpage()
ComplexHeatmap::draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE,
  newpage = FALSE
)


#Nasal/serum correlation scatter plots

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
  select(marker_id, Nasal, Serum)

matched_marker_pairs

format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))
}

spearman_label <- function(data, label) {
  cor_result <- suppressWarnings(
    cor.test(data$nasal_value, data$serum_value, method = "spearman")
  )

  paste0(
    label, ": rho = ", round(unname(cor_result$estimate), 2),
    ", p = ", format_p(cor_result$p.value),
    ", n = ", nrow(data)
  )
}

make_matched_spearman_scatter <- function(marker_id, Nasal, Serum, data) {
  
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

  subtitle_text <- paste(
    spearman_label(plot_data, "All volunteers"),
    spearman_label(filter(plot_data, Vaccination_status == "Vaccinated"), "Vaccinated only"),
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
    geom_point(size = 2.4, alpha = 0.85) +
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
      title = marker_id,
      subtitle = subtitle_text,
      x = paste("Serum", marker_id),
      y = paste("Nasal", marker_id),
      colour = "Participant group",
      shape = "Participant group"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
}

log_markers_with_vax <- df %>%
  select(
    Vaccinated,
    matches("BL"),
    -matches("Alpha|Omicron|Delta|Beta|Gamma|OC43|229E|HKU|CoV_1|NL63")
  ) %>%
  mutate(across(-Vaccinated, ~ log10(as.numeric(.x))))

matched_scatter_plots_coloured <- pmap(
  matched_marker_pairs,
  ~ make_matched_spearman_scatter(
    marker_id = ..1,
    Nasal = ..2,
    Serum = ..3,
    data = log_markers_with_vax
  )
)

names(matched_scatter_plots_coloured) <- matched_marker_pairs$marker_id

wrap_plots(matched_scatter_plots_coloured, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

