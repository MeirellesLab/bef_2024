#!/usr/bin/env R interactive
#' @author Leonardo Brait, Felipe Alexandre, Mariana Faustino
#' @date May 2024
library(tidyverse)
library(vegan)
library(dplyr)
library(ggpubr)
library(extrafont)
library(cowplot)
# Arial font
loadfonts(device = c("all"), quiet = TRUE)
if (Sys.info()[['sysname']] == "Linux") {
  font_import(pattern = "Arial", prompt = FALSE)
  arial_font <- "Arial"
} else {
  font_import(pattern = "arial", prompt = FALSE)
  arial_font <- "arial"
}
source("src/ggsave2.R")

################################# Functions ####################################

draw_linear_regression <- function(
  df,
  response,
  predictor,
  title,
  xlab,
  ylab,
  x_breaks = NULL,
  x_breaks_labels = NULL,
  y_breaks = NULL,
  y_breaks_labels = NULL
) {
  plot_y_range <- max(df[[response]]) - min(df[[response]])
  plot <-
    ggplot(df, aes(x = !!sym(predictor), y = !!sym(response))) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = title, x = xlab, y = ylab) +
    theme_pubr() +
    theme(
      axis.text = element_text(
        size = unit(20, "pt"), family = arial_font, hjust = 1
      ),
      axis.title = element_text(
        size = unit(20, "pt"), family = arial_font, face = "bold"
      ),
    ) +
    scale_x_continuous(breaks = x_breaks, labels = x_breaks_labels) +
    scale_y_continuous(breaks = y_breaks, labels = y_breaks_labels)


  # linear regression
  formula <- paste(response, "~", predictor, sep = " ")
  lm <- lm(formula, data = df) %>%
    summary()
    print(lm)
    
  # get p-value
  pvalue <- lm$coefficients[2,4]
  r_squared <- lm$r.squared

  #get coefficients
  coefficients <- coef(lm)
  intercept <- coefficients[1]
  slope <- coefficients[2]

  # add p-value and r squared to the plot
  plot <- plot +
    annotate(
      "text",
      x = max(df[[predictor]]) - max(df[[predictor]]) / 5,
      y = min(df[[response]]) + plot_y_range * 0.9,
      label = paste("p = ", round(pvalue, 3), sep = ""),
      size = 6, hjust = 0, vjust = 0, family = arial_font
    ) +
    annotate(
      "text", 
      x = max(df[[predictor]]) - max(df[[predictor]]) / 5,
      y = min(df[[response]]) + plot_y_range * 0.85,
      label = paste("R² = ", round(r_squared, 3), sep = ""),
      size = 6, hjust = 0, vjust = 0, family = arial_font
    ) +
    annotate(
      "text", 
      x = max(df[[predictor]]) - max(df[[predictor]]) / 4,
      y = min(df[[response]] + plot_y_range * 0.8),
      label = paste(
        "y = ", round(intercept, 3), "+", round(slope, 3), "x", sep = ""
      ),
      size = 6, hjust = 0, vjust = 0, family = "Arial"
    )

  return(plot)
}

##################################### Main #####################################
# Directories
exploratory_figures_dir <- "results/figures_exploratory/"
if (!dir.exists(exploratory_figures_dir)) {
  dir.create(exploratory_figures_dir, recursive = TRUE)
}
main_figures_dir <- "results/figures_main/"
if (!dir.exists(main_figures_dir)) {
  dir.create(main_figures_dir, recursive = TRUE)
}



# Paths to the csv files
csv_paths_2023 <- c(
  "data/taxa_richness_diversity/richness_shannon_phyla_relative_16s_2023.csv",
  "data/taxa_richness_diversity/richness_shannon_family_relative_16s_2023.csv",
  "data/taxa_richness_diversity/richness_shannon_genus_relative_16s_2023.csv"  
)


#Load metadata
metadata <- read_csv("data/metadata.csv") %>%
    select(sample, tree_richness, id_plot) %>%
    filter(!tree_richness == "forest") %>%
    filter(!tree_richness == "pasture") %>%
    mutate(
      tree_richness = case_when(
        tree_richness == "pasture" ~ "0",
        tree_richness == "mono" ~ "1",
        tree_richness == "tri" ~ "3",
        tree_richness == "hexa" ~ "6",
        tree_richness == "div" ~ "12",
        tree_richness == "alt_div" ~ "24",
        TRUE ~ tree_richness
      ),
      tree_richness = as.numeric(tree_richness)
    ) 

# load gases
gases_df <- read_csv("data/soil_abiotics/greenhouse_gases.csv")


################################################################################
############################### LOAD 2023 DATA #################################

# Merge metadata and gases data into the richness and shannon dataframes
phyla_df <- read_csv(csv_paths_2023[1]) %>%
  inner_join(metadata, by = "sample") %>%
  inner_join(gases_df, by = "sample") %>%
  rename(co2 = "co2(mg/m^2/day)", ch4 = "ch4(mg/m^2/day)") %>%
  mutate(richness = log(richness))
family_df <- read_csv(csv_paths_2023[2]) %>%
  inner_join(metadata, by = "sample") %>%
  inner_join(gases_df, by = "sample") %>%
  rename(co2 = "co2(mg/m^2/day)", ch4 = "ch4(mg/m^2/day)") %>%
  mutate(richness = log(richness))
genus_df <- read_csv(csv_paths_2023[3]) %>%
  inner_join(metadata, by = "sample") %>%
  inner_join(gases_df, by = "sample") %>%
  rename(co2 = "co2(mg/m^2/day)", ch4 = "ch4(mg/m^2/day)") %>%
  mutate(richness = log(richness))

################################# Exploratory ##################################
# Tree richness vs richness ----------------------------------------------------
phyla_richness_versus_trees_plot <- draw_linear_regression(
  df = phyla_df,
  response = "richness",
  predictor = "tree_richness",
  title = "phyla richness vs tree richness 2023",
  xlab = "Tree Richness",
  ylab = "Richness (log)"
)
family_richnees_versus_trees_plot <- draw_linear_regression(
  df = family_df,
  response = "richness",
  predictor = "tree_richness",
  title = "family richness vs tree richness 2023",
  xlab = "Tree Richness",
  ylab = "Richness (log)"
)
genus_richness_versus_trees_plot <- draw_linear_regression(
  df = genus_df,
  response = "richness",
  predictor = "tree_richness",
  title = "genus richness vs tree richness 2023",
  xlab = "Tree Richness",
  ylab = "Richness (log)"
)
panel_richness_versus_trees_2023 <- plot_grid(
  phyla_richness_versus_trees_plot,
  family_richnees_versus_trees_plot,
  genus_richness_versus_trees_plot,
  ncol = 3
)
# png
ggsave(
    filename = paste0(
      exploratory_figures_dir, 
      "panel_richness_versus_trees_2023.png"
    ),
    plot = panel_richness_versus_trees_2023,
  width = unit(21, "cm"),
  height = unit(7, "cm")
)
# svg
ggsave(
    filename = paste0(
      exploratory_figures_dir, 
      "panel_richness_versus_trees_2023.svg"
    ),
    plot = panel_richness_versus_trees_2023,
  width = unit(21, "cm"),
  height = unit(7, "cm")
)

# Co2 vs microbial Richness ----------------------------------------------------
co2_versus_richness_phyla_plot <- draw_linear_regression(
  df = phyla_df,
  response = "co2",
  predictor = "richness",
  title = "co2 vs microbial richness 2023",
  xlab = "Microbial Richness",
  ylab = "Co2 (mg/m²/day)"
)
co2_versus_richness_family_plot <- draw_linear_regression(
  df = family_df,
  response = "co2",
  predictor = "richness",
  title = "co2 vs microbial richness 2023",
  xlab = "Microbial Richness",
  ylab = "Co2 (mg/m²/day)"
)
co2_versus_richness_genus_plot <- draw_linear_regression(
  df = genus_df,
  response = "co2",
  predictor = "richness",
  title = "co2 vs microbial richness 2023",
  xlab = "Microbial Richness",
  ylab = "Co2 (mg/m²/day)"
)
panel_co2_versus_richness_2023 <- plot_grid(
  co2_versus_richness_phyla_plot,
  co2_versus_richness_family_plot,
  co2_versus_richness_genus_plot,
  ncol = 3
)

# png
ggsave(
    filename = paste0(
      exploratory_figures_dir, 
      "panel_co2_versus_richness_2023.png"
    ),
    plot = panel_co2_versus_richness_2023,
  width = unit(21, "cm"),
  height = unit(7, "cm")
)

# svg
ggsave(
    filename = paste0(
      exploratory_figures_dir, 
      "panel_co2_versus_richness_2023.svg"
    ),
    plot = panel_co2_versus_richness_2023,
  width = unit(21, "cm"),
  height = unit(7, "cm")
)

# Ch4 vs microbial Richness ----------------------------------------------------
ch4_versus_richness_phyla_plot <- draw_linear_regression(
  df = phyla_df,
  response = "ch4",
  predictor = "richness",
  title = "ch4 vs microbial richness 2023",
  xlab = "Microbial Richness",
  ylab = "Ch4 (mg/m²/day)"
)
ch4_versus_richness_family_plot <- draw_linear_regression(
  df = family_df,
  response = "ch4",
  predictor = "richness",
  title = "ch4 vs microbial richness 2023",
  xlab = "Microbial Richness",
  ylab = "Ch4 (mg/m²/day)"
)
ch4_versus_richness_genus_plot <- draw_linear_regression(
  df = genus_df,
  response = "ch4",
  predictor = "richness",
  title = "ch4 vs microbial richness 2023",
  xlab = "Microbial Richness",
  ylab = "Ch4 (mg/m²/day)"
)
panel_ch4_versus_richness_2023 <- plot_grid(
  ch4_versus_richness_phyla_plot,
  ch4_versus_richness_family_plot,
  ch4_versus_richness_genus_plot,
  ncol = 3
)

# png
ggsave(
    filename = paste0(
      exploratory_figures_dir, 
      "panel_ch4_versus_richness_2023.png"
    ),
    plot = panel_ch4_versus_richness_2023,
  width = unit(21, "cm"),
  height = unit(7, "cm")
)

# svg
ggsave(
    filename = paste0(
      exploratory_figures_dir, 
      "panel_ch4_versus_richness_2023.svg"
    ),
    plot = panel_ch4_versus_richness_2023,
  width = unit(21, "cm"),
  height = unit(7, "cm")
)

################################## Final Panels 2023 ###########################
# Panel 1 (Richness[phyla, genus] vs Trees)
phyla_richness_versus_trees_plot2 <- draw_linear_regression(
  df = phyla_df,
  response = "richness",
  predictor = "tree_richness",
  title = "",
  xlab = "Tree Richness",
  ylab = "Richness of Microbial Phyla (log)",
  x_breaks = c(0, 1, 3, 6, 12, 24),
  x_breaks_labels = c("0", "1", "3", "6", "12", "24"),
  y_breaks = c(3.05, 3.255, 3.48),
  y_breaks_labels = c("3.05", "3.255", "3.48")
)
genus_richness_versus_trees_plot <- draw_linear_regression(
  df = genus_df,
  response = "richness",
  predictor = "tree_richness",
  title = "",
  xlab = "Tree Richness",
  ylab = "Richness of Microbial Genus (log)",
  x_breaks = c(0, 1, 3, 6, 12, 24),
  x_breaks_labels = c("0", "1", "3", "6", "12", "24"),
  y_breaks = c(5.32, 5.55, 5.85),
  y_breaks_labels = c("5.32", "5.55", "5.85")
)
panel1 <- plot_grid(
  phyla_richness_versus_trees_plot2,
  genus_richness_versus_trees_plot,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 20,
  label_fontfamily = arial_font
)

ggsave2(
    filename = paste0(
      main_figures_dir, 
      "panel_richness_versus_trees_2023"
    ),
    plot = panel1,
  width = unit(21, "cm"),
  height = unit(7, "cm")
)

#panel 2 (Diversity[phyla, genus] vs Trees)
diversity_phyla_versus_trees_plot <- draw_linear_regression(
  df = phyla_df,
  response = "shannon",
  predictor = "tree_richness",
  title = "",
  xlab = "Tree Richness",
  ylab = "Shannon Diversity for Microbial Phyla",
  x_breaks = c(0, 1, 3, 6, 12, 24),
  x_breaks_labels = c("0", "1", "3", "6", "12", "24"),
  y_breaks = c(5.4, 6.5, 7.7),
  y_breaks_labels = c("5.4", "6.5", "7.7")
)
diversity_genus_versus_trees_plot <- draw_linear_regression(
  df = genus_df,
  response = "shannon",
  predictor = "tree_richness",
  title = "",
  xlab = "Tree Richness",
  ylab = "Shannon Diversity for Microbial Genus",
  x_breaks = c(0, 1, 3, 6, 12, 24),
  x_breaks_labels = c("0", "1", "3", "6", "12", "24"),
  y_breaks = c(42, 65, 85),
  y_breaks_labels = c("42", "65", "85")
)

panel2 <- plot_grid(
  diversity_phyla_versus_trees_plot,
  diversity_genus_versus_trees_plot,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 20,
  label_fontfamily = arial_font
)

ggsave2(
  filename = paste0(
      main_figures_dir,
      "panel_diversity_versus_trees_2023"
    ),
    plot = panel2,
  width = unit(21, "cm"),
  height = unit(7, "cm")
)

# Panel 3 (Co2 vs trees)
co2_versus_trees_plot <- draw_linear_regression(
  df = phyla_df,
  response = "co2",
  predictor = "tree_richness",
  title = "",
  xlab = "Tree Richness",
  ylab = "Co2 (mg/m²/day)",
  x_breaks = c(0, 1, 3, 6, 12, 24),
  x_breaks_labels = c("0", "1", "3", "6", "12", "24"),
  y_breaks = c(1000, 6500, 13000),
  y_breaks_labels = c("1000", "6500", "13000")
)

ch4_versus_trees_plot <- draw_linear_regression(
  df = phyla_df,
  response = "ch4",
  predictor = "tree_richness",
  title = "",
  xlab = "Tree Richness",
  ylab = "Ch4 (mg/m²/day)",
  x_breaks = c(0, 1, 3, 6, 12, 24),
  x_breaks_labels = c("0", "1", "3", "6", "12", "24"),
  y_breaks = c(-4.2, 2.5, 9.2),
  y_breaks_labels = c("-4.2", "2.5", "9.2")
)

panel3 <- plot_grid(
  co2_versus_trees_plot,
  ch4_versus_trees_plot,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 20,
  label_fontfamily = arial_font
)

ggsave2(
  filename = paste0(
      main_figures_dir,
      "panel_gases_versus_trees_2023"
    ),
    plot = panel3,
  width = unit(21, "cm"),
  height = unit(7, "cm")
)