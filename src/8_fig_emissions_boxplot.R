#!/usr/bin/env R Interactive
#' @Author Leonardo Brait, Felipe Alexandre, Mariana Faustino
#' @Date May 2024

library(tidyverse)
library(dplyr)
library(ggpubr)
library(cowplot)
library(extrafont)
library(ggplot2)
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

draw_boxplot <- function(
  df, x_var, y_var, fill_var, x_lab, y_lab, title, breaks, break_labels
) {
  # Function to draw a boxplot with jittered points.
  # Args:
  # *df: dataframe with the data.
  # *x_var: variable to plot in the x axis.
  # *y_var: variable to plot in the y axis.
  # *fill_var: variable to fill the boxes.
  # *x_lab: x axis label.
  # *y_lab: y axis label.
  # *title: plot title.
  # Returns:
  # *ggplot object.

  plot <-
    ggplot(
      df,
      aes(x = !!sym(x_var), y = !!sym(y_var), fill = !!sym(fill_var))
    ) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = rep("gray", 6)) +
    geom_jitter(width = 0.2, size = 1) +
    labs(x = x_lab, y = y_lab) +
    theme_pubr() +
    theme(
      axis.text.x = element_text(
        angle = 80, hjust = 1, family = arial_font, size = unit(20, "pt")
      ),
      axis.text.y = element_text(
        hjust = 1, family = arial_font, size = unit(20, "pt")
      ),
      axis.title = element_text(
        family = arial_font, size = unit(20, "pt"), face = "bold"
      ),
      legend.position = "none"
    ) +
    scale_y_continuous(breaks = breaks, labels = break_labels) +
    ggtitle(title)
  
  return(plot)
}

##################################### Main #####################################
# Directories
figures_dir <- "results/figures_main/"
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

metadata <- read_csv("data/metadata.csv") %>%
    select(sample, tree_richness, id_plot, year) %>%
    filter(!tree_richness == "florest") %>%
    filter(!tree_richness == "pasture") %>%
    mutate(
      tree_richness = case_when(
        tree_richness == "pasture" ~ "Pasture",
        tree_richness == "mono" ~ "1 sp.",
        tree_richness == "tri" ~ "3 spp.",
        tree_richness == "hexa" ~ "6 spp.",
        tree_richness == "div" ~ "12 spp.",
        tree_richness == "alt_div" ~ "24 spp.",
        TRUE ~ tree_richness
      ),
      tree_richness = factor(
        tree_richness, 
        levels = c(
          "Pasture", "1 sp.", "3 spp.", "6 spp.", "12 spp.", "24 spp."
        )
      )
    )

# load gases
gases_df <- read_csv("data/soil_abiotics/greenhouse_gases.csv") %>%
  rename(co2 = "co2(mg/m^2/day)", ch4 = "ch4(mg/m^2/day)") %>%
  inner_join(metadata, by = "sample")

#boxplot co2 ~ tree richness 2022
co2_2022_plot <- draw_boxplot(
  df = gases_df %>% filter(!year == 2022),
  x_var = "tree_richness",
  y_var = "co2",
  fill_var = "tree_richness",
  x_lab = "Tree richness",
  y_lab = "CO2 (mg / m² / day)",
  title = "2022",
  breaks = c(1000, 6500, 13000),
  break_labels = c("1000", "6,500", "13,000")
)

#boxplot ch4 ~ tree richness
ch4_2022_plot <- draw_boxplot(
  df = gases_df %>% filter(!year == 2022),
  x_var = "tree_richness",
  y_var = "ch4",
  fill_var = "tree_richness",
  x_lab = "Tree richness",
  y_lab = "CH4 (mg / m² / day)",
  title = "",
  breaks = c(-4, 2.4, 9),
  break_labels = c("-4", "2.4", "9")
)
panel <- plot_grid(co2_2022_plot, ch4_2022_plot, ncol = 2)

ggsave2(
  paste0(figures_dir, "emissions_boxplot2022"), 
  panel, 
  width = unit(7, "cm"), 
  height = unit(3.5, "cm")
)

#boxplot co2 ~ tree richness 2023
co2_2023_plot <- draw_boxplot(
  df = gases_df %>% filter(!year == 2023),
  x_var = "tree_richness",
  y_var = "co2",
  fill_var = "tree_richness",
  x_lab = "Tree richness",
  y_lab = "CO2 (mg / m² / day)",
  title = "2023",
  breaks = c(1000, 12000, 24000),
  break_labels = c("1000", "12,000", "24,000")
)

#boxplot ch4 ~ tree richness 2023
ch4__2023_plot <- draw_boxplot(
  df = gases_df %>% filter(!year == 2023),
  x_var = "tree_richness",
  y_var = "ch4",
  fill_var = "tree_richness",
  x_lab = "Tree richness",
  y_lab = "CH4 (mg / m² / day)",
  title = "",
  breaks = c(-2.3, 2.8, 8),
  break_labels = c("-2.3", "2.8", "8")
)

panel <- plot_grid(
  co2_2023_plot,
  ch4__2023_plot,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 20,
  label_fontfamily = arial_font
)
ggsave2(
  paste0(figures_dir, "emissions_boxplot2023"), 
  panel, 
  width = unit(7, "cm"), 
  height = unit(3.5, "cm")
)


