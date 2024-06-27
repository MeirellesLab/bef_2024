#!/usr/bin/env R Interactive
#' @Author Leonardo Brait, Felipe Alexandre, Mariana Faustino
#' @date May 2024

library(ggplot2)
library(tidyverse)
library(reshape2)
library(cowplot)
library(ggpubr)
library(extrafont)
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
rename_taxa <- function(df){
  colnames(df) <- gsub("\\.", " ", colnames(df))
  colnames(df) <- gsub("Unclassified Family", "Unc.", colnames(df))
  colnames(df) <- gsub("Unclassified Genus", "Unc.", colnames(df))
  return(df)
}

calculate_others <- function(df, threshold) {
  df_numeric <- df %>%
    select(-sample)

  abundants <- df_numeric %>%
    gather() %>%
    group_by(key) %>%
    summarise(taxa_mean = mean(value)) %>%
    filter(taxa_mean > threshold) %>%
    pull(key)

  others <- df %>% select(-sample, -abundants) %>% rowSums()
  df <- df %>% select(sample, abundants)
  df <- cbind(df, Others = others)

  return(df)
}


order_abundance <- function(df_melt) {
  taxa_order <- df_melt %>%
    group_by(variable) %>%
    summarise(mean_value = mean(value)) %>%
    arrange(desc(mean_value)) %>%
    pull(variable)
  df_melt$variable <- factor(df_melt$variable, levels = taxa_order)

  return(df_melt)
}


produce_abundance_df <- function(df, metadata) {
  df <- df %>%
    inner_join(metadata, by = "sample")
  df_melt <- melt(df, id.vars = c("sample", "id_plot", "tree_richness")) %>%
    order_abundance() %>%
    group_by(id_plot, variable, tree_richness) %>%
    summarise(value = mean(value)) %>%
    ungroup() %>%
    mutate(value = as.numeric(value)) %>%
    mutate(value = value * 100) %>%
    mutate(id_plot = factor(id_plot))
  return(df_melt)
}


draw_stackedbar <- function(df) {
  plot <-
    ggplot(df, aes(x = id_plot, y = value, fill = variable)) +
    geom_bar(stat = "identity",   position = "fill", width = 1) +
    theme_pubr() +
    facet_grid(cols = vars(tree_richness), space = "free", scales = "free") +

    theme(
      strip.background = element_blank(),
      strip.text = element_text(
        size = unit(18, "pt"), face = "bold", family = arial_font
      ),
      axis.text = element_text(
        size = unit(18, "pt"), family = arial_font, hjust = 1
      ),
      axis.title = element_text(
        size = unit(18, "pt"), family = arial_font, face = "bold"
      ),
      legend.text = element_text(
        size = unit(17.7, "pt"), family = arial_font
      ),
      legend.title = element_text(
        size = unit(18, "pt"), family = arial_font, face = "bold"
      ),
      panel.spacing = unit(0.5, "lines")
    ) +
    labs(
      x = "Plot ID",
      y = "Relative abundance ",
      fill = "Taxa"
    ) +
    scale_fill_manual(values = c(stacked_color))
}

#################################### Main ######################################
# Color palette
stacked_color <- c(
  "#B8CCC6",
  "#75E6DA",
  "#D4F1F4",
  "#27678B",
  "#81B622",
  "#0C2D48",
  "#2E8BC0",
  "#B1D4E0",
  "#D5CE58",
  "#ECF87F",
  "#AB351A",
  "#FF6666",
  "#FFFFBF",
  "#ABD9E9",
  "#104210",
  "#215B21",
  "#327E32",
  "#449F44",
  "#55BF55"
)


metadata_test <- read_csv("data/metadata.csv")
# Paths to the csv files
csv_paths_2022 <- c(
  "data/taxonomic_annotation/phyla_relative_16s_2022.csv",
  "data/taxonomic_annotation/family_relative_16s_2022.csv",
  "data/taxonomic_annotation/genus_relative_16s_2022.csv"
)
csv_paths_2023 <- c(
  "data/taxonomic_annotation/phyla_relative_16s_2023.csv",
  "data/taxonomic_annotation/family_relative_16s_2023.csv",
  "data/taxonomic_annotation/genus_relative_16s_2023.csv"
)

#Load metadata
metadata <- read_csv("data/metadata.csv") %>%
    select(sample, tree_richness, id_plot) %>%
    mutate(
      tree_richness = case_when(
        tree_richness == "pasture" ~ "Pasture",
        tree_richness == "mono" ~ "1 sp.",
        tree_richness == "tri" ~ "3 spp.",
        tree_richness == "hexa" ~ "6 spp.",
        tree_richness == "div" ~ "12 spp.",
        tree_richness == "alt_div" ~ "24 spp.",
        tree_richness == "forest" ~ "Forest",
        TRUE ~ tree_richness
      ),
      tree_richness = factor(
        tree_richness, 
        levels = c(
          "Pasture", "1 sp.", "3 spp.", "6 spp.", "12 spp.", "24 spp.", "Forest"
        )
      )
    ) 

# 2022 -------------------------------------------------------------------------
# phyla 22 
df_melt <- read_csv(csv_paths_2022[1]) %>%
  rename_taxa() %>%
  calculate_others(0.02) %>%
  produce_abundance_df(metadata %>% filter(!tree_richness == "Forest"))

phyla22_plot <- draw_stackedbar(df_melt) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  labs(fill = "Phyla") +
  ggtitle("2022")

# family 22
df_melt <- read_csv(csv_paths_2022[2]) %>%
  rename_taxa() %>%
  calculate_others(0.02) %>%
  produce_abundance_df(metadata %>% filter(!tree_richness == "Forest"))

family22_plot <- draw_stackedbar(df_melt) +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank(), 
    axis.title.x = element_blank()
  ) +
  labs(fill = "Family")

# genus 22
df_melt <- read_csv(csv_paths_2022[3]) %>%
  rename_taxa() %>%
  calculate_others(0.02) %>%
  produce_abundance_df(metadata %>% filter(!tree_richness == "Forest"))

genus22_plot <- draw_stackedbar(df_melt) +
  theme(strip.text = element_blank()) +
  labs(fill = "Genus")

# Panel
panel_2022 <- plot_grid(
  phyla22_plot, 
  family22_plot, 
  genus22_plot, 
  ncol = 1, 
  align = "v",
  labels = c("A", "B", "C"),
  label_size = 20,
  label_fontfamily = arial_font
)

ggsave2(
  "results/figures_main/stacked_bar_2022",
  panel_2022,
  width = unit(17, "cm"),
  height = unit(15, "cm"),
)

# 2023 -------------------------------------------------------------------------

# phyla 23
df_melt <- read_csv(csv_paths_2023[1]) %>%
  rename_taxa() %>%
  calculate_others(0.02) %>%
  produce_abundance_df(metadata)

phyla23_plot <- draw_stackedbar(df_melt) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  labs(fill = "Phyla") +
  ggtitle("2023")

# family 23
df_melt <- read_csv(csv_paths_2023[2]) %>%
  rename_taxa() %>%
  calculate_others(0.02) %>%
  produce_abundance_df(metadata )

family23_plot <- draw_stackedbar(df_melt) +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank(), 
    axis.title.x = element_blank()
  ) +
  labs(fill = "Family")

# genus 23
df_melt <- read_csv(csv_paths_2023[3]) %>%
  rename_taxa() %>%
  calculate_others(0.02) %>%
  produce_abundance_df(metadata)

genus23_plot <- draw_stackedbar(df_melt) +
  theme(strip.text = element_blank()) +
  labs(fill = "Genus")

# Panel
panel_2023 <- plot_grid(
  phyla23_plot,
  family23_plot,
  genus23_plot,
  ncol = 1,
  align = "v",
  labels = c("A", "B", "C"),
  label_size = 20,
  label_fontfamily = arial_font
) +
  ggtitle("2023")

ggsave2(
  "results/figures_main/stacked_bar_2023",
  panel_2023,
  width = unit(17, "cm"),
  height = unit(15, "cm"),
)

