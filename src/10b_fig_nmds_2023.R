#!/usr/bin/env R Interactive
#' @Author Leonardo Brait, Felipe Alexandre, Mariana Faustino
#' @date May 2024
library(vegan)
library(tidyverse)
library(extrafont)
library(ggpubr)
library(cowplot)
library(ggplot2)
source("src/ggsave2.R")

# Arial font
loadfonts(device = c("all"), quiet = TRUE)
if (Sys.info()[['sysname']] == "Linux") {
  font_import(pattern = "Arial", prompt = FALSE)
  arial_font <- "Arial"
} else {
  font_import(pattern = "arial", prompt = FALSE)
  arial_font <- "arial"
}


#################################### Functions #################################
draw_nmds <- function(nmds_df, stress, title, anosim_data) {
  plot <- ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = Category)) +
    geom_point(size = 4) +
    theme_bw() +
    labs(title = title, color = "Tree Richness")+
    theme(
      plot.title = element_text(
        face = "bold", size = unit(30, "pt"), family = "Arial"
      ),
      plot.tag = element_text(face = "bold", size = 30),
      legend.position = "right",
      legend.title = element_text(
        size = unit(20, "pt"), face = "bold", family = "Arial"
      ),
      legend.key = element_rect(fill = "white"),
      legend.text = element_text(size = unit(20, "pt"), family = "Arial"),
      axis.title = element_text(
        size = unit(20, "pt"), face = "bold", family = "Arial"
      ),
      axis.text = element_text(size = unit(20, "pt"), family = "Arial"),
      panel.border = element_rect(colour = "black", size = 1),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_colour_manual(values = c(
      "#43CD80",
      "#9400D3",
      "#c2c212",
      "#3eb5bf",
      "#7e3c03",
      "#0f0476",
      "#f72d02"
    )) + 
    annotate(
      "text", x = 0.025, y = 0.1, label = paste("Stress: ", stress), 
      hjust = 0, vjust = 1, size = 5
    ) +
    annotate(
      "text", x = 0.025, y = 0.09, 
      label = paste("ANOSIM R = ", round(anosim_data$statistic, 2)),
      hjust = 0, vjust = 1, size = 5
    ) +
    annotate(
      "text", x = 0.025, y = 0.08, 
      label = paste("p-value = ", round(anosim_data$signif, 2)),
      hjust = 0, vjust = 1, size = 5
    )

  return(plot)
}

#################################### Main ######################################
# General metadata ---------------------
metadata_df <- read_csv("data/metadata.csv") %>%
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

# Phyla ------------------------------------------------------------------------
phyla_df <- read_csv("data/taxonomic_annotation/phyla_relative_16s_2023.csv")
merged_df <- phyla_df %>%
  inner_join(metadata_df, by = "sample")
data_num <- merged_df %>%
  select(-sample, -tree_richness, -id_plot)


nmds_data <- metaMDS(data_num, try = 4999, distance = "euclidean")
nmds_df <- 
  cbind(
    sample = merged_df$sample,
    tree_richness = merged_df$tree_richness,
    data.frame(nmds_data$points)
  ) %>%
  rename( Category = tree_richness, Samples = sample)
stress <- round(nmds_data$stress, 2)
anosim <- anosim(
  data_num, 
  merged_df$tree_richness,
  permutations = 999,
  distance = "euclidean"
)

phyla_plot <- draw_nmds(nmds_df, stress, "Phyla", anosim)

# Family -----------------------------------------------------------------------
family_df <- read_csv("data/taxonomic_annotation/family_relative_16s_2023.csv")
merged_df <- family_df %>%
  inner_join(metadata_df, by = "sample")
data_num <- merged_df %>%
  select(-sample, -tree_richness, -id_plot)

nmds_data <- metaMDS(data_num, try = 4999, distance = "euclidean")
nmds_df <- 
  cbind(
    sample = merged_df$sample,
    tree_richness = merged_df$tree_richness,
    data.frame(nmds_data$points)
  ) %>%
  rename( Category = tree_richness, Samples = sample)
stress <- round(nmds_data$stress, 2)
anosim <- anosim(
  data_num, 
  merged_df$tree_richness,
  permutations = 4999,
  distance = "euclidean"
)

family_plot <- draw_nmds(nmds_df, stress, "Family", anosim)

# Genus ------------------------------------------------------------------------
genus_df <- read_csv("data/taxonomic_annotation/genus_relative_16s_2023.csv")
merged_df <- genus_df %>%
  inner_join(metadata_df, by = "sample")
data_num <- merged_df %>%
  select(-sample, -tree_richness, -id_plot)

nmds_data <- metaMDS(data_num, try = 4999, distance = "euclidean")
nmds_df <- 
  cbind(
    sample = merged_df$sample,
    tree_richness = merged_df$tree_richness,
    data.frame(nmds_data$points)
  ) %>%
  rename( Category = tree_richness, Samples = sample)
stress <- round(nmds_data$stress, 2)
anosim <- anosim(
  data_num, 
  merged_df$tree_richness,
  permutations = 4999,
  distance = "euclidean"
)

genus_plot <- draw_nmds(nmds_df, stress, "Genus", anosim)


# panel
panel <- ggarrange(
  phyla_plot,
  genus_plot,
  nrow = 1,
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom",
  labels = "AUTO",
  font.label = list(size = 20, family = arial_font)
)
ggsave2("results/figures_main/nmds_panel_2023", panel, width = 15, height = 7.5)
