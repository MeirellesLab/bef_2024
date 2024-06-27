#!/usr/bin/env R Interactive
#' @Author Leonardo Brait, Felipe Alexandre, Mariana Faustino
#' @Date May 2024

library(tidyverse)
library(vegan)
library(ggplot2)
library(gridExtra)
library(extrafont)
library(ggpubr)
# Arial font
loadfonts(device = c("all"), quiet = TRUE)
if (Sys.info()[['sysname']] == "Linux") {
  font_import(pattern = "Arial", prompt = FALSE)
  arial_font <- "Arial"
} else {
  font_import(pattern = "arial", prompt = FALSE)
  arial_font <- "arial"
}

################################## Functions ###################################
do_pca <- function(merged_df_numeric) {
  pca <- rda(merged_df_numeric, scale = TRUE)
  pca_sum <- summary(pca)
  exp <- pca_sum$cont
  exp <- data.frame(exp$importance)

  ord_df <- data.frame(pca_sum$sites, row.names = NULL)
  pca_arrows <- pca_sum$species
  pca_arrows <- data.frame(pca_arrows)
  pca_arrows <- data.frame(
    pc1 = pca_arrows$PC1,
    pc2 = pca_arrows$PC2,
    Species = rownames(pca_arrows)
  )
  pca_arrows <- pca_arrows[1:7, ]

  return(list(ord_df, pca_arrows, exp))

}

draw_pca <- function(pca_result, category) {

  # Get pca results
  ord_df <- pca_result[[1]]
  pca_arrows <- pca_result[[2]]
  exp <- pca_result[[3]]
  ord_df$category <- factor(category)

  # draw the plot
  plot <- ggplot(ord_df) +
    theme_bw() +

    geom_segment(
      data = pca_arrows,
      aes(x = 0, xend = pc1, y = 0, yend = pc2),
      arrow = arrow(length = unit(0.5, "cm")), size = 1, colour = "grey"
    ) +
    geom_point(
      mapping = aes(x = PC1, y = PC2, color = category), 
      size = 8,
      alpha = .8) +
    geom_text(
      data = pca_arrows,
      aes(x = pc1 * 1, y = pc2 * 1, label = Species),
      size = unit(8, "pt"),
      family = "Arial"
    ) +

    xlim(c(-1.75, NA)) +
    geom_hline(yintercept = 0, colour = "gray70") +
    geom_vline(xintercept = 0, colour = "gray70") +

    theme(
      panel.border = element_rect(colour = "white", size = 1),
      legend.key = element_rect(fill = "white"),
      axis.title = element_text(
        face = "bold", size = unit(30, "pt"), family = arial_font
      ),
      axis.text = element_text(
        size = unit(20, "pt"), family = arial_font
      ),
      legend.text = element_text(
        size = unit(20, "pt"), family = arial_font
      ),
      legend.title = element_text(
        size = unit(20, "pt"), family = arial_font, face = "bold"
      ),
      legend.position = "right",
      axis.line.y = element_blank(),
      axis.line.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()

    ) +
    guides(color = guide_legend("Tree Richness")) +
    labs(
      title = NULL,
      x = paste("PC1 (", round(exp[2, 1] * 100, digits = 2), "%)"), 
      y = paste("PC2 (", round(exp[2, 2] * 100, digits = 2), "%)")
    )


  return(plot)
}

##################################### Main #####################################

# General tables -----------------------
metadata <- read_csv("data/metadata.csv") %>%
    select(sample, tree_richness, id_plot, year) %>%
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
      tree_richness = as.factor(tree_richness)
    ) %>%
    filter(tree_richness != "Pasture", tree_richness != "Forest")
gases_df <- read_csv("data/soil_abiotics/greenhouse_gases.csv") %>%
  rename(
    "CO2" = "co2(mg/m^2/day)", 
    "CH4" = "ch4(mg/m^2/day)",
    "Moisture" = "moisture(g)",
    "Conductivity" = "conductivity(ms/cm)",
    "Organic Matter" = "organic_matter(g)",
    "pH" = "ph",
  )


# Produce plots ------------------------
merged_df <- gases_df %>%
  inner_join(metadata, by = "sample") %>%
  filter(year == 2022) %>%
  select(-year)
category <- merged_df$tree_richness

# Co2 plot
merged_df_numeric_co2 <- merged_df %>%
  select(-sample, -tree_richness, -id_plot, -CH4)
pca_co2 <- do_pca(merged_df_numeric_co2)
plot_co2 <- draw_pca(pca_co2, category)


# Ch4 plot
merged_df_numeric_ch4 <- merged_df %>%
  select(-sample, -tree_richness, -id_plot, -CO2)
pca_ch4 <- do_pca(merged_df_numeric_ch4)
plot_ch4 <- draw_pca(pca_ch4, category)

# panel
panel <- ggarrange(
  plot_co2,
  plot_ch4,
  nrow = 2,
  common.legend = TRUE,
  align = "hv",
  labels = "AUTO",
  font.label = list(size = 20, family = arial_font)
)
ggsave2(
  "results/figures_main/panel_pca_gases_2022",
  panel,
  width = unit(19, "cm"),
  height = unit(10, "cm"),
  dpi = 300
)

# Ch4 + Co2 plot
merged_df_numeric <- merged_df %>%
  select(-sample, -tree_richness, -id_plot)
pca <- do_pca(merged_df_numeric)
plot <- draw_pca(pca, category)
ggsave2(
  "results/figures_main/pca_gases_2022",
  plot,
  width = unit(19, "cm"),
  height = unit(10, "cm"),
  dpi = 300
)
