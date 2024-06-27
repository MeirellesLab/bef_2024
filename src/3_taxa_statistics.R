#!/usr/bin/env R Interactive
#' @title Produce taxa descriptive statistics
#' @description This script reads the tables with the relative abundance of
#'   taxa and calculates the mean, median, standard deviation and standard error
#'   for each taxa.
#' @Author Leonardo Brait
#' @date April 2024

library(tidyverse)

produce_taxa_stats <- function(table_name) {
  # This function reads a wide table with the abundance of taxa and
  # calculates the mean, median, standard deviation and standard error for each
  # taxa.
    df <- read_csv(paste0("data/taxonomic_annotation/", table_name))

    # transform in long format
    df_long <- df %>%
      pivot_longer(
        cols = -c(sample),
        names_to = "taxa",
        values_to = "relative_abundance"
      )
    
    # get statistics for taxa (mean, median, sd, se)
    taxa_descriptive_statistics <- df_long %>%
      group_by(taxa) %>%
      summarise(
        mean = mean(relative_abundance),
        median = median(relative_abundance),
        max = max(relative_abundance),
        min = min(relative_abundance),
        sd = sd(relative_abundance),
        se = sd(relative_abundance) / sqrt(n())
      )
    
    return(taxa_descriptive_statistics)
}

tables_names <- c(
    "class_relative_16s_2022",
    "family_relative_16s_2022",
    "genus_relative_16s_2022",
    "order_relative_16s_2022",
    "phyla_relative_16s_2022",
    
    "class_relative_16s_2023",
    "family_relative_16s_2023",
    "genus_relative_16s_2023",
    "order_relative_16s_2023",
    "phyla_relative_16s_2023"
)

if (!dir.exists("results/stats_descriptive_taxa")) {
    dir.create("results/stats_descriptive_taxa")
}

for (table_name in tables_names) {
    taxa_descriptive_statistics <- produce_taxa_stats(
      paste0(table_name, ".csv")
    )
    write_csv(
      taxa_descriptive_statistics, 
      paste0("results/stats_descriptive_taxa/", table_name, ".csv")
    )
}