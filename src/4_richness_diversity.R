#!/bin/R Interactive
#' @title Diversity Means
#' @description This script calculates the richness and diversity 
#' of taxonomic groups.
#' @Params:
#' *input_path: path to the directory with the tables.
#' *results_path: path to the directory to save the results.
#' @Author Leonardo Brait, Felipe Alexandre, Mariana Faustino
#' @Revisor Leonardo Brait
#' @date April 2024

library(dplyr)
library(tidyverse)
library(vegan)


calculate_richness_shannon <- function(taxa_abundance_path) {
  #' @Title Calculate Richness
  #' @Description: Function to calculate the richness and 
  #' EXPONENTIAL of Shanon diversity of taxonomic groups in a given dataframe.
  #' @Args:
  #' *taxa_abundance_path: path to the dataframe with the abundance of taxa.
  #' @Returns:
  #'  *richness_diversity_df: dataframe with the richness of taxonomic groups 
  #' and corresponding factors.
 
  print(paste0("Calculating richness for ", taxa_abundance_path))
  taxa_abundance_df <- read_csv(taxa_abundance_path) %>%
    mutate(sample = as.factor(sample))
  
  taxa_abundance_num <- taxa_abundance_df %>%
    select_if(is.numeric)
  
  richness_diversity_df <- taxa_abundance_df %>%
    mutate(
      richness = rowSums(taxa_abundance_num > 0),
      shannon = exp(diversity(taxa_abundance_num)) 
    ) %>%
    select(sample, richness, shannon)
  return(richness_diversity_df)
}

#################################### main ######################################
input_path <- "data/taxonomic_annotation/"
results_path <- "data/taxa_richness_diversity/"
if (!dir.exists(results_path)) {
  dir.create(results_path, recursive = TRUE)
}

table_paths <- list.files(input_path, full.names = TRUE)
for (i in 1:length(table_paths)) {
  result_path <- paste0(
    results_path,
    "richness_shannon_",
    basename(table_paths[i])
  )

  taxa_rich_shanon_by_treerichness <- calculate_richness_shannon(
    table_paths[i]
  )
  write_csv(taxa_rich_shanon_by_treerichness, result_path)
}
