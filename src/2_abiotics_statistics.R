#!/usr/bin/ R Interactive
#' @title Produce abiotics descriptive statistics
#' @description This function reads a table with abiotics data and produces 
#' descriptive statistics for each abiotic variable.
#' @param table_paths Path to the directory with abiotics data.
#' @param results_dir Path to the directory to save the results.
#' @param metadata_path Path to the metadata file.
#' @example Rscript src/2_abiotics_statistics.R data/soil_abiotics/ results/stats_descriptive_abiotics/ data/metadata.csv
#' @author Leonardo Brait
#' @date April 2024

library(tidyverse)

merge_dataframes <- function(df, metadata) {
    df <- df
    metadata <- metadata %>% select(sample, year)
    merged_df <- left_join(df, metadata, by = "sample")
    return(merged_df)
}

produce_abiotics_stats <- function(df) {
    df <- df

    # transform in long format
    df_long <- df %>%
      pivot_longer(
        cols = -c(sample),
        names_to = "abiotic",
        values_to = "measure"
      )
    
    abiotics_statistics <- df_long %>%
      group_by(abiotic) %>%
      summarise(
        mean = mean(measure),
        median = median(measure),
        max = max(measure),
        min = min(measure),
        sd = sd(measure),
        se = sd(measure) / sqrt(n())
      )
    
    return(abiotics_statistics)
}

##################################### main #####################################
args <- commandArgs(trailingOnly = TRUE)
input_dir <- "data/soil_abiotics/"
results_dir <- "results/stats_descriptive_abiotics/"
metadata_path <-  "data/metadata.csv"

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

table_paths <- list.files(input_dir, full.names = TRUE)
for (table_path in table_paths) {
  merged_df <- merge_dataframes(read_csv(table_path), read_csv(metadata_path))
  
  # Check if there are more than one year in the dataset
  # If so, split the dataset by year and save the statistics for each year
  years <- unique(merged_df$year)
  if (length(years) >= 2) {
    # df list gets a list of dataframes, one for each year
    df_list <- split(merged_df, merged_df$year)
    for (year in names(df_list)) {
      df_list[[year]] <- df_list[[year]] %>%
        select(-year)
      abiotics_statistics <- produce_abiotics_stats(df_list[[year]])
      write_csv(
        abiotics_statistics,
        paste0(results_dir, year, "_", basename(table_path))
      )
    }
  } else { # If there is only one year, save the statistics for the whole
    merged_df <- merged_df %>%
      select(-year)
    abiotics_statistics <- produce_abiotics_stats(merged_df)
    write_csv(
      abiotics_statistics,
      paste0(results_dir, basename(table_path))
    )
  }
}
