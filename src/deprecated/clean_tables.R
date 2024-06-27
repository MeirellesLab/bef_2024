#' @title Clean tables
#' @description Correct incongruences in the tables. Replace some
#'   samples collected in 2022 displaced in 2023.
#'   in wrong place.
#' @observation we renamed the "Richness" and "Riqueza" columns 
#'   to "tree-richness" before running this script.
#' @author Leonardo Brait
#' @date April 2024

library(tidyverse)

clean_tables <- function(table_2022, table_2023, displaced_samples) {
  
  # load tables, re-type and rename columns
  table_2022_data <- read_csv(paste0(
    "data/taxonomy_annotation/", table_2022
    )) %>%
    mutate(sample = as.character(sample))
  table_2023_data <- read_csv(paste0(
    "data/taxonomy_annotation/", table_2023)) %>%
     mutate(id_plot = as.character(id_plot))

  # take the samples that were displaced,
  # remove zero columns and bind the tables
  table_2022_plus <- table_2023_data %>%
    filter(sample %in% displaced_samples) %>%
    select_if(~ !all(. == 0)) %>%
    bind_rows(table_2022_data) %>%
    replace_na(list(id_plot = "missing_data")) %>%
    mutate_all(~ replace_na(., 0)) %>%
    distinct(sample, .keep_all = TRUE)
  
  table_2023_clean <- table_2023_data %>%
    filter(!sample %in% displaced_samples)
  
  # write the new tables
  write_csv(table_2022_plus, paste0("data/taxonomy_annotation/", table_2022))
  write_csv(table_2023_clean, paste0("data/taxonomy_annotation/", table_2023))
}

displaced_samples <- c("I07", "I22", "P66")
correspondent_tables <- c(
   "class_relative_16s_2022.csv" = "class_relative_16s_2023.csv",
   "family_relative_16s_2022.csv" = "family_relative_16s_2023.csv",
   "genus_relative_16s_2022.csv" = "genus_relative_16s_2023.csv",
   "order_relative_16s_2022.csv" = "order_relative_16s_2023.csv",
   "phyla_relative_16s_2022.csv" = "phyla_relative_16s_2023.csv"
)

for (table_2022 in names(correspondent_tables)) {
  clean_tables(
    table_2022,
    correspondent_tables[[table_2022]],
    displaced_samples
  )
}
