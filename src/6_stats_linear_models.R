#!/usr/bin/env R interactive
#' @author Leonardo Brait, Felipe Alexandre, Mariana Faustino
#' @date May 2024

library(tidyverse)
library(vegan)
library(dplyr)
library(gvlma)



################################### Functions ##################################

do_linear_model <- function(
  df, formula, table_name, stats_dir = "results/stats_linear_model/") {
  # Function to perform a linear model and save the results.
  # Args:
  # *df: dataframe with the data.
  # *formula: formula for the linear model.
  # *table_name: name of the table.
  # *stats_dir: directory to save the results.

  if (!dir.exists(stats_dir)) {
    dir.create(stats_dir, recursive = TRUE)
  }
  
  linear_model <- lm(
    formula,
    data = df
  )
  
  summary(linear_model)$coefficients %>%
    as.data.frame() %>%
    mutate(test = rownames(.)) %>%
    select(test, everything()) %>%
    write_csv(paste0(stats_dir, "coefficients_", table_name, ".csv"))
  
  gvlma(linear_model) %>%
    summary() %>%
    as.data.frame() %>%
    mutate(test = rownames(.)) %>%
    select(test, everything()) %>%
    write_csv(paste0(stats_dir, "assumptions_", table_name, ".csv"))

}


##################################### Main #####################################




# Paths to the csv files
csv_paths_2022 <- c(
  "data/taxa_richness_diversity/richness_shannon_phyla_relative_16s_2022.csv",
  "data/taxa_richness_diversity/richness_shannon_family_relative_16s_2022.csv",
  "data/taxa_richness_diversity/richness_shannon_genus_relative_16s_2022.csv"  
)
csv_paths_2023 <- c(
  "data/taxa_richness_diversity/richness_shannon_phyla_relative_16s_2023.csv",
  "data/taxa_richness_diversity/richness_shannon_family_relative_16s_2023.csv",
  "data/taxa_richness_diversity/richness_shannon_genus_relative_16s_2023.csv"
)



#Load metadata
metadata <- read_csv("data/metadata.csv") %>%
    select(sample, tree_richness, id_plot) %>%
    filter(!tree_richness == "florest") %>%
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
############################### LOAD 2022 DATA #################################

# Merge metadata and gases data into the richness and shannon dataframes
phyla_df <- read_csv(csv_paths_2022[1]) %>%
  inner_join(metadata, by = "sample") %>%
  inner_join(gases_df, by = "sample") %>%
  rename(co2 = "co2(mg/m^2/day)", ch4 = "ch4(mg/m^2/day)") %>%
  mutate(richness = log(richness))
family_df <- read_csv(csv_paths_2022[2]) %>%
  inner_join(metadata, by = "sample") %>%
  inner_join(gases_df, by = "sample") %>%
  rename(co2 = "co2(mg/m^2/day)", ch4 = "ch4(mg/m^2/day)") %>%
  mutate(richness = log(richness))
genus_df <- read_csv(csv_paths_2022[3]) %>%
  inner_join(metadata, by = "sample") %>%
  inner_join(gases_df, by = "sample") %>%
  rename(co2 = "co2(mg/m^2/day)", ch4 = "ch4(mg/m^2/day)") %>%
  mutate(richness = log(richness))

################################## Statistics ##################################
# Statstics Richness ~ tree_richness -------------------------------------------
do_linear_model(
  phyla_df, 
  as.formula("richness ~ tree_richness"), 
  "2022_richness_phyla_versus_trees"
)
do_linear_model(
  family_df, 
  richness ~ tree_richness, 
  "2022_richness_family_versus_trees"
)
do_linear_model(
  genus_df, 
  richness ~ tree_richness, 
  "2022_richness_genus_versus_trees"
)
# Statstics Shannon ~ tree_richness --------------------------------------------
do_linear_model(
  phyla_df, 
  shannon ~ tree_richness, 
  "2022_shannon_phyla_versus_trees"
)
do_linear_model(
  family_df, 
  shannon ~ tree_richness, 
  "2022_shannon_family_versus_trees"
)
do_linear_model(
  genus_df, 
  shannon ~ tree_richness, 
  "2022_shannon_genus_versus_trees"
)
# Co2 ~ richness ---------------------------------------------------------------
do_linear_model(
  phyla_df, 
  "co2 ~ richness",
  "2022_co2_versus_phyla_richness"
)
do_linear_model(
  family_df, 
  "co2 ~ richness",
  "2022_co2_versus_family_richness"
)
do_linear_model(
  genus_df, 
  "co2 ~ richness",
  "2022_co2_versus_genus_richness"
)
# Co2 ~ shannon ---------------------------------------------------------------
do_linear_model(
  phyla_df, 
  "co2 ~ shannon",
  "2022_co2_versus_phyla_shannon"
)
do_linear_model(
  family_df, 
  "co2 ~ shannon",
  "2022_co2_versus_family_shannon"
)
do_linear_model(
  genus_df, 
  "co2 ~ shannon",
  "2022_co2_versus_genus_shannon"
)
# CH4 ~ richness --------------------------------------------------------------
do_linear_model(
  phyla_df, 
  "ch4 ~ richness",
  "2022_ch4_versus_phyla_richness"
)
do_linear_model(
  family_df, 
  "ch4 ~ richness",
  "2022_ch4_versus_family_richness"
)
do_linear_model(
  genus_df, 
  "ch4 ~ richness",
  "2022_ch4_versus_genus_richness"
)
# CH4 ~ shannon --------------------------------------------------------------
do_linear_model(
  phyla_df, 
  "ch4 ~ shannon",
  "2022_ch4_versus_phyla_shannon"
)
do_linear_model(
  family_df, 
  "ch4 ~ shannon",
  "2022_ch4_versus_family_shannon"
)
do_linear_model(
  genus_df, 
  "ch4 ~ shannon",
  "2022_ch4_versus_genus_shannon"
)
################################################################################
############################### LOAD 2023 DATA #################################

# Merge metadata and gases data into the richness and shannon dataframes
phyla_df <- read_csv(csv_paths_2023[1]) %>%
  inner_join(metadata, by = "sample") %>%
  inner_join(gases_df, by = "sample") %>%
  rename(co2 = "co2(mg/m^2/day)", ch4 = "ch4(mg/m^2/day)")
family_df <- read_csv(csv_paths_2023[2]) %>%
  inner_join(metadata, by = "sample") %>%
  inner_join(gases_df, by = "sample") %>%
  rename(co2 = "co2(mg/m^2/day)", ch4 = "ch4(mg/m^2/day)")
genus_df <- read_csv(csv_paths_2023[3]) %>%
  inner_join(metadata, by = "sample") %>%
  inner_join(gases_df, by = "sample") %>%
  rename(co2 = "co2(mg/m^2/day)", ch4 = "ch4(mg/m^2/day)")

################################## Statistics ##################################
# Statstics Richness ~ tree_richness -------------------------------------------
do_linear_model(
  phyla_df, 
  as.formula("richness ~ tree_richness"), 
  "2023_richness_phyla_versus_trees"
)
do_linear_model(
  family_df, 
  richness ~ tree_richness, 
  "2023_richness_family_versus_trees"
)
do_linear_model(
  genus_df, 
  richness ~ tree_richness, 
  "2023_richness_genus_versus_trees"
)
# Statstics Shannon ~ tree_richness --------------------------------------------
do_linear_model(
  phyla_df, 
  shannon ~ tree_richness, 
  "2023_shannon_phyla_versus_trees"
)
do_linear_model(
  family_df, 
  shannon ~ tree_richness, 
  "2023_shannon_family_versus_trees"
)
do_linear_model(
  genus_df, 
  shannon ~ tree_richness, 
  "2023_shannon_genus_versus_trees"
)
# Co2 ~ richness ---------------------------------------------------------------
do_linear_model(
  phyla_df, 
  "co2 ~ richness",
  "2023_co2_versus_phyla_richness"
)
do_linear_model(
  family_df, 
  "co2 ~ richness",
  "2023_co2_versus_family_richness"
)
do_linear_model(
  genus_df, 
  "co2 ~ richness",
  "2023_co2_versus_genus_richness"
)
# Co2 ~ shannon ---------------------------------------------------------------
do_linear_model(
  phyla_df, 
  "co2 ~ shannon",
  "2023_co2_versus_phyla_shannon"
)
do_linear_model(
  family_df, 
  "co2 ~ shannon",
  "2023_co2_versus_family_shannon"
)
do_linear_model(
  genus_df, 
  "co2 ~ shannon",
  "2023_co2_versus_genus_shannon"
)
# CH4 ~ richness --------------------------------------------------------------
do_linear_model(
  phyla_df, 
  "ch4 ~ richness",
  "2023_ch4_versus_phyla_richness"
)
do_linear_model(
  family_df, 
  "ch4 ~ richness",
  "2023_ch4_versus_family_richness"
)
do_linear_model(
  genus_df, 
  "ch4 ~ richness",
  "2023_ch4_versus_genus_richness"
)
# CH4 ~ shannon --------------------------------------------------------------
do_linear_model(
  phyla_df, 
  "ch4 ~ shannon",
  "2023_ch4_versus_phyla_shannon"
)
do_linear_model(
  family_df, 
  "ch4 ~ shannon",
  "2023_ch4_versus_family_shannon"
)
do_linear_model(
  genus_df, 
  "ch4 ~ shannon",
  "2023_ch4_versus_genus_shannon"
)



