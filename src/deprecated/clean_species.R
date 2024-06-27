library(tidyverse)
df <- read_csv("data/taxonomy_annotation/species_relative_16s_2022.csv")
df <- df %>%
  mutate("id_plot" = "missing_data") %>%
  rename("tree_richness" = "tree-richness") %>%
    write_csv("data/taxonomy_annotation/species_relative_16s_2022.csv")
