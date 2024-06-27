library(tidyverse)
df_paths <- list.files(path = "data/taxonomic_annotation/", full.names = TRUE)

for (df_path in df_paths) {
  df <- read_csv(df_path)
  df <- df %>%
    select(-c("tree_richness","id_plot"))
  write_csv(df, df_path)
}