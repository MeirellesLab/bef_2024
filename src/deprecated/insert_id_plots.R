
library(tidyverse)

insert_id_plot <- function(df_path, reference_df_path) {
  df <- read_csv(df_path)
  reference_df <- read_csv(reference_df_path)

  df <- df %>%
    mutate(id_plot = reference_df$id_plot[match(sample, reference_df$sample)])

  return(df)
}

args <- commandArgs(trailingOnly = TRUE)
df_path <- args[1]
reference_df_path <- args[2]

list_of_df_paths <- list.files(
    path = df_path, pattern = "2022.csv", full.names = TRUE)
for (df_path in list_of_df_paths) {
  adjusted_df <- insert_id_plot(df_path, reference_df_path)
  write_csv(adjusted_df, df_path)
}