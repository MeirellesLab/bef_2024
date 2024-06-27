library(tidyverse)

df_2022 <- read_csv("data/gases_22.csv") %>%
 select(-year)
df_2023 <- read_csv("data/gases_23.csv") %>%
 rename(
    "ph" = "ph (6.0)",
    "co2" = "flux_co2",
    "ch4" = "flux_ch4",
    "conductivity" = "condutividade (mS/cm)",
    "moisture" = "umidade",
    "organic_matter" = "materia_organica",
    "inorganic_matter" = "materia_inorganica"
    )

total_df <- bind_rows(df_2022, df_2023)

total_df <- total_df %>%
  rename(
    "co2(mg/m^2/day)" = "co2",
    "ch4(mg/m^2/day)" = "ch4",
    "conductivity(ms/cm)" = "conductivity",
    "moisture(g)" = "moisture",
    "organic_matter(g)" = "organic_matter",
    "inorganic_matter(g)" = "inorganic_matter"

  )

write_csv(total_df, "data/gases.csv")
