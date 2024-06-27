#Running permutational ANOVA on the data
#Load libraries
library(vegan)
library(tidyverse)
library(EcolUtils)
library(ggpubr)
library(ggplot2)

#Load data
taxonomy_22_num <- read.csv("data/taxonomic_annotation/genus_relative_16s_2022.csv")
taxonomy_23_num <- read.csv("data/taxonomic_annotation/genus_relative_16s_2023.csv")
metadata <-  read.csv("data/metadata.csv")


#Filter taxonomy data to only include samples that are in the metadata
taxonomy_22_num <- taxonomy_22_num %>% filter(sample %in% metadata$sample)
taxonomy_23_num <- taxonomy_23_num %>% filter(sample %in% metadata$sample)

#Check which groups are common to the two years
common_groups <- intersect(colnames(taxonomy_22), colnames(taxonomy_23))

#Merge data
taxonomy_22 <- merge(metadata, taxonomy_22_num, by = "sample")

taxonomy_23 <- merge(metadata, taxonomy_23_num, by = "sample")

#Put column in the rownames
taxonomy_22_num <- taxonomy_22_num %>% column_to_rownames(var = "sample")
taxonomy_23_num <- taxonomy_23_num %>% column_to_rownames(var = "sample")

#Calculate diversity and richness for each year
richness_diversity_df_22 <- taxonomy_22_num %>%
    mutate(
    richness = rowSums(taxonomy_22_num > 0),
    shannon = exp(diversity(taxonomy_22_num))
    ) %>%
    select(richness, shannon)

richness_diversity_df_23 <- taxonomy_23_num %>%
    mutate(
    richness = rowSums(taxonomy_23_num > 0),
    shannon = exp(diversity(taxonomy_23_num))
    ) %>%
    select(richness, shannon)

#Create the factor variable
tree_richness_22 <- metadata %>%
    filter(year == 2022) %>%
        select(sample, tree_richness) %>%
            column_to_rownames(var = "sample") %>%
                as.factor()

#Add to the respective diversity and richness dataframes
richness_diversity_df_22 <- cbind(richness_diversity_df_22, tree_richness = metadata$tree_richness[metadata$year == 2022])

tree_richness_23 <- metadata %>%
    filter(year == 2023) %>%
        select(tree_richness) %>%
                as.factor()

richness_diversity_df_23 <- cbind(richness_diversity_df_23, tree_richness = metadata$tree_richness[metadata$year == 2023])


#Combining the two years creating a new column for the year
richness_diversity_df <- rbind(richness_diversity_df_22, richness_diversity_df_23) %>%
    mutate(year = c(rep(2022, nrow(richness_diversity_df_22)), rep(2023, nrow(richness_diversity_df_23))))

#Running ANOVA to check if richness and diversity are significantly different between the two years
anova_rich <- adonis2(richness_diversity_df$richness ~ year, data = richness_diversity_df, permutations = 999)
anova_rich

anova_div <- adonis2(richness_diversity_df$shannon ~ year, data = richness_diversity_df, permutations = 999)
anova_div

#Create a boxplot to visualize the differences
richness_diversity_df %>%
    ggplot(aes(x = as.factor(year), y = richness)) +
    geom_boxplot(fill = "gray") +
    theme_pubr() +
    labs(title = "Richness by year", x = "Year", y = "Richness")

#Save the plot
ggsave("results/richness_boxplot.png")
