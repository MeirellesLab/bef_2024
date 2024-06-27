#Remove dummy data from tables and save

library(tidyverse)

#Load data
phyla <- read.csv("data/taxonomic_annotation/phyla_relative_16s_2022.csv")
class <- read.csv("data/taxonomic_annotation/class_relative_16s_2022.csv")
order <- read.csv("data/taxonomic_annotation/order_relative_16s_2022.csv")
family <- read.csv("data/taxonomic_annotation/family_relative_16s_2022.csv")
genus <- read.csv("data/taxonomic_annotation/genus_relative_16s_2022.csv")

#Remove dummy samples named "9996", "9997", "9998", "9999"
phyla <- phyla %>% filter(sample != "9996" & sample != "9997" & sample != "9998" & sample != "9999")

class <- class %>% filter(sample != "9996" & sample != "9997" & sample != "9998" & sample != "9999")

order <- order %>% filter(sample != "9996" & sample != "9997" & sample != "9998" & sample != "9999")

family <- family %>% filter(sample != "9996" & sample != "9997" & sample != "9998" & sample != "9999")

genus <- genus %>% filter(sample != "9996" & sample != "9997" & sample != "9998" & sample != "9999")

#Save data
write.csv(phyla, "data/taxonomic_annotation/phyla_relative_16s_2022.csv", row.names = FALSE)
write.csv(class, "data/taxonomic_annotation/class_relative_16s_2022.csv", row.names = FALSE)
write.csv(order, "data/taxonomic_annotation/order_relative_16s_2022.csv", row.names = FALSE)
write.csv(family, "data/taxonomic_annotation/family_relative_16s_2022.csv", row.names = FALSE)
write.csv(genus, "data/taxonomic_annotation/genus_relative_16s_2022.csv", row.names = FALSE)