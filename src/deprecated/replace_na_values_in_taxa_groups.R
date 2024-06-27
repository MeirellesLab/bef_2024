#Script to replace NA values with the information about upper taxa
#This script basically performs a string replacement to remove NA names and put
#the information on the column before
library(stringr)

#Read data
OTU_TAX_table <- read.csv("outputs/OTU_TAX_table-10-03-2024.csv")

#Creating a vector with only taxon ranks info
rank_info <- OTU_TAX_table %>% select(Kingdom, Phylum, Class, Order, Family, Genus)

#Class
#Identify rows number that have NA's
class_na_info <- which(is.na(rank_info$Class))
#Based on the rows number wich have NA, replace Na with the name of the column 
# before
rank_info$Class[is.na(rank_info$Class)] <- paste("Unclassified Class in", 
                                                 rank_info[class_na_info, 2])

#Order
order_na_info <- which(is.na(rank_info$Order))

rank_info$Order[is.na(rank_info$Order)] <- paste("Unclassified Order in", 
                                                 rank_info[order_na_info, 3])
#If you have NA in an Order and in a Class, you'll have some repeated information
#So, this line just removes repeated "unclassified"
rank_info$Order <- rank_info$Order %>% str_replace(" Unclassified Class in", "")

#Now you can repeat the same steps to bottom taxa
#Family
family_na_info <- which(is.na(rank_info$Family))

rank_info$Family[is.na(rank_info$Family)] <- paste("Unclassified Family in", 
                                                 rank_info[family_na_info, 4])
#Removing repeated "unclassified" information
rank_info$Family <- rank_info$Family %>% str_replace(" Unclassified Order in", "")
#Genus
genus_na_info <- which(is.na(rank_info$Genus))

rank_info$Genus[is.na(rank_info$Genus)] <- paste("Unclassified Genus in", 
                                                 rank_info[genus_na_info, 5])
#Removing repeated "unclassified" information
rank_info$Genus <- rank_info$Genus %>% str_replace(" Unclassified Family in", "")

OTU_TAX_table <- OTU_TAX_table %>% select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus)

OTU_TAX_table <- cbind(OTU_TAX_table, rank_info)

write.csv(OTU_TAX_table, "inputs/OTU_TAX_table_without_na-02-01-2024.csv", row.names = FALSE)
