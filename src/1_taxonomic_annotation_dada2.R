#Load packages
print("Loading packages")
library(dada2)
library(tidyverse)
library(funrar)
library(stringr)
library(DECIPHER)

d <- Sys.Date() %>% str_replace_all("-", "_")

#Creating a vector with the name of directory where the samples are
path <- "inputs/inputs_sequences"

#Listing all the samples, forward and reverse
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))

print(fnFs[1:5])
#Creating a vector with samples names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

print("Plotting quality profiles")
#Plotting average quality numbers in each basepairs of samples
plotQualityProfile(fnFs[1:2])

#Creating a directory to put our filtered sequences
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Apply the filter in the sequences. You can modify parameters to remove primers
#from the start of the sequences, or remove basepairs based in the quality plots

#Remember that the first cut is the "trimLeft" option, and 
#then you have to select truncLen number based on 
#the number of [basepairs - trimLeft]
print("Filtering sequences")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220, 220),
  maxN = 0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
  compress = TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)


#Revealing errors in sequences
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

#Removing sequencing errors
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

#Merge pair-end sequences
print("Merging sequences")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 4, verbose=TRUE)
head(mergers[[1]])

#Make a table with the sequences merged
print("Making a table with sequences merged")
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

table(nchar(getSequences(seqtab)))

#Removing chimeras from the samples merged
print("Removing chimeras")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#Calculating the percentage of samples that remaining after removing chimeras
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))

#Creating a table with metadata information
print("Creating a table with metadata information")
track <- cbind(out, sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               rowSums(seqtab.nochim))

#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchimF")
rownames(track) <- sample.names

head(track)

#Taxonomic annotation using silva v138 train set
print("Taxonomic annotation using silva v138 train set")
taxa <- assignTaxonomy(seqtab.nochim, "inputs/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
#Adding species to taxonomic annotation using the same version of silva database
print("Adding species to taxonomic annotation using the same version of silva database")
taxa <- addSpecies(taxa, "inputs/silva_species_assignment_v138.1.fa.gz")
####carry on using the IdTaxa assignments
#taxa <- taxid
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#Save workspace
print("Saving workspace")
save.image(file = paste("workspace",d,".RData", sep = ""))

##### Alternative to assignTaxonomy from dada2 package #####

##The recently developed IdTaxa taxonomic classification method is 
#also available via the DECIPHER Bioconductor package.
#The paper introducing the IDTAXA algorithm reports classification 
#performance that is better than the long-time standard set by the naive 
#Bayesian classifier. 
#Here we include a code block that allows you to use IdTaxa as a drop-in 
#replacement for assignTaxonomy (and it's faster as well!). 
#Trained classifiers are available from http://DECIPHER.codes/Downloads.html.
#Download the SILVA SSU r132 (modified) file to follow along.


#To continue with this alternative, use the sequences after removing chimeras 
dna <- DNAStringSet(getSequences(seqtab.nochim))

#loading and training the dataset for DECIPHER annotation
load("SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)


#Then I recommend you to save your workspace, so you can have both alternatives
#to taxonomic annotation saved, and can use whatever and whenever you want

#Loading the workspace saved after taxonomic annotation
load("workspace2024_03_01.RData")

#Loading packages
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

#Reading metadata table
samdf <- read.csv("inputs/metadata_23.csv", sep = ";")
colnames(samdf)[1] <- "sampleID"
#samdf$sampleID <- toupper(samdf$sampleID)
row.names(samdf) <- samdf[,1]

#Using phyloseq to obtain an OTU table
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


#ASV abundance
asv_otu <- t(seqtab.nochim)

#ASV taxonomy
asv_tax <- taxa

#merging abundance and tax table
OTU_TAX_table <- merge(asv_otu, asv_tax, by=0)
#removing NA's in kingom and phylum, which does not provide any information
OTU_TAX_table <- OTU_TAX_table %>% drop_na(Kingdom)
OTU_TAX_table <- OTU_TAX_table %>% drop_na(Phylum)
#removing eukaryotes groups
OTU_TAX_table <- OTU_TAX_table[50 != "Eukaryota",]

#writting output files
write.table(OTU_TAX_table, "outputs/OTU_TAX_table-10-03-2024.csv", sep= ",", quote = F, row.names = F)

#After save, use the other script called "replace_na_values_in_taxa_groups"
#to remove NA values and replace with upper taxa names
#### Attention!! You have to change the name (or date) of the table in the script before running
source("src/replace_na_values_in_taxa_groups.R")

##### Transforming output from dada2 into abundance matrix ####
#Phyla
OTU_TAX_table <- read.csv("inputs/OTU_TAX_table_without_na-02-01-2024.csv")
#Leave just the taxonomic level you want
OTU_TAX_table_phyla <- OTU_TAX_table %>% select(-Kingdom, Phylum, -Class, -Order, -Family, -Genus)
#Creating a table with only numeric values
OTU_TAX_table_phyla_num <- OTU_TAX_table_phyla %>% select_if(is.numeric)
#Aggregate the values in samples by each phylum
data_summ <- aggregate(OTU_TAX_table_phyla_num, by=list(OTU_TAX_table_phyla$Phylum), FUN=sum)

#Transpose data, give colnames and removing the line
phyla_matrix <- t(data_summ)
colnames(phyla_matrix) <- phyla_matrix[1,]
phyla_matrix <- phyla_matrix[-1,]

#Binding numeric values with samples information
phyla_matrix <- cbind(sample=samdf$sampleID, Richness=samdf$richess,
  id_plot=samdf$id_plot, phyla_matrix)
phyla_matrix <- as.data.frame(phyla_matrix)

#Write the matrix. This is the absolute matrix, with the number of sequences in a 
# sample classified in each phylum
d <- Sys.Date() %>% str_replace_all("-", "_")
write.csv(phyla_matrix, file = paste("inputs/2nd_run/16s_phyla_absolute.csv", sep = ""), row.names = F)


#Read that table again
#Make sure to use the matrix you created on the steps above
phyla_absolute <- read.csv("inputs/2nd_run/16s_phyla_absolute.csv")
#Only numeric values
phyla_absolute_num <- phyla_absolute %>% select(-sample, -Richness, -id_plot)

phyla_absolute_num <- as.matrix(phyla_absolute_num)
#Converting values to relative abundance.
#This functions divide the value in each cell by the rowsums of a sample
phyla_relative <- make_relative(phyla_absolute_num)
phyla_relative <- as.data.frame(phyla_relative)

#The sum of rows should be 1. To certificate, run:
rowSums(phyla_relative)

#Binding with sample information again
phyla_relative <- cbind(phyla_absolute[, c("sample", "Richness", "id_plot")], phyla_relative)

#Write tables with relative abundance
write.csv(phyla_relative, file = paste("inputs/2nd_run/16s_phyla_relative.csv", sep = ""), row.names = F)


#The steps below just repeat the same that you already did, but for different
#taxonomic levels

#Class
OTU_TAX_table <- read.csv("inputs/OTU_TAX_table_without_na-02-01-2024.csv", sep = ",")

OTU_TAX_table_class <- OTU_TAX_table %>% select(-Kingdom, -Phylum, -Order, -Family, -Genus)

OTU_TAX_table_class_num <- OTU_TAX_table_class %>% select_if(is.numeric)

data_summ <- aggregate(OTU_TAX_table_class_num, by=list(OTU_TAX_table_class$Class), FUN=sum)

class_matrix <- t(data_summ)
colnames(class_matrix) <- class_matrix[1,]
class_matrix <- class_matrix[-1,]
class_matrix <- cbind(sample=samdf$sampleID, Richness=samdf$richess, id_plot=samdf$id_plot, class_matrix)
class_matrix <- as.data.frame(class_matrix)

write.csv(class_matrix, file = paste("inputs/2nd_run/16s_class_absolute.csv", sep = ""), row.names = F)

class_absolute <- read.csv("inputs/2nd_run/16s_class_absolute.csv")

class_absolute_num <- class_absolute %>% select(-sample, -Richness, -id_plot)

class_absolute_num <- as.matrix(class_absolute_num)

class_relative <- make_relative(class_absolute_num)
class_relative <- as.data.frame(class_relative)
rowSums(class_relative)

class_relative <- cbind(class_absolute[,c("sample", "Richness", "id_plot")], class_relative)

write.csv(class_relative, file = paste("inputs/2nd_run/16s_class_relative.csv", sep = ""), row.names = F)

#Order
OTU_TAX_table <- read.csv("inputs/OTU_TAX_table_without_na-02-01-2024.csv", sep = ",")

OTU_TAX_table_order <- OTU_TAX_table %>% select(-Kingdom, -Phylum, -Class, -Family, -Genus)

OTU_TAX_table_order_num <- OTU_TAX_table_order %>% select_if(is.numeric)

data_summ <- aggregate(OTU_TAX_table_order_num, by=list(OTU_TAX_table_order$Order), FUN=sum)

order_matrix <- t(data_summ)
colnames(order_matrix) <- order_matrix[1,]
order_matrix <- order_matrix[-1,]
order_matrix <- cbind(sample=samdf$sampleID, Richness=samdf$richess, id_plot=samdf$id_plot, order_matrix)
order_matrix <- as.data.frame(order_matrix)

write.csv(order_matrix, file = paste("inputs/2nd_run/16s_order_absolute.csv", sep = ""), row.names = F)

order_absolute <- read.csv("inputs/2nd_run/16s_order_absolute.csv")

order_absolute_num <- order_absolute %>% select(-sample, -Richness, -id_plot)

order_absolute_num <- as.matrix(order_absolute_num)

order_relative <- make_relative(order_absolute_num)
order_relative <- as.data.frame(order_relative)
rowSums(order_relative)

order_relative <- cbind(order_absolute[, c("sample", "Richness", "id_plot")], order_relative)

write.csv(order_relative, file = paste("inputs/2nd_run/16s_order_relative.csv", sep = ""), row.names = F)

#Family
OTU_TAX_table <- read.csv("inputs/OTU_TAX_table_without_na-02-01-2024.csv", sep = ",")

OTU_TAX_table_family <- OTU_TAX_table %>% select(-Kingdom, -Phylum, -Class, -Order, -Genus)

OTU_TAX_table_family_num <- OTU_TAX_table_family %>% select_if(is.numeric)

data_summ <- aggregate(OTU_TAX_table_family_num, by=list(OTU_TAX_table_family$Family), FUN=sum)

family_matrix <- t(data_summ)
colnames(family_matrix) <- family_matrix[1,]
family_matrix <- family_matrix[-1,]
family_matrix <- cbind(sample=samdf$sampleID, Richness=samdf$richess, id_plot=samdf$id_plot, family_matrix)
family_matrix <- as.data.frame(family_matrix)

write.csv(family_matrix, file = paste("inputs/2nd_run/16s_family_absolute.csv", sep = ""), row.names = F)

family_absolute <- read.csv("inputs/2nd_run/16s_family_absolute.csv")

family_absolute_num <- family_absolute %>% select(-sample, -Richness, -id_plot)

family_absolute_num <- as.matrix(family_absolute_num)

family_relative <- make_relative(family_absolute_num)
family_relative <- as.data.frame(family_relative)
rowSums(family_relative)

family_relative <- cbind(family_absolute[,c("sample", "Richness", "id_plot")], family_relative)

write.csv(family_relative, file = paste("inputs/2nd_run/16s_family_relative.csv", sep = ""), row.names = F)

#Genus
OTU_TAX_table <- read.csv("inputs/OTU_TAX_table_without_na-02-01-2024.csv", sep = ",")

OTU_TAX_table_genus <- OTU_TAX_table %>% select(-Kingdom, -Phylum, -Class, -Order, -Family)

OTU_TAX_table_genus_num <- OTU_TAX_table_genus %>% select_if(is.numeric)

data_summ <- aggregate(OTU_TAX_table_genus_num, by=list(OTU_TAX_table_genus$Genus), FUN=sum)

genus_matrix <- t(data_summ)

colnames(genus_matrix) <- genus_matrix[1,]

genus_matrix <- genus_matrix[-1,]

genus_matrix <- cbind(sample=samdf$sampleID,
  Richness=samdf$richess,
  id_plot=samdf$id_plot,
  genus_matrix)
genus_matrix <- as.data.frame(genus_matrix)

write.csv(genus_matrix,
  file = paste("inputs/2nd_run/16s_genus_absolute.csv", sep = ""),
    row.names = F)

genus_absolute <- read.csv("inputs/2nd_run/16s_genus_absolute.csv")

genus_absolute_num <- genus_absolute %>%
  select(-sample, -Richness, -id_plot)

genus_absolute_num <- as.matrix(genus_absolute_num)

genus_relative <- make_relative(genus_absolute_num)
genus_relative <- as.data.frame(genus_relative)
rowSums(genus_relative)

genus_relative <- cbind(genus_absolute[, c("sample", "Richness", "id_plot")],
  genus_relative)

write.csv(genus_relative,
  file = paste("inputs/2nd_run/16s_genus_relative.csv", sep = ""),
  row.names = F)
