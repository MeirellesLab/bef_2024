#This is the main script to perform all analysis of the BEF microbial community
#The script are not ready yet, but will be improved as we proceed with the analysis

#Check libraries and install if necessary
if (!require("vegan")) install.packages("vegan")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggfortify")) install.packages("ggfortify")
if (!require("plotly")) install.packages("plotly")

#Load input data
phyla <- read.csv("inputs//bef_phyla_matrix_relative_2023_06_07.csv", header = TRUE, sep = ",")

metadata <- read.delim("inputs//metadados_micro.csv", header = TRUE, sep = ",")

#Rename id_micro column to "Sample"
metadata <- metadata %>% rename(Sample = id_micro)

#Merge metadata and phyla data
phyla_metadata <- merge(phyla, metadata, by = "Sample")

#Now, perform a PCA analysis
source("src//pca_script.R")
