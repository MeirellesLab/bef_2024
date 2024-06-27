################################################################
# Script to perform PCA analysis with BEF emission gases
# Mariana Faustino | 24/08/23
################################################################

# Load necessary libraries
library(dplyr)
library(vegan)
library(ggplot2)

setwd("C:/Users/maria/OneDrive/Área de Trabalho/BioInfo")

dir()

### 01. SEND DATA

## Load data from CSV files
#input_dados <- read.csv("bef_genus_matrix_relative_2023_06_07.csv")%>%
#  rename(category = Riqueza)

#input_dados[1:5,1:5]

variable_data <- read.csv2("metadados.csv") %>%
  dplyr::rename(Sample = id_micro)
  
#variable_data <- variable_data  %>% 
#  filter(!grepl('pasto', riqueza))

## Removal of unnecessary columns
variable_data <- variable_data[, -which(names(variable_data) %in% c("id_block", "id_plot", "Sample", "CO2"))]

variable_data <- variable_data %>%
  dplyr::rename(MO = materia_organica, MI = materia_inorganica)

### 02. PREPARE DATA 

# Assuming both files have a column named "Sample" to join them
#merged_data <- inner_join(input_dados, variable_data, by = "Sample")

# Select only numeric variables
variable_numeric <- variable_data %>% select_if(is.numeric)

# Removing Sample column
#variable_numeric <- variable_numeric %>% select(-Sample)

# Standardize the data
variable_numeric <- scale(variable_numeric)

# Perform PCA
pca <- rda(variable_numeric, scale = TRUE)



#### 03. LIPE SCRIPT

#Get summary
pca_sum<-summary(pca)

# Get explained variance
exp <- pca_sum$cont$importance

# Access specific values
exp_pc1 <- exp[2, 1]
exp_pc2 <- exp[2, 2]

exp_pc1
exp_pc2

#Combining PCA results with metadata
ord_df <- cbind(category = variable_data$riqueza, data.frame(pca_sum$sites,row.names=NULL))


ord_df$category <- factor(ord_df$category)


#Get PCA arrows
pca_arrows <-pca_sum$species
pca_arrows <- data.frame(pca_arrows)
pca_arrows<-data.frame(pc1=pca_arrows$PC1, pc2=pca_arrows$PC2, Species = rownames(pca_arrows))
pca_arrows

#Getting only the first five pca_arrows
pca_arrows <- pca_arrows[1:5,]



#Plot PCA

ord_df <- ord_df %>%
  mutate(
    category = case_when(
      category == "mono" ~ "1sp",
      category == "tri" ~ "3spp",
      category == "hexa" ~ "6spp",
      category == "div" ~ "12spp",
      category == "alt_div" ~ "24spp")) %>%
  mutate(category = as.factor(category)) %>%
  mutate(category = fct_relevel(
    category, "1sp", "3spp", "6spp", "12spp", "24spp"))

pca_plot <- ggplot(ord_df) +
  geom_point(mapping = aes(x=PC1, y=PC2, color=category),size = 8)+
  geom_segment(data = pca_arrows,
               aes(x = 0, xend = pc1, y = 0, yend = pc2),
               arrow = arrow(length = unit(0.5, "cm")), size=1, colour = "grey") +
  geom_text(data = pca_arrows, 
            aes(x = pc1*1.1, y = pc2*1.1, label = Species),
            size = 8) +
  xlim(c(-1.75, NA))+
  geom_hline(yintercept = 0, colour = "gray70") +
  geom_vline(xintercept = 0, colour = "gray70") +
  # scale_colour_manual(values=color_plot)+ 
  theme_bw()+
  theme(panel.border = element_rect(colour = "white", size=1))+
  theme(legend.key=element_rect(fill='white'))+
  theme(legend.key = element_rect(colour = "white"))+
  labs(title = "Principal Component Analysis in the Context of Environmental Variables and CH4", x = paste("PC1 (",round(exp[2,1]*100,digits=2),"%)"), y = paste("PC2 (",round(exp[2,2]*100,digits=2),"%)")) +
  theme(axis.title.x = element_text(face="bold", size=18))+
  theme(axis.title.y = element_text(face="bold", size=18))+
  theme(axis.text.x = element_text(size=10,color="black"))+
  theme(axis.text.y = element_text(size=10,color="black"))+
  theme(legend.title = element_text(color="black", size=20, face = "bold"),
        legend.text = element_text(size = 15))+
  theme(axis.line.x=element_blank())+
  theme(axis.line.y=element_blank())+
  guides(color=guide_legend("Treatments"))+
  theme( # remove the vertical grid lines
    panel.grid.major.x = element_blank() ,
    panel.grid.minor.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank())

pca_plot

#Save PCA plot

ggsave("results//pca_plot_ch4_pasto_09.png", width = 10, height = 10, dpi = 300)


