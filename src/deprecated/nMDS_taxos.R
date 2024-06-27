################################################
# ANALYSI metaMDS TAXOS
# Mariana Faustino | 15/06/23
################################################

library(vegan)

setwd("C:/Users/maria/OneDrive/Área de Trabalho/BioInfo")

dir()


## CARREGANDO OS DADOS
input_dados <- read.csv("bef_genus_matrix_relative_2023_06_07.csv") %>%
    dplyr::rename(category = Riqueza)

input_dados <-input_dados  %>% 
  filter(!grepl('floresta|pasto', category))

input_dados[1:5,1:5]

library(dplyr)

data_num <-select_if((input_dados), is.numeric)
data_num <-data_num[, -1]

data_num[1:5,1:5]
rowSums(data_num)



############## metaMDS

library(vegan)

nmds_phyla <- metaMDS(data_num)

#A distancia varia conforma o comportamento dos seus dados, no caso de relativo se aconselha a utilizacao da dist. "euclidean", diferente da dist. "jaccard" que utiliza para dados binarios.

nmds_phyla <- metaMDS(data_num, distance = "euclidean") 

samples <- input_dados$Sample
category <- input_dados$category

nmds_phyla_df <- cbind(samples,category,data.frame(nmds_phyla$points))

colnames(nmds_phyla_df) <- c("Samples", "Category", "MDS1","MDS2")

attach(nmds_phyla_df)

##### MUDANCA DO NOME DAS CLASSES

nmds_phyla_df <- nmds_phyla_df %>%
mutate(
  Category = case_when(
    Category == "mono" ~ "1sp",
    Category == "tri" ~ "3spp",
    Category == "hexa" ~ "6spp",
    Category == "div" ~ "12spp",
    Category == "alt_div" ~ "24spp",
    TRUE ~ category)) %>%
  mutate(Category = as.factor(Category)) %>%
  mutate(Category = fct_relevel(
    Category, "1sp", "3spp", "6spp", "12spp", "24spp"))


##### PLOTANDO O GRAFICO DE NMDS
plot_tax <- ggplot(nmds_phyla_df, aes(x=MDS1, y=MDS2, color=Category)) + 
  geom_point(size = 4) +
  theme_bw()+
  labs(title = "Order")+
  theme(plot.title = element_text(face = "bold", size = 28), plot.tag = element_text(face = "bold", size = 30))+
  theme(legend.position = "right", legend.title = element_text(size = 28, face = "bold"))+
  theme(panel.border = element_rect(colour = "black", size=1))+
  theme(legend.key=element_rect(fill='white'))+
  theme(legend.key = element_rect(colour = "white"))+
  scale_colour_manual(values=c("#43CD80","#9400D3","#c2c212", "#3eb5bf", "#7e3c03", "#0f0476","#f72d02"))+ 
  theme(legend.text=element_text(size=18))+
  theme(axis.title.x = element_text(face="bold", size=18))+
  theme(axis.title.y = element_text(face="bold", size=18))+
  theme(axis.text.x = element_text(size=18,color="black"))+
  theme(axis.text.y = element_text(size=18,color="black"))+
  theme(axis.line.x=element_blank())+
  theme(axis.line.y=element_blank())+
  annotate("text", x = 0.1, y = 0.1, label = paste("Stress=",round(nmds_phyla$stress, digits=3)), size = 3)+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank() ,
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

plot_tax

##### MUDANCA DO NOME DAS CLASSES

#ord_df <- ord_df %>% mutate(Category = str_replace_all(Category, "alt_div", "alta diversidade") %>% str_replace_all("hexa", "hexacultura"))
#ord_df <- ord_df %>% mutate(Category = str_replace_all(Category, "mono", "monocultura") %>% str_replace_all("tri", "tricultura"))


levels(as.factor(nmds_phyla_df$Category))

category <- nmds_phyla_df$Category

################# THE MORE 

table(factor(Category))
### Descobrimos aqui que as sequencias q estao ausentes sao referentes a riqueza de hexacultura.


###### VALIDACAO ESTATISTICA DO NMDS | anosim (ANALISE DE SIMILARIDADE)

danosim <- anosim(data_num,Category)
danosim

## Não foi significativo os valores da diferenca entre os grupos a nivel de FILO, logo seguimos para analise das demais classes taxonomicas | ANOSIM statistic R: -0.006368  Significance: 0.506 



############# PLOT PHYLA final  


plot_tax <- ggplot(nmds_phyla_df, aes(x=MDS1, y=MDS2, color=Category)) + 
  geom_point(size = 4) +
  theme_bw()+
  labs(title = "Genus")+
  theme(plot.title = element_text(face = "bold", size = 19))+
  theme(legend.position = "right", legend.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_rect(colour = "black", size=1))+
  theme(legend.key=element_rect(fill='white'))+
  theme(legend.key = element_rect(colour = "white"))+
  scale_colour_manual(values=c("#43CD80","#9400D3","#c2c212", "#3eb5bf", "#7e3c03", "#0f0476","#f72d02"))+ 
  theme(legend.text=element_text(size=16))+
  theme(axis.title.x = element_text(face="bold", size=16))+
  theme(axis.title.y = element_text(face="bold", size=16))+
  theme(axis.text.x = element_text(size=16,color="black"))+
  theme(axis.text.y = element_text(size=16,color="black"))+
  theme(axis.line.x=element_blank())+
  theme(axis.line.y=element_blank())+
  annotate("text", x = 0.05, y = 0.1, label = paste("Stress=",round(nmds_phyla$stress, digits=3)), size = 4)+
  annotate("text", x = 0.05, y = 0.13, label = paste("R²=",round(danosim$statistic, digits = 3)), size = 4)+
  annotate("text", x=0.05, y=0.15, label = paste("pvalue=",round(danosim$signif, digits = 4)), size = 4)+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank() ,
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

plot_tax

d <- Sys.Date() %>% str_replace_all("-","_")

ggsave(plot = plot_tax, filename = paste("outputs/genus_nmds", d, ".png", sep = ""))



############# PLOT GENUS final  

## Não foi significativo os valores da diferenca entre os grupos a nivel de GENERO, logo seguimos para analise das demais classes taxonomicas | ANOSIM statistic R: 0.01237  Significance: 0.423 


############# PLOT CLASS final  

## Não foi significativo os valores da diferenca entre os grupos a nivel de CLASSE, logo seguimos para analise das demais classes taxonomicas | ANOSIM statistic R: -0.006368  Significance: 0.5  


############# PLOT FAMILY final  

## Não foi significativo os valores da diferenca entre os grupos a nivel de FILO, logo seguimos para analise das demais classes taxonomicas | ANOSIM statistic R: 0.009181  Significance: 0.413 


############# PLOT ORDER final  

## Não foi significativo os valores da diferenca entre os grupos a nivel de FILO, logo seguimos para analise das demais classes taxonomicas | ANOSIM statistic R: -0.003795  Significance: 0.522 

############# PLOT PHYLA final  

## Não foi significativo os valores da diferenca entre os grupos a nivel de FILO, logo seguimos para analise das demais classes taxonomicas | ANOSIM statistic R: -0.003795  Significance: 0.521 