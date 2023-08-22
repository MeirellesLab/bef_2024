#Script to perform PCA analysis with BEF microbial community

#Select only numeric variables
phyla_metadata_numeric <- phyla_metadata %>% select_if(is.numeric)

#Removing Sample column
phyla_metadata_numeric <- phyla_metadata_numeric %>% select(-Sample)
#Perform PCA
pca <- rda(phyla_metadata_numeric, scale = TRUE)

#Get summary
pca_sum<-summary(pca)

#Get explained variance
exp<-pca_sum$cont
exp<-data.frame(exp$importance)
exp[2,2]
exp[2,1]

#Combining PCA results with metadata
ord_df <- cbind(category = phyla_metadata$Riqueza, data.frame(pca_sum$sites,row.names=NULL))

ord_df$category <- factor(ord_df$category)

#Get PCA arrows
pca_arrows <-pca_sum$species
pca_arrows <- data.frame(pca_arrows)
pca_arrows<-data.frame(pc1=pca_arrows$PC1, pc2=pca_arrows$PC2, Species = rownames(pca_arrows))
pca_arrows

#Getting only the first five pca_arrows
pca_arrows <- pca_arrows[1:5,]

#Plot PCA
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
  labs(title = NULL, x = paste("PC1 (",round(exp[2,1]*100,digits=2),"%)"), y = paste("PC2 (",round(exp[2,2]*100,digits=2),"%)")) +
  theme(axis.title.x = element_text(face="bold", size=18))+
  theme(axis.title.y = element_text(face="bold", size=18))+
  theme(axis.text.x = element_text(size=16,color="black"))+
  theme(axis.text.y = element_text(size=16,color="black"))+
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
ggsave("results//pca_plot.png", width = 10, height = 10, dpi = 300)
