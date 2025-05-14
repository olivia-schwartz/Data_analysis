#calling the necessary packages:
library(ggplot2)
library(stringr)
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(reshape)     
library(gplots)
library(pheatmap)
library(heatmap3)
library(magrittr) 
library(tibble) 
library(RColorBrewer)
library(paletteer)
library(clr)
library(tidyverse) #used for data science. The eight core packages inside this library are: ggplot2 (data visualisation), dplyr (data manipulation), tidyr, readr, purrr, tibble, stringr, and forcats
library(KODAMA) # to use the normalisation function
library(ggrepel) #mainly used to repel overlapping text labels in ggplots
library(vegan) #popular library for analysing ecological diversity and for multivariate analysis of community data. Here, we use it for PCoA
library(svglite) #to save the plots in support vector graphics (svg) format
library(factoextra) #for extracting and visualizing outputs of multivariate analyses such as PCA, k-means
library(ggsci) #provides color palettes for ggplot2 that can be used for scientific journals
library(matrixStats) #contains highly optimized functions to perform statistics on matrix data
library(cowplot) #efficient functions to arrange several plots

#Files Input:
# (1) Clr-transported feature table as .csv
# (2) metadata as .csv
# where both tables have file names in their first column

############INPUTS############

#Files
setwd("C:/Users/Olivia.Schwartz/OneDrive - University of Denver/Projects/Alzheimers NAU/202504_HILIC_IntStd_Runs/Frontal Cortex")
Imp_clr <- read.csv("processed/IMP_clr.csv", header = T, check.names = F, sep = ",")
metadata <- read.csv("202504_HILIC_FC_md_16O.csv")

feat_test <- "Treatment" #Feature being tested for dissimilarity 
plot_title <- "HILIC FC"

#Data filters
  #Set unused filters to NA
  #Attribute = Column name
  #Condition = value to keep after filtering

attribute_1 <- "Isotope"
condition_1 <- "16O"

attribute_2 <- "Sample_Type"
condition_2 <- "Frontal Cortex"

attribute_3 <- "Strain"
condition_3 <- "WT"

############INPUTS############


colnames(metadata)<-gsub("ATTRIBUTE_","",colnames(metadata)) #Remove "ATTRIBUTE_"

if (is.na(attribute_1)) { 
} else if (is.na(attribute_2)) {
  # Filter 1 condition
  metadata <- metadata %>% dplyr::filter(metadata[[attribute_1]] != "NA" & metadata[[attribute_1]] == condition_1)
  plot_title <- paste0(plot_title,'_',attribute_1,'_',condition_1,'_',feat_test)
} else if (is.na(attribute_3)) {
  # Filter 2 conditions
  metadata <- metadata %>% dplyr::filter(metadata[[attribute_1]] == condition_1 & metadata[[attribute_2]] == condition_2)
  plot_title <- paste0(plot_title,'_',condition_1,'_',condition_2,'_',feat_test)
} else {
  # Filter 3 conditions
  metadata <- metadata %>% dplyr::filter(metadata[[attribute_1]] == condition_1 & metadata[[attribute_2]] == condition_2 & metadata[[attribute_3]] == condition_3)
  plot_title <- paste0(plot_title,'_',condition_1,'_',condition_2,'_',condition_3,'_',feat_test)
}

colnames(metadata)[1] <- "filename"
colnames(Imp_clr)[1] <- "filename"


#merge metadata column with clr-transformed feature table
data_merge <- metadata %>% dplyr::select("filename",feat_test) %>% 
  left_join(Imp_clr) %>% 
  column_to_rownames("filename")

#make scree plot 
res.pca <- data_merge %>% dplyr::filter(data_merge[[feat_test]] != "NA") %>%
  dplyr::select(-feat_test) %>% 
  prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)

#plot pca plot
pca <- fviz_pca_ind(res.pca, col.ind = (data_merge %>% dplyr::filter(data_merge[[feat_test]] != "NA"))[[feat_test]], addEllipses = TRUE, label = "none",title =  paste0('PCA_',plot_title))



#output loadings plot
pc_loadings <- res.pca$rotation
pc_loadings <- as_tibble(pc_loadings, rownames = "feature")
write.csv(pc_loadings, paste0('PCA_all_',plot_title,'s.csv'),row.names = TRUE)

#plot driving features
pca_plot <- fviz_pca_var(res.pca, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, title =  paste0('PCA_',plot_title)     # Avoid text overlapping
) 
pca_plot
ggsave(file= paste0('Driving_',plot_title,'.svg') , plot=pca_plot, width=5, height=5)

pca_feat_table <- as.data.frame(pca_plot[["data"]][["name"]])
write.csv(pca_feat_table, paste0('PCA_feat_',plot_title,'.csv'))
# PERMANOVA
metadata <- data_merge[, 1]
metadata_df <- data.frame(metadata = metadata)
metabolites <- data_merge[, 2:ncol(data_merge)]
dist_metabolites <- vegdist(metabolites, method = "euclidean", na.rm = TRUE)
permanova_all <- adonis2(dist_metabolites ~ metadata_df$metadata, metadata_df, na.action = na.omit)
write.csv(permanova_all, paste0('Permanova',plot_title,'.csv'),row.names = TRUE)

pca_pvalue <- permanova_all$`Pr(>F)`[1]

# Add text in the top-left corner



pca <- pca + annotate("text", x = -20, y = 13, label = paste0('Pr(>F):',pca_pvalue), size = 5)
pca
ggsave(file= paste0('PCA_',plot_title,'.svg'), plot=pca, width=10, height=8)