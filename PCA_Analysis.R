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

setwd("/Users/liv/Library/CloudStorage/OneDrive-UniversityofDenver/Projects/Proteus Proteobactin/mzMINE/pos")
Imp_clr <- read.csv("CLR_Scaled_pos.csv", header = T, check.names = F, sep = ",")
metadata <- read.csv("md_pbt_nrp.csv")
feat_test <- "ATTRIBUTE_group" #Feature being tested for dissimilarity 

#SPECIFIC CONDITION# -> Comment out if not applicable 
# Test dissimilarity of this feature under a specific condition (ex. specific timepoint or phenotype)
# feature = "week" #Feature under which the condition falls
# feat_condition = "52" #Value for condition
# metadata <- metadata %>% dplyr::filter(metadata[[feature]] != "NA" & metadata[[feature]] == feat_condition)
#SPECIFIC CONDIITON#

plot_title <- "Proteobactin_Pbt_vs_Nrp"#Name attached to plots and file names

############INPUTS############

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
pca
ggsave(file= paste0('PCA_',plot_title,'.svg'), plot=pca, width=10, height=8)

#output loadings plot
pc_loadings <- res.pca$rotation
pc_loadings <- as_tibble(pc_loadings, rownames = "feature")
write.csv(pc_loadings, paste0('PCA_all_',feat_test,'s.csv'),row.names = TRUE)

#plot driving features
pca_plot <- fviz_pca_var(res.pca, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, title =  paste0('PCA_',plot_title)     # Avoid text overlapping
)
pca_plot
ggsave(file= paste0('Driving_',plot_title,'.svg') , plot=pca_plot, width=10, height=8)

# PERMANOVA
# metadata <- data_merge[, 1]
# metadata_df <- data.frame(metadata = metadata)
# metabolites <- data_merge[, 2:ncol(data_merge)]
# dist_metabolites <- vegdist(metabolites, method = "euclidean", na.rm = TRUE)
# permanova_all <- adonis2(dist_metabolites ~ metadata_df$metadata, metadata_df, na.action = na.omit)

