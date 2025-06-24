library(tidyverse)

#This is for inputting a set of features found in one sample type and calculating the coefficient of variation (CV)

#####Inputs####

setwd("C:/Users/Olivia.Schwartz/OneDrive - University of Denver/Projects/Alzheimers NAU/20250527_HILIC_BATCH3/20250527_HILIC_BATCH3_FC")

# (1) quant table (.csv) as formatted by mzmine
# (2) metadata table (.csv), "filename"

quant_table <- read.csv("mzmine/20250527_HILIC_NAU_FC_Batch3_iimn_gnps_quant.csv", header = T, check.names = F, sep = ",")
metadata <- read.csv("C:/Users/Olivia.Schwartz/OneDrive - University of Denver/Projects/Alzheimers NAU/20250527_HILIC_BATCH3/20250527_HILIC_BATCH3_FC/20250527_HILIC_BATCH3_FC_md.csv", header = T, check.names = F, sep = ",")
outfolder <- "C:/Users/Olivia.Schwartz/OneDrive - University of Denver/Projects/Alzheimers NAU/20250527_HILIC_BATCH3/20250527_HILIC_BATCH3_FC/"

ft_CV <- c("FT5210","FT6644","FT7037","FT7799","FT7815")
ft_names <- c("Tyramine","Choline","Methionine","Proline","Tyrosine")
category <- "Sample_Type"
QC_name <- "QC"


###############

#Metadata
colnames(metadata)<-gsub("ATTRIBUTE_","",colnames(metadata)) #Remove "ATTRIBUTE_"
colnames(metadata)[1] <- "filename"
metadata <- as.data.frame(metadata)
metadata <- metadata %>% dplyr::filter(metadata[[category]] == paste(QC_name))
#Feature Table
colnames(quant_table)<-gsub(".Peak.area","",colnames(quant_table)) #Remove ".Peak.area." 
colnames(quant_table)[1] <- "ID"
quant_table <- quant_table[-c(4:13)] #Retain only ID, mz, RT
colnames(quant_table)[2] <- "mz"
colnames(quant_table)[3] <- "rt"
quant_ft <- t(quant_table)
quant_ft <- as.data.frame(quant_ft)
colnames(quant_ft) <- quant_ft[1,]

filename <- row.names(quant_ft)
quant_ft <- cbind(filename, quant_ft)
quant_ft <- quant_ft[-c(1:3), ]   # notice the -
colnames(quant_ft)[2:ncol(quant_ft)] <- paste("FT", colnames(quant_ft)[2:ncol(quant_ft)], sep = "")
quant_ft <- quant_ft[-nrow(quant_ft), ]

category_ID_vec <- metadata$filename[metadata[[category]] == QC_name]
filt_quant <- quant_ft[quant_ft$filename %in% category_ID_vec, ]
filt_quant <- as.data.frame(filt_quant)

cv_table <- as.data.frame(matrix(NA, nrow = 4, ncol = length(ft_CV)))
cv_table[1, ] <- ft_names
colnames(cv_table) <- ft_CV
rownames(cv_table) <- c("Compound","mean","SD","CV")

for (i in 1:length(ft_CV)) {
  
  ft_test <- ft_CV[i]
  mean_value <- mean(filt_quant[[ft_test]])
  sd_value <- sd(filt_quant[[ft_test]])
  cv_value <- (sd_value/mean_value)*100
  cv_table[2,i] <- mean_value
  cv_table[3,i] <- sd_value
  cv_table[4,i] <- cv_value
  
}

cv_table <- t(cv_table)

# Barplot
ggplot(cv_table, aes(x=Compound, y=CV)) + 
  geom_bar(stat = "identity") +
  coord_flip() + ggtitle("QC Mix Compounds CV") +
  theme(text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18))
ggsave(paste0(outfolder,"QC_Mix_Compounds_CV.svg"),plot=last_plot())