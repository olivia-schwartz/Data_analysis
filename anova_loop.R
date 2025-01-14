library(tidyverse)

# (1) Quant table (.csv)
#     As produced by mzmine
# (2) Metadata (.csv)
#     1st col: file names
setwd("C:/Users/Olivia.Schwartz/OneDrive - University of Denver/DU/Aron Lab/Experiments/20250106_NAU_HILIC/mzmine/Hippocampus_Run")

metadata <- read.csv("20250106_NAU_HILIC_Hip_16O_md.csv")
quant_table <- read.csv("20250113_Hip_16O_quant.csv", sep = ",")

## Transform quant file ##
colnames(quant_table)[1] <- "ID"
colnames(quant_table)[2] <- "mz"
colnames(quant_table)[3] <- "rt"

#remove cols 4:13
colnames(quant_table)<-gsub("..Peak.area","",colnames(quant_table)) #Remove ".Peak.area."  
quant_table <- quant_table[, -c(2:13)]

rownames(quant_table) <- quant_table$ID
quant_table <- quant_table[, -1]  


ft <- as.data.frame(t(quant_table))
ft$ID <- rownames(ft)
ft <- ft[, c(ncol(ft), 1:(ncol(ft) - 1))]
rownames(ft) <- NULL
colnames(ft)[1] <- "filename"
colnames(ft)[2:ncol(ft)] <- paste0("FT", colnames(ft)[2:ncol(ft)])

## Transform metadata ##
colnames(metadata)<-gsub("ATTRIBUTE_","",colnames(metadata)) #Remove "_ATTRIBUTE" 

## Merge ##
data_merge <- left_join(metadata, ft)
data_merge <- data_merge %>% dplyr::filter(data_merge[["Sample_Type"]] != "Standard" & data_merge[["Sample_Type"]] != "Blank")

## Anova ##

#Empty vector and list for anova results to be appended into
anova_p_values <- c()
anova_results <- list()

#Loop through each feature to calculate anovas based on Strain_Treatment, output the p-values
for (i in names(data_merge)[-c(1:14)]) {
  anova_results[[i]] <- summary(aov(get(i) ~Strain_Treatment, data = data_merge)) 
  anova_p_values[i] <- print(anova_results[[i]][[1]][["Pr(>F)"]][[1]])

}

anova_p_values <- as.data.frame(anova_p_values)
colnames(anova_p_values)[1] <- "p_value"

write.csv(anova_p_values, file = "Hip_16O_anova_pval.csv")
