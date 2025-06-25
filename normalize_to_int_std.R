library(tidyverse)
library(vegan) #rclr transform

setwd("C:/Users/Olivia.Schwartz/OneDrive - University of Denver/Projects/Alzheimers NAU/20250527_HILIC_BATCH3/20250527_HILIC_BATCH3_FC")

#Inputs:
# (1) quant table (.csv)
# (2) metadata table (.csv), "filename"

ft <- read.csv("mzmine/20250530_HILIC_BATCH3_FC_batchcorrect_quant.csv", header = T, check.names = F, sep = ",")
outfile <- "HILIC_FC_batchcorrected_normtostd" #outfile name (no extension)

new_ft <- ft
#Removing Peak area extensions
colnames(new_ft) <- gsub(' Peak area','',colnames(new_ft))

#Changing the row names of the files
#new_ft <- data.frame(new_ft)
rownames(new_ft) <- paste(new_ft$'row ID',round(new_ft$'row m/z',digits = 3),round(new_ft$'row retention time',digits = 3), sep = '_')

#Picking only the files with column names containing 'mzML'
new_ft <- new_ft[,grep('mzML',colnames(new_ft))]
new_ft <- as.data.frame(t(new_ft))
norm_ft <- new_ft

#Normalize to internal standard
#FT ID = 7817
#for every feature of applicable samples, divide by the value given in 7817
#Applicable samples: Blank, GABA_Std, Frontal Cortex

std_feat <- "7817_184.094_2.289"
std_index <- grep(paste(std_feat), colnames(new_ft))
r <- 1 #row index
c <- 1 #column index

  for (r in 1:nrow(norm_ft)) {
    for (c in 1:ncol(norm_ft)) {
      norm_ft[r,c] <- norm_ft[r,c]/norm_ft[r,std_index]
}
  }

#RCLR (robut centered log ratio) transform
RCLR_t <- new_ft %>% decostand(method = "rclr")
write.csv(RCLR_t, paste(outfile,"_RCLR.csv"),row.names = TRUE)

#Reformat to match mzmine output formatting
norm_ft<-as.data.frame(t(norm_ft))
colnames(norm_ft) <- paste0(colnames(norm_ft), " Peak area")
norm_ft$X <- rownames(norm_ft)
norm_ft[c('row ID', 'row m/z', "row retention time")] <- str_split_fixed(norm_ft$X, '_', 3)
norm_ft <- norm_ft[, !names(norm_ft) %in% "X"]
norm_ft <- norm_ft[, c((ncol(norm_ft)-2):(ncol(norm_ft)),1:(ncol(norm_ft)-3))]
rownames(norm_ft) <- NULL

#Recombining with the first 13 columns output from mzmine
ft_13col <- ft[,c(1,4:13)]
ft_13col$`row ID` <- as.character(ft_13col$`row ID`)

output_norm <- inner_join(ft_13col, norm_ft)

write.csv(output_norm, file=paste(outfile,".csv"))

