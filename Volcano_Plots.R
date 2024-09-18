library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(dplyr)

#Files Input:
# (1) Quant table as .csv
# (2) metadata as .csv
# where quant table column 1 is the feature IDs
# metadata table column 1 is the file names

############INPUTS############
setwd("/Users/liv/Library/CloudStorage/OneDrive-UniversityofDenver/Projects/Proteus Proteobactin/mzMINE/pos")
dset<-read.csv("Normalised_Quant_table_pos.csv", header = TRUE, check.names = FALSE)
metadata <- read.csv("md_new(1).csv")

feature <- "ATTRIBUTE_group"
ko_condition <- "Nrp_Pbt_mutant" #condiiton 1
wt_condition <- "Nrp_mutant" #condition 2

plot_title <- "Proteus_DoubleKO_Nrp_mutant" #Name attached to files and plot titles
############INPUTS############

data<- as.matrix(dset[,2:ncol(dset)])
featureIDs <- dset[,1]
# Check the loaded dataset
# Dimension of the dataset
dim(data) 
head(data)

#Grab the file names for specified conditions
colnames(metadata[1]) <- "filename"
metadata_ko <- metadata %>% dplyr::filter(metadata[[feature]] == ko_condition)
metadata_wt <- metadata %>% dplyr::filter(metadata[[feature]] == wt_condition)
ko_names <- metadata_ko$filename
wt_names <- metadata_wt$filename

#Separate data matrix according to respective file names for each condition
ko <- data[, (colnames(data) %in% ko_names)] # "Condition 1"
wt <- data[, (colnames(data) %in% wt_names)] # "Condition 2"

# Compute the means of the samples of each condition
wt.mean = apply(wt, 1, mean)
ko.mean = apply(ko, 1, mean)

# Just get the maximum of all the means for plotting
limit = max(wt.mean, ko.mean)
means <- cbind.data.frame(wt.mean, ko.mean)

#Make sure full dataset only contains tested conditions
data <- cbind(as.data.frame(ko), as.data.frame(wt))
###############################################################
######### WT VS KO ANALYSIS FROM THIS POINT ONWARD ###########
DEres.df <- as.data.frame(data)
DEres.df <- cbind.data.frame(featureIDs, DEres.df)
#compute fold
fold = wt.mean / ko.mean
log2fold <- log2(fold)
#add to results dataframe
DEres.df$logFC <- log2fold

# Generate histogram of the fold differences
foldhist <- ggplot(DEres.df, aes(x=logFC)) +
  geom_histogram( binwidth=0.01, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle(paste0('Fold_Change_',plot_title)) +
  labs(x= "Log fold change in abundance between groups", y = "Count") +
  theme_minimal() 
foldhist

#compute stat sig - t-test
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics
for(i in 1 : nrow(data)) { # For each gene : 
  x = wt[i,] # WT of gene number i
  y = ko[i,] # KO of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y, na.rm=TRUE, paired=FALSE, exact=FALSE)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}
# add results to results dataframe
DEres.df$pvalue <- pvalue

# Generate -log10(pvalues)
pvallogged <- -log10(pvalue)
# add results to results dataframe
DEres.df$logpval <- pvallogged
# Generate histogram of -log10(pvalues)
pval_hist <- ggplot(DEres.df, aes(x=logpval)) +
  geom_histogram( binwidth=0.05, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle(paste0('-log10(pvalue)_',plot_title)) +
  labs(x= "Significance levels occuring in D.E. results", y = "Count") +
  xlim(0,4) +
  theme_minimal() 
pval_hist

# Volcano Plots:
# put the biological significance (fold changes)
# and statistical significance (p-value) in one plot
# Generate the volcano plot
volcanoPlot <- function(results.df, fold_cutoff=0.5,
                        pvalue_cutoff=0.01, title = 
                          paste0("Volcano_",plot_title)){
  # Inputs: ##############################################
  
  # results.df: dataframe containing columns: 
  #             1) "featureIDs": feature IDs
  #             2) "logFC": foldchange values from DE analysis
  #             3) "logpval": log-transformed p-values from DE analysis
  # fold cutoff: significance threshold to render visually on plot;
  #               denotes fold difference between mutant and wildtype
  #               also referred to as "biologial signal"
  # p_value_cutoff: significance threshold to render visually on plot;
  #                 denotes statistical significance
  ########################################################
  # create factor denoting differential expression status
  stats.df <- results.df %>% mutate(diffexpressed = 
                                      ifelse(logFC < -(fold_cutoff) & 
                                               logpval >= -log10(pvalue_cutoff),
                                             ko_condition, 
                                             ifelse(logFC > fold_cutoff & 
                                                      logpval >= -log10(pvalue_cutoff), 
                                                    wt_condition, "Insignificant")))
  
  stats.df$diffexpressed <- as.factor(stats.df$diffexpressed)
  # Base plot
  plot <- ggplot(stats.df, aes(x = logFC, y = logpval, col = diffexpressed, text = featureIDs)) +
    geom_point() +
    ggtitle(title) + 
    ylab("-log10(pvalue)") +
    geom_vline(xintercept=c(-fold_cutoff, fold_cutoff), col="red") +
    geom_hline(yintercept=-log10(pvalue_cutoff), col="darkturquoise") +
    scale_color_manual(values=c("black", "blue", "red"), 
                       name = "Differential expression\nstatus") +
    theme_minimal()
  return(plot)
  
}

library("plotly")
# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
plot <- volcanoPlot(DEres.df, fold_cutoff=0.5, pvalue_cutoff=0.01)
ggplotly(plot)

install.packages('svglite')
ggsave(file=paste0("Volcano",plot_title,".svg"), plot=plot, width=10, height=8)

write.csv(DEres.df, paste0('R_Statistics_',plot_title,".csv"), row.names=FALSE)

# 
# ################TEST FOR NORMALITY#############################
# library("dplyr")
# library("tidyverse")
# #my_data3 <- read_csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797.csv") %>% column_to_rownames("group")
# my_data2 <- read.csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797_2.csv")
# my_data2 <- my_data2 %>% column_to_rownames(., var = 'group')
# sapply(my_data2, class) 
# my_data2$norm_peak_area <- as.numeric(as.character(my_data2$norm_peak_area))
# my_matrix2 <- data.matrix(my_data2, rownames.force = NA)
# hist(my_matrix2)
# shapiro.test(my_matrix2)
# 
# #A7
# my_data_A7 <- read.csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797_A7.csv")
# my_data_A7 <- my_data_A7 %>% column_to_rownames(., var = 'group')
# my_matrixA7 <- data.matrix(my_data_A7, rownames.force = NA)
# hist(my_matrixA7)
# shapiro.test(my_matrixA7)
# 
# #KO
# my_data_A7 <- read.csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797_KO.csv")
# my_data_A7 <- my_data_A7 %>% column_to_rownames(., var = 'group')
# my_matrixA7 <- data.matrix(my_data_A7, rownames.force = NA)
# hist(my_matrixA7)
# shapiro.test(my_matrixA7)
# 
# #WT
# my_data_A7 <- read.csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797_WT.csv")
# my_data_A7 <- my_data_A7 %>% column_to_rownames(., var = 'group')
# my_matrixA7 <- data.matrix(my_data_A7, rownames.force = NA)
# hist(my_matrixA7)
# shapiro.test(my_matrixA7)
# 
# #WTlim
# my_data_A7 <- read.csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797_WT2.csv")
# my_data_A7 <- my_data_A7 %>% column_to_rownames(., var = 'group')
# my_matrixA7 <- data.matrix(my_data_A7, rownames.force = NA)
# hist(my_matrixA7)
# shapiro.test(my_matrixA7)
# 
# 
# ############################################################################
# #####box plot with wilcoxon and krustal-wallis stats for feature m/z 797####
# library("ggpubr")
# library("datatable")
# my_data2 <- read.csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797.csv")
# p <- ggboxplot(my_data2,x="group", y="norm_peak_area", color="group", add = "jitter", shape = "group", yscale="log10")
# p
# 
# my_comparisons <- list( c("A7", "KO"), c("KO", "WT"), c("A7", "WT"), c("WT", "WT_limited") )
# fig <- p + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
#   stat_compare_means(label.y = .0012)
# fig
# ggsave(file="wilcox_boxplot.svg", plot=fig, width=10, height=8)
# 
# ############################################################################
# #####box plot with anova and stats for feature m/z 797####
# 
# library("ggpubr")
# library("data.table")
# my_data2 <- read.csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797.csv")
# p <- ggboxplot(my_data2,x="group", y="norm_peak_area", color="group", add = "jitter", shape = "group", yscale="log10")
# p
# 
# my_comparisons <- list( c("A7", "KO"), c("KO", "WT"), c("A7", "WT"), c("WT", "WT_limited") )
# fig <- p + stat_compare_means(comparisons = my_comparisons, method = "t.test") + stat_compare_means(method = "anova", label.y = .0012)
# fig
# 
# ggsave(file="anova_boxplot.svg", plot=fig, width=10, height=8)
