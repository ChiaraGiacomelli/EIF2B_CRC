#load libraries
library(DESeq2)
library(tidyverse)
library(vsn)
library(viridis)
library(apeglm)
library(ggrepel)
library(conflicted)
library(ashr)

# directories
parent_dir <- '~/data/CGIACOME/AAA_eIF2B/AAA_Collaborative_stuff/Wurzburg_collab/AAA_Riboseq/Analysis/Ribo-seq-allreps'

plot_dir <- paste0(parent_dir,"/plots/PDF_for_paper/")
tables_dir <- paste0(parent_dir,"/Analysis/Tables_for_paper/")

#read in the most abundant transcripts per gene csv file----
most_abundant_transcripts <- read_csv(file = file.path(parent_dir, "Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

# NTC corrected ------

#create shared treatment variable
treatment <- 'plus'

#___202303___2B1 + NTC-----

#read in common variables
source("common_variables_2B1.R")

#read in data----
#the following for loop reads in each final CDS counts file and renames the counts column by the sample name and saves each data frame to a list
data_list <- list()
for (sample in RPF_sample_names) {
  df <- read_csv(file = file.path(parent_dir, "Analysis/UTR5_counts", paste0(sample, "_pc_final_UTR5_counts.csv")), col_names = T)
  colnames(df) <- c("transcript", sample)
  data_list[[sample]] <- df
}

#merge all data within the above list using reduce
#Using full join retains all transcripts, but the NAs need to be replaced with 0 as these are transcipts that had 0 counts in that sample
#DESeq2 needs the transcripts to be as rownames, not as a column
data_list %>%
  purrr::reduce(full_join, by = "transcript") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames("transcript") -> RPF_counts

#create a data frame with the condition/batch information----
#you need to make sure this data frame is correct for your samples
sample_info <- data.frame(row.names = RPF_sample_names,
                          shRNA = factor(c(rep('NTC', 8), rep(shRNA,8))),
                          condition = factor(c(rep("min", 4), rep(treatment, 4), rep("min", 4), rep(treatment, 4))),
                          batch = factor(c(rep(1:4, 4))))

#print the data frame to visually check it has been made as expected
sample_info

#only way ton convince DEseq2 to use NTC as a reference is to add it BEFORE it goes into the analyis - this way the final comparison is shRNA2B1:conditionplus
sample_info$shRNA <- relevel(sample_info$shRNA, ref = "NTC")

#make a DESeq data set from imported data----
DESeq2data <- DESeqDataSetFromMatrix(countData = RPF_counts,
                                     colData = sample_info,
                                     design = ~ batch + shRNA + condition + shRNA:condition)
DESeq2data$condition <- relevel(DESeq2data$condition, ref = "min")

model.matrix(~ batch + shRNA + condition + shRNA:condition, data=sample_info)

#pre-filter to remove genes with less than an average of 10 counts across all samples----
keep <- rowMeans(counts(DESeq2data)) >= 10
table(keep)
DESeq2data <- DESeq2data[keep,]

#run DESeq on DESeq data set----
dds <- DESeq(DESeq2data)

#extract results for each comparison----
resultsNames(dds) #need the name or coef for further analyses
results(dds) # can ready check that at the moment the results is for your last comparison, which is the integrated info of shRNA 2B1 + treatment
res <- results(dds)

#summarise results----
summary(res)

#apply LFC shrinkage for each comparison----
lfc_shrink <- lfcShrink(dds, coef="shRNA2B1.conditionplus", type="apeglm")

### plot MA on the spot to check for problemns -----
# if anything weird happens with apeglm shrinkage, consider other shrinking methods
FDR = 0.1
DESeq2::plotMA(res, alpha = FDR, main = paste0("Unshrunken MA-Plot (FDR = ", FDR, ")"), ylim = c(-5,5))
DESeq2::plotMA(lfc_shrink, alpha = FDR, main = paste0("apeglm MA-Plot (FDR = ", FDR, ")"), ylim = c(-5,5))

#write reslts to csv----
as.data.frame(lfc_shrink[order(lfc_shrink$padj),]) %>%
  rownames_to_column("transcript") %>%
  inner_join(most_abundant_transcripts, by = "transcript") %>%
  relocate("gene", "gene_sym", .after = "transcript")-> DEseq2_output
write_csv(DEseq2_output, file = paste0(tables_dir, paste0("RPFs_5UTRs_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))

#extract normalised counts and plot SD vs mean----
ntd <- normTransform(dds) #this gives log2(n + 1)
vsd <- vst(dds, blind=FALSE) #Variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) #Regularized log transformation

meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

#Regularized log transformation looks preferable for this data. Check for your own data and select the appropriate one
#write out normalised counts data----
as.data.frame(assay(rld)) %>%
  rownames_to_column("transcript") %>%
  inner_join(most_abundant_transcripts, by = "transcript") %>%
  relocate("gene", "gene_sym", .after = "transcript") -> normalised_counts
write_csv(normalised_counts, file = paste0(tables_dir, paste0("RPFs_5UTRs_", shRNA, "_NTCcorrected_normalised_counts.csv")))

#plot PCA----
pcaData <- plotPCA(rld, intgroup=c("condition", "batch","shRNA"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file = paste0(plot_dir, paste0(shRNA, "_NTCcorrected_RPFs_5UTRs_PCA.pdf")), width = 8, height = 7.5)

ggplot(pcaData, aes(PC1, PC2, color=shRNA, shape=condition, label = batch)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  scale_color_manual(values = PCA_cols) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  ggtitle(paste(shRNA, "RPFs 5'UTR")) +
  geom_text_repel(aes(label=batch), colour = 'black',size = 6)
dev.off()

#apply batch correct and re-plot heatmap and PCA----
mat <- assay(rld)
mat <- limma::removeBatchEffect(mat, rld$batch)
assay(rld) <- mat

#PCA
pcaData <- plotPCA(rld, intgroup=c("condition", "batch","shRNA"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file = paste0(plot_dir, paste0(shRNA, "_NTCcorrected_RPFs_5UTRs_batchCorrect_PCA.pdf")), width = 8, height = 7.5)
ggplot(pcaData, aes(PC1, PC2, color=shRNA, shape=condition, label = batch)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  scale_color_manual(values = PCA_cols) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  ggtitle(paste(shRNA, "RPFs 5'UTR batch corrected")) +
  geom_text_repel(aes(label=batch), colour = 'black',size = 6)
dev.off()

#___202303___2B4 + NTC-----

#read in common variables
source("common_variables_2B4.R")

#read in data----
#the following for loop reads in each final CDS counts file and renames the counts column by the sample name and saves each data frame to a list
data_list <- list()
for (sample in RPF_sample_names) {
  df <- read_csv(file = file.path(parent_dir, "Analysis/UTR5_counts", paste0(sample, "_pc_final_UTR5_counts.csv")), col_names = T)
  colnames(df) <- c("transcript", sample)
  data_list[[sample]] <- df
}

#merge all data within the above list using reduce
#Using full join retains all transcripts, but the NAs need to be replaced with 0 as these are transcipts that had 0 counts in that sample
#DESeq2 needs the transcripts to be as rownames, not as a column
data_list %>%
  purrr::reduce(full_join, by = "transcript") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames("transcript") -> RPF_counts

#create a data frame with the condition/batch information----
#you need to make sure this data frame is correct for your samples
sample_info <- data.frame(row.names = RPF_sample_names,
                          shRNA = factor(c(rep('NTC', 8), rep(shRNA,8))),
                          condition = factor(c(rep("min", 4), rep(treatment, 4), rep("min", 4), rep(treatment, 4))),
                          batch = factor(c(rep(1:4, 4))))

#print the data frame to visually check it has been made as expected
sample_info

#only way ton convince DEseq2 to use NTC as a reference is to add it BEFORE it goes into the analyis - this way the final comparison is shRNA2B1:conditionplus
sample_info$shRNA <- relevel(sample_info$shRNA, ref = "NTC")

#make a DESeq data set from imported data----
DESeq2data <- DESeqDataSetFromMatrix(countData = RPF_counts,
                                     colData = sample_info,
                                     design = ~ batch + shRNA + condition + shRNA:condition)
DESeq2data$condition <- relevel(DESeq2data$condition, ref = "min")

model.matrix(~ batch + shRNA + condition + shRNA:condition, data=sample_info)

#pre-filter to remove genes with less than an average of 10 counts across all samples----
keep <- rowMeans(counts(DESeq2data)) >= 10
table(keep)
DESeq2data <- DESeq2data[keep,]

#run DESeq on DESeq data set----
dds <- DESeq(DESeq2data)

#extract results for each comparison----
resultsNames(dds) #need the name or coef for further analyses
results(dds) # can ready check that at the moment the results is for your last comparison, which is the integrated info of shRNA 2B1 + treatment
res <- results(dds)

#summarise results----
summary(res)

#apply LFC shrinkage for each comparison----
lfc_shrink <- lfcShrink(dds, coef="shRNA2B4.conditionplus", type="apeglm")

### plot MA on the spot to check for problemns -----
# if anything weird happens with apeglm shrinkage, consider other shrinking methods
FDR = 0.1
DESeq2::plotMA(res, alpha = FDR, main = paste0("Unshrunken MA-Plot (FDR = ", FDR, ")"), ylim = c(-5,5))
DESeq2::plotMA(lfc_shrink, alpha = FDR, main = paste0("apeglm MA-Plot (FDR = ", FDR, ")"), ylim = c(-5,5))

#write reslts to csv----
as.data.frame(lfc_shrink[order(lfc_shrink$padj),]) %>%
  rownames_to_column("transcript") %>%
  inner_join(most_abundant_transcripts, by = "transcript") %>%
  relocate("gene", "gene_sym", .after = "transcript")-> DEseq2_output
write_csv(DEseq2_output, file = paste0(tables_dir, paste0("RPFs_5UTRs_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))

#extract normalised counts and plot SD vs mean----
ntd <- normTransform(dds) #this gives log2(n + 1)
vsd <- vst(dds, blind=FALSE) #Variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) #Regularized log transformation

meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

#Regularized log transformation looks preferable for this data. Check for your own data and select the appropriate one
#write out normalised counts data----
as.data.frame(assay(rld)) %>%
  rownames_to_column("transcript") %>%
  inner_join(most_abundant_transcripts, by = "transcript") %>%
  relocate("gene", "gene_sym", .after = "transcript") -> normalised_counts
write_csv(normalised_counts, file = paste0(tables_dir, paste0("RPFs_5UTRs_", shRNA, "_NTCcorrected_normalised_counts.csv")))

#plot PCA----
pcaData <- plotPCA(rld, intgroup=c("condition", "batch","shRNA"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file = paste0(plot_dir, paste0(shRNA, "_NTCcorrected_RPFs_5UTRs_PCA.pdf")), width = 8, height = 7.5)

ggplot(pcaData, aes(PC1, PC2, color=shRNA, shape=condition, label = batch)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  scale_color_manual(values = PCA_cols) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  ggtitle(paste(shRNA, "RPFs 5'UTR")) +
  geom_text_repel(aes(label=batch), colour = 'black',size = 6)
dev.off()

#apply batch correct and re-plot heatmap and PCA----
mat <- assay(rld)
mat <- limma::removeBatchEffect(mat, rld$batch)
assay(rld) <- mat

#PCA
pcaData <- plotPCA(rld, intgroup=c("condition", "batch","shRNA"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file = paste0(plot_dir, paste0(shRNA, "_NTCcorrected_RPFs_5UTRs_batchCorrect_PCA.pdf")), width = 8, height = 7.5)
ggplot(pcaData, aes(PC1, PC2, color=shRNA, shape=condition, label = batch)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  scale_color_manual(values = PCA_cols) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  ggtitle(paste(shRNA, "RPFs 5'UTR batch corrected")) +
  geom_text_repel(aes(label=batch), colour = 'black',size = 6)
dev.off()

#___202303___2B5 + NTC-----

#read in common variables
source("common_variables_2B5.R")

#read in data----
#the following for loop reads in each final CDS counts file and renames the counts column by the sample name and saves each data frame to a list
data_list <- list()
for (sample in RPF_sample_names) {
  df <- read_csv(file = file.path(parent_dir, "Analysis/UTR5_counts", paste0(sample, "_pc_final_UTR5_counts.csv")), col_names = T)
  colnames(df) <- c("transcript", sample)
  data_list[[sample]] <- df
}

#merge all data within the above list using reduce
#Using full join retains all transcripts, but the NAs need to be replaced with 0 as these are transcipts that had 0 counts in that sample
#DESeq2 needs the transcripts to be as rownames, not as a column
data_list %>%
  purrr::reduce(full_join, by = "transcript") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames("transcript") -> RPF_counts

#create a data frame with the condition/batch information----
#you need to make sure this data frame is correct for your samples
sample_info <- data.frame(row.names = RPF_sample_names,
                          shRNA = factor(c(rep('NTC', 8), rep(shRNA,8))),
                          condition = factor(c(rep("min", 4), rep(treatment, 4), rep("min", 4), rep(treatment, 4))),
                          batch = factor(c(rep(1:4, 4))))

#print the data frame to visually check it has been made as expected
sample_info

#only way ton convince DEseq2 to use NTC as a reference is to add it BEFORE it goes into the analyis - this way the final comparison is shRNA2B1:conditionplus
sample_info$shRNA <- relevel(sample_info$shRNA, ref = "NTC")

#make a DESeq data set from imported data----
DESeq2data <- DESeqDataSetFromMatrix(countData = RPF_counts,
                                     colData = sample_info,
                                     design = ~ batch + shRNA + condition + shRNA:condition)
DESeq2data$condition <- relevel(DESeq2data$condition, ref = "min")

model.matrix(~ batch + shRNA + condition + shRNA:condition, data=sample_info)

#pre-filter to remove genes with less than an average of 10 counts across all samples----
keep <- rowMeans(counts(DESeq2data)) >= 10
table(keep)
DESeq2data <- DESeq2data[keep,]

#run DESeq on DESeq data set----
dds <- DESeq(DESeq2data)

#extract results for each comparison----
resultsNames(dds) #need the name or coef for further analyses
results(dds) # can ready check that at the moment the results is for your last comparison, which is the integrated info of shRNA 2B1 + treatment
res <- results(dds)

#summarise results----
summary(res)

#apply LFC shrinkage for each comparison----
lfc_shrink <- lfcShrink(dds, coef="shRNA2B5.conditionplus", type="apeglm")

### plot MA on the spot to check for problemns -----
# if anything weird happens with apeglm shrinkage, consider other shrinking methods
FDR = 0.1
DESeq2::plotMA(res, alpha = FDR, main = paste0("Unshrunken MA-Plot (FDR = ", FDR, ")"), ylim = c(-5,5))
DESeq2::plotMA(lfc_shrink, alpha = FDR, main = paste0("apeglm MA-Plot (FDR = ", FDR, ")"), ylim = c(-5,5))

#write reslts to csv----
as.data.frame(lfc_shrink[order(lfc_shrink$padj),]) %>%
  rownames_to_column("transcript") %>%
  inner_join(most_abundant_transcripts, by = "transcript") %>%
  relocate("gene", "gene_sym", .after = "transcript")-> DEseq2_output
write_csv(DEseq2_output, file = paste0(tables_dir, paste0("RPFs_5UTRs_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))

#extract normalised counts and plot SD vs mean----
ntd <- normTransform(dds) #this gives log2(n + 1)
vsd <- vst(dds, blind=FALSE) #Variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) #Regularized log transformation

meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

#Regularized log transformation looks preferable for this data. Check for your own data and select the appropriate one
#write out normalised counts data----
as.data.frame(assay(rld)) %>%
  rownames_to_column("transcript") %>%
  inner_join(most_abundant_transcripts, by = "transcript") %>%
  relocate("gene", "gene_sym", .after = "transcript") -> normalised_counts
write_csv(normalised_counts, file = paste0(tables_dir, paste0("RPFs_5UTRs_", shRNA, "_NTCcorrected_normalised_counts.csv")))

#plot PCA----
pcaData <- plotPCA(rld, intgroup=c("condition", "batch","shRNA"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file = paste0(plot_dir, paste0(shRNA, "_NTCcorrected_RPFs_5UTRs_PCA.pdf")), width = 8, height = 7.5)

ggplot(pcaData, aes(PC1, PC2, color=shRNA, shape=condition, label = batch)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  scale_color_manual(values = PCA_cols) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  ggtitle(paste(shRNA, "RPFs 5'UTR")) +
  geom_text_repel(aes(label=batch), colour = 'black',size = 6)
dev.off()

#apply batch correct and re-plot heatmap and PCA----
mat <- assay(rld)
mat <- limma::removeBatchEffect(mat, rld$batch)
assay(rld) <- mat

#PCA
pcaData <- plotPCA(rld, intgroup=c("condition", "batch","shRNA"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file = paste0(plot_dir, paste0(shRNA, "_NTCcorrected_RPFs_5UTRs_batchCorrect_PCA.pdf")), width = 8, height = 7.5)
ggplot(pcaData, aes(PC1, PC2, color=shRNA, shape=condition, label = batch)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  scale_color_manual(values = PCA_cols) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
  ggtitle(paste(shRNA, "RPFs 5'UTR batch corrected")) +
  geom_text_repel(aes(label=batch), colour = 'black',size = 6)
dev.off()

