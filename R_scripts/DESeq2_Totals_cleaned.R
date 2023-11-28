#load libraries
library(DESeq2)
library(tidyverse)
library(tximport)
library(vsn)
library(viridis)
library(ashr)
library(ggrepel)
#library(conflicted)

# NTC corrected ------

#create shared treatment variable
treatment <- 'plus'

# directories
parent_dir <- '~/data/CGIACOME/AAA_eIF2B/AAA_Collaborative_stuff/Wurzburg_collab/AAA_Riboseq/Analysis/Ribo-seq-allreps'

plot_dir <- paste0(parent_dir,"/plots/PDF_for_paper/")
tables_dir <- paste0(parent_dir,"/Analysis/Tables_for_paper/")

#read in gene to transcript IDs map and rename and select ENSTM and ENSGM columns----
#this is used by DESeq2 and needs to be in this structure
tx2gene <- read_tsv(file = "~/data/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_gene_IDs.txt", col_names = F)
tx2gene %>%
  dplyr::rename(GENEID = X1,
                TXNAME = X2) %>%
  select(TXNAME, GENEID) -> tx2gene

#read in the most abundant transcripts per gene csv file----
most_abundant_transcripts <- read_csv(file = file.path(parent_dir, "Analysis/most_abundant_transcripts/most_abundant_transcripts_IDs.csv"))

#__202310____2B1 NTC corrected__-----
#read in common variables
source("common_variables_2B1.R")

#import rsem data----
#set directory where rsem output is located
rsem_dir <- file.path(parent_dir, 'rsem')

#create a named vector of files (with path)
files <- file.path(rsem_dir, paste0(Total_sample_names, ".isoforms.results"))
names(files) <- Total_sample_names

#import data with txi
txi <- tximport(files, type="rsem", tx2gene=tx2gene)

#create a data frame with the condition/batch information----
#you need to make sure this data frame is correct for your samples, the below creates one for a n=3 with EFT226 treatment.

sample_info <- data.frame(row.names = Total_sample_names,
                          shRNA = factor(c(rep('NTC', 8), rep(shRNA,8))),
                          condition = factor(c(rep("min", 4), rep(treatment, 4), rep("min", 4), rep(treatment, 4))),
                          batch = factor(c(rep(1:4, 4))))

sample_info$shRNA <- relevel(sample_info$shRNA, ref = "NTC")

#print the data frame to visually check it has been made as expected
sample_info

#make a DESeq data set from imported data----
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sample_info,
                                   design = ~ batch + shRNA + condition + shRNA:condition)

model.matrix(~ batch + shRNA + condition + shRNA:condition, data=sample_info)

#make sure levels are set appropriately so that Ctrl is "untreated"
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "min")

#pre-filter to remove genes with less than an average of 10 counts across all samples----
keep <- rowMeans(counts(ddsTxi)) >= 10
table(keep)
ddsTxi <- ddsTxi[keep,]

#run DESeq on DESeq data set----
dds <- DESeq(ddsTxi)

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
DESeq2::plotMA(res, alpha = FDR, main = paste0(shRNA, " unshrunken MA-Plot (FDR = ", FDR, ")"), ylim = c(-5,5))
DESeq2::plotMA(lfc_shrink, alpha = FDR, main = paste0(shRNA, " apeglm MA-Plot (FDR = ", FDR, ")"), ylim = c(-5,5))

#write reslts to csv----
as.data.frame(lfc_shrink[order(lfc_shrink$padj),]) %>%
  rownames_to_column("gene") %>%
  inner_join(most_abundant_transcripts, by = "gene") %>%
  relocate("transcript") %>%
  relocate("gene", "gene_sym", .after = "transcript")-> DEseq2_output
write_csv(DEseq2_output, file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))

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
  rownames_to_column("gene") %>%
  inner_join(most_abundant_transcripts, by = "gene") %>%
  relocate("transcript") %>%
  relocate("gene", "gene_sym", .after = "transcript") -> normalised_counts
write_csv(normalised_counts, file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_normalised_counts.csv")))

#plot PCA----
pcaData <- plotPCA(rld, intgroup=c("shRNA","condition", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file = paste0(plot_dir, shRNA, "_NTCcorrected_Totals_PCA.pdf"), width = 8, height = 7.5)
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
  ggtitle(paste(shRNA, "Totals")) +
  geom_text_repel(aes(label=batch), colour = 'black',size = 6)
dev.off()

#apply batch correct and re-plot heatmap and PCA----
mat <- assay(rld)
mat <- limma::removeBatchEffect(mat, rld$batch)
assay(rld) <- mat

#PCA
pcaData <- plotPCA(rld, intgroup=c("condition", "batch","shRNA"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file = paste0(plot_dir, shRNA, "_NTCcorrected_Totals_batchCorrect_PCA.pdf"), width = 8, height = 7.5)
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
  ggtitle(paste(shRNA, "Totals batch corrected")) +
  geom_text_repel(aes(label=batch), colour = 'black',size = 6)
dev.off()

#__202310____2B4 NTC corrected__-----
#read in common variables
source("common_variables_2B4.R")

#import rsem data----
#set directory where rsem output is located
rsem_dir <- file.path(parent_dir, 'rsem')

#create a named vector of files (with path)
files <- file.path(rsem_dir, paste0(Total_sample_names, ".isoforms.results"))
names(files) <- Total_sample_names

#import data with txi
txi <- tximport(files, type="rsem", tx2gene=tx2gene)

#create a data frame with the condition/batch information----
#you need to make sure this data frame is correct for your samples, the below creates one for a n=3 with EFT226 treatment.

sample_info <- data.frame(row.names = Total_sample_names,
                          shRNA = factor(c(rep('NTC', 8), rep(shRNA,8))),
                          condition = factor(c(rep("min", 4), rep(treatment, 4), rep("min", 4), rep(treatment, 4))),
                          batch = factor(c(rep(1:4, 4))))

sample_info$shRNA <- relevel(sample_info$shRNA, ref = "NTC")

#print the data frame to visually check it has been made as expected
sample_info

#make a DESeq data set from imported data----
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sample_info,
                                   design = ~ batch + shRNA + condition + shRNA:condition)

model.matrix(~ batch + shRNA + condition + shRNA:condition, data=sample_info)

#make sure levels are set appropriately so that Ctrl is "untreated"
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "min")

#pre-filter to remove genes with less than an average of 10 counts across all samples----
keep <- rowMeans(counts(ddsTxi)) >= 10
table(keep)
ddsTxi <- ddsTxi[keep,]

#run DESeq on DESeq data set----
dds <- DESeq(ddsTxi)

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
DESeq2::plotMA(res, alpha = FDR, main = paste0(shRNA, " unshrunken MA-Plot (FDR = ", FDR, ")"), ylim = c(-5,5))
DESeq2::plotMA(lfc_shrink, alpha = FDR, main = paste0(shRNA, " apeglm MA-Plot (FDR = ", FDR, ")"), ylim = c(-5,5))

#write reslts to csv----
as.data.frame(lfc_shrink[order(lfc_shrink$padj),]) %>%
  rownames_to_column("gene") %>%
  inner_join(most_abundant_transcripts, by = "gene") %>%
  relocate("transcript") %>%
  relocate("gene", "gene_sym", .after = "transcript")-> DEseq2_output
write_csv(DEseq2_output, file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))

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
  rownames_to_column("gene") %>%
  inner_join(most_abundant_transcripts, by = "gene") %>%
  relocate("transcript") %>%
  relocate("gene", "gene_sym", .after = "transcript") -> normalised_counts
write_csv(normalised_counts, file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_normalised_counts.csv")))


#plot PCA----
pcaData <- plotPCA(rld, intgroup=c("shRNA","condition", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file = paste0(plot_dir, shRNA, "_NTCcorrected_Totals_PCA.pdf"), width = 8, height = 7.5)
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
  ggtitle(paste(shRNA, "Totals")) +
  geom_text_repel(aes(label=batch), colour = 'black',size = 6)
dev.off()

#apply batch correct and re-plot heatmap and PCA----
mat <- assay(rld)
mat <- limma::removeBatchEffect(mat, rld$batch)
assay(rld) <- mat

#PCA
pcaData <- plotPCA(rld, intgroup=c("condition", "batch","shRNA"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file = paste0(plot_dir, shRNA, "_NTCcorrected_Totals_batchCorrect_PCA.pdf"), width = 8, height = 7.5)
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
  ggtitle(paste(shRNA, "Totals batch corrected")) +
  geom_text_repel(aes(label=batch), colour = 'black',size = 6)
dev.off()

#__202310____2B5 NTC corrected__-----
#read in common variables
source("common_variables_2B5.R")

#import rsem data----
#set directory where rsem output is located
rsem_dir <- file.path(parent_dir, 'rsem')

#create a named vector of files (with path)
files <- file.path(rsem_dir, paste0(Total_sample_names, ".isoforms.results"))
names(files) <- Total_sample_names

#import data with txi
txi <- tximport(files, type="rsem", tx2gene=tx2gene)

#create a data frame with the condition/batch information----
#you need to make sure this data frame is correct for your samples, the below creates one for a n=3 with EFT226 treatment.

sample_info <- data.frame(row.names = Total_sample_names,
                          shRNA = factor(c(rep('NTC', 8), rep(shRNA,8))),
                          condition = factor(c(rep("min", 4), rep(treatment, 4), rep("min", 4), rep(treatment, 4))),
                          batch = factor(c(rep(1:4, 4))))

sample_info$shRNA <- relevel(sample_info$shRNA, ref = "NTC")

#print the data frame to visually check it has been made as expected
sample_info

#make a DESeq data set from imported data----
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sample_info,
                                   design = ~ batch + shRNA + condition + shRNA:condition)

model.matrix(~ batch + shRNA + condition + shRNA:condition, data=sample_info)

#make sure levels are set appropriately so that Ctrl is "untreated"
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "min")

#pre-filter to remove genes with less than an average of 10 counts across all samples----
keep <- rowMeans(counts(ddsTxi)) >= 10
table(keep)
ddsTxi <- ddsTxi[keep,]

#run DESeq on DESeq data set----
dds <- DESeq(ddsTxi)

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
DESeq2::plotMA(res, alpha = FDR, main = paste0(shRNA, " unshrunken MA-Plot (FDR = ", FDR, ")"), ylim = c(-5,5))
DESeq2::plotMA(lfc_shrink, alpha = FDR, main = paste0(shRNA, " apeglm MA-Plot (FDR = ", FDR, ")"), ylim = c(-5,5))

#write reslts to csv----
as.data.frame(lfc_shrink[order(lfc_shrink$padj),]) %>%
  rownames_to_column("gene") %>%
  inner_join(most_abundant_transcripts, by = "gene") %>%
  relocate("transcript") %>%
  relocate("gene", "gene_sym", .after = "transcript")-> DEseq2_output
write_csv(DEseq2_output, file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))

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
  rownames_to_column("gene") %>%
  inner_join(most_abundant_transcripts, by = "gene") %>%
  relocate("transcript") %>%
  relocate("gene", "gene_sym", .after = "transcript") -> normalised_counts
write_csv(normalised_counts, file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_normalised_counts.csv")))

#plot PCA----
pcaData <- plotPCA(rld, intgroup=c("shRNA","condition", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file = paste0(plot_dir, shRNA, "_NTCcorrected_Totals_PCA.pdf"), width = 8, height = 7.5)
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
  ggtitle(paste(shRNA, "Totals")) +
  geom_text_repel(aes(label=batch), colour = 'black',size = 6)
dev.off()

#apply batch correct and re-plot heatmap and PCA----
mat <- assay(rld)
mat <- limma::removeBatchEffect(mat, rld$batch)
assay(rld) <- mat

#PCA
pcaData <- plotPCA(rld, intgroup=c("condition", "batch","shRNA"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file = paste0(plot_dir, shRNA, "_NTCcorrected_Totals_batchCorrect_PCA.pdf"), width = 8, height = 7.5)
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
  ggtitle(paste(shRNA, "Totals batch corrected")) +
  geom_text_repel(aes(label=batch), colour = 'black',size = 6)
dev.off()
