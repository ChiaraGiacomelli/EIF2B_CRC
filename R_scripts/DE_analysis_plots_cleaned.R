# this script is used to make the scatter plots of RPF vs total RNA changes

# load libraries
library(tidyverse)
library(viridis)
library(conflicted)

# themes----
mytheme <- theme_classic()+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))

viridis_colors <- c("#000004","#2c115f", "#721f81", "#b73779","#f1605d", "#feb078")

# create a variable for what the treatment is----
treatment <- "plus"

parent_dir <- '~/data/CGIACOME/AAA_eIF2B/AAA_Collaborative_stuff/Wurzburg_collab/AAA_Riboseq/Analysis/Ribo-seq-allreps'
plot_dir <- paste0(parent_dir,"/plots/PDF_for_paper/")
tables_dir <- paste0(parent_dir,"/Analysis/Tables_for_paper/")

# 202310___2B1 NTCcorrected --------

#read in common variables
source("common_variables_2B1.R")

## write down the thresholds for significance ----
padj_thr = 0.1
log2FC_thr = 0

## read in DESeq2 output----
totals <- read_csv(file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))
RPFs <- read_csv(file = file.path(tables_dir, paste0("RPFs_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))

# ## plot volcanos----
# RPFs |>
#   dplyr::filter(!(is.na(padj))) |>
#   mutate(sig = factor(case_when(padj < padj_thr ~ "*",
#                                 padj >= padj ~ "NS"))) |>
#   ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sig))+
#   geom_point(size=3,alpha = 0.5)+
#   mytheme+
#   xlab("log2FC")+
#   ylab("-log10(padj)")+
#   ggtitle(paste(shRNA, "RPFs"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr)) +
#   geom_hline(yintercept = 1, linetype='dashed', linewidth=0.5)-> RPFs_volcano
# 
# pdf(file = paste0(plot_dir, "RPFs_", shRNA, "_NTCcorrected_volcano_fixedlimitswGrey.pdf"), width = 7, height = 7)
# print(RPFs_volcano)
# dev.off()
# 
# totals |>
#   dplyr::filter(!(is.na(padj))) |>
#   mutate(sig = factor(case_when(padj < padj_thr ~ "*",
#                                 padj >= padj_thr ~ "NS"))) |>
#   ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sig))+
#   geom_point(size=3, alpha = 0.5)+
#   mytheme+
#   xlab("log2FC")+
#   ylab("-log10(padj)")+
#   ggtitle(paste(shRNA, "Totals"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr))+
#   geom_hline(yintercept = 1, linetype='dashed', linewidth=0.5)-> totals_volcano
# 
# pdf(file = paste0(plot_dir, "Totals_", shRNA, "_NTCcorrected_volcano_fixedlimitswGrey.pdf"), width = 7, height = 7)
# print(totals_volcano)
# dev.off()
# 
# ## plot MAs----
# RPFs |>
#   dplyr::filter(!(is.na(padj))) |>
#   mutate(sig = factor(case_when(padj < padj_thr ~ "*",
#                                 padj >= padj_thr ~ "NS"))) |>
#   ggplot(aes(x = log10(baseMean), y = log2FoldChange, colour = sig))+
#   geom_point(alpha = 0.5)+
#   mytheme+
#   xlab("log10(mean expression)")+
#   ylab("log2FC")+
#   ggtitle(paste(shRNA, "RPFs"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr, "\napeglm shrinking")) -> RPFs_MA
# 
# pdf(file = paste0(plot_dir, paste0("RPFs_", shRNA, "_NTCcorrected_MA_fixedlimitswGrey.pdf")), width = 5, height = 5)
# print(RPFs_MA)
# dev.off()
# 
# totals |>
#   dplyr::filter(!(is.na(padj))) |>
#   mutate(sig = factor(case_when(padj < 0.01 ~ "*",
#                                 padj >= 0.01 ~ "NS"))) |>
#   ggplot(aes(x = log10(baseMean), y = log2FoldChange, colour = sig))+
#   geom_point(alpha = 0.5)+
#   mytheme+
#   xlab("log10(mean expression)")+
#   ylab("log2FC")+
#   ggtitle(paste(shRNA, "Totals"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr, "\napeglm shrinking")) -> totals_MA
# 
# pdf(file = paste0(plot_dir, paste0("Totals_", shRNA, "_NTCcorrected_MA_fixedlimitswGrey.pdf")), width = 5, height = 5)
# print(totals_MA)
# dev.off()

## merge RPF with totals data----
RPFs |>
  dplyr::select(gene, log2FoldChange, padj,gene_sym,transcript) |>
  dplyr::rename(RPFs_log2FC = log2FoldChange,
                RPFs_padj = padj) |>
  inner_join(totals[,c("gene", "log2FoldChange", "padj")], by = "gene") |>
  dplyr::rename(totals_log2FC = log2FoldChange,
                totals_padj = padj)|>
  dplyr::relocate(gene_sym,transcript,gene) |>
  mutate(group = factor(case_when(RPFs_padj < padj_thr & RPFs_log2FC < log2FC_thr & (totals_padj >= padj_thr | is.na(totals_padj)) ~ "RPFs down",
                                  RPFs_padj < padj_thr & RPFs_log2FC > log2FC_thr & (totals_padj >= padj_thr | is.na(totals_padj)) ~ "RPFs up",
                                  (RPFs_padj >= padj_thr | is.na(RPFs_padj)) & totals_padj < padj_thr & totals_log2FC < log2FC_thr ~"Totals down",
                                  (RPFs_padj >= padj_thr | is.na(RPFs_padj)) & totals_padj < padj_thr & totals_log2FC > log2FC_thr ~"Totals up",
                                  RPFs_padj < padj_thr & totals_padj < padj_thr & RPFs_log2FC < log2FC_thr & totals_log2FC < log2FC_thr ~ "both down",
                                  RPFs_padj < padj_thr & totals_padj < padj_thr & RPFs_log2FC > log2FC_thr & totals_log2FC > log2FC_thr ~ "both up")),
         alpha_score = case_when(is.na(group) ~ 0.01,
                                 !(is.na(group)) ~ 0.5)) -> merged_data
summary(merged_data)

## Write CSVs out ------
write_csv(merged_data, file = file.path(tables_dir, paste0(shRNA, "_NTCcorrected_TE_genes_CDSvsTotals.csv")))

### plot TE scatter----
lims <- c(-2, 5.8)

merged_data |>
  arrange(alpha_score == 0.5) |>
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = group, alpha = alpha_score))+
  geom_point()+
  scale_alpha(guide = "none")+
  mytheme+
  scale_colour_manual(values = viridis_colors, na.value = "grey")+
  xlab("Total RNA log2FC")+
  ylab("RPFs log2FC")+
  ggtitle(paste(shRNA, "TE scatter"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr))+
  xlim(lims)+
  ylim(lims)+
  geom_abline(lty=1)+
  geom_hline(yintercept = 0, lty=1)+
  geom_hline(yintercept = 1, lty=2)+
  geom_hline(yintercept = -1, lty=2)+
  geom_vline(xintercept = 0, lty=1)+
  geom_vline(xintercept = 1, lty=2)+
  geom_vline(xintercept = -1, lty=2) -> TE_scatter_plot

pdf(file = paste0(plot_dir, shRNA, "_NTCcorrected_TE_scatter_fixedlimitswGrey.pdf"), width = 7.7, height = 7)
print(TE_scatter_plot)
dev.off()

#write out group sizes
write.table(file = file.path(tables_dir, paste0(shRNA, "_NTCcorrected_TE_scatter_groups.txt")), summary(merged_data$group), col.names = F, quote = F)

## tighter thresholding----
padj_thr = 0.05
log2FC_thr = 0.5

RPFs |>
  dplyr::select(gene, log2FoldChange, padj,gene_sym,transcript) |>
  dplyr::rename(RPFs_log2FC = log2FoldChange,
                RPFs_padj = padj) |>
  inner_join(totals[,c("gene", "log2FoldChange", "padj")], by = "gene") |>
  dplyr::rename(totals_log2FC = log2FoldChange,
                totals_padj = padj)|>
  dplyr::relocate(gene_sym,transcript,gene) |>
  mutate(group = factor(case_when(RPFs_padj < padj_thr & RPFs_log2FC < -log2FC_thr & (totals_padj >= padj_thr | is.na(totals_padj)) ~ "RPFs down",
                                  RPFs_padj < padj_thr & RPFs_log2FC > log2FC_thr & (totals_padj >= padj_thr | is.na(totals_padj)) ~ "RPFs up",
                                  (RPFs_padj >= padj_thr | is.na(RPFs_padj)) & totals_padj < padj_thr & totals_log2FC < -log2FC_thr ~"Totals down",
                                  (RPFs_padj >= padj_thr | is.na(RPFs_padj)) & totals_padj < padj_thr & totals_log2FC > log2FC_thr ~"Totals up",
                                  RPFs_padj < padj_thr & totals_padj < padj_thr & RPFs_log2FC < -log2FC_thr & totals_log2FC < -log2FC_thr ~ "both down",
                                  RPFs_padj < padj_thr & totals_padj < padj_thr & RPFs_log2FC > log2FC_thr & totals_log2FC > log2FC_thr ~ "both up"))) -> merged_data
summary(merged_data)

## Write CSVs out ------
write_csv(merged_data, file = file.path(tables_dir, paste0(shRNA, "_NTCcorrected_tighter_TE_genes_CDSvsTotals.csv")))

### plot TE scatter----
lims <- c(-2, 5.8)

merged_data |>
  mutate(alpha_score = case_when(is.na(group) ~ 0.01,
                                 !(is.na(group)) ~ 0.5)) |>
  arrange(alpha_score == 0.5) |>
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = group, alpha = alpha_score))+
  geom_point()+
  scale_alpha(guide = "none")+
  mytheme+
  scale_colour_manual(values = viridis_colors, na.value = "grey")+
  xlab("Total RNA log2FC")+
  ylab("RPFs log2FC")+
  ggtitle(paste(shRNA, "TE scatter"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr))+
  xlim(lims)+
  ylim(lims)+
  geom_abline(lty=1)+
  geom_hline(yintercept = 0, lty=1)+
  geom_hline(yintercept = 1, lty=2)+
  geom_hline(yintercept = -1, lty=2)+
  geom_vline(xintercept = 0, lty=1)+
  geom_vline(xintercept = 1, lty=2)+
  geom_vline(xintercept = -1, lty=2) -> TE_scatter_plot

pdf(file = paste0(plot_dir, shRNA, "_NTCcorrected_tighter_TE_scatter_fixedlimitswGrey.pdf"), width = 7.7, height = 7)
print(TE_scatter_plot)
dev.off()

#write out group sizes
write.table(file = file.path(tables_dir, paste0(shRNA, "_NTCcorrected_tighter_TE_scatter_groups.txt")), summary(merged_data$group), col.names = F, quote = F)

# 202310___2B4 NTCcorrected --------

#read in common variables
source("common_variables_2B4.R")

## write down the thresholds for significance ----
padj_thr = 0.1
log2FC_thr = 0

## read in DESeq2 output----
totals <- read_csv(file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))
RPFs <- read_csv(file = file.path(tables_dir, paste0("RPFs_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))

# ## plot volcanos----
# RPFs |>
#   dplyr::filter(!(is.na(padj))) |>
#   mutate(sig = factor(case_when(padj < padj_thr ~ "*",
#                                 padj >= padj ~ "NS"))) |>
#   ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sig))+
#   geom_point(size=3,alpha = 0.5)+
#   mytheme+
#   xlab("log2FC")+
#   ylab("-log10(padj)")+
#   ggtitle(paste(shRNA, "RPFs"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr)) +
#   geom_hline(yintercept = 1, linetype='dashed', linewidth=0.5)-> RPFs_volcano
# 
# pdf(file = paste0(plot_dir, "RPFs_", shRNA, "_NTCcorrected_volcano_fixedlimitswGrey.pdf"), width = 7, height = 7)
# print(RPFs_volcano)
# dev.off()
# 
# totals |>
#   dplyr::filter(!(is.na(padj))) |>
#   mutate(sig = factor(case_when(padj < padj_thr ~ "*",
#                                 padj >= padj_thr ~ "NS"))) |>
#   ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sig))+
#   geom_point(size=3, alpha = 0.5)+
#   mytheme+
#   xlab("log2FC")+
#   ylab("-log10(padj)")+
#   ggtitle(paste(shRNA, "Totals"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr))+
#   geom_hline(yintercept = 1, linetype='dashed', linewidth=0.5)-> totals_volcano
# 
# pdf(file = paste0(plot_dir, "Totals_", shRNA, "_NTCcorrected_volcano_fixedlimitswGrey.pdf"), width = 7, height = 7)
# print(totals_volcano)
# dev.off()
# 
# ## plot MAs----
# RPFs |>
#   dplyr::filter(!(is.na(padj))) |>
#   mutate(sig = factor(case_when(padj < padj_thr ~ "*",
#                                 padj >= padj_thr ~ "NS"))) |>
#   ggplot(aes(x = log10(baseMean), y = log2FoldChange, colour = sig))+
#   geom_point(alpha = 0.5)+
#   mytheme+
#   xlab("log10(mean expression)")+
#   ylab("log2FC")+
#   ggtitle(paste(shRNA, "RPFs"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr, "\napeglm shrinking")) -> RPFs_MA
# 
# pdf(file = paste0(plot_dir, paste0("RPFs_", shRNA, "_NTCcorrected_MA_fixedlimitswGrey.pdf")), width = 5, height = 5)
# print(RPFs_MA)
# dev.off()
# 
# totals |>
#   dplyr::filter(!(is.na(padj))) |>
#   mutate(sig = factor(case_when(padj < 0.01 ~ "*",
#                                 padj >= 0.01 ~ "NS"))) |>
#   ggplot(aes(x = log10(baseMean), y = log2FoldChange, colour = sig))+
#   geom_point(alpha = 0.5)+
#   mytheme+
#   xlab("log10(mean expression)")+
#   ylab("log2FC")+
#   ggtitle(paste(shRNA, "Totals"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr, "\napeglm shrinking")) -> totals_MA
# 
# pdf(file = paste0(plot_dir, paste0("Totals_", shRNA, "_NTCcorrected_MA_fixedlimitswGrey.pdf")), width = 5, height = 5)
# print(totals_MA)
# dev.off()

## merge RPF with totals data----
RPFs |>
  dplyr::select(gene, log2FoldChange, padj,gene_sym,transcript) |>
  dplyr::rename(RPFs_log2FC = log2FoldChange,
                RPFs_padj = padj) |>
  inner_join(totals[,c("gene", "log2FoldChange", "padj")], by = "gene") |>
  dplyr::rename(totals_log2FC = log2FoldChange,
                totals_padj = padj)|>
  dplyr::relocate(gene_sym,transcript,gene) |>
  mutate(group = factor(case_when(RPFs_padj < padj_thr & RPFs_log2FC < log2FC_thr & (totals_padj >= padj_thr | is.na(totals_padj)) ~ "RPFs down",
                                  RPFs_padj < padj_thr & RPFs_log2FC > log2FC_thr & (totals_padj >= padj_thr | is.na(totals_padj)) ~ "RPFs up",
                                  (RPFs_padj >= padj_thr | is.na(RPFs_padj)) & totals_padj < padj_thr & totals_log2FC < log2FC_thr ~"Totals down",
                                  (RPFs_padj >= padj_thr | is.na(RPFs_padj)) & totals_padj < padj_thr & totals_log2FC > log2FC_thr ~"Totals up",
                                  RPFs_padj < padj_thr & totals_padj < padj_thr & RPFs_log2FC < log2FC_thr & totals_log2FC < log2FC_thr ~ "both down",
                                  RPFs_padj < padj_thr & totals_padj < padj_thr & RPFs_log2FC > log2FC_thr & totals_log2FC > log2FC_thr ~ "both up")),
         alpha_score = case_when(is.na(group) ~ 0.01,
                                 !(is.na(group)) ~ 0.5)) -> merged_data
summary(merged_data)

## Write CSVs out ------
write_csv(merged_data, file = file.path(tables_dir, paste0(shRNA, "_NTCcorrected_TE_genes_CDSvsTotals.csv")))

### plot TE scatter----
lims <- c(-2, 5.8)

merged_data |>
  arrange(alpha_score == 0.5) |>
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = group, alpha = alpha_score))+
  geom_point()+
  scale_alpha(guide = "none")+
  mytheme+
  scale_colour_manual(values = viridis_colors, na.value = "grey")+
  xlab("Total RNA log2FC")+
  ylab("RPFs log2FC")+
  ggtitle(paste(shRNA, "TE scatter"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr))+
  xlim(lims)+
  ylim(lims)+
  geom_abline(lty=1)+
  geom_hline(yintercept = 0, lty=1)+
  geom_hline(yintercept = 1, lty=2)+
  geom_hline(yintercept = -1, lty=2)+
  geom_vline(xintercept = 0, lty=1)+
  geom_vline(xintercept = 1, lty=2)+
  geom_vline(xintercept = -1, lty=2) -> TE_scatter_plot

pdf(file = paste0(plot_dir, shRNA, "_NTCcorrected_TE_scatter_fixedlimitswGrey.pdf"), width = 7.7, height = 7)
print(TE_scatter_plot)
dev.off()

#write out group sizes
write.table(file = file.path(tables_dir, paste0(shRNA, "_NTCcorrected_TE_scatter_groups.txt")), summary(merged_data$group), col.names = F, quote = F)

## tighter thresholding----
padj_thr = 0.05
log2FC_thr = 0.5

RPFs |>
  dplyr::select(gene, log2FoldChange, padj,gene_sym,transcript) |>
  dplyr::rename(RPFs_log2FC = log2FoldChange,
                RPFs_padj = padj) |>
  inner_join(totals[,c("gene", "log2FoldChange", "padj")], by = "gene") |>
  dplyr::rename(totals_log2FC = log2FoldChange,
                totals_padj = padj)|>
  dplyr::relocate(gene_sym,transcript,gene) |>
  mutate(group = factor(case_when(RPFs_padj < padj_thr & RPFs_log2FC < -log2FC_thr & (totals_padj >= padj_thr | is.na(totals_padj)) ~ "RPFs down",
                                  RPFs_padj < padj_thr & RPFs_log2FC > log2FC_thr & (totals_padj >= padj_thr | is.na(totals_padj)) ~ "RPFs up",
                                  (RPFs_padj >= padj_thr | is.na(RPFs_padj)) & totals_padj < padj_thr & totals_log2FC < -log2FC_thr ~"Totals down",
                                  (RPFs_padj >= padj_thr | is.na(RPFs_padj)) & totals_padj < padj_thr & totals_log2FC > log2FC_thr ~"Totals up",
                                  RPFs_padj < padj_thr & totals_padj < padj_thr & RPFs_log2FC < -log2FC_thr & totals_log2FC < -log2FC_thr ~ "both down",
                                  RPFs_padj < padj_thr & totals_padj < padj_thr & RPFs_log2FC > log2FC_thr & totals_log2FC > log2FC_thr ~ "both up")),
         alpha_score = case_when(is.na(group) ~ 0.01,
                                 !(is.na(group)) ~ 0.5)) -> merged_data
summary(merged_data)

## Write CSVs out ------
write_csv(merged_data, file = file.path(tables_dir, paste0(shRNA, "_NTCcorrected_tighter_TE_genes_CDSvsTotals.csv")))

### plot TE scatter----
lims <- c(-2, 5.8)

merged_data |>
  arrange(alpha_score == 0.5) |>
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = group, alpha = alpha_score))+
  geom_point()+
  scale_alpha(guide = "none")+
  mytheme+
  scale_colour_manual(values = viridis_colors, na.value = "grey")+
  xlab("Total RNA log2FC")+
  ylab("RPFs log2FC")+
  ggtitle(paste(shRNA, "TE scatter"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr))+
  xlim(lims)+
  ylim(lims)+
  geom_abline(lty=1)+
  geom_hline(yintercept = 0, lty=1)+
  geom_hline(yintercept = 1, lty=2)+
  geom_hline(yintercept = -1, lty=2)+
  geom_vline(xintercept = 0, lty=1)+
  geom_vline(xintercept = 1, lty=2)+
  geom_vline(xintercept = -1, lty=2) -> TE_scatter_plot

pdf(file = paste0(plot_dir, shRNA, "_NTCcorrected_tighter_TE_scatter_fixedlimitswGrey.pdf"), width = 7.7, height = 7)
print(TE_scatter_plot)
dev.off()

#write out group sizes
write.table(file = file.path(tables_dir, paste0(shRNA, "_NTCcorrected_tighter_TE_scatter_groups.txt")), summary(merged_data$group), col.names = F, quote = F)

# 202310___2B5 NTCcorrected --------

#read in common variables
source("common_variables_2B5.R")

## write down the thresholds for significance ----
padj_thr = 0.1
log2FC_thr = 0

## read in DESeq2 output----
totals <- read_csv(file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))
RPFs <- read_csv(file = file.path(tables_dir, paste0("RPFs_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))

# ## plot volcanos----
# RPFs |>
#   dplyr::filter(!(is.na(padj))) |>
#   mutate(sig = factor(case_when(padj < padj_thr ~ "*",
#                                 padj >= padj ~ "NS"))) |>
#   ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sig))+
#   geom_point(size=3,alpha = 0.5)+
#   mytheme+
#   xlab("log2FC")+
#   ylab("-log10(padj)")+
#   ggtitle(paste(shRNA, "RPFs"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr)) +
#   geom_hline(yintercept = 1, linetype='dashed', linewidth=0.5)-> RPFs_volcano
# 
# pdf(file = paste0(plot_dir, "RPFs_", shRNA, "_NTCcorrected_volcano_fixedlimitswGrey.pdf"), width = 7, height = 7)
# print(RPFs_volcano)
# dev.off()
# 
# totals |>
#   dplyr::filter(!(is.na(padj))) |>
#   mutate(sig = factor(case_when(padj < padj_thr ~ "*",
#                                 padj >= padj_thr ~ "NS"))) |>
#   ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = sig))+
#   geom_point(size=3, alpha = 0.5)+
#   mytheme+
#   xlab("log2FC")+
#   ylab("-log10(padj)")+
#   ggtitle(paste(shRNA, "Totals"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr))+
#   geom_hline(yintercept = 1, linetype='dashed', linewidth=0.5)-> totals_volcano
# 
# pdf(file = paste0(plot_dir, "Totals_", shRNA, "_NTCcorrected_volcano_fixedlimitswGrey.pdf"), width = 7, height = 7)
# print(totals_volcano)
# dev.off()
# 
# ## plot MAs----
# RPFs |>
#   dplyr::filter(!(is.na(padj))) |>
#   mutate(sig = factor(case_when(padj < padj_thr ~ "*",
#                                 padj >= padj_thr ~ "NS"))) |>
#   ggplot(aes(x = log10(baseMean), y = log2FoldChange, colour = sig))+
#   geom_point(alpha = 0.5)+
#   mytheme+
#   xlab("log10(mean expression)")+
#   ylab("log2FC")+
#   ggtitle(paste(shRNA, "RPFs"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr, "\napeglm shrinking")) -> RPFs_MA
# 
# pdf(file = paste0(plot_dir, paste0("RPFs_", shRNA, "_NTCcorrected_MA_fixedlimitswGrey.pdf")), width = 5, height = 5)
# print(RPFs_MA)
# dev.off()
# 
# totals |>
#   dplyr::filter(!(is.na(padj))) |>
#   mutate(sig = factor(case_when(padj < 0.01 ~ "*",
#                                 padj >= 0.01 ~ "NS"))) |>
#   ggplot(aes(x = log10(baseMean), y = log2FoldChange, colour = sig))+
#   geom_point(alpha = 0.5)+
#   mytheme+
#   xlab("log10(mean expression)")+
#   ylab("log2FC")+
#   ggtitle(paste(shRNA, "Totals"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr, "\napeglm shrinking")) -> totals_MA
# 
# pdf(file = paste0(plot_dir, paste0("Totals_", shRNA, "_NTCcorrected_MA_fixedlimitswGrey.pdf")), width = 5, height = 5)
# print(totals_MA)
# dev.off()

## merge RPF with totals data----
RPFs |>
  dplyr::select(gene, log2FoldChange, padj,gene_sym,transcript) |>
  dplyr::rename(RPFs_log2FC = log2FoldChange,
                RPFs_padj = padj) |>
  inner_join(totals[,c("gene", "log2FoldChange", "padj")], by = "gene") |>
  dplyr::rename(totals_log2FC = log2FoldChange,
                totals_padj = padj)|>
  dplyr::relocate(gene_sym,transcript,gene) |>
  mutate(group = factor(case_when(RPFs_padj < padj_thr & RPFs_log2FC < log2FC_thr & (totals_padj >= padj_thr | is.na(totals_padj)) ~ "RPFs down",
                                  RPFs_padj < padj_thr & RPFs_log2FC > log2FC_thr & (totals_padj >= padj_thr | is.na(totals_padj)) ~ "RPFs up",
                                  (RPFs_padj >= padj_thr | is.na(RPFs_padj)) & totals_padj < padj_thr & totals_log2FC < log2FC_thr ~"Totals down",
                                  (RPFs_padj >= padj_thr | is.na(RPFs_padj)) & totals_padj < padj_thr & totals_log2FC > log2FC_thr ~"Totals up",
                                  RPFs_padj < padj_thr & totals_padj < padj_thr & RPFs_log2FC < log2FC_thr & totals_log2FC < log2FC_thr ~ "both down",
                                  RPFs_padj < padj_thr & totals_padj < padj_thr & RPFs_log2FC > log2FC_thr & totals_log2FC > log2FC_thr ~ "both up")),
         alpha_score = case_when(is.na(group) ~ 0.01,
                                 !(is.na(group)) ~ 0.5)) -> merged_data
summary(merged_data)

## Write CSVs out ------
write_csv(merged_data, file = file.path(tables_dir, paste0(shRNA, "_NTCcorrected_TE_genes_CDSvsTotals.csv")))

### plot TE scatter----
lims <- c(-2, 5.8)

merged_data |>
  arrange(alpha_score == 0.5) |>
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = group, alpha = alpha_score))+
  geom_point()+
  scale_alpha(guide = "none")+
  mytheme+
  scale_colour_manual(values = viridis_colors, na.value = "grey")+
  xlab("Total RNA log2FC")+
  ylab("RPFs log2FC")+
  ggtitle(paste(shRNA, "TE scatter"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr))+
  xlim(lims)+
  ylim(lims)+
  geom_abline(lty=1)+
  geom_hline(yintercept = 0, lty=1)+
  geom_hline(yintercept = 1, lty=2)+
  geom_hline(yintercept = -1, lty=2)+
  geom_vline(xintercept = 0, lty=1)+
  geom_vline(xintercept = 1, lty=2)+
  geom_vline(xintercept = -1, lty=2) -> TE_scatter_plot

pdf(file = paste0(plot_dir, shRNA, "_NTCcorrected_TE_scatter_fixedlimitswGrey.pdf"), width = 7.7, height = 7)
print(TE_scatter_plot)
dev.off()

#write out group sizes
write.table(file = file.path(tables_dir, paste0(shRNA, "_NTCcorrected_TE_scatter_groups.txt")), summary(merged_data$group), col.names = F, quote = F)

## tighter thresholding----
padj_thr = 0.05
log2FC_thr = 0.5

RPFs |>
  dplyr::select(gene, log2FoldChange, padj,gene_sym,transcript) |>
  dplyr::rename(RPFs_log2FC = log2FoldChange,
                RPFs_padj = padj) |>
  inner_join(totals[,c("gene", "log2FoldChange", "padj")], by = "gene") |>
  dplyr::rename(totals_log2FC = log2FoldChange,
                totals_padj = padj)|>
  dplyr::relocate(gene_sym,transcript,gene) |>
  mutate(group = factor(case_when(RPFs_padj < padj_thr & RPFs_log2FC < -log2FC_thr & (totals_padj >= padj_thr | is.na(totals_padj)) ~ "RPFs down",
                                  RPFs_padj < padj_thr & RPFs_log2FC > log2FC_thr & (totals_padj >= padj_thr | is.na(totals_padj)) ~ "RPFs up",
                                  (RPFs_padj >= padj_thr | is.na(RPFs_padj)) & totals_padj < padj_thr & totals_log2FC < -log2FC_thr ~"Totals down",
                                  (RPFs_padj >= padj_thr | is.na(RPFs_padj)) & totals_padj < padj_thr & totals_log2FC > log2FC_thr ~"Totals up",
                                  RPFs_padj < padj_thr & totals_padj < padj_thr & RPFs_log2FC < -log2FC_thr & totals_log2FC < -log2FC_thr ~ "both down",
                                  RPFs_padj < padj_thr & totals_padj < padj_thr & RPFs_log2FC > log2FC_thr & totals_log2FC > log2FC_thr ~ "both up")),
         alpha_score = case_when(is.na(group) ~ 0.01,
                                 !(is.na(group)) ~ 0.5)) -> merged_data
summary(merged_data)

## Write CSVs out ------
write_csv(merged_data, file = file.path(tables_dir, paste0(shRNA, "_NTCcorrected_tighter_TE_genes_CDSvsTotals.csv")))

### plot TE scatter----
lims <- c(-2, 5.8)

merged_data |>
  arrange(alpha_score == 0.5) |>
  ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = group, alpha = alpha_score))+
  geom_point()+
  scale_alpha(guide = "none")+
  mytheme+
  scale_colour_manual(values = viridis_colors, na.value = "grey")+
  xlab("Total RNA log2FC")+
  ylab("RPFs log2FC")+
  ggtitle(paste(shRNA, "TE scatter"), subtitle = paste("p-value thr", padj_thr,"\nlog2FC thr +/-", log2FC_thr))+
  xlim(lims)+
  ylim(lims)+
  geom_abline(lty=1)+
  geom_hline(yintercept = 0, lty=1)+
  geom_hline(yintercept = 1, lty=2)+
  geom_hline(yintercept = -1, lty=2)+
  geom_vline(xintercept = 0, lty=1)+
  geom_vline(xintercept = 1, lty=2)+
  geom_vline(xintercept = -1, lty=2) -> TE_scatter_plot

pdf(file = paste0(plot_dir, shRNA, "_NTCcorrected_tighter_TE_scatter_fixedlimitswGrey.pdf"), width = 7.7, height = 7)
print(TE_scatter_plot)
dev.off()

#write out group sizes
write.table(file = file.path(tables_dir, paste0(shRNA, "_NTCcorrected_tighter_TE_scatter_groups.txt")), summary(merged_data$group), col.names = F, quote = F)

