# get all my data in a single table - merging all normalised and all log2FC######
# as the venn diagram had less numbers in some categories, I will try full_joins instead of left_joins
# these data then goes into tables_dir and plot_dirs which are for publication
# print in PDF for high resolution and vector graphics

library(tidyverse)
library(tximport)
library(vsn)
library(viridis)
library(conflicted)
library(UpSetR)
library(VennDiagram)
library(vennplot)
library(ggplot2)

## themes----
mytheme <- theme_minimal()+
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))


## read in common variables & files ----
source("common_variables.R")
shRNAs <- c('2B1','2B4','2B5')

# From Mac parent directory
# parent_dir <- '/Volumes/data-1/CGIACOME/AAA_eIF2B/AAA_Collaborative_stuff/Wurzburg_collab/AAA_Riboseq/Analysis/Ribo-seq-allreps'
# plot_dir <- paste0(parent_dir,"/plots/PDF_for_paper/")
# tables_dir <- paste0(parent_dir,"/Analysis/tables_for_paper")

## DEseq2 results for totals, RPFs, and 5UTRs ------
totals_2B1 <- read_csv(file = file.path(tables_dir, "Totals_2B1_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv"))
RPFs_2B1 <- read_csv(file = file.path(tables_dir, "RPFs_2B1_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv"))

totals_2B4 <- read_csv(file = file.path(tables_dir, "Totals_2B4_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv"))
RPFs_2B4 <- read_csv(file = file.path(tables_dir, "RPFs_2B4_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv"))

## Rld normalised counts ------
totals_nc_2B1 <- read_csv(file = file.path(tables_dir, "Totals_2B1_NTCcorrected_normalised_counts.csv"))
RPFs_nc_2B1 <- read_csv(file = file.path(tables_dir, "RPFs_2B1_NTCcorrected_normalised_counts.csv"))

totals_nc_2B4 <- read_csv(file = file.path(tables_dir, "Totals_2B4_NTCcorrected_normalised_counts.csv"))
RPFs_nc_2B4 <- read_csv(file = file.path(tables_dir, "RPFs_2B4_NTCcorrected_normalised_counts.csv"))

## TE data ------
TE_NTCcor_2B1 <- read_csv(file = file.path(tables_dir, "TE_2B1_NTCcorrected_4fgsea.csv"))
TE_NTCcor_2B4 <- read_csv(file = file.path(tables_dir, "TE_2B4_NTCcorrected_4fgsea.csv"))

## TE data including groups
groups_2B1 <- read_csv(file = file.path(tables_dir, "2B1_NTCcorrected_TE_genes_CDSvsTotals.csv"))
groups_2B4 <- read_csv(file = file.path(tables_dir, "2B4_NTCcorrected_TE_genes_CDSvsTotals.csv"))

# merge all normalised counts-----
## might be useful if one wants to plot genes from different conditions

RPFs_nc_2B1 |>
  inner_join(RPFs_nc_2B4[,c(2, 12:19)], by = 'gene') |>
  relocate(transcript,gene,gene_sym)-> RPFs_nc_all
  
totals_nc_2B1 |>
  inner_join(totals_nc_2B4[,c(2, 12:19)], by = 'gene') |>
  relocate(transcript,gene,gene_sym)-> totals_nc_all

## Write CSVs out ------
write_csv(RPFs_nc_all, file = file.path(tables_dir, "RPFs_2B1_2B4_NTCcorrected_normalised_counts.csv"))
write_csv(totals_nc_all, file = file.path(tables_dir, "Totals_2B1_2B4_NTCcorrected_normalised_counts.csv"))

# merge all basemeans,log2FC,adjpvals -----
RPFs_2B1 |>
  dplyr::select(transcript, gene, gene_sym, baseMean, log2FoldChange, padj) |>
  full_join(RPFs_2B4[,c('gene', 'baseMean', 'log2FoldChange', 'padj')], by = 'gene', suffix = c('_2B1','_2B4')) |>
  relocate(gene)-> RPFs_all

totals_2B1 |>
  dplyr::select(transcript, gene, gene_sym, baseMean, log2FoldChange, padj) |>
  full_join(totals_2B4[,c('gene', 'baseMean', 'log2FoldChange', 'padj')], by = 'gene', suffix = c('_2B1','_2B4')) |>
  relocate(gene)-> totals_all

## Write CSVs out ------
write_csv(RPFs_all, file = file.path(tables_dir, "RPFs_2B1_2B4_NTCcorrected_DEGs.csv"))
write_csv(totals_all, file = file.path(tables_dir, "Totals_2B1_2B4_NTCcorrected_DEGs.csv"))

# baseline TE analysis -----
# merge normalised counts for RPFs and totals within same shRNA KD
# this is to calculate baseline TE in uninduced conditions & plot just to look how it looks

RPFs_nc_2B1 |>
  inner_join(totals_nc_2B1[,c(2, 12:19)], by = 'gene') |>
  dplyr::select(-c(4:11)) |>                               # remove NTC data
  dplyr::select(-c(8:11, 16:19)) |>                        # remove induced data
  relocate(transcript,gene,gene_sym) |>
  rowwise() |>
  mutate(ave_RPFs = base::mean(c_across(c(4:7)), na.rm = TRUE)) |>
  mutate(ave_tots = base::mean(c_across(c(8:11)), na.rm = TRUE)) |>
  mutate(Baseline_TE = ave_RPFs-ave_tots) -> B1_nc_uninduced

RPFs_nc_2B4 |>
  inner_join(totals_nc_2B4[,c(2, 12:19)], by = 'gene') |>
  dplyr::select(-c(4:11)) |>                               # remove NTC data
  dplyr::select(-c(8:11, 16:19)) |>                        # remove induced data
  relocate(transcript,gene,gene_sym) |>
  rowwise() |>
  mutate(ave_RPFs = base::mean(c_across(c(4:7)), na.rm = TRUE)) |>
  mutate(ave_tots = base::mean(c_across(c(8:11)), na.rm = TRUE)) |>
  mutate(Baseline_TE = ave_RPFs-ave_tots) -> B4_nc_uninduced

B1_nc_uninduced |>
  ggplot(aes(x = ave_tots, y = ave_RPFs))+
  stat_binhex(bins = 40)+
  mytheme+
  scale_fill_viridis(discrete=FALSE,option='G', begin = 0, end = 0.8)+
  xlab("normalised total RNA (log2)")+
  ylab("normalised RPFs (log2)")+
  ggtitle("Baseline Translation 2B1") -> TE_scatter_plot

pdf(file = paste0(plot_dir, "2B1_baseline_translation.pdf"), width = 5.5, height = 5)
print(TE_scatter_plot)
dev.off()

B4_nc_uninduced |>
  ggplot(aes(x = ave_tots, y = ave_RPFs))+
  stat_binhex(bins = 40)+
  mytheme+
  scale_fill_viridis(discrete=FALSE,option='G', begin = 0, end = 0.8)+
  xlab("normalised total RNA (log2)")+
  ylab("normalised RPFs (log2)")+
  ggtitle("Baseline Translation 2B4") -> TE_scatter_plot

pdf(file = paste0(plot_dir, "2B4_baseline_translation.pdf"), width = 5.5, height = 5)
print(TE_scatter_plot)
dev.off()

# redo also by delta TE
B1_nc_uninduced |>
  inner_join(TE_NTCcor_2B1[,c(1,5,8)]) |>
  mutate(alpha_score = case_when(RPFs_padj <= 0.1 ~ 1, RPFs_padj >0.1 ~ 0.1)) |>
  ggplot(aes(x = ave_tots, y = ave_RPFs, color = TE, alpha = alpha_score))+
  geom_point()+
  scale_alpha(guide = "none")+
  mytheme+
  scale_color_viridis(discrete=FALSE,option='G', begin = 0, end = 0.8)+
  xlab("normalised total RNA (log2)")+
  ylab("normalised RPFs (log2)")+
  labs(color = "log2\nTE")+
  ggtitle("Baseline Translation 2B1") -> TE_scatter_plot

pdf(file = paste0(plot_dir, "2B1_baseline_translation_TEcolor.pdf"), width = 5.5, height = 5)
print(TE_scatter_plot)
dev.off()

B4_nc_uninduced |>
  inner_join(TE_NTCcor_2B4[,c(1,5,8)]) |>
  mutate(alpha_score = case_when(RPFs_padj <= 0.1 ~ 1, RPFs_padj >0.1 ~ 0.1)) |>
  ggplot(aes(x = ave_tots, y = ave_RPFs, color = TE, alpha = alpha_score))+
  geom_point()+
  scale_alpha(guide = "none")+
  mytheme+
  scale_color_viridis(discrete=FALSE,option='G', begin = 0, end = 0.8)+
  xlab("normalised total RNA (log2)")+
  ylab("normalised RPFs (log2)")+
  labs(color = "log2\nTE")+
  ggtitle("Baseline Translation 2B4") -> TE_scatter_plot

pdf(file = paste0(plot_dir, "2B4_baseline_translation_TEcolor.pdf"), width = 5.5, height = 50)
print(TE_scatter_plot)
dev.off()

# redo by coloring subgroups from DE analyses
viridis_colors <- c("#000004","#2c115f", "#721f81", "#b73779","#f1605d", "#feb078")

B1_nc_uninduced |>
  inner_join(groups_2B1[,c(2,8,9)]) |>
  ggplot(aes(x = ave_tots, y = ave_RPFs, color = group, alpha = alpha_score))+
  scale_alpha(guide = "none")+
  geom_point()+
  mytheme+
  scale_colour_manual(values = viridis_colors, na.value = NA)+
  xlab("normalised total RNA (log2)")+
  ylab("normalised RPFs (log2)")+
  ggtitle("Baseline Translation 2B1") -> TE_scatter_plot

pdf(file = paste0(plot_dir, "2B1_baseline_translation_group_color.pdf"), width = 6, height = 5)
print(TE_scatter_plot)
dev.off()

B4_nc_uninduced |>
  inner_join(groups_2B4[,c(2,8,9)]) |>
  ggplot(aes(x = ave_tots, y = ave_RPFs, color = group, alpha = alpha_score))+
  scale_alpha(guide = "none")+
  geom_point()+
  mytheme+
  scale_colour_manual(values = viridis_colors, na.value = NA)+
  xlab("normalised total RNA (log2)")+
  ylab("normalised RPFs (log2)")+
  ggtitle("Baseline Translation 2B4") -> TE_scatter_plot

pdf(file = paste0(plot_dir, "2B4_baseline_translation_group_color.pdf"), width = 6, height = 5)
print(TE_scatter_plot)
dev.off()

# plote baseline translation vs delta TE correlation w/ R
myR <- function(x) {
  r <- base::round(base::as.numeric(x$estimate), digits = 2)
  p <- base::as.numeric(x$p.value)
  if (p > 0.001) {
    rounded_p <- base::round(p, digits = 3)
    p_label <- paste('P =', rounded_p)
  }else{
    if (p<2.2E-16) {
      p_label <- 'P < 2.2e-16'
    } else {
      rounded_p <- formatC(p, format = "e", digits = 2)
      p_label <- paste('P =', rounded_p)
    }
  }
  return(paste0('r = ', r, '\n', p_label))
}

B1_nc_uninduced |>
  inner_join(TE_NTCcor_2B1[,c(1,5,8)]) -> TE_BL_vs_Delta

r <- cor.test(x = TE_BL_vs_Delta$Baseline_TE, y = TE_BL_vs_Delta$TE)
r_label <- myR(x = r)

TE_BL_vs_Delta |>
  ggplot(aes(x = Baseline_TE, y = TE)) +
  stat_binhex(bins = 40)+
  scale_fill_viridis()+
  ylim(c(-2.5,2.5))+
  mytheme+
  xlab("Baseline TE")+
  ylab("TE log2FC upon KD")+
  ggtitle("2B1 Baseline TE vs Delta", subtitle = r_label) -> scatter

pdf(file = paste0(plot_dir, "2B1_baseline_translation_deltaTE.pdf"), width = 5.5, height = 5)
print(scatter)
dev.off()

B4_nc_uninduced |>
  inner_join(TE_NTCcor_2B4[,c(1,5,8)]) -> TE_BL_vs_Delta

r <- cor.test(x = TE_BL_vs_Delta$Baseline_TE, y = TE_BL_vs_Delta$TE)
r_label <- myR(x = r)

TE_BL_vs_Delta |>
  ggplot(aes(x = Baseline_TE, y = TE)) +
  stat_binhex(bins = 40)+
  scale_fill_viridis()+
  ylim(c(-2.5,2.5))+
  mytheme+
  xlab("Baseline TE")+
  ylab("TE log2FC upon KD")+
  ggtitle("2B4 Baseline TE vs Delta", subtitle = r_label) -> scatter

pdf(file = paste0(plot_dir, "2B4_baseline_translation_deltaTE.pdf"), width = 5.5, height = 5)
print(scatter)
dev.off()

## Write out the baseline TEs -----
write_csv(B1_nc_uninduced, file = file.path(tables_dir, "BaselineTE_2B1.csv"))
write_csv(B4_nc_uninduced, file = file.path(tables_dir, "BaselineTE_2B4.csv"))

# shRNA-specific and shared effects ----
# Merge TE data
# check same transcript?
TE_NTCcor_2B1 |>
  inner_join(TE_NTCcor_2B4, by = 'gene', suffix = c('_2B1','_2B4')) |>
  mutate(same_transcript = factor(case_when(transcript_2B1 == transcript_2B4 ~ "yes",
                                            transcript_2B1 != transcript_2B4 ~ "no"))) |>
  relocate(same_transcript)-> RPFs_nc_all

summary(RPFs_nc_all)

# compare the RPF only and both categories from DE analyses ------
groups_2B1 |>
  dplyr::select(-alpha_score) |>
  dplyr::mutate(TE = RPFs_log2FC - totals_log2FC) -> TE_2B1

groups_2B4 |>
  dplyr::select(-alpha_score) |>
  dplyr::mutate(TE = RPFs_log2FC - totals_log2FC) -> TE_2B4

# TE comparison for CDS vs Totals ------
# tested first to categorise all combos, but some groups are empty, so here below are the ones that have at least 1 gene inside

TE_2B1 |>
  full_join(TE_2B4[,c(1,4:9)], by = 'gene_sym', suffix = c('_2B1','_2B4')) |>
  mutate(cumulative_group = factor(case_when(group_2B1 == "both up" & group_2B4 == "both up" ~"both up 2B1 & 2B4",
                                             group_2B1 == "both down" & group_2B4 == "both down" ~"both down 2B1 & 2B4",
                                             group_2B1 == "RPFs up" & group_2B4 == "RPFs up" ~"RPFs up 2B1 & 2B4",
                                             group_2B1 == "RPFs down" & group_2B4 == "RPFs down" ~"RPFs down 2B1 & 2B4",
                                             group_2B1 == "both up" & (is.na(group_2B4) | group_2B4 == "Totals up" | group_2B4 == "Totals down") ~"both up 2B1 only",
                                             group_2B1 == "both down" & (is.na(group_2B4) | group_2B4 == "Totals up" | group_2B4 == "Totals down") ~"both down 2B1 only",
                                             (is.na(group_2B1) | group_2B1 == "Totals up" | group_2B1 == "Totals down") & group_2B4 == "both up" ~"both up 2B4 only",
                                             (is.na(group_2B1) | group_2B1 == "Totals up" | group_2B1 == "Totals down") & group_2B4 == "both down" ~"both down 2B4 only",
                                             group_2B1 == "RPFs up" & (is.na(group_2B4) | group_2B4 == "Totals up" | group_2B4 == "Totals down") ~"RPFs up 2B1 only",
                                             group_2B1 == "RPFs down" & (is.na(group_2B4) | group_2B4 == "Totals up" | group_2B4 == "Totals down") ~"RPFs down 2B1 only",
                                             (is.na(group_2B1) | group_2B1 == "Totals up" | group_2B1 == "Totals down") & group_2B4 == "RPFs up" ~"RPFs up 2B4 only",
                                             (is.na(group_2B1) | group_2B1 == "Totals up" | group_2B1 == "Totals down") & group_2B4 == "RPFs down" ~"RPFs down 2B4 only",
                                             group_2B1 == "both up" & group_2B4 == "RPFs up" ~"both up 2B1 & RPFs up 2B4",
                                             group_2B1 == "RPFs up" & group_2B4 == "both up" ~"RPFs up 2B1, both up 2B4",
                                             group_2B1 == "both down" & group_2B4 == "RPFs down" ~"both down 2B1, RPFs down 2B4",
                                             group_2B1 == "RPFs down" & group_2B4 == "both down" ~"RPFs down 2B1, both down 2B4"))) |>
  relocate(c(cumulative_group, group_2B1, group_2B4, TE_2B1, TE_2B4), .after = "gene") -> merged_TE

summary(merged_TE)

## write out size of groups ----
write.table(file = file.path(tables_dir, paste0('Merged_2B1_2B4', "_NTCcorrected_cumulative_groups.txt")), summary(merged_TE$cumulative_group), col.names = F, quote = F)

## plot scattes for the 2B1 TE vs 2B4 TE ------
# caluclate axis limits
B1_lims <- max(abs(merged_TE$TE_2B1))
B4_lims <- max(abs(merged_TE$TE_2B4))
lim <- max(c(B1_lims, B4_lims))
lims <- c(-lim, lim)

merged_TE |>
  ggplot(aes(x = TE_2B1, y = TE_2B4))+
  geom_point()+
  mytheme+
  scale_color_viridis(discrete=TRUE,option='G', begin = 0, end = 0.8)+
  xlab("TE EIF2B1 KD \n(RPFs log2FC - Total log2FC)")+
  ylab("TE EIF2B4 KD \n(RPFs log2FC - Total log2FC)")+
  ggtitle("TE 2B1 vs 2B4")+
  xlim(lims)+
  ylim(lims)+
  geom_abline(lty=1)+
  geom_hline(yintercept = 0, lty=1)+
  geom_hline(yintercept = 1, lty=2)+
  geom_hline(yintercept = -1, lty=2)+
  geom_vline(xintercept = 0, lty=1)+
  geom_vline(xintercept = 1, lty=2)+
  geom_vline(xintercept = -1, lty=2) -> TE_scatter_plot

pdf(file = paste0(plot_dir, "2B1_vs_2B4_TE.pdf"), width = 5.5, height = 5)
print(TE_scatter_plot)
dev.off()

merged_TE |>
  mutate(group = case_when((cumulative_group == "RPFs up 2B1 only" | cumulative_group == "both up 2B1 only") ~"2B1 up \nRPFs + both",
                           (cumulative_group == "RPFs down 2B1 only" | cumulative_group == "both down 2B1 only") ~"2B1 down \nRPFs + both",
                           (cumulative_group == "RPFs up 2B4 only" | cumulative_group == "both up 2B4 only") ~"2B4 up \nRPFs + both",
                           (cumulative_group == "RPFs down 2B4 only" | cumulative_group == "both down 2B4 only") ~"2B4 down \nRPFs + both")) |>
  ggplot(aes(x = TE_2B1, y = TE_2B4, colour = group))+
  geom_point()+
  mytheme+
  scale_color_viridis(discrete=TRUE,option='G', begin = 0, end = 0.8)+
  xlab("TE EIF2B1 KD \n(RPFs log2FC - Total log2FC)")+
  ylab("TE EIF2B4 KD \n(RPFs log2FC - Total log2FC)")+
  ggtitle("TE 2B1 vs 2B4")+
  xlim(lims)+
  ylim(lims)+
  geom_abline(lty=1)+
  geom_hline(yintercept = 0, lty=1)+
  geom_hline(yintercept = 1, lty=2)+
  geom_hline(yintercept = -1, lty=2)+
  geom_vline(xintercept = 0, lty=1)+
  geom_vline(xintercept = 1, lty=2)+
  geom_vline(xintercept = -1, lty=2) -> TE_scatter_plot

pdf(file = paste0(plot_dir, "2B1_vs_2B4_TE_groups_exclusive.pdf"), width = 5, height = 4)
print(TE_scatter_plot)
dev.off()

merged_TE |>
  mutate(group = case_when(cumulative_group == "RPFs up 2B1 only" ~"2B1 up",
                           cumulative_group == "RPFs down 2B1 only" ~"2B1 down",
                           cumulative_group == "RPFs up 2B4 only" ~"2B4 up",
                           cumulative_group == "RPFs down 2B4 only" ~"2B4 down")) |>
  ggplot(aes(x = TE_2B1, y = TE_2B4, colour = group))+
  geom_point()+
  mytheme+
  scale_color_viridis(discrete=TRUE,option='G', begin = 0, end = 0.8)+
  xlab("TE EIF2B1 KD \n(RPFs log2FC - Total log2FC)")+
  ylab("TE EIF2B4 KD \n(RPFs log2FC - Total log2FC)")+
  ggtitle("TE 2B1 vs 2B4")+
  xlim(lims)+
  ylim(lims)+
  geom_abline(lty=1)+
  geom_hline(yintercept = 0, lty=1)+
  geom_hline(yintercept = 1, lty=2)+
  geom_hline(yintercept = -1, lty=2)+
  geom_vline(xintercept = 0, lty=1)+
  geom_vline(xintercept = 1, lty=2)+
  geom_vline(xintercept = -1, lty=2) -> TE_scatter_plot

pdf(file = paste0(plot_dir, "2B1_vs_2B4_TE_groups_RPFs_exclusive.pdf"), width = 5, height = 4)
print(TE_scatter_plot)
dev.off()

merged_TE |>
  mutate(group = case_when(cumulative_group == "both up 2B1 only" ~"2B1 up",
                           cumulative_group == "both down 2B1 only" ~"2B1 down",
                           cumulative_group == "both up 2B4 only" ~"2B4 up",
                           cumulative_group == "both down 2B4 only" ~"2B4 down")) |>
  ggplot(aes(x = TE_2B1, y = TE_2B4, colour = group))+
  geom_point()+
  mytheme+
  scale_color_viridis(discrete=TRUE,option='G', begin = 0, end = 0.8)+
  xlab("TE EIF2B1 KD \n(RPFs log2FC - Total log2FC)")+
  ylab("TE EIF2B4 KD \n(RPFs log2FC - Total log2FC)")+
  ggtitle("TE 2B1 vs 2B4")+
  xlim(lims)+
  ylim(lims)+
  geom_abline(lty=1)+
  geom_hline(yintercept = 0, lty=1)+
  geom_hline(yintercept = 1, lty=2)+
  geom_hline(yintercept = -1, lty=2)+
  geom_vline(xintercept = 0, lty=1)+
  geom_vline(xintercept = 1, lty=2)+
  geom_vline(xintercept = -1, lty=2) -> TE_scatter_plot

pdf(file = paste0(plot_dir, "2B1_vs_2B4_TE_groups_both_exclusive.pdf"), width = 5, height = 4)
print(TE_scatter_plot)
dev.off()

#  VENN diagrams ------
# # get group sizes
# merged_TE |>
#   count(cumulative_group) |>
#   mutate(cumulative_group = fct_na_value_to_level(cumulative_group, level = "ns")) |>
#   column_to_rownames(var = "cumulative_group")-> cum_group_counts
# 
# cum_group_counts

# get values analysis only the distribution of both only (RPF + totals) going up or down
TE_2B1 |>
  full_join(TE_2B4[,c(1,4:9)], by = 'gene_sym', suffix = c('_2B1','_2B4')) |>
  mutate(group_2B1 = fct_na_value_to_level(group_2B1, level = "ns"),
         group_2B4 = fct_na_value_to_level(group_2B4, level = "ns")) |>
  mutate(cumulative_group = factor(case_when(group_2B1 == "both up" & group_2B4 == "both up" ~"both up 2B1 & 2B4",
                                             group_2B1 == "both down" & group_2B4 == "both down" ~"both down 2B1 & 2B4",
                                             group_2B1 == "both up" & (group_2B4 != "both up") ~"both up 2B1 only",
                                             group_2B4 == "both up" & (group_2B1 != "both up") ~"both up 2B4 only",
                                             group_2B1 == "both down" & (group_2B4 != "both down") ~"both down 2B1 only",
                                             group_2B4 == "both down" & (group_2B1 != "both down") ~"both down 2B4 only")
  )) |>
  count(cumulative_group) |>
  mutate(cumulative_group = fct_na_value_to_level(cumulative_group, level = "ns")) |>
  column_to_rownames(var = "cumulative_group")-> cum_group_counts_both_only

# get values analysis only the distribution of RPFs only going up or down

TE_2B1 |>
  full_join(TE_2B4[,c(1,4:9)], by = 'gene_sym', suffix = c('_2B1','_2B4')) |>
  mutate(group_2B1 = fct_na_value_to_level(group_2B1, level = "ns"),
         group_2B4 = fct_na_value_to_level(group_2B4, level = "ns")) |>
  mutate(cumulative_group = factor(case_when(group_2B1 == "RPFs up" & group_2B4 == "RPFs up" ~"RPFs up 2B1 & 2B4",
                                             group_2B1 == "RPFs down" & group_2B4 == "RPFs down" ~"RPFs down 2B1 & 2B4",
                                             group_2B1 == "RPFs up" & (group_2B4 != "RPFs up" | (is.na(group_2B4))) ~"RPFs up 2B1 only",
                                             group_2B4 == "RPFs up" & (group_2B1 != "RPFs up"| (is.na(group_2B1))) ~"RPFs up 2B4 only",
                                             group_2B1 == "RPFs down" & (group_2B4 != "RPFs down"| (is.na(group_2B4))) ~"RPFs down 2B1 only",
                                             group_2B4 == "RPFs down" & (group_2B1 != "RPFs down"| (is.na(group_2B1))) ~"RPFs down 2B4 only"))) |>
  count(cumulative_group) |>
  mutate(cumulative_group = fct_na_value_to_level(cumulative_group, level = "ns")) |>
  column_to_rownames(var = "cumulative_group") -> cum_group_counts_RPFs_only

#last three rows need to be removed if one wants to write out the cumulative groups for ORA

# write file out to do ORA
#write_csv(cum_group_counts_RPFs_only, file = file.path(tables_dir, "20231101_RPFs_groups_for_ORA.csv"))

## both up 2B1 vs both up 2B4 -----

A1 <- (cum_group_counts_both_only["both up 2B1 only",])+(cum_group_counts_both_only["both up 2B1 & 2B4",])
A2 <- (cum_group_counts_both_only["both up 2B4 only",])+(cum_group_counts_both_only["both up 2B1 & 2B4",])
A3 <- cum_group_counts_both_only["both up 2B1 & 2B4",]

pdf(file = paste0(plot_dir, "Venn_2B1_vs_2B4_both_up.pdf"), width = 4, height = 4)
print(draw.pairwise.venn(area1 = A1,
                         area2 = A2,
                         cross.area = A3,
                         fill = c("#440154FF", "#33638DFF"),
                         lwd = rep(1,2),
                         category = c("2B1 \nboth up", "2B4 \nboth up"),
                         cat.pos = c(0, 0),
                         cat.dist = 0.05,
                         cat.fontfamily = "sans",
                         fontfamily = "sans"
                         ))
dev.off()

## both down 2B1 vs both down 2B4 -----

A1 <- (cum_group_counts_both_only["both down 2B1 only",])+(cum_group_counts_both_only["both down 2B1 & 2B4",])
A2 <- (cum_group_counts_both_only["both down 2B4 only",])+(cum_group_counts_both_only["both down 2B1 & 2B4",])
A3 <- cum_group_counts_both_only["both down 2B1 & 2B4",]

pdf(file = paste0(plot_dir, "Venn_2B1_vs_2B4_both_down.pdf"), width = 5, height = 5)
print(draw.pairwise.venn(area1 = A1,
                         area2 = A2,
                         cross.area = A3,
                         fill = c("#440154FF", "#33638DFF"),
                         lwd = rep(1,2),
                         category = c("2B1 \nboth down", "2B4 \nboth down"),
                         cat.pos = c(0,0),
                         cat.dist = 0.05,
                         cat.fontfamily = "sans",
                         fontfamily = "sans"))
dev.off()

## RPFs up 2B1 vs RPFs up 2B4 -----

A1 <- (cum_group_counts_RPFs_only["RPFs up 2B1 only",])+(cum_group_counts_RPFs_only["RPFs up 2B1 & 2B4",])
A2 <- (cum_group_counts_RPFs_only["RPFs up 2B4 only",])+(cum_group_counts_RPFs_only["RPFs up 2B1 & 2B4",])
A3 <- cum_group_counts_RPFs_only["RPFs up 2B1 & 2B4",]

pdf(file = paste0(plot_dir, "Venn_2B1_vs_2B4_RPFs_up.pdf"), width = 5, height = 5)
print(draw.pairwise.venn(area1 = A1,
                         area2 = A2,
                         cross.area = A3,
                         fill = c("#440154FF", "#33638DFF"),
                         lwd = rep(1,2),
                         category = c("2B1 \nRPFs up", "2B4 \nRPFs up"),
                         cat.pos = c(0, 0),
                         cat.dist = 0.05,
                         cat.fontfamily = "sans",
                         fontfamily = "sans"))
dev.off()

## RPFs down 2B1 vs RPFs down 2B4 -----

A1 <- (cum_group_counts_RPFs_only["RPFs down 2B1 only",])+(cum_group_counts_RPFs_only["RPFs down 2B1 & 2B4",])
A2 <- (cum_group_counts_RPFs_only["RPFs down 2B4 only",])+(cum_group_counts_RPFs_only["RPFs down 2B1 & 2B4",])
A3 <- cum_group_counts_RPFs_only["RPFs down 2B1 & 2B4",]

pdf(file = paste0(plot_dir, "Venn_2B1_vs_2B4_RPFs_down.pdf"), width = 5, height = 5)
print(draw.pairwise.venn(area1 = A1,
                         area2 = A2,
                         cross.area = A3,
                         fill = c("#440154FF", "#33638DFF"),
                         lwd = rep(1,2),
                         category = c("2B1 \nRPFs down", "2B4 \nRPFs down"),
                         cat.pos = c(0, 0),
                         cat.dist = 0.05,
                         cat.fontfamily = "sans",
                         fontfamily = "sans"))
dev.off()

# UpSet plots -------
## prepare CSV for upset plot ----

TE_2B1 |>
  dplyr::select(c(1,8,9)) |>
  full_join(TE_2B4[,c(1,8,9)], by = 'gene_sym', suffix = c('_2B1','_2B4')) |>
  mutate(RPF_up_2B1 = as.integer(case_when(group_2B1 == "RPFs up" ~"1", .default = "0")),
         RPF_down_2B1 = as.integer(case_when(group_2B1 == "RPFs down" ~"1", .default = "0")),
         both_up_2B1 = as.integer(case_when(group_2B1 == "both up" ~"1", .default = "0")),
         both_down_2B1 = as.integer(case_when(group_2B1 == "both down" ~"1", .default = "0")),
         RPFs_up_2B4 = as.integer(case_when(group_2B4 == "RPFs up" ~"1", .default = "0")),
         RPF_down_2B4 = as.integer(case_when(group_2B4 == "RPFs down" ~"1", .default = "0")),
         both_up_2B4 = as.integer(case_when(group_2B4 == "both up" ~"1", .default = "0")),
         both_down_2B4 = as.integer(case_when(group_2B4 == "both down" ~"1", .default = "0"))) -> TE_4_upset

write_csv(TE_4_upset, file = file.path(tables_dir, "2B1_vs_2B4_TE_forUpset.csv"))

TE_4_upset |>
  as.data.frame() -> TE_4_upset_df

pdf(file = paste0(plot_dir, "2B1_2B4_upset.pdf"))
upset(TE_4_upset_df, nsets = 18, mb.ratio = c(0.5, 0.5),
      order.by = c("freq"), decreasing = c(TRUE),
      boxplot.summary = c("TE_2B1", "TE_2B4"))
dev.off()

# 5' UTR vs CDS data
TE_UTR5_2B1 |>
  dplyr::select(c(1,7,8)) |>
  inner_join(TE_UTR5_2B4[,c(1,7,8)], by = 'gene_sym', suffix = c('_2B1','_2B4')) |>
  mutate(CDS_up_2B1 = as.integer(case_when(group_2B1 == "CDS up" ~"1", .default = "0")),
         CDS_down_2B1 = as.integer(case_when(group_2B1 == "CDS down" ~"1", .default = "0")),
         both_up_2B1 = as.integer(case_when(group_2B1 == "both up" ~"1", .default = "0")),
         both_down_2B1 = as.integer(case_when(group_2B1 == "both down" ~"1", .default = "0")),
         CDS_up_2B4 = as.integer(case_when(group_2B4 == "CDS up" ~"1", .default = "0")),
         CDS_down_2B4 = as.integer(case_when(group_2B4 == "CDS down" ~"1", .default = "0")),
         both_up_2B4 = as.integer(case_when(group_2B4 == "both up" ~"1", .default = "0")),
         both_down_2B4 = as.integer(case_when(group_2B4 == "both down" ~"1", .default = "0")),
         UTR5_up_2B4 = as.integer(case_when(group_2B4 == "5' up" ~"1", .default = "0")),
         UTR5_down_2B4 = as.integer(case_when(group_2B4 == "5' down" ~"1", .default = "0")),
         UTR5_up_2B1 = as.integer(case_when(group_2B1 == "5' up" ~"1", .default = "0")),
         UTR5_down_2B1 = as.integer(case_when(group_2B1 == "5' down" ~"1", .default = "0"))) -> TE_4_upset

write_csv(TE_4_upset, file = file.path(parent_dir, "Analysis/DESeq2_output", "2B1_vs_2B4_TE_UTR5_forUpset.csv"))

TE_4_upset |>
  as.data.frame() -> TE_4_upset_df

pdf(file = paste0(plot_dir, "2B1_2B4_UTR5_upset.pdf"))
upset(TE_4_upset_df, nsets = 18, mb.ratio = c(0.5, 0.5),
      order.by = c("freq"), decreasing = c(TRUE),
      boxplot.summary = c("TE_2B1", "TE_2B4"))
dev.off()
