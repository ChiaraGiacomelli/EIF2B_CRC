# get all my data in a single table - merging all normalised and all log2FC
# Required for Venn Diagram making
# Calculate Baseline TE in uninduced conditions for Gradient Boosting 

# version from 20231101

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

plot_dir <- paste0(parent_dir,"/plots/PDF_for_paper/")
tables_dir <- paste0(parent_dir,"/Analysis/tables_for_paper")

## DEseq2 results for totals & RPF ------
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
# merge normalised counts for RPFs and totals within same shRNA KD &
# calculate baseline TE in UNINDUCED condition - to be used in gradient boosting 

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

## Write out the baseline TEs -----
write_csv(B1_nc_uninduced, file = file.path(tables_dir, "BaselineTE_2B1.csv"))
write_csv(B4_nc_uninduced, file = file.path(tables_dir, "BaselineTE_2B4.csv"))

#  VENN diagrams ------

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

#last three rows need to be removed if one wants to write out the cumulative groups for ORA - here below
# write file out to do ORA
#write_csv(cum_group_counts_RPFs_only, file = file.path(tables_dir, "20231101_RPFs_groups_for_ORA.csv"))

## both up 2B1 vs both up 2B4 (both means RNA and RPF changing concordantly) ----- 
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

## both down 2B1 vs both down 2B4 (both means RNA and RPF changing concordantly) -----
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

## RPFs up 2B1 vs RPFs up 2B4 (both means RNA and RPF changing concordantly) -----
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

## RPFs down 2B1 vs RPFs down 2B4 (both means RNA and RPF changing concordantly) -----
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
