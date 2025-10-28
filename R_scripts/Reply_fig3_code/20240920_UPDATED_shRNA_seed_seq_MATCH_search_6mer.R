# This script searches for the shRNA seed sequences in 3' UTRs and CDSs of transcripts quantified
# As we used only 1x shRNA per condition in Riboseq this is supposed to check for potential off-targets effects
library(tidyverse)
library(viridis)
library(conflicted)
library(ggplot2)
library(ggrepel)
library(seqinr)
library(spgs)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(base::intersect)
conflicts_prefer(base::setdiff)
options(ggrepel.max.overlaps = 50)

# shared info & variables ------
## themes----
mytheme <- theme_minimal()+
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

violin_theme <- theme_bw()+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        plot.title = element_text(hjust = 1, vjust = 0, size=18, face="bold"),
        legend.title = element_blank(),
        legend.text = element_blank())

source("common_variables.R")

# FASTA directory
fasta_dir <- ("your/fasta/folder")

# directories for saving
plot_dir <- paste0(parent_dir,"/plots/PDF_for_paper/")
tables_dir <- paste0(parent_dir,"/Analysis/tables_for_paper")

# sequences of Dox-inducible shRNAs as described in Reagents and Tools table
shCTRL <- "cttactctcgcccaagcgagag"
shEIF2B1 <- "ttatacttaaacttatctggga"
shEIF2B4 <- "tagactagattcaacaaccgta"

# seed-matches as defined in fig 1A of PMID: 17612493
# 6mer = nt 2 to 6
seed_CTRL <- "ttactc"
seed_2B1 <- "tatact"
seed_2B4 <- "agacta"

match_CTRL <- spgs::reverseComplement(seed_CTRL)
match_2B1 <- spgs::reverseComplement(seed_2B1)
match_2B4 <- spgs::reverseComplement(seed_2B4)

# read in DE analyses results -------
## DEseq2 results for totals, RPFs, and 5UTRs ------
totals_2B1 <- read_csv(file = file.path(tables_dir, "Totals_2B1_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv"))
RPFs_2B1 <- read_csv(file = file.path(tables_dir, "RPFs_2B1_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv"))

totals_2B4 <- read_csv(file = file.path(tables_dir, "Totals_2B4_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv"))
RPFs_2B4 <- read_csv(file = file.path(tables_dir, "RPFs_2B4_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv"))

## TE data ------
TE_NTCcor_2B1 <- read_csv(file = file.path(tables_dir, "TE_2B1_NTCcorrected_4fgsea.csv"))
TE_NTCcor_2B4 <- read_csv(file = file.path(tables_dir, "TE_2B4_NTCcorrected_4fgsea.csv"))

## TE data including groups
groups_2B1 <- read_csv(file = file.path(tables_dir, "2B1_NTCcorrected_TE_genes_CDSvsTotals.csv"))
groups_2B4 <- read_csv(file = file.path(tables_dir, "2B4_NTCcorrected_TE_genes_CDSvsTotals.csv"))

## merge all basemeans,log2FC,adjpvals -----
RPFs_2B1 |>
  dplyr::select(transcript, gene, gene_sym, baseMean, log2FoldChange, padj) |>
  full_join(RPFs_2B4[,c('gene', 'baseMean', 'log2FoldChange', 'padj')], by = 'gene', suffix = c('_2B1','_2B4')) |>
  relocate(gene)-> RPFs_all

totals_2B1 |>
  dplyr::select(transcript, gene, gene_sym, baseMean, log2FoldChange, padj) |>
  full_join(totals_2B4[,c('gene', 'baseMean', 'log2FoldChange', 'padj')], by = 'gene', suffix = c('_2B1','_2B4')) |>
  relocate(gene)-> totals_all

# Merge TE data
TE_NTCcor_2B1 |>
  inner_join(TE_NTCcor_2B4, by = 'gene', suffix = c('_2B1','_2B4')) %>% 
  select(-c(9, 10)) %>% 
  rename(transcript = transcript_2B1, gene_sym = gene_sym_2B1) %>% 
  inner_join(totals_all[c(1:4,7)]) %>% 
  inner_join(RPFs_all[c(1:4,7)], by = join_by(transcript, gene, gene_sym), suffix = c("_tots", "_RPFs")) -> Master_table

## compare the RPF only and both categories from DE analyses ------
groups_2B1 |>
  dplyr::select(-alpha_score) |>
  dplyr::mutate(TE = RPFs_log2FC - totals_log2FC) -> TE_2B1

groups_2B4 |>
  dplyr::select(-alpha_score) |>
  dplyr::mutate(TE = RPFs_log2FC - totals_log2FC) -> TE_2B4

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
  relocate(c(cumulative_group, group_2B1, group_2B4, TE_2B1, TE_2B4), .after = "gene") %>% 
  inner_join(Master_table[c(1:3, 14:17)])-> merged_TE


# CDS =======
## read in the fasta & count sites -------
CDS <- read.fasta(paste0(fasta_dir, "/gencode.v38.pc_transcripts_filtered_CDS.fa"))
CDS <- lapply(CDS, str_c, collapse = "")

counts_site <- sapply(CDS, str_count, pattern = match_CTRL)
counts_site %>% 
  map_dfr(~.x %>% 
            as_tibble(), .id = "transcript_id") -> CDS_ctrl_seed_match

counts_site <- sapply(CDS, str_count, pattern = match_2B1)
counts_site %>% 
  map_dfr(~.x %>% 
            as_tibble(), .id = "transcript_id") -> CDS_2B1_seed_match

counts_site <- sapply(CDS, str_count, pattern = match_2B4)
counts_site %>% 
  map_dfr(~.x %>% 
            as_tibble(), .id = "transcript_id") -> CDS_2B4_seed_match

counts_site %>% 
  as_tibble

## merge with TE data for plots -------
merged_TE %>% 
  inner_join(CDS_ctrl_seed_match, by = join_by(transcript == transcript_id)) %>%
  rename(CDS_seed_matchmatches_shCTRL = value) %>% 
  inner_join(CDS_2B1_seed_match,  by = join_by(transcript == transcript_id)) %>% 
  rename(CDS_seed_matchmatches_sh2B1 = value) %>% 
  inner_join(CDS_2B4_seed_match, by = join_by(transcript == transcript_id)) %>% 
  rename(CDS_seed_matchmatches_sh2B4 = value) %>% 
  mutate(CDS_seed_matchmatches = factor(case_when(CDS_seed_matchmatches_shCTRL > 0 & CDS_seed_matchmatches_sh2B1 > 0 & CDS_seed_matchmatches_sh2B4 > 0 ~ "CDS seed matches\nall shRNAs",
                                                  CDS_seed_matchmatches_shCTRL == 0 & CDS_seed_matchmatches_sh2B1 > 0 & CDS_seed_matchmatches_sh2B4 > 0 ~ "CDS seed matches\nsh2B1 & sh2B4",
                                                  CDS_seed_matchmatches_shCTRL == 0 & CDS_seed_matchmatches_sh2B1 == 0 & CDS_seed_matchmatches_sh2B4 > 0 ~ "CDS seed matches\nsh2B4 only",
                                                  CDS_seed_matchmatches_shCTRL == 0 & CDS_seed_matchmatches_sh2B1 > 0 & CDS_seed_matchmatches_sh2B4 == 0 ~ "CDS seed matches\nsh2B1 only",
                                                  CDS_seed_matchmatches_shCTRL > 0 & CDS_seed_matchmatches_sh2B1 == 0 & CDS_seed_matchmatches_sh2B4 == 0 ~ "CDS seed matches\nshCTRL only",
                                                  CDS_seed_matchmatches_shCTRL > 0 & CDS_seed_matchmatches_sh2B1 > 0 & CDS_seed_matchmatches_sh2B4 == 0 ~ "CDS seed matches\nshCTRL & sh2B1",
                                                  CDS_seed_matchmatches_shCTRL > 0 & CDS_seed_matchmatches_sh2B1 == 0 & CDS_seed_matchmatches_sh2B4 > 0 ~ "CDS seed matches\nshCTRL & sh2B4"))) %>% 
  relocate(CDS_seed_matchmatches, .after = "gene") -> merged_TE_CDS_seed_matches

summary(merged_TE_CDS_seed_matches)

## Cumulative Fractions ----
### 2B1 shRNA in 2B1 data ========
merged_TE_CDS_seed_matches %>%
  filter(gene_sym != "EIF2B1") %>% 
  mutate(CDS_seed_matchmatches_sh2B1 = as_factor(CDS_seed_matchmatches_sh2B1)) %>% 
  select(CDS_seed_matchmatches_sh2B1, TE_2B1) %>%
  ggplot(aes(x = TE_2B1, colour = CDS_seed_matchmatches_sh2B1)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "sh2B1 6mer\nMatches in CDS", direction = -1,
                      labels = merged_TE_CDS_seed_matches %>%
                        group_by(CDS_seed_matchmatches_sh2B1) %>%
                        summarize(num = n(), med = base::round(mean(TE_2B1), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", CDS_seed_matchmatches_sh2B1, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("Translational Efficiency (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Translation Efficiency") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_CDS_seed_matches_2B1_TE_wo_2B1.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_CDS_seed_matches %>% 
  filter(gene_sym != "EIF2B1") %>% 
  mutate(CDS_seed_matchmatches_sh2B1 = as_factor(CDS_seed_matchmatches_sh2B1)) %>% 
  select(CDS_seed_matchmatches_sh2B1, RPFs_log2FC_2B1) %>% 
  ggplot(aes(x = RPFs_log2FC_2B1, colour = CDS_seed_matchmatches_sh2B1)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "sh2B1 6mer\nMatches in CDS", direction = -1,
                      labels = merged_TE_CDS_seed_matches %>%
                        group_by(CDS_seed_matchmatches_sh2B1) %>%
                        summarize(num = n(), med = base::round(mean(RPFs_log2FC_2B1), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", CDS_seed_matchmatches_sh2B1, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("RPFs (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Ribosome Occupancy") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_CDS_seed_matches_2B1_RPF_wo_2B1.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_CDS_seed_matches %>% 
  filter(gene_sym != "EIF2B1") %>% 
  mutate(CDS_seed_matchmatches_sh2B1 = as_factor(CDS_seed_matchmatches_sh2B1)) %>% 
  select(CDS_seed_matchmatches_sh2B1, totals_log2FC_2B1) %>% 
  ggplot(aes(x = totals_log2FC_2B1, colour = CDS_seed_matchmatches_sh2B1)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "sh2B1 6mer\nMatches in CDS", direction = -1,
                      labels = merged_TE_CDS_seed_matches %>%
                        group_by(CDS_seed_matchmatches_sh2B1) %>%
                        summarize(num = n(), med = base::round(mean(totals_log2FC_2B1), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", CDS_seed_matchmatches_sh2B1, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("total RNA (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Cytoplasmic Tot RNA") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_CDS_seed_matches_2B1_tots_wo_2B1.pdf"), width = 8, height = 6)
print(plot)
dev.off()

### 2B4 shRNA in 2B4 data ========
merged_TE_CDS_seed_matches %>% 
  filter(gene_sym != "EIF2B4") %>% 
  mutate(CDS_seed_matchmatches_sh2B4 = as_factor(CDS_seed_matchmatches_sh2B4)) %>% 
  select(CDS_seed_matchmatches_sh2B4, TE_2B4) %>%
  ggplot(aes(x = TE_2B4, colour = CDS_seed_matchmatches_sh2B4)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "sh2B4 6mer\nMatches in CDS", direction = -1,
                      labels = merged_TE_CDS_seed_matches %>%
                        group_by(CDS_seed_matchmatches_sh2B4) %>%
                        summarize(num = n(), med = base::round(mean(TE_2B4), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", CDS_seed_matchmatches_sh2B4, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("Translational Efficiency (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Translation Efficiency") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_CDS_seed_matches_2B4_TE_wo_2B4.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_CDS_seed_matches %>% 
  filter(gene_sym != "EIF2B4") %>% 
  mutate(CDS_seed_matchmatches_sh2B4 = as_factor(CDS_seed_matchmatches_sh2B4)) %>% 
  select(CDS_seed_matchmatches_sh2B4, RPFs_log2FC_2B4) %>% 
  ggplot(aes(x = RPFs_log2FC_2B4, colour = CDS_seed_matchmatches_sh2B4)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "sh2B4 6mer\nMatches in CDS", direction = -1,
                      labels = merged_TE_CDS_seed_matches %>%
                        group_by(CDS_seed_matchmatches_sh2B4) %>%
                        summarize(num = n(), med = base::round(mean(RPFs_log2FC_2B4), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", CDS_seed_matchmatches_sh2B4, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("RPFs (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Ribosome Occupancy") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_CDS_seed_matches_2B4_RPF_wo_2B4.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_CDS_seed_matches %>% 
  filter(gene_sym != "EIF2B4") %>% 
  mutate(CDS_seed_matchmatches_sh2B4 = as_factor(CDS_seed_matchmatches_sh2B4)) %>% 
  select(CDS_seed_matchmatches_sh2B4, totals_log2FC_2B4) %>% 
  ggplot(aes(x = totals_log2FC_2B4, colour = CDS_seed_matchmatches_sh2B4)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "sh2B4 6mer\nMatches in CDS", direction = -1,
                      labels = merged_TE_CDS_seed_matches %>%
                        group_by(CDS_seed_matchmatches_sh2B4) %>%
                        summarize(num = n(), med = base::round(mean(totals_log2FC_2B4), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", CDS_seed_matchmatches_sh2B4, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("total RNA (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Cytoplasmic Tot RNA") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_CDS_seed_matches_2B4_tots_wo_2B4.pdf"), width = 8, height = 6)
print(plot)
dev.off()

### CTRL shRNA in 2B1 data ========
merged_TE_CDS_seed_matches %>% 
  mutate(CDS_seed_matchmatches_shCTRL = as_factor(CDS_seed_matchmatches_shCTRL)) %>% 
  select(CDS_seed_matchmatches_shCTRL, TE_2B1) %>%
  ggplot(aes(x = TE_2B1, colour = CDS_seed_matchmatches_shCTRL)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "shCTRL 6mer\nMatches in CDS", direction = -1,
                      labels = merged_TE_CDS_seed_matches %>%
                        group_by(CDS_seed_matchmatches_shCTRL) %>%
                        summarize(num = n(), med = base::round(mean(TE_2B1), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", CDS_seed_matchmatches_shCTRL, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("Translational Efficiency (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Translation Efficiency") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_CDS_seed_matches_CTRLin2B1_TE.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_CDS_seed_matches %>% 
  mutate(CDS_seed_matchmatches_shCTRL = as_factor(CDS_seed_matchmatches_shCTRL)) %>% 
  select(CDS_seed_matchmatches_shCTRL, RPFs_log2FC_2B1) %>% 
  ggplot(aes(x = RPFs_log2FC_2B1, colour = CDS_seed_matchmatches_shCTRL)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "shCTRL 6mer\nMatches in CDS", direction = -1,
                      labels = merged_TE_CDS_seed_matches %>%
                        group_by(CDS_seed_matchmatches_shCTRL) %>%
                        summarize(num = n(), med = base::round(mean(RPFs_log2FC_2B1), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", CDS_seed_matchmatches_shCTRL, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("RPFs (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Ribosome Occupancy") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_CDS_seed_matches_CTRLin2B1_RPF.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_CDS_seed_matches %>% 
  mutate(CDS_seed_matchmatches_shCTRL = as_factor(CDS_seed_matchmatches_shCTRL)) %>% 
  select(CDS_seed_matchmatches_shCTRL, totals_log2FC_2B1) %>% 
  ggplot(aes(x = totals_log2FC_2B1, colour = CDS_seed_matchmatches_shCTRL)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "shCTRL 6mer\nMatches in CDS", direction = -1,
                      labels = merged_TE_CDS_seed_matches %>%
                        group_by(CDS_seed_matchmatches_shCTRL) %>%
                        summarize(num = n(), med = base::round(mean(totals_log2FC_2B1), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", CDS_seed_matchmatches_shCTRL, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("total RNA (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Cytoplasmic Tot RNA") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_CDS_seed_matches_CTRLin2B1_tots.pdf"), width = 8, height = 6)
print(plot)
dev.off()

### CTRL shRNA in 2B4 data ========
merged_TE_CDS_seed_matches %>% 
  mutate(CDS_seed_matchmatches_shCTRL = as_factor(CDS_seed_matchmatches_shCTRL)) %>% 
  select(CDS_seed_matchmatches_shCTRL, TE_2B4) %>%
  ggplot(aes(x = TE_2B4, colour = CDS_seed_matchmatches_shCTRL)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "shCTRL 6mer\nMatches in CDS", direction = -1,
                      labels = merged_TE_CDS_seed_matches %>%
                        group_by(CDS_seed_matchmatches_shCTRL) %>%
                        summarize(num = n(), med = base::round(mean(TE_2B4), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", CDS_seed_matchmatches_shCTRL, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("Translational Efficiency (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Translation Efficiency") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_CDS_seed_matches_CTRLin2B4_TE.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_CDS_seed_matches %>% 
  mutate(CDS_seed_matchmatches_shCTRL = as_factor(CDS_seed_matchmatches_shCTRL)) %>% 
  select(CDS_seed_matchmatches_shCTRL, RPFs_log2FC_2B4) %>% 
  ggplot(aes(x = RPFs_log2FC_2B4, colour = CDS_seed_matchmatches_shCTRL)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "shCTRL 6mer\nMatches in CDS", direction = -1,
                      labels = merged_TE_CDS_seed_matches %>%
                        group_by(CDS_seed_matchmatches_shCTRL) %>%
                        summarize(num = n(), med = base::round(mean(RPFs_log2FC_2B4), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", CDS_seed_matchmatches_shCTRL, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("RPFs (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Ribosome Occupancy") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_CDS_seed_matches_CTRLin2B4_RPF.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_CDS_seed_matches %>% 
  mutate(CDS_seed_matchmatches_shCTRL = as_factor(CDS_seed_matchmatches_shCTRL)) %>% 
  select(CDS_seed_matchmatches_shCTRL, totals_log2FC_2B4) %>% 
  ggplot(aes(x = totals_log2FC_2B4, colour = CDS_seed_matchmatches_shCTRL)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "shCTRL 6mer\nMatches in CDS", direction = -1,
                      labels = merged_TE_CDS_seed_matches %>%
                        group_by(CDS_seed_matchmatches_shCTRL) %>%
                        summarize(num = n(), med = base::round(mean(totals_log2FC_2B4), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", CDS_seed_matchmatches_shCTRL, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("total RNA (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Cytoplasmic Tot RNA") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_CDS_seed_matches_CTRLin2B4_tots.pdf"), width = 8, height = 6)
print(plot)
dev.off()

# 3' UTR ----------
## read in the fasta & count sites -------
UTR3 <- read.fasta(paste0(fasta_dir, "/gencode.v38.pc_transcripts_filtered_UTR3.fa"))
UTR3 <- lapply(UTR3, str_c, collapse = "")

counts_site <- sapply(UTR3, str_count, pattern = match_CTRL)
counts_site %>% 
  map_dfr(~.x %>% 
            as_tibble(), .id = "transcript_id") -> UTR3_ctrl_seed_match

counts_site <- sapply(UTR3, str_count, pattern = match_2B1)
counts_site %>% 
  map_dfr(~.x %>% 
            as_tibble(), .id = "transcript_id") -> UTR3_2B1_seed_match

counts_site <- sapply(UTR3, str_count, pattern = match_2B4)
counts_site %>% 
  map_dfr(~.x %>% 
            as_tibble(), .id = "transcript_id") -> UTR3_2B4_seed_match

## merge with TE data for plots -------
merged_TE %>% 
  inner_join(UTR3_ctrl_seed_match, by = join_by(transcript == transcript_id)) %>%
  rename(UTR3_seed_matchmatches_shCTRL = value) %>% 
  inner_join(UTR3_2B1_seed_match,  by = join_by(transcript == transcript_id)) %>% 
  rename(UTR3_seed_matchmatches_sh2B1 = value) %>% 
  inner_join(UTR3_2B4_seed_match, by = join_by(transcript == transcript_id)) %>% 
  rename(UTR3_seed_matchmatches_sh2B4 = value) %>% 
  mutate(UTR3_seed_matchmatches = factor(case_when(UTR3_seed_matchmatches_shCTRL > 0 & UTR3_seed_matchmatches_sh2B1 > 0 & UTR3_seed_matchmatches_sh2B4 > 0 ~ "UTR3 seed matches\nall shRNAs",
                                                   UTR3_seed_matchmatches_shCTRL == 0 & UTR3_seed_matchmatches_sh2B1 > 0 & UTR3_seed_matchmatches_sh2B4 > 0 ~ "UTR3 seed matches\nsh2B1 & sh2B4",
                                                   UTR3_seed_matchmatches_shCTRL == 0 & UTR3_seed_matchmatches_sh2B1 == 0 & UTR3_seed_matchmatches_sh2B4 > 0 ~ "UTR3 seed matches\nsh2B4 only",
                                                   UTR3_seed_matchmatches_shCTRL == 0 & UTR3_seed_matchmatches_sh2B1 > 0 & UTR3_seed_matchmatches_sh2B4 == 0 ~ "UTR3 seed matches\nsh2B1 only",
                                                   UTR3_seed_matchmatches_shCTRL > 0 & UTR3_seed_matchmatches_sh2B1 == 0 & UTR3_seed_matchmatches_sh2B4 == 0 ~ "UTR3 seed matches\nshCTRL only",
                                                   UTR3_seed_matchmatches_shCTRL > 0 & UTR3_seed_matchmatches_sh2B1 > 0 & UTR3_seed_matchmatches_sh2B4 == 0 ~ "UTR3 seed matches\nshCTRL & sh2B1",
                                                   UTR3_seed_matchmatches_shCTRL > 0 & UTR3_seed_matchmatches_sh2B1 == 0 & UTR3_seed_matchmatches_sh2B4 > 0 ~ "UTR3 seed matches\nshCTRL & sh2B4"))) %>% 
  relocate(UTR3_seed_matchmatches, .after = "gene") -> merged_TE_UTR3_seed_matches

summary(merged_TE_UTR3_seed_matches)

## Cumulative Fractions ----
### 2B1 shRNA in 2B1 data ========
merged_TE_UTR3_seed_matches %>% 
  filter(gene_sym != "EIF2B1") %>% 
  mutate(UTR3_seed_matchmatches_sh2B1 = as_factor(UTR3_seed_matchmatches_sh2B1)) %>% 
  select(UTR3_seed_matchmatches_sh2B1, TE_2B1) %>%
  ggplot(aes(x = TE_2B1, colour = UTR3_seed_matchmatches_sh2B1)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "sh2B1 6mer\nMatches in 3' UTR", direction = -1,
                      labels = merged_TE_UTR3_seed_matches %>%
                        group_by(UTR3_seed_matchmatches_sh2B1) %>%
                        summarize(num = n(), med = base::round(mean(TE_2B1), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", UTR3_seed_matchmatches_sh2B1, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("Translational Efficiency (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Translation Efficiency") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_UTR3_seed_matches_2B1_TE_wo_2B1.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_UTR3_seed_matches %>% 
  filter(gene_sym != "EIF2B1") %>% 
  mutate(UTR3_seed_matchmatches_sh2B1 = as_factor(UTR3_seed_matchmatches_sh2B1)) %>% 
  select(UTR3_seed_matchmatches_sh2B1, RPFs_log2FC_2B1) %>% 
  ggplot(aes(x = RPFs_log2FC_2B1, colour = UTR3_seed_matchmatches_sh2B1)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "sh2B1 6mer\nMatches in 3' UTR", direction = -1,
                      labels = merged_TE_UTR3_seed_matches %>%
                        group_by(UTR3_seed_matchmatches_sh2B1) %>%
                        summarize(num = n(), med = base::round(mean(RPFs_log2FC_2B1), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", UTR3_seed_matchmatches_sh2B1, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("RPFs (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Ribosome Occupancy") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_UTR3_seed_matches_2B1_RPF_wo_2B1.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_UTR3_seed_matches %>% 
  filter(gene_sym != "EIF2B1") %>% 
  mutate(UTR3_seed_matchmatches_sh2B1 = as_factor(UTR3_seed_matchmatches_sh2B1)) %>% 
  select(UTR3_seed_matchmatches_sh2B1, totals_log2FC_2B1) %>% 
  ggplot(aes(x = totals_log2FC_2B1, colour = UTR3_seed_matchmatches_sh2B1)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "sh2B1 6mer\nMatches in 3' UTR", direction = -1,
                      labels = merged_TE_UTR3_seed_matches %>%
                        group_by(UTR3_seed_matchmatches_sh2B1) %>%
                        summarize(num = n(), med = base::round(mean(totals_log2FC_2B1), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", UTR3_seed_matchmatches_sh2B1, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("total RNA (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Cytoplasmic Tot RNA") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_UTR3_seed_matches_2B1_tots_wo_2B1.pdf"), width = 8, height = 6)
print(plot)
dev.off()

### 2B4 shRNA in 2B4 data ========
merged_TE_UTR3_seed_matches %>% 
  filter(gene_sym != "EIF2B4") %>% 
  mutate(UTR3_seed_matchmatches_sh2B4 = as_factor(UTR3_seed_matchmatches_sh2B4)) %>% 
  select(UTR3_seed_matchmatches_sh2B4, TE_2B4) %>%
  ggplot(aes(x = TE_2B4, colour = UTR3_seed_matchmatches_sh2B4)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "sh2B4 6mer\nMatches in 3' UTR", direction = -1,
                      labels = merged_TE_UTR3_seed_matches %>%
                        group_by(UTR3_seed_matchmatches_sh2B4) %>%
                        summarize(num = n(), med = base::round(mean(TE_2B4), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", UTR3_seed_matchmatches_sh2B4, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("Translational Efficiency (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Translation Efficiency") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_UTR3_seed_matches_2B4_TE_wo_2B4.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_UTR3_seed_matches %>% 
  filter(gene_sym != "EIF2B4") %>% 
  mutate(UTR3_seed_matchmatches_sh2B4 = as_factor(UTR3_seed_matchmatches_sh2B4)) %>% 
  select(UTR3_seed_matchmatches_sh2B4, RPFs_log2FC_2B4) %>% 
  ggplot(aes(x = RPFs_log2FC_2B4, colour = UTR3_seed_matchmatches_sh2B4)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "sh2B4 6mer\nMatches in 3' UTR", direction = -1,
                      labels = merged_TE_UTR3_seed_matches %>%
                        group_by(UTR3_seed_matchmatches_sh2B4) %>%
                        summarize(num = n(), med = base::round(mean(RPFs_log2FC_2B4), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", UTR3_seed_matchmatches_sh2B4, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("RPFs (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Ribosome Occupancy") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_UTR3_seed_matches_2B4_RPF_wo_2B4.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_UTR3_seed_matches %>% 
  filter(gene_sym != "EIF2B4") %>% 
  mutate(UTR3_seed_matchmatches_sh2B4 = as_factor(UTR3_seed_matchmatches_sh2B4)) %>% 
  select(UTR3_seed_matchmatches_sh2B4, totals_log2FC_2B4) %>% 
  ggplot(aes(x = totals_log2FC_2B4, colour = UTR3_seed_matchmatches_sh2B4)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "sh2B4 6mer\nMatches in 3' UTR", direction = -1,
                      labels = merged_TE_UTR3_seed_matches %>%
                        group_by(UTR3_seed_matchmatches_sh2B4) %>%
                        summarize(num = n(), med = base::round(mean(totals_log2FC_2B4), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", UTR3_seed_matchmatches_sh2B4, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("total RNA (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Cytoplasmic Tot RNA") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_UTR3_seed_matches_2B4_tots_wo_2B4.pdf"), width = 8, height = 6)
print(plot)
dev.off()

### CTRL shRNA in 2B1 data ========
merged_TE_UTR3_seed_matches %>% 
  mutate(UTR3_seed_matchmatches_shCTRL = as_factor(UTR3_seed_matchmatches_shCTRL)) %>% 
  select(UTR3_seed_matchmatches_shCTRL, TE_2B1) %>%
  ggplot(aes(x = TE_2B1, colour = UTR3_seed_matchmatches_shCTRL)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "shCTRL 6mer\nMatches in 3' UTR", direction = -1,
                      labels = merged_TE_UTR3_seed_matches %>%
                        group_by(UTR3_seed_matchmatches_shCTRL) %>%
                        summarize(num = n(), med = base::round(mean(TE_2B1), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", UTR3_seed_matchmatches_shCTRL, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("Translational Efficiency (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Translation Efficiency") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_UTR3_seed_matches_CTRLin2B1_TE.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_UTR3_seed_matches %>% 
  mutate(UTR3_seed_matchmatches_shCTRL = as_factor(UTR3_seed_matchmatches_shCTRL)) %>% 
  select(UTR3_seed_matchmatches_shCTRL, RPFs_log2FC_2B1) %>% 
  ggplot(aes(x = RPFs_log2FC_2B1, colour = UTR3_seed_matchmatches_shCTRL)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "shCTRL 6mer\nMatches in 3' UTR", direction = -1,
                      labels = merged_TE_UTR3_seed_matches %>%
                        group_by(UTR3_seed_matchmatches_shCTRL) %>%
                        summarize(num = n(), med = base::round(mean(RPFs_log2FC_2B1), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", UTR3_seed_matchmatches_shCTRL, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("RPFs (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Ribosome Occupancy") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_UTR3_seed_matches_CTRLin2B1_RPF.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_UTR3_seed_matches %>% 
  mutate(UTR3_seed_matchmatches_shCTRL = as_factor(UTR3_seed_matchmatches_shCTRL)) %>% 
  select(UTR3_seed_matchmatches_shCTRL, totals_log2FC_2B1) %>% 
  ggplot(aes(x = totals_log2FC_2B1, colour = UTR3_seed_matchmatches_shCTRL)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "shCTRL 6mer\nMatches in 3' UTR", direction = -1,
                      labels = merged_TE_UTR3_seed_matches %>%
                        group_by(UTR3_seed_matchmatches_shCTRL) %>%
                        summarize(num = n(), med = base::round(mean(totals_log2FC_2B1), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", UTR3_seed_matchmatches_shCTRL, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("total RNA (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Cytoplasmic Tot RNA") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_UTR3_seed_matches_CTRLin2B1_tots.pdf"), width = 8, height = 6)
print(plot)
dev.off()

### CTRL shRNA in 2B4 data ========
merged_TE_UTR3_seed_matches %>% 
  mutate(UTR3_seed_matchmatches_shCTRL = as_factor(UTR3_seed_matchmatches_shCTRL)) %>% 
  select(UTR3_seed_matchmatches_shCTRL, TE_2B4) %>%
  ggplot(aes(x = TE_2B4, colour = UTR3_seed_matchmatches_shCTRL)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "shCTRL 6mer\nMatches in 3' UTR", direction = -1,
                      labels = merged_TE_UTR3_seed_matches %>%
                        group_by(UTR3_seed_matchmatches_shCTRL) %>%
                        summarize(num = n(), med = base::round(mean(TE_2B4), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", UTR3_seed_matchmatches_shCTRL, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("Translational Efficiency (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Translation Efficiency") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_UTR3_seed_matches_CTRLin2B4_TE.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_UTR3_seed_matches %>% 
  mutate(UTR3_seed_matchmatches_shCTRL = as_factor(UTR3_seed_matchmatches_shCTRL)) %>% 
  select(UTR3_seed_matchmatches_shCTRL, RPFs_log2FC_2B4) %>% 
  ggplot(aes(x = RPFs_log2FC_2B4, colour = UTR3_seed_matchmatches_shCTRL)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "shCTRL 6mer\nMatches in 3' UTR", direction = -1,
                      labels = merged_TE_UTR3_seed_matches %>%
                        group_by(UTR3_seed_matchmatches_shCTRL) %>%
                        summarize(num = n(), med = base::round(mean(RPFs_log2FC_2B4), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", UTR3_seed_matchmatches_shCTRL, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("RPFs (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Ribosome Occupancy") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_UTR3_seed_matches_CTRLin2B4_RPF.pdf"), width = 8, height = 6)
print(plot)
dev.off()

merged_TE_UTR3_seed_matches %>% 
  mutate(UTR3_seed_matchmatches_shCTRL = as_factor(UTR3_seed_matchmatches_shCTRL)) %>% 
  select(UTR3_seed_matchmatches_shCTRL, totals_log2FC_2B4) %>% 
  ggplot(aes(x = totals_log2FC_2B4, colour = UTR3_seed_matchmatches_shCTRL)) +
  stat_ecdf(linewidth = 1) +
  mytheme +
  scale_color_viridis(discrete = TRUE, name = "shCTRL 6mer\nMatches in 3' UTR", direction = -1,
                      labels = merged_TE_UTR3_seed_matches %>%
                        group_by(UTR3_seed_matchmatches_shCTRL) %>%
                        summarize(num = n(), med = base::round(mean(totals_log2FC_2B4), digits = 2)) %>%
                        mutate(lab = paste0("\nmatches = ", UTR3_seed_matchmatches_shCTRL, "\nn = ", num, "\nmean = ", med)) %>%
                        pull(lab)) +
  xlab("total RNA (log2)") +
  ylab("Cumulative Fraction") +
  ggtitle("CDF of Cytoplasmic Tot RNA") -> plot

pdf(file = paste0(plot_dir, "seed_matches/", "ECDF_UTR3_seed_matches_CTRLin2B4_tots.pdf"), width = 8, height = 6)
print(plot)
dev.off()

# write data out ----
merged_TE_CDS_seed_matches %>% 
  inner_join(merged_TE_UTR3_seed_matches[,c(2,4,22:24)]) -> merged_TE_all_seed_matches

write_csv(merged_TE_all_seed_matches, file = file.path(tables_dir, "shRNA_seed_match_sequence_matches_6mers.csv"))
