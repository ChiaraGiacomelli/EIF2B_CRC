# load packages----
library(tidyverse)
library(ggrepel)
library(viridis)

# read in common variables----
source("common_variables.R")

# theme----
myTheme <- theme_classic()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.4))

# load RSCU data list----
load("~/data/R11/bioinformatics_resources/FASTAs/human/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_RSCU.Rdata")

#read in DESeq2 output----
merged_DEseq_data <- read_csv(file = file.path(tables_dir, "20231101_RPFs_groups_for_ORA.csv"))

summary(merged_DEseq_data)

# extract IDs----
merged_DEseq_data %>%
  dplyr::filter(cumulative_group == "RPFs down 2B1 only") %>%
  pull(transcript) -> RPFs_dn_2B1

merged_DEseq_data %>%
  dplyr::filter(cumulative_group == "RPFs up 2B1 only") %>%
  pull(transcript) -> RPFs_up_2B1

merged_DEseq_data %>%
  dplyr::filter(cumulative_group == "RPFs down 2B4 only") %>%
  pull(transcript) -> RPFs_dn_2B4

merged_DEseq_data %>%
  dplyr::filter(cumulative_group == "RPFs up 2B4 only") %>%
  pull(transcript) -> RPFs_up_2B4

# subset list----
RPFs_down_2B1_rscu_list <- RSCU_list[RPFs_dn_2B1]
RPFs_down_2B1_rscu <- do.call("rbind", RPFs_down_2B1_rscu_list)

RPFs_up_2B1_rscu_list <- RSCU_list[RPFs_up_2B1]
RPFs_up_2B1_rscu <- do.call("rbind", RPFs_up_2B1_rscu_list)

RPFs_down_2B4_rscu_list <- RSCU_list[RPFs_dn_2B4]
RPFs_down_2B4_rscu <- do.call("rbind", RPFs_down_2B4_rscu_list)

RPFs_up_2B4_rscu_list <- RSCU_list[RPFs_up_2B4]
RPFs_up_2B4_rscu <- do.call("rbind", RPFs_up_2B4_rscu_list)

# summarise---- # removes stop codons, start codon, and Tryptophan
RPFs_down_2B1_rscu %>%
  group_by(codon, AA) %>%
  summarise(RPFs_2B1_down_rscu = mean(RSCU, na.rm = T)) %>%
  dplyr::filter(!(codon %in% c("UAA", "UAG", "UGA", "AUG", "UGG"))) -> RPFs_down_2B1_summarised

RPFs_up_2B1_rscu %>%
  group_by(codon, AA) %>%
  summarise(RPFs_up_2B1_rscu = mean(RSCU, na.rm = T)) %>%
  dplyr::filter(!(codon %in% c("UAA", "UAG", "UGA", "AUG", "UGG"))) -> RPFs_up_2B1_summarised

RPFs_down_2B4_rscu %>%
  group_by(codon, AA) %>%
  summarise(RPFs_2B4_down_rscu = mean(RSCU, na.rm = T)) %>%
  dplyr::filter(!(codon %in% c("UAA", "UAG", "UGA", "AUG", "UGG"))) -> RPFs_down_2B4_summarised

RPFs_up_2B4_rscu %>%
  group_by(codon, AA) %>%
  summarise(RPFs_up_2B4_rscu = mean(RSCU, na.rm = T)) %>%
  dplyr::filter(!(codon %in% c("UAA", "UAG", "UGA", "AUG", "UGG"))) -> RPFs_up_2B4_summarised

# plot----
## 2B1 up vs down ------
RPFs_down_2B1_summarised %>%
  inner_join(RPFs_up_2B1_summarised[c(1,3)], by = "codon") %>%
  mutate(wobble = factor(str_sub(codon, 3,3), levels = c("A", "U", "G", "C"))) %>%
  ggplot(aes(x = RPFs_2B1_down_rscu, y = RPFs_up_2B1_rscu, colour = wobble, label = AA))+
  scale_color_viridis(discrete = TRUE, option = "H", begin = 0.8, end = 0.2) +
  geom_point(size = 2.5, alpha = 0.8)+
  geom_text_repel(hjust = "inward", nudge_x = 0.1, color = "black", max.overlaps = 15) +
  xlab("RPFs down")+
  ylab("RPFs up")+
  ggtitle("Relative Synonymous Codon Usage", subtitle = "RPFs down vs up upon 2B1 KD")+
  xlim(c(0,3.5))+
  ylim(c(0,3.5))+
  geom_abline(lty=2)+
  myTheme -> RPFs_down_plot

pdf(file = paste0(plot_dir, "Suppl_5_RSCU/RPFs_down_vs_up_rscu_2B1_exclusive.pdf"), width = 8, height = 7)
print(RPFs_down_plot)
dev.off()

## 2B4 up vs down ------
RPFs_down_2B4_summarised %>%
  inner_join(RPFs_up_2B4_summarised[c(1,3)], by = "codon") %>%
  mutate(wobble = factor(str_sub(codon, 3,3), levels = c("A", "U", "G", "C"))) %>%
  ggplot(aes(x = RPFs_2B4_down_rscu, y = RPFs_up_2B4_rscu, colour = wobble, label = AA))+
  scale_color_viridis(discrete = TRUE, option = "H", begin = 0.8, end = 0.2) +
  geom_point(size = 2.5, alpha = 0.8)+
  geom_text_repel(hjust = "inward", nudge_x = 0.1, color = "black", max.overlaps = 15) +
  xlab("RPFs down")+
  ylab("RPFs up")+
  ggtitle("Relative Synonymous Codon Usage", subtitle = "RPFs down vs up upon 2B4 KD")+
  xlim(c(0,3.5))+
  ylim(c(0,3.5))+
  geom_abline(lty=2)+
  myTheme -> RPFs_down_plot

pdf(file = paste0(plot_dir, "Suppl_5_RSCU/RPFs_down_vs_up_rscu_2B4_exclusive.pdf"), width = 8, height = 7)
print(RPFs_down_plot)
dev.off()

## 2B1 up vs down NO LABELS ------
RPFs_down_2B1_summarised %>%
  inner_join(RPFs_up_2B1_summarised, by = "codon") %>%
  mutate(wobble = factor(str_sub(codon, 3,3), levels = c("A", "U", "G", "C"))) %>%
  ggplot(aes(x = RPFs_2B1_down_rscu, y = RPFs_up_2B1_rscu, colour = wobble))+
  scale_color_viridis(discrete = TRUE, option = "H", begin = 0.8, end = 0.2) +
  geom_point(size = 3.5, alpha = 0.7)+
  xlab("RPFs down")+
  ylab("RPFs up")+
  ggtitle("Relative Synonymous Codon Usage", subtitle = "RPFs down vs up upon 2B1 KD")+
  xlim(c(0,3.5))+
  ylim(c(0,3.5))+
  geom_abline(lty=2)+
  myTheme -> RPFs_down_plot

pdf(file = paste0(plot_dir, "Suppl_5_RSCU/RPFs_down_vs_up_rscu_2B1_noLabel_exclusive.pdf"), width = 8, height = 7)
print(RPFs_down_plot)
dev.off()

## 2B4 up vs down NO LABELS ------
RPFs_down_2B4_summarised %>%
  inner_join(RPFs_up_2B4_summarised, by = "codon") %>%
  mutate(wobble = factor(str_sub(codon, 3,3), levels = c("A", "U", "G", "C"))) %>%
  ggplot(aes(x = RPFs_2B4_down_rscu, y = RPFs_up_2B4_rscu, colour = wobble))+
  scale_color_viridis(discrete = TRUE, option = "H", begin = 0.8, end = 0.2) +
  geom_point(size = 3.5, alpha = 0.7)+
  xlab("RPFs down")+
  ylab("RPFs up")+
  ggtitle("Relative Synonymous Codon Usage", subtitle = "RPFs down vs up upon 2B4 KD")+
  xlim(c(0,3.5))+
  ylim(c(0,3.5))+
  geom_abline(lty=2)+
  myTheme -> RPFs_down_plot

pdf(file = paste0(plot_dir, "Suppl_5_RSCU/RPFs_down_vs_up_rscu_2B4_noLabel_exclusive.pdf"), width = 8, height = 7)
print(RPFs_down_plot)
dev.off()

