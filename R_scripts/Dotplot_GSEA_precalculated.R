# Wurzburg sent over some GSEA results performed from RNAsequencing of SW480 cells with restored APC expression (keep in mind still KRAS mut)
# Steffi wrote that for the files gsea_report_for_shRNA_EtOH_1539157605992 and gsea_report_for_ETOH_1506082050131 the DE data imputed in was the Dox vs EtOH
# where DOX induces the re-expression of APC - therefore the pathways go down when APC is reintroduced 

library(tidyverse)
library(viridis)

# folder with data and where to put plots
#home <- 'N:/CGIACOME/AAA_eIF2B/AAA_Collaborative_stuff/Wurzburg_collab/20230925_GSEA_data_transfer'
home <- '~/data/CGIACOME/AAA_eIF2B/AAA_Collaborative_stuff/Wurzburg_collab/20230925_GSEA_data_transfer'

mytheme <- theme_classic()+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# read in tables
RNA_seq_only <- read_csv(file = file.path(home, 'gsea_report_for_ETOH_1506082050131.csv'))
From_RiboSeq_dataset <- read_csv(file = file.path(home, 'gsea_report_for_shRNA_EtOH_1539157605992.csv'))

RNA_seq_only |>
  select('NAME', 'NES', 'FDR q-val') |>
  rename(FDR_qval = 'FDR q-val') |>
  filter(FDR_qval <= 0.25) |>
  mutate(NAME = str_remove(NAME, "HALLMARK_")) |>
  mutate(NAME = str_replace(NAME, "_", " ")) |>
  mutate(NAME = str_replace(NAME, "_", " ")) |>
  arrange(-NES) |>
  mutate(NAME = factor(NAME, levels = NAME)) |>
  ggplot(aes(x = NAME, y = NES)) +
  geom_point(aes(size = -NES, color = FDR_qval)) +
  coord_flip() +
  mytheme +
  scale_color_viridis(option = 'B', begin = 0.8, end = 0 ) +
  ggtitle("GSEA upon APC re-expression in SW480 cells") +
  xlab("Hallmark pathways") -> plot

png(filename = file.path(home, "gsea_RNAseq.png"), width = 650, height = 500)
print(plot)
dev.off()

pdf("gsea_RNAseq.pdf", width = 9, height = 6)
print(plot)
dev.off()

From_RiboSeq_dataset |>
  select('NAME', 'NES', 'FDR q-val') |>
  rename(FDR_qval = 'FDR q-val') |>
  filter(FDR_qval <= 0.25) |>
  mutate(NAME = str_remove(NAME, "HALLMARK_")) |>
  mutate(NAME = str_replace(NAME, "_", " ")) |>
  mutate(NAME = str_replace(NAME, "_", " ")) |>
  arrange(-NES) |>
  mutate(NAME = factor(NAME, levels = NAME)) |>
  ggplot(aes(x = NAME, y = NES)) +
  geom_point(aes(size = -NES, color = FDR_qval)) +
  coord_flip() +
  mytheme +
  scale_color_viridis(option = 'B', begin = 0.8, end = 0 ) +
  ggtitle("GSEA upon APC re-expression in SW480 cells") +
  xlab("Hallmark pathways") -> plot

png(filename = file.path(home, "gsea_RNAseq_fromRiboseq_experiment.png"), width = 650, height = 500)
print(plot)
dev.off()

pdf("gsea_RNAseq_fromRiboseq_experiment.pdf", width = 9, height = 6)
print(plot)
dev.off()
