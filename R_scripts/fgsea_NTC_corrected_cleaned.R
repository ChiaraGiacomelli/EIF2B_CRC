#load libraries----
library(tidyverse)
library(fgsea)
library(Glimma)
library(data.table)
library(viridis)
library(conflicted)
library(msigdbr)


#variables----
treatment <- 'plus'
ctrl <- 'min'
source("common_variables.R")
source("read_human_GSEA_pathways.R") # from personal folder

#themes----
mytheme <- theme_classic()+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))

#functions----
run_fgsea <- function(named_vector, pathway) {
  results <- fgsea(pathways = pathway,
                   stats=named_vector,
                   minSize = 20,
                   maxSize = 1000,
                   nPermSimple = 10000)
  return(results)
}

collapse_pathways <- function(fgsea_results, named_vector, fgsea_pathway, padj_threshold) {
  collapsedPathways <- collapsePathways(fgsea_results[order(pval)][padj < padj_threshold], 
                                        fgsea_pathway, named_vector)
  mainPathways <- fgsea_results[pathway %in% collapsedPathways$mainPathways][
    order(-NES), pathway]
  
  return(mainPathways)
}

extract_pathways <- function(fgsea_results, named_vectors, gsea_set, padj) {
  RPF_collapsed_pathways <- collapse_pathways(fgsea_results = fgsea_results[[1]],
                                              named_vector = named_vectors[[1]],
                                              gsea_set,
                                              padj_threshold = padj)
  
  totals_collapsed_pathways <- collapse_pathways(fgsea_results = fgsea_results[[2]],
                                                 named_vector = named_vectors[[2]],
                                                 gsea_set,
                                                 padj_threshold = padj)
  
  TE_collapsed_pathways <- collapse_pathways(fgsea_results = fgsea_results[[3]],
                                             named_vector = named_vectors[[3]],
                                             gsea_set,
                                             padj_threshold = padj)
  
  
  all_pathways <- unique(c(RPF_collapsed_pathways, totals_collapsed_pathways, TE_collapsed_pathways))
  
  return(all_pathways)
}

make_plot <- function(fgsea_result, padj_threshold, title) {
  plot <- ggplot(data = fgsea_result[fgsea_result$padj < padj_threshold], aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = padj)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=title) + 
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=12))
  return(plot)
}

plot_scatters <- function(df, gsea_set, pathway, dir) {
  gene_names <- gsea_set[[pathway]]
  
  df %>%
    mutate(group = factor(gene_sym %in% gene_names),
           alpha_score = case_when(group == T ~ 1,
                                   group == F ~ 0.01)) %>%
    arrange(group) %>%
    ggplot(aes(x = totals_log2FC, y = RPFs_log2FC, colour = group, alpha = alpha_score))+
    geom_point()+
    scale_colour_manual(values=c("grey", shRNA_color))+
    scale_alpha(guide = "none")+
    geom_abline(lty = 2)+
    geom_hline(yintercept = 0, lty = 2)+
    geom_vline(xintercept = 0, lty = 2)+
    ylim(c(-2.5,2.5))+
    xlim(c(-2.5,2.5))+
    mytheme+
    xlab("Total RNA log2FC")+
    ylab("RPFs log2FC")+
    labs(colour = "Genes \nin set")+
    ggtitle(paste(shRNA, str_replace_all(pathway, "_", " "), sep = "\n")) -> scatter_plot
  
  pdf(file = paste0(plot_dir, "scatters/", dir, paste(shRNA, pathway, "NTCcorrected_TE_scatter_plot.pdf", sep = "_")), width = 6, height = 6)
  print(scatter_plot)
  dev.off()
}

############ 2B1  ####
## read in common variables----
source("common_variables_2B1.R")

## read in DESeq2 output----
totals <- read_csv(file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))
RPFs <- read_csv(file = file.path(tables_dir, paste0("RPFs_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))
totals_norm_counts <- read_csv(file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_normalised_counts.csv")))
RPFs_norm_counts <- read_csv(file = file.path(tables_dir, paste0("RPFs_", shRNA, "_NTCcorrected_normalised_counts.csv")))

## select just the Ctrl and treatment normalised counts----
RPFs_norm_counts %>%
  column_to_rownames("gene") %>%
  dplyr::select(matches(c(ctrl, treatment))) %>%
  rownames_to_column("GeneID") -> filtered_RPFs_norm_counts

totals_norm_counts %>%
  column_to_rownames("gene") %>%
  dplyr::select(matches(c(ctrl, treatment))) %>%
  rownames_to_column("GeneID") -> filtered_totals_norm_counts

## merge data from deseq2 output----
RPFs %>%
  dplyr::select(transcript, gene, gene_sym, log2FoldChange, padj) %>%
  dplyr::rename(RPFs_log2FC = log2FoldChange,
                RPFs_padj = padj) %>%
  inner_join(totals[,c("transcript", "log2FoldChange", "padj")], by = "transcript") %>%
  dplyr::rename(totals_log2FC = log2FoldChange,
                totals_padj = padj) %>%
  mutate(TE = RPFs_log2FC - totals_log2FC) %>%
  dplyr::filter(!(is.na(RPFs_padj)) & !(is.na(totals_padj))) -> merged_data

# write for me the list of genes
write_csv(merged_data, file = file.path(tables_dir, paste0("TE_", shRNA, "_NTCcorrected_4fgsea.csv")))

## make named vectors----
merged_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(RPFs_log2FC)) %>%
  deframe() -> RPFs_named_vector

merged_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(totals_log2FC)) %>%
  deframe() -> totals_named_vector

merged_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(TE)) %>%
  deframe() -> TE_named_vector

named_vectors <- list(RPFs_named_vector, totals_named_vector, TE_named_vector)

### hallmark----
#carry out fgsea
hallmark_results <- lapply(named_vectors, run_fgsea, pathway = pathways.hallmark)

#set adjusted p-value
padj <- 0.01

#plot enriched pathways
pdf(file = paste0(plot_dir, shRNA, "_RPFs_hallmark.pdf"), width = 7, height = 10)
make_plot(fgsea_result = hallmark_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Hallmark gene sets"))
dev.off()

pdf(file = paste0(plot_dir, shRNA, "_totals_hallmark.pdf"), width = 7, height = 10)
make_plot(fgsea_result = hallmark_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Hallmark gene sets"))
dev.off()

pdf(file = paste0(plot_dir, shRNA, "_TEs_hallmark.pdf"), width = 7, height = 10)
make_plot(fgsea_result = hallmark_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Hallmark gene sets"))
dev.off()

# extract pathways
# # can loose here the threshold if you want
# padj <- 0.05
# all_pathways <- extract_pathways(fgsea_results = hallmark_results, named_vectors = named_vectors, gsea_set = pathways.hallmark, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/hallmark"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.hallmark, dir = "hallmark")

# #write out table
# df_RPF<-as.data.frame(hallmark_results[[1]])
# fwrite(df_RPF,file=paste(shRNA,'RPF NTC corr HM.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_Tot<-as.data.frame(hallmark_results[[2]])
# fwrite(df_Tot,file=paste(shRNA,'Tot NTC corr HM.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_TE<-as.data.frame(hallmark_results[[3]])
# fwrite(df_TE,file=paste(shRNA,'TE NTC corr HM.tsv',sep = ''),sep='\t',sep2=c('',' ',''))

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.hallmark, df = merged_data, dir = "hallmark")

### biological processes----
#carry out fgsea
bio_processes_results <- lapply(named_vectors, run_fgsea, pathway = pathways.bio_processes)

#set adjusted p-value
padj <- 0.001

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_bio_processes.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = bio_processes_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Biological Processes gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_bio_processes.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = bio_processes_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Biological Processes gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_bio_processes.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = bio_processes_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Biological Processes gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.05
# all_pathways <- extract_pathways(fgsea_results = bio_processes_results, named_vectors = named_vectors, gsea_set = pathways.bio_processes, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/bio_processes"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.bio_processes, dir = "bio_processes")

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.bio_processes, df = merged_data, dir = "bio_processes")

### molecular functions----
#carry out fgsea
mol_funs_results <- lapply(named_vectors, run_fgsea, pathway = pathways.mol_funs)

#set adjusted p-value
padj <- 0.005

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_mol_funs.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = mol_funs_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Molecular Functions gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_mol_funs.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = mol_funs_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Molecular Functions gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_mol_funs.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = mol_funs_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Molecular Functions gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.05
# all_pathways <- extract_pathways(fgsea_results = mol_funs_results, named_vectors = named_vectors, gsea_set = pathways.mol_funs, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/mol_funs"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.mol_funs, dir = "mol_funs")

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.mol_funs, df = merged_data, dir = "mol_funs")

### cellular component----
#carry out fgsea
cell_comp_results <- lapply(named_vectors, run_fgsea, pathway = pathways.cell_comp)

#set adjusted p-value
padj <- 0.01

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_cell_comp.pdf", sep = "_")), width = 10, height = 12)
make_plot(fgsea_result = cell_comp_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Cellular Component gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_cell_comp.pdf", sep = "_")), width = 10, height = 12)
make_plot(fgsea_result = cell_comp_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Cellular Component gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_cell_comp.pdf", sep = "_")), width = 10, height = 12)
make_plot(fgsea_result = cell_comp_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Cellular Component gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.05
# all_pathways <- extract_pathways(fgsea_results = cell_comp_results, named_vectors = named_vectors, gsea_set = pathways.cell_comp, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/cell_comp"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.cell_comp, dir = "cell_comp")
# 
# #write out table
# df_RPF<-as.data.frame(cell_comp_results[[1]])
# fwrite(df_RPF,file=paste(shRNA,'RPF NTC corr cell component.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_Tot<-as.data.frame(cell_comp_results[[2]])
# fwrite(df_Tot,file=paste(shRNA,'Tot NTC corr cell component.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_TE<-as.data.frame(cell_comp_results[[3]])
# fwrite(df_TE,file=paste(shRNA,'TE NTC corr cell component.tsv',sep = ''),sep='\t',sep2=c('',' ',''))

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.cell_comp, df = merged_data, dir = "cell_comp")

### onco----
#carry out fgsea
onco_results <- lapply(named_vectors, run_fgsea, pathway = pathways.onco)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_onco.pdf", sep = "_")), width = 10, height = 7)
make_plot(fgsea_result = onco_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Oncogenic Signature gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_onco.pdf", sep = "_")), width = 10, height = 7)
make_plot(fgsea_result = onco_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Oncogenic Signature gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_onco.pdf", sep = "_")), width = 10, height = 5)
make_plot(fgsea_result = onco_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Oncogenic Signature gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.1
# all_pathways <- extract_pathways(fgsea_results = onco_results, named_vectors = named_vectors, gsea_set = pathways.onco, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/onco"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.onco, dir = "onco")

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.onco, df = merged_data, dir = "onco")

### kegg----
#carry out fgsea
kegg_results <- lapply(named_vectors, run_fgsea, pathway = pathways.kegg)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_kegg.pdf", sep = "_")), width = 7, height = 7)
make_plot(fgsea_result = kegg_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Kegg gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_kegg.pdf", sep = "_")), width = 7, height = 7)
make_plot(fgsea_result = kegg_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Kegg gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_kegg.pdf", sep = "_")), width = 7, height = 7)
make_plot(fgsea_result = kegg_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Kegg gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.1
# all_pathways <- extract_pathways(fgsea_results = kegg_results, named_vectors = named_vectors, gsea_set = pathways.kegg, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/kegg"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.kegg, dir = "kegg")

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.kegg, df = merged_data, dir = "kegg")

#write out table
# df_RPF<-as.data.frame(kegg_results[[1]])
# fwrite(df_RPF,file=paste(shRNA,'RPF NTC corr KEGG.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_Tot<-as.data.frame(kegg_results[[2]])
# fwrite(df_Tot,file=paste(shRNA,'Tot NTC corr KEGG.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_TE<-as.data.frame(kegg_results[[3]])
# fwrite(df_TE,file=paste(shRNA,'TE NTC corr KEGG.tsv',sep = ''),sep='\t',sep2=c('',' ',''))

## plot enrichments as people seem to like them -----
# Integrated stress response is part of GOBP
# pathways <- gmtPathways("~/data/CGIACOME/Protocols_and_Documentation/Sequencing/Gene_signatures/c5.go.bp.v2023.1.Hs.symbols.gmt")
# ISR <- pathways[["GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_ISR.png")), width = 600, height = 400)
# print(plotEnrichment(ISR, TE_named_vector)+labs(title="ISR TE 2B1", subtitle = "GO:BP Integrated Stress Response"))
# dev.off()
# 
# res <- fgsea(pathways = pathways, stats = TE_named_vector, minSize = 15, maxSize = 1000)
# res %>%
#   dplyr::filter(pathway == "GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING") -> ISR_res
# ISR_res <- as_tibble(ISR_res)
# 
# pathways <- gmtPathways("~/data/CGIACOME/Protocols_and_Documentation/Sequencing/Gene_signatures/c2.cgp.v2023.1.Hs.symbols.gmt")
# ATF4_UP <- pathways[["IGARASHI_ATF4_TARGETS_UP"]]
# ATF4_DN <- pathways[["IGARASHI_ATF4_TARGETS_DN"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_ATF4_UP.png")), width = 600, height = 400)
# print(plotEnrichment(ATF4_UP, TE_named_vector)+labs(title="ATF4 targets UP TE 2B1", subtitle = "genes upregulated after ATF4 KD in A549"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_ATF4_DN.png")), width = 600, height = 400)
# print(plotEnrichment(ATF4_DN, TE_named_vector)+labs(title="ATF4 targets DOWN TE 2B1", subtitle = "genes downregulated after ATF4 KD in A549"))
# dev.off()
# 
# # sansom signatures - also part of the C2 CGP set
# OS_1 <- pathways[["SANSOM_APC_MYC_TARGETS"]]
# OS_2 <- pathways[["SANSOM_APC_TARGETS"]]
# OS_3 <- pathways[["SANSOM_APC_TARGETS_DN"]]
# OS_4 <- pathways[["SANSOM_APC_TARGETS_REQUIRE_MYC"]]
# OS_5 <- pathways[["SANSOM_APC_TARGETS_UP"]]
# OS_6 <- pathways[["SANSOM_WNT_PATHWAY_REQUIRE_MYC"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_1.png")), width = 600, height = 400)
# print(plotEnrichment(OS_1, TE_named_vector)+labs(title="SANSOM APC MYC TARGETS 2B1 TE", subtitle = "genes downregulated after APC and MYC KO"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_2.png")), width = 600, height = 400)
# print(plotEnrichment(OS_2, TE_named_vector)+labs(title="SANSOM APC TARGETS 2B1 TE", subtitle = "genes upregulated after APC KO"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_3.png")), width = 600, height = 400)
# print(plotEnrichment(OS_3, TE_named_vector)+labs(title="SANSOM APC TARGETS DOWN 2B1 TE", subtitle = "top genes downregulated after APC KO (5 days pi)"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_4.png")), width = 600, height = 400)
# print(plotEnrichment(OS_4, TE_named_vector)+labs(title="SANSOM APC TARGETS REQUIRING MYC 2B1 TE", subtitle = "genes upregulated after APC KO requiring functional MYC"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_5.png")), width = 600, height = 400)
# print(plotEnrichment(OS_5, TE_named_vector)+labs(title="SANSOM APC TARGETS UP 2B1 TE", subtitle = "top genes upregulated after APC KO (5 days pi)"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_6.png")), width = 600, height = 400)
# print(plotEnrichment(OS_6, TE_named_vector)+labs(title="SANSOM WNT PATHWAY REQUIRING MYC 2B1 TE", subtitle = "Wnt-target genes upregulated after APC and requiring functional MYC"))
# dev.off()
# 
# LIN_APC <- pathways[["LIN_APC_TARGETS"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_Lin_APC.png")), width = 600, height = 400)
# print(plotEnrichment(LIN_APC, TE_named_vector)+labs(title="LIN APC TARGETS 2B1 TE", subtitle = "Genes up-regulated by forced APC expression in SW480"))
# dev.off()
# 
# # get results out (NES/adjp/LeadingEdge)
# res <- fgsea(pathways = pathways, stats = TE_named_vector, minSize = 3, maxSize = 1000)
# res %>%
#   dplyr::filter(pathway == "IGARASHI_ATF4_TARGETS_UP") -> AU_res
# AU_res <- as_tibble(AU_res)
# 
# res %>%
#   dplyr::filter(pathway == "IGARASHI_ATF4_TARGETS_DN") -> AD_res
# AD_res <- as_tibble(AD_res)
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_MYC_TARGETS") -> r1
#                 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_TARGETS") -> r2
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_TARGETS_DN") -> r3
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_TARGETS_REQUIRE_MYC") -> r4
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_TARGETS_UP") -> r5
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_WNT_PATHWAY_REQUIRE_MYC") -> r6
# 
# res %>%
#   dplyr::filter(pathway == "LIN_APC_TARGETS") -> r_lin
# 
# # 2004 OS signature
# pathways <- gmtPathways("~/data/CGIACOME/Protocols_and_Documentation/Sequencing/Gene_signatures/20230821_Sansom_Apc.Hs.symbols.gmt")
# OS_7 <- pathways[["SANSOM_2004_APC_KO_UP"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_Sansom_2004_APC.png")), width = 600, height = 400)
# print(plotEnrichment(OS_7, TE_named_vector)+labs(title="SANSOM 2004 2B1 TE", subtitle = "Genes upregulated upon APC KO"))
# dev.off()
# # Gives almost same distribution as OS_5 signature
# 
# # get results out (NES/adjp/LeadingEdge)
# res <- fgsea(pathways = pathways, stats = TE_named_vector, minSize = 3, maxSize = 1000)
# res %>%
#   dplyr::filter(pathway == "SANSOM_2004_APC_KO_UP") -> r7
# 
# # get all results out as a data table
# ISR_res %>%
#   bind_rows(AU_res, AD_res, r1, r2, r3, r4, r5, r6, r_lin, r7) %>%
#   unnest(leadingEdge)-> full_res
# 
# # write it out
# write_csv(full_res, file = file.path(parent_dir, "Analysis/fgsea", paste0("TE_", shRNA, "_fgsea_curated_ofInterest.csv")))
# 
# 2B4  ####
## read in common variables----
source("common_variables_2B4.R")

## read in DESeq2 output----
totals <- read_csv(file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))
RPFs <- read_csv(file = file.path(tables_dir, paste0("RPFs_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))
totals_norm_counts <- read_csv(file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_normalised_counts.csv")))
RPFs_norm_counts <- read_csv(file = file.path(tables_dir, paste0("RPFs_", shRNA, "_NTCcorrected_normalised_counts.csv")))

## select just the Ctrl and treatment normalised counts----
RPFs_norm_counts %>%
  column_to_rownames("gene") %>%
  dplyr::select(matches(c(ctrl, treatment))) %>%
  rownames_to_column("GeneID") -> filtered_RPFs_norm_counts

totals_norm_counts %>%
  column_to_rownames("gene") %>%
  dplyr::select(matches(c(ctrl, treatment))) %>%
  rownames_to_column("GeneID") -> filtered_totals_norm_counts

## merge data from deseq2 output----
RPFs %>%
  dplyr::select(transcript, gene, gene_sym, log2FoldChange, padj) %>%
  dplyr::rename(RPFs_log2FC = log2FoldChange,
                RPFs_padj = padj) %>%
  inner_join(totals[,c("transcript", "log2FoldChange", "padj")], by = "transcript") %>%
  dplyr::rename(totals_log2FC = log2FoldChange,
                totals_padj = padj) %>%
  mutate(TE = RPFs_log2FC - totals_log2FC) %>%
  dplyr::filter(!(is.na(RPFs_padj)) & !(is.na(totals_padj))) -> merged_data

# write for me the list of genes
write_csv(merged_data, file = file.path(tables_dir, paste0("TE_", shRNA, "_NTCcorrected_4fgsea.csv")))

## make named vectors----
merged_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(RPFs_log2FC)) %>%
  deframe() -> RPFs_named_vector

merged_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(totals_log2FC)) %>%
  deframe() -> totals_named_vector

merged_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(TE)) %>%
  deframe() -> TE_named_vector

named_vectors <- list(RPFs_named_vector, totals_named_vector, TE_named_vector)

### hallmark----
#carry out fgsea
hallmark_results <- lapply(named_vectors, run_fgsea, pathway = pathways.hallmark)

#set adjusted p-value
padj <- 0.01

#plot enriched pathways
pdf(file = paste0(plot_dir, shRNA, "_RPFs_hallmark.pdf"), width = 7, height = 10)
make_plot(fgsea_result = hallmark_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Hallmark gene sets"))
dev.off()

pdf(file = paste0(plot_dir, shRNA, "_totals_hallmark.pdf"), width = 7, height = 10)
make_plot(fgsea_result = hallmark_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Hallmark gene sets"))
dev.off()

pdf(file = paste0(plot_dir, shRNA, "_TEs_hallmark.pdf"), width = 7, height = 10)
make_plot(fgsea_result = hallmark_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Hallmark gene sets"))
dev.off()

# extract pathways
# # can loose here the threshold if you want
# padj <- 0.05
# all_pathways <- extract_pathways(fgsea_results = hallmark_results, named_vectors = named_vectors, gsea_set = pathways.hallmark, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/hallmark"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.hallmark, dir = "hallmark")

# #write out table
# df_RPF<-as.data.frame(hallmark_results[[1]])
# fwrite(df_RPF,file=paste(shRNA,'RPF NTC corr HM.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_Tot<-as.data.frame(hallmark_results[[2]])
# fwrite(df_Tot,file=paste(shRNA,'Tot NTC corr HM.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_TE<-as.data.frame(hallmark_results[[3]])
# fwrite(df_TE,file=paste(shRNA,'TE NTC corr HM.tsv',sep = ''),sep='\t',sep2=c('',' ',''))

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.hallmark, df = merged_data, dir = "hallmark")

### biological processes----
#carry out fgsea
bio_processes_results <- lapply(named_vectors, run_fgsea, pathway = pathways.bio_processes)

#set adjusted p-value
padj <- 0.001

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_bio_processes.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = bio_processes_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Biological Processes gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_bio_processes.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = bio_processes_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Biological Processes gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_bio_processes.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = bio_processes_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Biological Processes gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.05
# all_pathways <- extract_pathways(fgsea_results = bio_processes_results, named_vectors = named_vectors, gsea_set = pathways.bio_processes, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/bio_processes"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.bio_processes, dir = "bio_processes")

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.bio_processes, df = merged_data, dir = "bio_processes")

### molecular functions----
#carry out fgsea
mol_funs_results <- lapply(named_vectors, run_fgsea, pathway = pathways.mol_funs)

#set adjusted p-value
padj <- 0.005

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_mol_funs.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = mol_funs_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Molecular Functions gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_mol_funs.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = mol_funs_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Molecular Functions gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_mol_funs.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = mol_funs_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Molecular Functions gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.05
# all_pathways <- extract_pathways(fgsea_results = mol_funs_results, named_vectors = named_vectors, gsea_set = pathways.mol_funs, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/mol_funs"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.mol_funs, dir = "mol_funs")

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.mol_funs, df = merged_data, dir = "mol_funs")

### cellular component----
#carry out fgsea
cell_comp_results <- lapply(named_vectors, run_fgsea, pathway = pathways.cell_comp)

#set adjusted p-value
padj <- 0.01

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_cell_comp.pdf", sep = "_")), width = 10, height = 12)
make_plot(fgsea_result = cell_comp_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Cellular Component gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_cell_comp.pdf", sep = "_")), width = 10, height = 12)
make_plot(fgsea_result = cell_comp_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Cellular Component gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_cell_comp.pdf", sep = "_")), width = 10, height = 12)
make_plot(fgsea_result = cell_comp_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Cellular Component gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.05
# all_pathways <- extract_pathways(fgsea_results = cell_comp_results, named_vectors = named_vectors, gsea_set = pathways.cell_comp, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/cell_comp"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.cell_comp, dir = "cell_comp")
# 
# #write out table
# df_RPF<-as.data.frame(cell_comp_results[[1]])
# fwrite(df_RPF,file=paste(shRNA,'RPF NTC corr cell component.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_Tot<-as.data.frame(cell_comp_results[[2]])
# fwrite(df_Tot,file=paste(shRNA,'Tot NTC corr cell component.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_TE<-as.data.frame(cell_comp_results[[3]])
# fwrite(df_TE,file=paste(shRNA,'TE NTC corr cell component.tsv',sep = ''),sep='\t',sep2=c('',' ',''))

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.cell_comp, df = merged_data, dir = "cell_comp")

### onco----
#carry out fgsea
onco_results <- lapply(named_vectors, run_fgsea, pathway = pathways.onco)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_onco.pdf", sep = "_")), width = 10, height = 7)
make_plot(fgsea_result = onco_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Oncogenic Signature gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_onco.pdf", sep = "_")), width = 10, height = 7)
make_plot(fgsea_result = onco_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Oncogenic Signature gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_onco.pdf", sep = "_")), width = 10, height = 5)
make_plot(fgsea_result = onco_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Oncogenic Signature gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.1
# all_pathways <- extract_pathways(fgsea_results = onco_results, named_vectors = named_vectors, gsea_set = pathways.onco, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/onco"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.onco, dir = "onco")

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.onco, df = merged_data, dir = "onco")

### kegg----
#carry out fgsea
kegg_results <- lapply(named_vectors, run_fgsea, pathway = pathways.kegg)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_kegg.pdf", sep = "_")), width = 7, height = 7)
make_plot(fgsea_result = kegg_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Kegg gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_kegg.pdf", sep = "_")), width = 7, height = 7)
make_plot(fgsea_result = kegg_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Kegg gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_kegg.pdf", sep = "_")), width = 7, height = 7)
make_plot(fgsea_result = kegg_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Kegg gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.1
# all_pathways <- extract_pathways(fgsea_results = kegg_results, named_vectors = named_vectors, gsea_set = pathways.kegg, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/kegg"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.kegg, dir = "kegg")

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.kegg, df = merged_data, dir = "kegg")

#write out table
# df_RPF<-as.data.frame(kegg_results[[1]])
# fwrite(df_RPF,file=paste(shRNA,'RPF NTC corr KEGG.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_Tot<-as.data.frame(kegg_results[[2]])
# fwrite(df_Tot,file=paste(shRNA,'Tot NTC corr KEGG.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_TE<-as.data.frame(kegg_results[[3]])
# fwrite(df_TE,file=paste(shRNA,'TE NTC corr KEGG.tsv',sep = ''),sep='\t',sep2=c('',' ',''))

## plot enrichments as people seem to like them -----
# Integrated stress response is part of GOBP
# pathways <- gmtPathways("~/data/CGIACOME/Protocols_and_Documentation/Sequencing/Gene_signatures/c5.go.bp.v2023.1.Hs.symbols.gmt")
# ISR <- pathways[["GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_ISR.png")), width = 600, height = 400)
# print(plotEnrichment(ISR, TE_named_vector)+labs(title="ISR TE 2B1", subtitle = "GO:BP Integrated Stress Response"))
# dev.off()
# 
# res <- fgsea(pathways = pathways, stats = TE_named_vector, minSize = 15, maxSize = 1000)
# res %>%
#   dplyr::filter(pathway == "GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING") -> ISR_res
# ISR_res <- as_tibble(ISR_res)
# 
# pathways <- gmtPathways("~/data/CGIACOME/Protocols_and_Documentation/Sequencing/Gene_signatures/c2.cgp.v2023.1.Hs.symbols.gmt")
# ATF4_UP <- pathways[["IGARASHI_ATF4_TARGETS_UP"]]
# ATF4_DN <- pathways[["IGARASHI_ATF4_TARGETS_DN"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_ATF4_UP.png")), width = 600, height = 400)
# print(plotEnrichment(ATF4_UP, TE_named_vector)+labs(title="ATF4 targets UP TE 2B1", subtitle = "genes upregulated after ATF4 KD in A549"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_ATF4_DN.png")), width = 600, height = 400)
# print(plotEnrichment(ATF4_DN, TE_named_vector)+labs(title="ATF4 targets DOWN TE 2B1", subtitle = "genes downregulated after ATF4 KD in A549"))
# dev.off()
# 
# # sansom signatures - also part of the C2 CGP set
# OS_1 <- pathways[["SANSOM_APC_MYC_TARGETS"]]
# OS_2 <- pathways[["SANSOM_APC_TARGETS"]]
# OS_3 <- pathways[["SANSOM_APC_TARGETS_DN"]]
# OS_4 <- pathways[["SANSOM_APC_TARGETS_REQUIRE_MYC"]]
# OS_5 <- pathways[["SANSOM_APC_TARGETS_UP"]]
# OS_6 <- pathways[["SANSOM_WNT_PATHWAY_REQUIRE_MYC"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_1.png")), width = 600, height = 400)
# print(plotEnrichment(OS_1, TE_named_vector)+labs(title="SANSOM APC MYC TARGETS 2B1 TE", subtitle = "genes downregulated after APC and MYC KO"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_2.png")), width = 600, height = 400)
# print(plotEnrichment(OS_2, TE_named_vector)+labs(title="SANSOM APC TARGETS 2B1 TE", subtitle = "genes upregulated after APC KO"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_3.png")), width = 600, height = 400)
# print(plotEnrichment(OS_3, TE_named_vector)+labs(title="SANSOM APC TARGETS DOWN 2B1 TE", subtitle = "top genes downregulated after APC KO (5 days pi)"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_4.png")), width = 600, height = 400)
# print(plotEnrichment(OS_4, TE_named_vector)+labs(title="SANSOM APC TARGETS REQUIRING MYC 2B1 TE", subtitle = "genes upregulated after APC KO requiring functional MYC"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_5.png")), width = 600, height = 400)
# print(plotEnrichment(OS_5, TE_named_vector)+labs(title="SANSOM APC TARGETS UP 2B1 TE", subtitle = "top genes upregulated after APC KO (5 days pi)"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_6.png")), width = 600, height = 400)
# print(plotEnrichment(OS_6, TE_named_vector)+labs(title="SANSOM WNT PATHWAY REQUIRING MYC 2B1 TE", subtitle = "Wnt-target genes upregulated after APC and requiring functional MYC"))
# dev.off()
# 
# LIN_APC <- pathways[["LIN_APC_TARGETS"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_Lin_APC.png")), width = 600, height = 400)
# print(plotEnrichment(LIN_APC, TE_named_vector)+labs(title="LIN APC TARGETS 2B1 TE", subtitle = "Genes up-regulated by forced APC expression in SW480"))
# dev.off()
# 
# # get results out (NES/adjp/LeadingEdge)
# res <- fgsea(pathways = pathways, stats = TE_named_vector, minSize = 3, maxSize = 1000)
# res %>%
#   dplyr::filter(pathway == "IGARASHI_ATF4_TARGETS_UP") -> AU_res
# AU_res <- as_tibble(AU_res)
# 
# res %>%
#   dplyr::filter(pathway == "IGARASHI_ATF4_TARGETS_DN") -> AD_res
# AD_res <- as_tibble(AD_res)
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_MYC_TARGETS") -> r1
#                 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_TARGETS") -> r2
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_TARGETS_DN") -> r3
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_TARGETS_REQUIRE_MYC") -> r4
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_TARGETS_UP") -> r5
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_WNT_PATHWAY_REQUIRE_MYC") -> r6
# 
# res %>%
#   dplyr::filter(pathway == "LIN_APC_TARGETS") -> r_lin
# 
# # 2004 OS signature
# pathways <- gmtPathways("~/data/CGIACOME/Protocols_and_Documentation/Sequencing/Gene_signatures/20230821_Sansom_Apc.Hs.symbols.gmt")
# OS_7 <- pathways[["SANSOM_2004_APC_KO_UP"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_Sansom_2004_APC.png")), width = 600, height = 400)
# print(plotEnrichment(OS_7, TE_named_vector)+labs(title="SANSOM 2004 2B1 TE", subtitle = "Genes upregulated upon APC KO"))
# dev.off()
# # Gives almost same distribution as OS_5 signature
# 
# # get results out (NES/adjp/LeadingEdge)
# res <- fgsea(pathways = pathways, stats = TE_named_vector, minSize = 3, maxSize = 1000)
# res %>%
#   dplyr::filter(pathway == "SANSOM_2004_APC_KO_UP") -> r7
# 
# # get all results out as a data table
# ISR_res %>%
#   bind_rows(AU_res, AD_res, r1, r2, r3, r4, r5, r6, r_lin, r7) %>%
#   unnest(leadingEdge)-> full_res
# 
# # write it out
# write_csv(full_res, file = file.path(parent_dir, "Analysis/fgsea", paste0("TE_", shRNA, "_fgsea_curated_ofInterest.csv")))
# 
# 2B5  ####
## read in common variables----
source("common_variables_2B5.R")

## read in DESeq2 output----
totals <- read_csv(file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))
RPFs <- read_csv(file = file.path(tables_dir, paste0("RPFs_", shRNA, "_NTCcorrected_DEseq2_apeglm_LFC_shrinkage.csv")))
totals_norm_counts <- read_csv(file = file.path(tables_dir, paste0("Totals_", shRNA, "_NTCcorrected_normalised_counts.csv")))
RPFs_norm_counts <- read_csv(file = file.path(tables_dir, paste0("RPFs_", shRNA, "_NTCcorrected_normalised_counts.csv")))

## select just the Ctrl and treatment normalised counts----
RPFs_norm_counts %>%
  column_to_rownames("gene") %>%
  dplyr::select(matches(c(ctrl, treatment))) %>%
  rownames_to_column("GeneID") -> filtered_RPFs_norm_counts

totals_norm_counts %>%
  column_to_rownames("gene") %>%
  dplyr::select(matches(c(ctrl, treatment))) %>%
  rownames_to_column("GeneID") -> filtered_totals_norm_counts

## merge data from deseq2 output----
RPFs %>%
  dplyr::select(transcript, gene, gene_sym, log2FoldChange, padj) %>%
  dplyr::rename(RPFs_log2FC = log2FoldChange,
                RPFs_padj = padj) %>%
  inner_join(totals[,c("transcript", "log2FoldChange", "padj")], by = "transcript") %>%
  dplyr::rename(totals_log2FC = log2FoldChange,
                totals_padj = padj) %>%
  mutate(TE = RPFs_log2FC - totals_log2FC) %>%
  dplyr::filter(!(is.na(RPFs_padj)) & !(is.na(totals_padj))) -> merged_data

# write for me the list of genes
write_csv(merged_data, file = file.path(tables_dir, paste0("TE_", shRNA, "_NTCcorrected_4fgsea.csv")))

## make named vectors----
merged_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(RPFs_log2FC)) %>%
  deframe() -> RPFs_named_vector

merged_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(totals_log2FC)) %>%
  deframe() -> totals_named_vector

merged_data %>%
  group_by(gene_sym) %>%
  summarise(stat = mean(TE)) %>%
  deframe() -> TE_named_vector

named_vectors <- list(RPFs_named_vector, totals_named_vector, TE_named_vector)

### hallmark----
#carry out fgsea
hallmark_results <- lapply(named_vectors, run_fgsea, pathway = pathways.hallmark)

#set adjusted p-value
padj <- 0.01

#plot enriched pathways
pdf(file = paste0(plot_dir, shRNA, "_RPFs_hallmark.pdf"), width = 7, height = 10)
make_plot(fgsea_result = hallmark_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Hallmark gene sets"))
dev.off()

pdf(file = paste0(plot_dir, shRNA, "_totals_hallmark.pdf"), width = 7, height = 10)
make_plot(fgsea_result = hallmark_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Hallmark gene sets"))
dev.off()

pdf(file = paste0(plot_dir, shRNA, "_TEs_hallmark.pdf"), width = 7, height = 10)
make_plot(fgsea_result = hallmark_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Hallmark gene sets"))
dev.off()

# extract pathways
# # can loose here the threshold if you want
# padj <- 0.05
# all_pathways <- extract_pathways(fgsea_results = hallmark_results, named_vectors = named_vectors, gsea_set = pathways.hallmark, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/hallmark"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.hallmark, dir = "hallmark")

# #write out table
# df_RPF<-as.data.frame(hallmark_results[[1]])
# fwrite(df_RPF,file=paste(shRNA,'RPF NTC corr HM.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_Tot<-as.data.frame(hallmark_results[[2]])
# fwrite(df_Tot,file=paste(shRNA,'Tot NTC corr HM.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_TE<-as.data.frame(hallmark_results[[3]])
# fwrite(df_TE,file=paste(shRNA,'TE NTC corr HM.tsv',sep = ''),sep='\t',sep2=c('',' ',''))

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.hallmark, df = merged_data, dir = "hallmark")

### biological processes----
#carry out fgsea
bio_processes_results <- lapply(named_vectors, run_fgsea, pathway = pathways.bio_processes)

#set adjusted p-value
padj <- 0.001

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_bio_processes.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = bio_processes_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Biological Processes gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_bio_processes.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = bio_processes_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Biological Processes gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_bio_processes.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = bio_processes_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Biological Processes gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.05
# all_pathways <- extract_pathways(fgsea_results = bio_processes_results, named_vectors = named_vectors, gsea_set = pathways.bio_processes, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/bio_processes"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.bio_processes, dir = "bio_processes")

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.bio_processes, df = merged_data, dir = "bio_processes")

### molecular functions----
#carry out fgsea
mol_funs_results <- lapply(named_vectors, run_fgsea, pathway = pathways.mol_funs)

#set adjusted p-value
padj <- 0.005

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_mol_funs.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = mol_funs_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Molecular Functions gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_mol_funs.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = mol_funs_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Molecular Functions gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_mol_funs.pdf", sep = "_")), width = 10, height = 10)
make_plot(fgsea_result = mol_funs_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Molecular Functions gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.05
# all_pathways <- extract_pathways(fgsea_results = mol_funs_results, named_vectors = named_vectors, gsea_set = pathways.mol_funs, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/mol_funs"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.mol_funs, dir = "mol_funs")

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.mol_funs, df = merged_data, dir = "mol_funs")

### cellular component----
#carry out fgsea
cell_comp_results <- lapply(named_vectors, run_fgsea, pathway = pathways.cell_comp)

#set adjusted p-value
padj <- 0.01

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_cell_comp.pdf", sep = "_")), width = 10, height = 12)
make_plot(fgsea_result = cell_comp_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Cellular Component gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_cell_comp.pdf", sep = "_")), width = 10, height = 12)
make_plot(fgsea_result = cell_comp_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Cellular Component gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_cell_comp.pdf", sep = "_")), width = 10, height = 12)
make_plot(fgsea_result = cell_comp_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Cellular Component gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.05
# all_pathways <- extract_pathways(fgsea_results = cell_comp_results, named_vectors = named_vectors, gsea_set = pathways.cell_comp, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/cell_comp"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.cell_comp, dir = "cell_comp")
# 
# #write out table
# df_RPF<-as.data.frame(cell_comp_results[[1]])
# fwrite(df_RPF,file=paste(shRNA,'RPF NTC corr cell component.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_Tot<-as.data.frame(cell_comp_results[[2]])
# fwrite(df_Tot,file=paste(shRNA,'Tot NTC corr cell component.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_TE<-as.data.frame(cell_comp_results[[3]])
# fwrite(df_TE,file=paste(shRNA,'TE NTC corr cell component.tsv',sep = ''),sep='\t',sep2=c('',' ',''))

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.cell_comp, df = merged_data, dir = "cell_comp")

### onco----
#carry out fgsea
onco_results <- lapply(named_vectors, run_fgsea, pathway = pathways.onco)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_onco.pdf", sep = "_")), width = 10, height = 7)
make_plot(fgsea_result = onco_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Oncogenic Signature gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_onco.pdf", sep = "_")), width = 10, height = 7)
make_plot(fgsea_result = onco_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Oncogenic Signature gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_onco.pdf", sep = "_")), width = 10, height = 5)
make_plot(fgsea_result = onco_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Oncogenic Signature gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.1
# all_pathways <- extract_pathways(fgsea_results = onco_results, named_vectors = named_vectors, gsea_set = pathways.onco, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/onco"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.onco, dir = "onco")

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.onco, df = merged_data, dir = "onco")

### kegg----
#carry out fgsea
kegg_results <- lapply(named_vectors, run_fgsea, pathway = pathways.kegg)

#set adjusted p-value
padj <- 0.05

#plot enriched pathways
pdf(file = paste0(plot_dir, paste(shRNA, "RPFs_kegg.pdf", sep = "_")), width = 7, height = 7)
make_plot(fgsea_result = kegg_results[[1]], padj_threshold = padj, title = paste(shRNA, "RPFs\nGSEA Kegg gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "totals_kegg.pdf", sep = "_")), width = 7, height = 7)
make_plot(fgsea_result = kegg_results[[2]], padj_threshold = padj, title = paste(shRNA, "Total RNA\nGSEA Kegg gene sets"))
dev.off()

pdf(file = paste0(plot_dir, paste(shRNA, "TE_kegg.pdf", sep = "_")), width = 7, height = 7)
make_plot(fgsea_result = kegg_results[[3]], padj_threshold = padj, title = paste(shRNA, "TE\nGSEA Kegg gene sets"))
dev.off()

#extract pathways
# can loose here the threshold if you want
# padj <- 0.1
# all_pathways <- extract_pathways(fgsea_results = kegg_results, named_vectors = named_vectors, gsea_set = pathways.kegg, padj = padj)
# 
# #plot overlaid scatters
# dir.create(file.path(plot_dir, "scatters/kegg"))
# lapply(all_pathways, plot_scatters, df = merged_data, gsea_set = pathways.kegg, dir = "kegg")

#plot interactive scatters
#lapply(all_pathways, make_interactive_scatter, gsea_set = pathways.kegg, df = merged_data, dir = "kegg")

#write out table
# df_RPF<-as.data.frame(kegg_results[[1]])
# fwrite(df_RPF,file=paste(shRNA,'RPF NTC corr KEGG.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_Tot<-as.data.frame(kegg_results[[2]])
# fwrite(df_Tot,file=paste(shRNA,'Tot NTC corr KEGG.tsv',sep = ''),sep='\t',sep2=c('',' ',''))
# 
# df_TE<-as.data.frame(kegg_results[[3]])
# fwrite(df_TE,file=paste(shRNA,'TE NTC corr KEGG.tsv',sep = ''),sep='\t',sep2=c('',' ',''))

## plot enrichments as people seem to like them -----
# Integrated stress response is part of GOBP
# pathways <- gmtPathways("~/data/CGIACOME/Protocols_and_Documentation/Sequencing/Gene_signatures/c5.go.bp.v2023.1.Hs.symbols.gmt")
# ISR <- pathways[["GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_ISR.png")), width = 600, height = 400)
# print(plotEnrichment(ISR, TE_named_vector)+labs(title="ISR TE 2B1", subtitle = "GO:BP Integrated Stress Response"))
# dev.off()
# 
# res <- fgsea(pathways = pathways, stats = TE_named_vector, minSize = 15, maxSize = 1000)
# res %>%
#   dplyr::filter(pathway == "GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING") -> ISR_res
# ISR_res <- as_tibble(ISR_res)
# 
# pathways <- gmtPathways("~/data/CGIACOME/Protocols_and_Documentation/Sequencing/Gene_signatures/c2.cgp.v2023.1.Hs.symbols.gmt")
# ATF4_UP <- pathways[["IGARASHI_ATF4_TARGETS_UP"]]
# ATF4_DN <- pathways[["IGARASHI_ATF4_TARGETS_DN"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_ATF4_UP.png")), width = 600, height = 400)
# print(plotEnrichment(ATF4_UP, TE_named_vector)+labs(title="ATF4 targets UP TE 2B1", subtitle = "genes upregulated after ATF4 KD in A549"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_ATF4_DN.png")), width = 600, height = 400)
# print(plotEnrichment(ATF4_DN, TE_named_vector)+labs(title="ATF4 targets DOWN TE 2B1", subtitle = "genes downregulated after ATF4 KD in A549"))
# dev.off()
# 
# # sansom signatures - also part of the C2 CGP set
# OS_1 <- pathways[["SANSOM_APC_MYC_TARGETS"]]
# OS_2 <- pathways[["SANSOM_APC_TARGETS"]]
# OS_3 <- pathways[["SANSOM_APC_TARGETS_DN"]]
# OS_4 <- pathways[["SANSOM_APC_TARGETS_REQUIRE_MYC"]]
# OS_5 <- pathways[["SANSOM_APC_TARGETS_UP"]]
# OS_6 <- pathways[["SANSOM_WNT_PATHWAY_REQUIRE_MYC"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_1.png")), width = 600, height = 400)
# print(plotEnrichment(OS_1, TE_named_vector)+labs(title="SANSOM APC MYC TARGETS 2B1 TE", subtitle = "genes downregulated after APC and MYC KO"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_2.png")), width = 600, height = 400)
# print(plotEnrichment(OS_2, TE_named_vector)+labs(title="SANSOM APC TARGETS 2B1 TE", subtitle = "genes upregulated after APC KO"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_3.png")), width = 600, height = 400)
# print(plotEnrichment(OS_3, TE_named_vector)+labs(title="SANSOM APC TARGETS DOWN 2B1 TE", subtitle = "top genes downregulated after APC KO (5 days pi)"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_4.png")), width = 600, height = 400)
# print(plotEnrichment(OS_4, TE_named_vector)+labs(title="SANSOM APC TARGETS REQUIRING MYC 2B1 TE", subtitle = "genes upregulated after APC KO requiring functional MYC"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_5.png")), width = 600, height = 400)
# print(plotEnrichment(OS_5, TE_named_vector)+labs(title="SANSOM APC TARGETS UP 2B1 TE", subtitle = "top genes upregulated after APC KO (5 days pi)"))
# dev.off()
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_OS_6.png")), width = 600, height = 400)
# print(plotEnrichment(OS_6, TE_named_vector)+labs(title="SANSOM WNT PATHWAY REQUIRING MYC 2B1 TE", subtitle = "Wnt-target genes upregulated after APC and requiring functional MYC"))
# dev.off()
# 
# LIN_APC <- pathways[["LIN_APC_TARGETS"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_Lin_APC.png")), width = 600, height = 400)
# print(plotEnrichment(LIN_APC, TE_named_vector)+labs(title="LIN APC TARGETS 2B1 TE", subtitle = "Genes up-regulated by forced APC expression in SW480"))
# dev.off()
# 
# # get results out (NES/adjp/LeadingEdge)
# res <- fgsea(pathways = pathways, stats = TE_named_vector, minSize = 3, maxSize = 1000)
# res %>%
#   dplyr::filter(pathway == "IGARASHI_ATF4_TARGETS_UP") -> AU_res
# AU_res <- as_tibble(AU_res)
# 
# res %>%
#   dplyr::filter(pathway == "IGARASHI_ATF4_TARGETS_DN") -> AD_res
# AD_res <- as_tibble(AD_res)
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_MYC_TARGETS") -> r1
#                 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_TARGETS") -> r2
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_TARGETS_DN") -> r3
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_TARGETS_REQUIRE_MYC") -> r4
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_APC_TARGETS_UP") -> r5
# 
# res %>%
#   dplyr::filter(pathway == "SANSOM_WNT_PATHWAY_REQUIRE_MYC") -> r6
# 
# res %>%
#   dplyr::filter(pathway == "LIN_APC_TARGETS") -> r_lin
# 
# # 2004 OS signature
# pathways <- gmtPathways("~/data/CGIACOME/Protocols_and_Documentation/Sequencing/Gene_signatures/20230821_Sansom_Apc.Hs.symbols.gmt")
# OS_7 <- pathways[["SANSOM_2004_APC_KO_UP"]]
# 
# png(file = paste0(parent_dir, "plots/fgsea", paste0(shRNA, "_Sansom_2004_APC.png")), width = 600, height = 400)
# print(plotEnrichment(OS_7, TE_named_vector)+labs(title="SANSOM 2004 2B1 TE", subtitle = "Genes upregulated upon APC KO"))
# dev.off()
# # Gives almost same distribution as OS_5 signature
# 
# # get results out (NES/adjp/LeadingEdge)
# res <- fgsea(pathways = pathways, stats = TE_named_vector, minSize = 3, maxSize = 1000)
# res %>%
#   dplyr::filter(pathway == "SANSOM_2004_APC_KO_UP") -> r7
# 
# # get all results out as a data table
# ISR_res %>%
#   bind_rows(AU_res, AD_res, r1, r2, r3, r4, r5, r6, r_lin, r7) %>%
#   unnest(leadingEdge)-> full_res
# 
# # write it out
# write_csv(full_res, file = file.path(parent_dir, "Analysis/fgsea", paste0("TE_", shRNA, "_fgsea_curated_ofInterest.csv")))
# 

######### keep interactives out of the way for now -----
#make_interactive_scatter <- function(gsea_set, pathway, df, dir) {
#  gene_names <- gsea_set[[pathway]]

#  df %>%
#    filter(Human_gene_name %in% gene_names) %>%
#    mutate(groupings = factor(case_when(RPFs_log2FC < 0 & RPFs_padj < 0.1 ~ -1,
#                                        RPFs_log2FC > 0 & RPFs_padj < 0.1 ~ 1,
#                                        RPFs_padj >= 0.1 ~ 0))) -> filtered_data 


#make gene ID/sym annotation table
#  filtered_data %>%
#    select(gene, gene_sym, transcript) %>%
#    dplyr::rename(GeneID = gene) %>%
#    as.data.frame() -> gene_anno
#  rownames(gene_anno) <- gene_anno$GeneID

#merge norm counts with gene annotation so that row names are in the same order and then select all counts columns in preferred order and make geneIDs row names
#  gene_anno %>%
#    inner_join(filtered_RPFs_norm_counts, by = "GeneID") %>%
#    inner_join(filtered_totals_norm_counts, by = "GeneID") %>%
#    column_to_rownames("GeneID") %>%
#    select(-c("transcript", "gene_sym")) -> merged_norm_counts
#  
#check row names match up
#  if (all(rownames(gene_anno) == filtered_data$gene) & all(rownames(gene_anno) == row.names(merged_norm_counts))) {
#make html files
#    glXYPlot(x = filtered_data$totals_log2FC,
#             y = filtered_data$RPFs_log2FC,
#             xlab = "RNA log2FC",
#             ylab = "RPFs log2FC",
#             main = paste(treatment, pathway),
#             status = filtered_data$groupings,
#             #cols = col_pal[1:3],
#             counts = merged_norm_counts,
#             side.xlab = "Sample",
#             side.ylab = "norm counts",
#             sample.cols = rep(c("#F8766D", "#00BA38", "#619CFF"),4), #these are the colours for each sample within the norm counts, needs to be the same length as groups
#             groups = factor(c(rep("Ctrl RPF", 3), rep(paste(treatment, "RPF"), 3),
#                               rep("Ctrl RNA", 3), rep(paste(treatment, "RNA"), 3)), 
#                             levels = c("Ctrl RPF",paste(treatment, "RPF"), "Ctrl RNA", paste(treatment, "RNA")), ordered = T),
#             anno = gene_anno,
#             path = file.path(parent_dir, "plots/fgsea/Interactive_scatters"),
#             folder = dir,
#             html = paste(treatment, pathway, sep = "_"))
#  }
#}

#make a mouse to human gene conversion table
#this is only applicable for mouse data as the gsea gene names are all human
#Mouse2HumanTable <- read_tsv(file = "\\\\data.beatson.gla.ac.uk/data/R11/bioinformatics_resources/useful_tables/mouse_to_human_gene_IDs.tsv")
#Mouse2HumanTable %>%
#  dplyr::select(Gene_stable_ID_version, Human_gene_name) %>%
#  filter(!(is.na(Human_gene_name))) %>%
#  group_by(Gene_stable_ID_version) %>%
#  sample_n(size = 1) %>%
#  dplyr::rename(gene = Gene_stable_ID_version) -> Mouse2HumanTable
