# load all the libraries needed
library(tidyverse)
#library(conflicted)
library(amap)
library(clusterProfiler)
library(reshape2)
library(pathviewr)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(ggnewscale)
library(ggridges)
library(scales)
library(tximport)
library(vsn)
library(viridis)
library(msigdbr)

# source common variables
source("common_variables.R")

# read in the VENN plots data -----
TE_4_venn <- read_csv(file = file.path(tables_dir, "20231101_RPFs_groups_for_ORA.csv"))

# bkg set are all the genes in the TE_4 venn list, which was already filtered for low reads during DEseq2
# needs also to be changed from symbol to entrez ID
TE_4_venn %>%
  pull(gene_sym) -> BKG

BKG_Entrez = bitr(BKG, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)

# make vectors - check  numbers as venn ----
# 2B1 only RPFs up 884 to 875
# 2B4 only RPFs up 56 to 54
# 2B1 only RPFs down 1353 to 1342
# 2B4 only RPFs down 78 to 78
# 2B1 and 2B4 RPFs down in both 130 to 130
# 2B1 and 2B4 RPFs up in both 103 to 103

TE_4_venn %>%
  dplyr::filter(cumulative_group == "RPFs up 2B1 only") %>%
  pull(gene_sym) %>%
  bitr(fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) %>%
  pull(ENTREZID)-> RPF_up_2B1
   
TE_4_venn %>%
  dplyr::filter(cumulative_group == "RPFs up 2B4 only") %>%
  pull(gene_sym) %>%
  bitr(fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) %>%
  pull(ENTREZID) -> RPF_up_2B4

TE_4_venn %>%
  dplyr::filter(cumulative_group == "RPFs down 2B1 only") %>%
  pull(gene_sym) %>%
  bitr(fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) %>%
  pull(ENTREZID) -> RPF_down_2B1

TE_4_venn %>%
  dplyr::filter(cumulative_group == "RPFs down 2B4 only") %>%
  pull(gene_sym) %>%
  bitr(fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) %>%
  pull(ENTREZID) -> RPF_down_2B4

TE_4_venn %>%
  dplyr::filter(cumulative_group == "RPFs down 2B1 & 2B4") %>%
  pull(gene_sym) %>%
  bitr(fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) %>%
  pull(ENTREZID) -> RPF_down_both

TE_4_venn %>%
  dplyr::filter(cumulative_group == "RPFs up 2B1 & 2B4") %>%
  pull(gene_sym) %>%
  bitr(fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) %>%
  pull(ENTREZID) -> RPF_up_both

# To run Ora with enrichGO - works only on Molecular Functions, Cell Components, and Biological processes
## Molecular Functions ------
enrichGO(gene = RPF_down_2B1, OrgDb = org.Hs.eg.db,readable = T, ont = "MF", 
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID) -> mf_ORA_2B1_dn

pdf(file = paste0(plot_dir,"RPF_down_2B1_MF_dot.pdf"), height = 5, width = 5)
dotplot(mf_ORA_2B1_dn, showCategory=20, title = "RPF_down_2B1_MF") + scale_y_discrete(labels = function(mf_ORA_2B1_dn) str_wrap(mf_ORA_2B1_dn, width = 80))
dev.off()

#save results
#save(file = file.path(parent_dir, "Analysis/fgsea/ORA_RPF_down_2B1_mol_funs_results_venn.Rdata"), mf_ORA_2B1_dn)

enrichGO(gene = RPF_down_2B4, OrgDb = org.Hs.eg.db,readable = T, ont = "MF",
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID) -> go_enrich

pdf(file = paste0(plot_dir,"RPF_down_2B4_MF_dot.pdf"), height = 5, width = 5)
dotplot(go_enrich, showCategory=20, title = "RPF_down_2B4_MF")+ scale_y_discrete (labels = function(go_enrich) str_wrap(go_enrich, width = 80))
dev.off()

enrichGO(gene = RPF_up_2B1, OrgDb = org.Hs.eg.db,readable = T, ont = "MF", 
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID) -> mf_ORA_2B1_up

pdf(file = paste0(plot_dir,"RPF_up_2B1_MF_dot.pdf"), height = 5, width = 5)
dotplot(mf_ORA_2B1_up, showCategory=20, title = "RPF_up_2B1_MF") + scale_y_discrete(labels = function(mf_ORA_2B1_up) str_wrap(mf_ORA_2B1_up, width = 80))
dev.off()

#save(file = file.path(parent_dir, "Analysis/fgsea/ORA_RPF_up_2B1_mol_funs_results_venn.Rdata"), mf_ORA_2B1_up)

enrichGO(gene = RPF_up_2B4, OrgDb = org.Hs.eg.db,readable = T, ont = "MF", 
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID) -> mf_ORA_2B4_up

pdf(file = paste0(plot_dir,"RPF_up_2B4_MF_dot.pdf"), height = 5, width = 5)
dotplot(mf_ORA_2B4_up, showCategory=20, title = "RPF_up_2B4_MF") + scale_y_discrete(labels = function(mf_ORA_2B4_up) str_wrap(mf_ORA_2B4_up, width = 80))
dev.off()

#save(file = file.path(parent_dir, "Analysis/fgsea/ORA_RPF_up_2B4_mol_funs_results_venn.Rdata"), mf_ORA_2B4_up)

## Cell components ------
enrichGO(gene = RPF_down_2B1, OrgDb = org.Hs.eg.db,readable = T, ont = "CC", 
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID) -> cc_ORA_2B1_dn

pdf(file = paste0(plot_dir,"RPF_down_2B1_CC_dot.pdf"), height = 5, width = 5)
dotplot(cc_ORA_2B1_dn, showCategory=20, title = "RPF_down_2B1_CC") + scale_y_discrete(labels = function(cc_ORA_2B1_dn) str_wrap(cc_ORA_2B1_dn, width = 80))
dev.off()

#save(file = file.path(parent_dir, "Analysis/fgsea/ORA_RPF_down_2B1_cell_comp_results_venn.Rdata"), cc_ORA_2B1_dn)

enrichGO(gene = RPF_down_2B4, OrgDb = org.Hs.eg.db,readable = T, ont = "CC",
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID) -> cc_ORA_2B4_dn

pdf(file = paste0(plot_dir,"RPF_down_2B4_CC_dot.pdf"), height = 5, width = 5)
dotplot(cc_ORA_2B4_dn, showCategory=20, title = "RPF_down_2B4_CC") + scale_y_discrete(labels = function(cc_ORA_2B4_dn) str_wrap(cc_ORA_2B4_dn, width = 80))
dev.off()

#save(file = file.path(parent_dir, "Analysis/fgsea/ORA_RPF_down_2B4_cell_comp_results_venn.Rdata"), cc_ORA_2B4_dn)

enrichGO(gene = RPF_up_2B1, OrgDb = org.Hs.eg.db,readable = T, ont = "CC", 
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID) -> cc_ORA_2B1_up

pdf(file = paste0(plot_dir,"RPF_up_2B1_CC_dot.pdf"), height = 5, width = 5)
dotplot(cc_ORA_2B1_up, showCategory=20, title = "RPF_up_2B1_CC") + scale_y_discrete(labels = function(cc_ORA_2B1_up) str_wrap(cc_ORA_2B1_up, width = 80))
dev.off()

#save(file = file.path(parent_dir, "Analysis/fgsea/ORA_RPF_up_2B1_cell_comp_results_venn.Rdata"), cc_ORA_2B1_up)

enrichGO(gene = RPF_up_2B4, OrgDb = org.Hs.eg.db,readable = T, ont = "CC", 
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID) -> cc_ORA_2B4_up

pdf(file = paste0(plot_dir,"RPF_up_2B4_CC_dot.pdf"), height = 5, width = 5)
dotplot(cc_ORA_2B4_up, showCategory=20, title = "RPF_up_2B4_CC") + scale_y_discrete(labels = function(cc_ORA_2B4_up) str_wrap(cc_ORA_2B4_up, width = 80))
dev.off()

#save(file = file.path(parent_dir, "Analysis/fgsea/ORA_RPF_up_2B4_cell_comp_results_venn.Rdata"), cc_ORA_2B4_up)

## Biological Processes ------
enrichGO(gene = RPF_down_2B1, OrgDb = org.Hs.eg.db,readable = T, ont = "BP", 
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID) -> bp_ORA_2B1_dn

pdf(file = paste0(plot_dir,"RPF_down_2B1_BP_dot.pdf"), height = 5, width = 5)
dotplot(bp_ORA_2B1_dn, showCategory=20, title = "RPF_down_2B1_BP") + scale_y_discrete(labels = function(bp_ORA_2B1_dn) str_wrap(bp_ORA_2B1_dn, width = 80))
dev.off()

#save(file = file.path(parent_dir, "Analysis/fgsea/ORA_RPF_down_2B1_biol_process_venn.Rdata"), bp_ORA_2B1_dn)

enrichGO(gene = RPF_down_2B4, OrgDb = org.Hs.eg.db,readable = T, ont = "BP",
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID) -> bp_ORA_2B4_dn

pdf(file = paste0(plot_dir,"RPF_down_2B4_BP_dot.pdf"), height = 5, width = 5)
dotplot(bp_ORA_2B4_dn, showCategory=20, title = "RPF_down_2B4_BP") + scale_y_discrete(labels = function(bp_ORA_2B4_dn) str_wrap(bp_ORA_2B4_dn, width = 80))
dev.off()

#save(file = file.path(parent_dir, "Analysis/fgsea/ORA_RPF_down_2B4_biol_process_venn.Rdata"), bp_ORA_2B4_dn)

enrichGO(gene = RPF_up_2B1, OrgDb = org.Hs.eg.db,readable = T, ont = "BP", 
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID) -> bp_ORA_2B1_up

pdf(file = paste0(plot_dir,"RPF_up_2B1_BP_dot.pdf"), height = 5, width = 5)
dotplot(bp_ORA_2B1_up, showCategory=20, title = "RPF_up_2B1_BP") + scale_y_discrete(labels = function(bp_ORA_2B1_up) str_wrap(bp_ORA_2B1_up, width = 80))
dev.off()

#save(file = file.path(parent_dir, "Analysis/fgsea/ORA_RPF_up_2B1_biol_process_venn.Rdata"), bp_ORA_2B1_up)

enrichGO(gene = RPF_up_2B4, OrgDb = org.Hs.eg.db,readable = T, ont = "BP",
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID) -> go_enrich

pdf(file = paste0(plot_dir,"RPF_up_2B4_BP_dot.pdf"), height = 5, width = 5)
dotplot(go_enrich, showCategory=20, title = "RPF_up_2B4_BP") + scale_y_discrete(labels = function(go_enrich) str_wrap(go_enrich, width = 80))
dev.off()
# 
# save(file = file.path(parent_dir, "Analysis/fgsea/ORA_RPF_up_2B4_biol_process_venn.Rdata"), go_enrich)

## hallmarks from MSigDB -----
msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::rename(ont = gs_name, gene = entrez_gene)

enricher(gene = RPF_down_2B1, 
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID, TERM2GENE = msig_h) -> go_enrich

pdf(file = paste0(plot_dir,"RPF_down_2B1_HM_dot.pdf"), height = 5, width = 8)
dotplot(go_enrich, showCategory=20, title = " RPF down 2B1 hallmarks") + scale_y_discrete(labels = function(go_enrich) str_wrap(go_enrich, width = 80))
dev.off()

enricher(gene = RPF_down_2B4,
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID, TERM2GENE = msig_h) -> go_enrich

pdf(file = paste0(plot_dir,"RPF_down_2B4_HM_dot.pdf"), height = 5, width = 8)
dotplot(go_enrich, showCategory=20, title = " RPF down 2B4 hallmarks")
dev.off()

enricher(gene = RPF_up_2B1,
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID, TERM2GENE = msig_h) -> go_enrich

pdf(file = paste0(plot_dir,"RPF_up_2B1_HM_dot.pdf"), height = 5, width = 8)
dotplot(go_enrich, showCategory=20, title = " RPF up 2B1 hallmarks")
dev.off()

enricher(gene = RPF_up_2B4,
         pvalueCutoff = 0.1, qvalueCutoff = 0.10, universe = BKG_Entrez$ENTREZID, TERM2GENE = msig_h) -> go_enrich

pdf(file = paste0(plot_dir,"RPF_up_2B4_HM_dot.pdf"), height = 5, width = 8)
dotplot(go_enrich, showCategory=20, title = " RPF up 2B4 hallmarks") + scale_y_discrete(labels = function(go_enrich) str_wrap(go_enrich, width = 80))
dev.off()
