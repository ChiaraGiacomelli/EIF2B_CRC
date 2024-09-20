# VM path
base_dir <- "~/data/CGIACOME/Protocols_and_Documentation/Sequencing/Gene_signatures/msigdb_v2023.1.Hs_GMTs/"
# VPN path
#base_dir <- "Volume/data-1/CGIACOME/Protocols_and_Documentation/Sequencing/Gene_signatures/msigdb_v2023.1.Hs_GMTs/"

pathways.hallmark <- gmtPathways(paste0(base_dir, "h.all.v2023.1.Hs.symbols.gmt"))
pathways.curated <- gmtPathways(paste0(base_dir, "c2.all.v2023.1.Hs.symbols.gmt"))
pathways.miRNA_targets <- gmtPathways(paste0(base_dir, "c3.mir.v2023.1.Hs.symbols.gmt"))
pathways.transcription_factors <- gmtPathways(paste0(base_dir, "c3.tft.v2023.1.Hs.symbols.gmt"))
pathways.mol_funs <- gmtPathways(paste0(base_dir, "c5.go.mf.v2023.1.Hs.symbols.gmt"))
pathways.bio_processes <- gmtPathways(paste0(base_dir, "c5.go.bp.v2023.1.Hs.symbols.gmt"))
pathways.cell_comp <- gmtPathways(paste0(base_dir, "c5.go.cc.v2023.1.Hs.symbols.gmt"))
pathways.tumour_phen_onto <- gmtPathways(paste0(base_dir, "c6.all.v2023.1.Hs.symbols.gmt"))
pathways.cell_type_sig <- gmtPathways(paste0(base_dir, "c8.all.v2023.1.Hs.symbols.gmt"))
pathways.biocarta <- gmtPathways(paste0(base_dir, "c2.cp.biocarta.v2023.1.Hs.symbols.gmt"))
pathways.reactome <- gmtPathways(paste0(base_dir, "c2.cp.reactome.v2023.1.Hs.symbols.gmt"))
pathways.wp <- gmtPathways(paste0(base_dir, "c2.cp.wikipathways.v2023.1.Hs.symbols.gmt"))
pathways.onco <- gmtPathways(paste0(base_dir, "c6.all.v2023.1.Hs.symbols.gmt"))
pathways.kegg <- gmtPathways(paste0(base_dir, "c2.cp.kegg.v2023.1.Hs.symbols.gmt")) #now is kegg_medicus in use
