#Set the parent directory (this should be the same directory as is set in the common_variables.sh script

# from VM
parent_dir <- '~/data/CGIACOME/AAA_eIF2B/AAA_Collaborative_stuff/Wurzburg_collab/AAA_Riboseq/Analysis/Ribo-seq-allreps'
# Windows desktop in lab
#parent_dir <- 'N:/CGIACOME/AAA_eIF2B/AAA_Collaborative_stuff/Wurzburg_collab/AAA_Riboseq/Analysis/Ribo-seq-allreps'
# Mac parent directory
#parent_dir <- '/Volumes/data-1/CGIACOME/AAA_eIF2B/AAA_Collaborative_stuff/Wurzburg_collab/AAA_Riboseq/Analysis/Ribo-seq-allreps'
#setwd(paste0(parent_dir, "/R_scripts"))
plot_dir <- paste0(parent_dir,"/plots/PDF_for_paper/")
tables_dir <- paste0(parent_dir,"/Analysis/tables_for_paper")

#PCA_cols <-c("NTC" = "#fde725", "2B5" = '#20A387FF', "2B1" = '#440154FF', "2B4" = '#33638DFF')

#set sample names - only 3_c and not 3
RPF_sample_names <- c('NTC_min_RPF_1', 'NTC_min_RPF_2', 'NTC_min_RPF_3_c', 'NTC_min_RPF_4',
                      '2B1_3A_min_RPF_1', '2B1_3A_min_RPF_2','2B1_3A_min_RPF_3_c', '2B1_3A_min_RPF_4',
                      '2B4_2_min_RPF_1', '2B4_2_min_RPF_2','2B4_2_min_RPF_3_c', '2B4_2_min_RPF_4',
                      '2B5_min_RPF_1', '2B5_min_RPF_2', '2B5_min_RPF_3_c', '2B5_min_RPF_4',
                      'NTC_plus_RPF_1', 'NTC_plus_RPF_2', 'NTC_plus_RPF_3_c', 'NTC_plus_RPF_4',
                      '2B1_3A_plus_RPF_1', '2B1_3A_plus_RPF_2','2B1_3A_plus_RPF_3_c', '2B1_3A_plus_RPF_4',
                      '2B4_2_plus_RPF_1', '2B4_2_plus_RPF_2','2B4_2_plus_RPF_3_c', '2B4_2_plus_RPF_4',
                      '2B5_plus_RPF_1', '2B5_plus_RPF_2', '2B5_plus_RPF_3_c', '2B5_plus_RPF_4')

Total_sample_names <- c('NTC_min_Tot_1', 'NTC_min_Tot_2', 'NTC_min_Tot_3_c', 'NTC_min_Tot_4',
                        '2B1_3A_min_Tot_1', '2B1_3A_min_Tot_2', '2B1_3A_min_Tot_3_c', '2B1_3A_min_Tot_4',
                        '2B4_2_min_Tot_1', '2B4_2_min_Tot_2','2B4_2_min_Tot_3_c', '2B4_2_min_Tot_4',
                        '2B5_min_Tot_1', '2B5_min_Tot_2', '2B5_min_Tot_3_c', '2B5_min_Tot_4',
                        'NTC_plus_Tot_1', 'NTC_plus_Tot_2', 'NTC_plus_Tot_3_c', 'NTC_plus_Tot_4',
                        '2B1_3A_plus_Tot_1', '2B1_3A_plus_Tot_2','2B1_3A_plus_Tot_3_c', '2B1_3A_plus_Tot_4',
                        '2B4_2_plus_Tot_1', '2B4_2_plus_Tot_2','2B4_2_plus_Tot_3_c', '2B4_2_plus_Tot_4',
                        '2B5_plus_Tot_1', '2B5_plus_Tot_2', '2B5_plus_Tot_3_c', '2B5_plus_Tot_4')

RPF_sample_info <- data.frame(sample = RPF_sample_names,
                              condition = factor(c(rep("min", 16), rep("plus", 16))),
                              replicate = factor(rep(c("1", "2", "3", "4"), 8)),
                              shRNA = c(rep("NTC",4),rep("2B1",4),rep("2B4",4),rep("2B5",4), rep("NTC",4),rep("2B1",4),rep("2B4",4),rep("2B5",4)))

#set read lengths used for library QC
lengths <- 25:35

#gene IDs of subunits & other interesting stuff
EIF2B5 <- 'ENSG00000145191.15'
EIF2B4 <- 'ENSG00000115211.16'
EIF2B3 <- 'ENSG00000070785.17'
EIF2B2 <- 'ENSG00000119718.11'
EIF2B1 <- 'ENSG00000111361.13'
ATF4 <- 'ENSG00000128272.17'
DDIT3_CHOP <- 'ENSG00000175197.13'
GADD34 <- 'ENSG00000087074.8'
DENR <- 'ENSG00000139726.11'
MCTS1 <- 'ENSG00000232119.8'
EIF2D <- 'ENSG00000143486.16'
MYC <- 'ENSG00000136997.21'

EIF2B_subunits = c(EIF2B1, EIF2B2, EIF2B3, EIF2B4, EIF2B5)
GOI = c(ATF4, DDIT3_CHOP, GADD34, DENR, MCTS1, EIF2D, MYC)

colors <- viridis(12)

RPF_NTC <- c('NTC_min_RPF_1', 'NTC_min_RPF_2', 'NTC_min_RPF_3_c', 'NTC_min_RPF_4',
             'NTC_plus_RPF_1', 'NTC_plus_RPF_2', 'NTC_plus_RPF_3_c', 'NTC_plus_RPF_4')

RPF_2B1 <- c('2B1_3A_min_RPF_1', '2B1_3A_min_RPF_2','2B1_3A_min_RPF_3_c', '2B1_3A_min_RPF_4',
             '2B1_3A_plus_RPF_1', '2B1_3A_plus_RPF_2','2B1_3A_plus_RPF_3_c', '2B1_3A_plus_RPF_4')

RPF_2B4 <- c('2B4_2_min_RPF_1', '2B4_2_min_RPF_2','2B4_2_min_RPF_3_c', '2B4_2_min_RPF_4',
             '2B4_2_plus_RPF_1', '2B4_2_plus_RPF_2','2B4_2_plus_RPF_3_c', '2B4_2_plus_RPF_4')


#set sample names all "5" replicates
#RPF_sample_names <- c('NTC_min_RPF_1', 'NTC_min_RPF_2', 'NTC_min_RPF_3','NTC_min_RPF_3_c', 'NTC_min_RPF_4',
                      #'2B1_3A_min_RPF_1', '2B1_3A_min_RPF_2', '2B1_3A_min_RPF_3','2B1_3A_min_RPF_3_c', '2B1_3A_min_RPF_4',
                      #'2B4_2_min_RPF_1', '2B4_2_min_RPF_2','2B4_2_min_RPF_3','2B4_2_min_RPF_3_c', '2B4_2_min_RPF_4',
                      #'2B5_min_RPF_1', '2B5_min_RPF_2', '2B5_min_RPF_3', '2B5_min_RPF_3_c', '2B5_min_RPF_4',
                      #'NTC_plus_RPF_1', 'NTC_plus_RPF_2', 'NTC_plus_RPF_3','NTC_plus_RPF_3_c', 'NTC_plus_RPF_4',
                      #'2B1_3A_plus_RPF_1', '2B1_3A_plus_RPF_2', '2B1_3A_plus_RPF_3','2B1_3A_plus_RPF_3_c', '2B1_3A_plus_RPF_4',
                      #'2B4_2_plus_RPF_1', '2B4_2_plus_RPF_2','2B4_2_plus_RPF_3','2B4_2_plus_RPF_3_c', '2B4_2_plus_RPF_4',
                      #'2B5_plus_RPF_1', '2B5_plus_RPF_2', '2B5_plus_RPF_3', '2B5_plus_RPF_3_c', '2B5_plus_RPF_4')

#Total_sample_names <- c('NTC_min_Tot_1', 'NTC_min_Tot_2', 'NTC_min_Tot_3','NTC_min_Tot_3_c', 'NTC_min_Tot_4',
                        #'2B1_3A_min_Tot_1', '2B1_3A_min_Tot_2', '2B1_3A_min_Tot_3','2B1_3A_min_Tot_3_c', '2B1_3A_min_Tot_4',
                        #'2B4_2_min_Tot_1', '2B4_2_min_Tot_2','2B4_2_min_Tot_3','2B4_2_min_Tot_3_c', '2B4_2_min_Tot_4',
                        #'2B5_min_Tot_1', '2B5_min_Tot_2', '2B5_min_Tot_3', '2B5_min_Tot_3_c', '2B5_min_Tot_4',
                        #'NTC_plus_Tot_1', 'NTC_plus_Tot_2', 'NTC_plus_Tot_3','NTC_plus_Tot_3_c', 'NTC_plus_Tot_4',
                        #'2B1_3A_plus_Tot_1', '2B1_3A_plus_Tot_2', '2B1_3A_plus_Tot_3','2B1_3A_plus_Tot_3_c', '2B1_3A_plus_Tot_4',
                        #'2B4_2_plus_Tot_1', '2B4_2_plus_Tot_2','2B4_2_plus_Tot_3','2B4_2_plus_Tot_3_c', '2B4_2_plus_Tot_4',
                        #'2B5_plus_Tot_1', '2B5_plus_Tot_2', '2B5_plus_Tot_3', '2B5_plus_Tot_3_c', '2B5_plus_Tot_4')
