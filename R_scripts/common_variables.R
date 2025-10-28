# Set the parent directory (this should be the same directory as is set in the CV_foureps.sh script

# from VM
parent_dir <- 'your/home/folder' # where all data is stored, same parent_dir as used in the shell scripts
setwd(paste0(parent_dir, "/R_scripts"))
plot_dir <- paste0(parent_dir,"/plots/PDF_for_paper/")
tables_dir <- paste0(parent_dir,"/Analysis/tables_for_paper")

PCA_cols <-c("NTC" = "#fde725", "2B5" = '#20A387FF', "2B1" = '#440154FF', "2B4" = '#33638DFF')

#set sample names - only 3_c and not 3
RPF_sample_names <- c('NTC_min_RPF_1', 'NTC_min_RPF_2', 'NTC_min_RPF_3_c', 'NTC_min_RPF_4',
                      '2B1_3A_min_RPF_1', '2B1_3A_min_RPF_2','2B1_3A_min_RPF_3_c', '2B1_3A_min_RPF_4',
                      '2B4_2_min_RPF_1', '2B4_2_min_RPF_2','2B4_2_min_RPF_3_c', '2B4_2_min_RPF_4',
                      'NTC_plus_RPF_1', 'NTC_plus_RPF_2', 'NTC_plus_RPF_3_c', 'NTC_plus_RPF_4',
                      '2B1_3A_plus_RPF_1', '2B1_3A_plus_RPF_2','2B1_3A_plus_RPF_3_c', '2B1_3A_plus_RPF_4',
                      '2B4_2_plus_RPF_1', '2B4_2_plus_RPF_2','2B4_2_plus_RPF_3_c', '2B4_2_plus_RPF_4')

Total_sample_names <- c('NTC_min_Tot_1', 'NTC_min_Tot_2', 'NTC_min_Tot_3_c', 'NTC_min_Tot_4',
                        '2B1_3A_min_Tot_1', '2B1_3A_min_Tot_2', '2B1_3A_min_Tot_3_c', '2B1_3A_min_Tot_4',
                        '2B4_2_min_Tot_1', '2B4_2_min_Tot_2','2B4_2_min_Tot_3_c', '2B4_2_min_Tot_4',
                        'NTC_plus_Tot_1', 'NTC_plus_Tot_2', 'NTC_plus_Tot_3_c', 'NTC_plus_Tot_4',
                        '2B1_3A_plus_Tot_1', '2B1_3A_plus_Tot_2','2B1_3A_plus_Tot_3_c', '2B1_3A_plus_Tot_4',
                        '2B4_2_plus_Tot_1', '2B4_2_plus_Tot_2','2B4_2_plus_Tot_3_c', '2B4_2_plus_Tot_4')

RPF_sample_info <- data.frame(sample = RPF_sample_names,
                              condition = factor(c(rep("min", 12), rep("plus", 12))),
                              replicate = factor(rep(c("1", "2", "3", "4"), 6)),
                              shRNA = c(rep("NTC",4),rep("2B1",4),rep("2B4",4), rep("NTC",4),rep("2B1",4),rep("2B4",4)))

#set read lengths used for library QC
lengths <- 25:35

EIF2B_subunits = c(EIF2B1, EIF2B2, EIF2B3, EIF2B4, EIF2B5)
GOI = c(ATF4, DDIT3_CHOP, GADD34, DENR, MCTS1, EIF2D, MYC)

colors <- viridis(12)

RPF_NTC <- c('NTC_min_RPF_1', 'NTC_min_RPF_2', 'NTC_min_RPF_3_c', 'NTC_min_RPF_4',
             'NTC_plus_RPF_1', 'NTC_plus_RPF_2', 'NTC_plus_RPF_3_c', 'NTC_plus_RPF_4')

RPF_2B1 <- c('2B1_3A_min_RPF_1', '2B1_3A_min_RPF_2','2B1_3A_min_RPF_3_c', '2B1_3A_min_RPF_4',
             '2B1_3A_plus_RPF_1', '2B1_3A_plus_RPF_2','2B1_3A_plus_RPF_3_c', '2B1_3A_plus_RPF_4')

RPF_2B4 <- c('2B4_2_min_RPF_1', '2B4_2_min_RPF_2','2B4_2_min_RPF_3_c', '2B4_2_min_RPF_4',
             '2B4_2_plus_RPF_1', '2B4_2_plus_RPF_2','2B4_2_plus_RPF_3_c', '2B4_2_plus_RPF_4')

