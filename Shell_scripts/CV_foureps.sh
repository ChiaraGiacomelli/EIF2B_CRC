#!/usr/bin/env bash

###filenames
#These are the filenames for all RPF and Total RNA-seq samples (without the <.fastq> or any alternative extension)

# Analyses E1, E2, E4, E3_c (appended seq runs for E3)
RPF_filenames='NTC_min_RPF_1 NTC_plus_RPF_1 2B1_3A_min_RPF_1 2B1_3A_plus_RPF_1 2B4_2_min_RPF_1 2B4_2_plus_RPF_1 NTC_min_RPF_2 NTC_plus_RPF_2 2B1_3A_min_RPF_2 2B1_3A_plus_RPF_2 2B4_2_min_RPF_2 2B4_2_plus_RPF_2 NTC_min_RPF_4 NTC_plus_RPF_4 2B1_3A_min_RPF_4 2B1_3A_plus_RPF_4 2B4_2_min_RPF_4 2B4_2_plus_RPF_4 NTC_min_RPF_3_c NTC_plus_RPF_3_c 2B1_3A_min_RPF_3_c 2B1_3A_plus_RPF_3_c 2B4_2_min_RPF_3_c 2B4_2_plus_RPF_3_c'
Totals_filenames='NTC_min_Tot_1 NTC_plus_Tot_1 2B1_3A_min_Tot_1 2B1_3A_plus_Tot_1 2B4_2_min_Tot_1 2B4_2_plus_Tot_1 NTC_min_Tot_2 NTC_plus_Tot_2 2B1_3A_min_Tot_2 2B1_3A_plus_Tot_2 2B4_2_min_Tot_2 2B4_2_plus_Tot_2 NTC_min_Tot_4 NTC_plus_Tot_4 2B1_3A_min_Tot_4 2B1_3A_plus_Tot_4 2B4_2_min_Tot_4 2B4_2_plus_Tot_4 NTC_min_Tot_3_c NTC_plus_Tot_3_c 2B1_3A_min_Tot_3_c 2B1_3A_plus_Tot_3_c 2B4_2_min_Tot_3_c 2B4_2_plus_Tot_3_c'

###adaptor
#This is the sequence of the 3' adaptor that was used in the library prep. Common sequences are below, unhash the correct one if present, or if not enter it as a variable

#RPF adaptors
RPF_adaptor='TGGAATTCTCGGGTGCCAAGG' #this is the adaptor used in the nextflex small RNA library kit

#Totals adaptors
Totals_adaptor='AGATCGGAAGAG' #this is the adaptor used in the LEXOGEN CORALL Total RNA-Seq Library Prep Kit

###paths
parent_dir='/home/local/BICR/cgiacome/data/CGIACOME/AAA_eIF2B/AAA_Collaborative_stuff/Wurzburg_collab/AAA_Riboseq/Analysis/Ribo-seq-allreps' #This is the path to the parent directory that contains all the data and where all the processed data will be saved

#The following directories are where all the processed data will be saved. These all need to be created prior to starting the analysis

fastq_dir=${parent_dir}/fastq_files
fastqc_dir=${parent_dir}/fastQC_files
SAM_dir=${parent_dir}/SAM_files
BAM_dir=${parent_dir}/BAM_files
log_dir=${parent_dir}/logs
counts_dir=${parent_dir}/Counts_files
csv_counts_dir=${parent_dir}/Counts_files/csv_files
STAR_dir=${parent_dir}/STAR_files

rsem_dir=${parent_dir}/rsem

#The following directories are where all the csv files that are used as input into R will be saved
analysis_dir=${parent_dir}/Analysis

region_counts_dir=${analysis_dir}/region_counts
spliced_counts_dir=${analysis_dir}/spliced_counts
periodicity_dir=${analysis_dir}/periodicity
cds_counts_dir=${analysis_dir}/CDS_counts
UTR5_counts_dir=$analysis_dir/UTR5_counts
codon_counts_dir=${analysis_dir}/codon_counts
most_abundant_transcripts_dir=${analysis_dir}/most_abundant_transcripts
DESeq2_dir=${analysis_dir}/DESeq2_output
reads_summary_dir=${analysis_dir}/reads_summary

#The following directories are where all the plots generated in R will be saved
plots_dir=${parent_dir}/plots

summed_counts_plots_dir=${plots_dir}/summed_counts
periodicity_plots_dir=${plots_dir}/periodicity
offset_plots_dir=${plots_dir}/offset
heatmaps_plots_dir=${plots_dir}/heatmaps
DE_analysis_dir=${plots_dir}/DE_analysis
PCA_dir=${plots_dir}/PCAs
Interactive_scatters_dir=${plots_dir}/Interactive_scatters
fgsea_dir=${plots_dir}/fgsea
fgsea_scatters_dir=${plots_dir}/fgsea/scatters
fgsea_interactive_scatters_dir=${plots_dir}/fgsea/Interactive_scatters
read_counts_summary_dir=${plots_dir}/read_counts_summary
binned_plots_dir=${plots_dir}/binned_plots

#Fastas
fasta_dir='/home/local/BICR/cgiacome/data/R11/bioinformatics_resources/FASTAs/human'

rRNA_fasta=${fasta_dir}/rRNA/human_rRNA.fa
tRNA_fasta=${fasta_dir}/tRNA/human_mature_tRNA.fa
mito_fasta=${fasta_dir}/GENCODE/v38/filtered/gencode.v38.mito_transcripts.fa
pc_fasta=${fasta_dir}/GENCODE/v38/filtered/gencode.v38.pc_transcripts_filtered.fa
rpf_fasta=${orf_dir}/PhaseI_ORFs_20210405_for_GENCODE.cds.fa
rsem_index=${fasta_dir}/GENCODE/v38/filtered/rsem_bowtie2_index/gencode.v38.pc_transcripts_filtered
STAR_index=${fasta_dir}/GENCODE/v38/original/STAR_index
STAR_GTF=${fasta_dir}/GENCODE/v38/original/gencode.v38.annotation.gtf
most_abundant_fasta=$most_abundant_transcripts_dir/most_abundant_transcripts.fa #this needs to be created for each specific project

###fasta info
#The below needs to point to a <.csv> file that contains the following information for all transcripts within the protein coding FASTA
#transcript_ID,5'UTR length,CDS length,3'UTR length
#Running the Filter_GENCODE_FASTA.py script will generate this file as one of its outputs

region_lengths=${fasta_dir}/GENCODE/v38/transcript_info/gencode.v38.pc_transcripts_region_lengths.csv
