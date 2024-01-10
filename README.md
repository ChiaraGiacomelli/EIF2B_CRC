# EIF2B CRC - exploring transcriptional and translational changes upon modulation of eIF2B activity in colorectal cancer

This repository contains the scripts used to analyse Riboseq and cytoplasmic RNA sequencing results associate with the manuscript "Sensing of p-eIF2alpha levels by the eIF2B complex is vital for colorectal cancer" (html here)
Fastq files are deposited with the code (XXX) in the Gene Expression Omnibus

The project utilizes codes developed primarily by Dr. Joseph Waldron at the Bushell lab (Cancer Research UK Scotland Institute) and available on github at the main project branch https://github.com/Bushell-lab/Ribo-seq.

For running shell scripts, two conda environments are used: one for total cytoplasmic RNA analysis (scripts A_ and B_), and one for the RPF analysis (scripts C_, D_, E_). Conda environments are generated as per https://github.com/Bushell-lab/Ribo-seq description.
Please note that the current pipeline on the Bushell-lab is slightly different than the one used for this project.
The main difference relies on the stage at which deduplication is performed - this pipeline performs deduplication on all data prior to alingment and removal of read for rRNA/tRNA/mitochondrial RNAs contaminants. The current pipeline at https://github.com/Bushell-lab/Ribo-seq runs deduplication after RNAs contaminats reads are removed.

Shell scripts analyses were run on a terminal in Ubuntu 20.04.6 LTS
R scripts were run in RStudio 2023.12.0, Build 369 using R version 4.3.2 (2023-10-31)  
