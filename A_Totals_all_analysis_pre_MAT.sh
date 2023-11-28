#!/usr/bin/env bash

#This script runs fastQC on all your raw fastq files and outputs them in the fastQC directory

#read in variables
source CV_Foureps.sh

#run fastQC
for filename in $Totals_filenames
do
fastqc $fastq_dir/${filename}.fastq --outdir=$fastqc_dir &
done
wait

#run cutadapt
for filename in $Totals_filenames
do
cutadapt $fastq_dir/${filename}.fastq -a $Totals_adaptor --nextseq-trim=20 -m 30 -o $fastq_dir/${filename}_cutadapt.fastq 1> $log_dir/${filename}_cutadapt_log.txt &
done
wait

#run fastqc on cutadapt output
for filename in $Totals_filenames
do
fastqc $fastq_dir/${filename}_cutadapt.fastq --outdir=$fastqc_dir &
done
wait

#read deduplication modified one sample at a time
for filename in $Totals_filenames
do
cd-hit-dup -i $fastq_dir/${filename}_cutadapt.fastq -o $fastq_dir/${filename}_cdhitdup.fastq -e 0
done

#run fastqc on cd-hit-dup output
for filename in $Totals_filenames
do
fastqc $fastq_dir/${filename}_cdhitdup.fastq --outdir=$fastqc_dir &
done
wait

#run cutadapt with -u 12 to remove 12nt from the 5' end of all reads
for filename in $Totals_filenames
do
cutadapt $fastq_dir/${filename}_cdhitdup.fastq -u 12 -o $fastq_dir/${filename}_UMIremoved.fastq &
done
wait

#run fastqc on UMIremoved output
for filename in $Totals_filenames
do
fastqc $fastq_dir/${filename}_UMIremoved.fastq --outdir=$fastqc_dir &
done
wait

#Align to rRNA
for filename in $Totals_filenames
do
bbmap.sh in=$fastq_dir/${filename}_UMIremoved.fastq out=$SAM_dir/${filename}_rRNA.sam ref=$rRNA_fasta outm=$fastq_dir/${filename}_rRNA.fastq outu=$fastq_dir/${filename}_non_rRNA.fastq ambiguous=best nodisk threads=30 2> $log_dir/${filename}_rRNA_log.txt
done

#Align to protein coding transcriptomes
for filename in $Totals_filenames
do
rsem-calculate-expression --strandedness forward --bowtie2 --fragment-length-mean 300 --fragment-length-sd 100 $fastq_dir/${filename}_non_rRNA.fastq $rsem_index $rsem_dir/${filename} &
done
wait
