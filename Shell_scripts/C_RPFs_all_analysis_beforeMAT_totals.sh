#!/usr/bin/env bash
#This script runs the RPF analysis until before the alignment - for which you need to have done Most Abundant Transcript Fasta file from your totals

#read in variables
source CV_Foureps.sh

#run fastQC
for filename in $RPF_filenames
do
fastqc $fastq_dir/${filename}.fastq --outdir=$fastqc_dir &
done
wait

#run cutadapt
for filename in $RPF_filenames
do
cutadapt $fastq_dir/${filename}.fastq -a $RPF_adaptor --nextseq-trim=20 -m 30 -M 50 -o $fastq_dir/${filename}_cutadapt.fastq 1> $log_dir/${filename}_cutadapt_log.txt &
done
wait

#run fastqc on cutadapt output
for filename in $RPF_filenames
do
fastqc $fastq_dir/${filename}_cutadapt.fastq --outdir=$fastqc_dir &
done
wait

#read deduplication one sample at a time
for filename in $RPF_filenames
do
cd-hit-dup -i $fastq_dir/${filename}_cutadapt.fastq -o $fastq_dir/${filename}_cdhitdup.fastq -e 0
done

#run fastqc on cd-hit-dup output
for filename in $RPF_filenames
do
fastqc $fastq_dir/${filename}_cdhitdup.fastq --outdir=$fastqc_dir &
done
wait

#run cutadapt with -u 4 -u -4 to remove 4nt from either end of all reads
for filename in $RPF_filenames
do
cutadapt $fastq_dir/${filename}_cdhitdup.fastq -u 4 -u -4 -o $fastq_dir/${filename}_UMIremoved.fastq &
done
wait

#run fastqc on UMIremoved output
for filename in $RPF_filenames
do
fastqc $fastq_dir/${filename}_UMIremoved.fastq --outdir=$fastqc_dir &
done
wait
