#!/usr/bin/env bash

#This script uses bbmap to align reads.
#in specificies the input <.fastq> file.
#out specificies the output <.SAM> file.
#ref specifies the <.fasta> file to use as a reference. bbmap will use this to make an index. As this is much quicker than other alignment programs, we use the nodisk option so that this isn't written to file
#outm and outu specificies filenames to write <.fastq> files for all reads that either align or do not align respectively
#ambigous specifies how to treat multimapped reads. We use best (keeps the highest scored alignment).
#2> stores the text that is printed to the screen as a log
#We first align to rRNA, then use everything that did not align to rRNA as input to align to tRNAs, then mitochondrial mRNAs and finally protein coding mRNAs. We then use fastQC on all the output <.fastq> files

#read in variables
source CV_foureps.sh

# Align to rRNA
for filename in $RPF_filenames
do
bbmap.sh in=$fastq_dir/${filename}_UMIremoved.fastq out=$SAM_dir/${filename}_rRNA.sam ref=$rRNA_fasta outm=$fastq_dir/${filename}_rRNA.fastq outu=$fastq_dir/${filename}_non_rRNA.fastq ambiguous=best nodisk threads=24 2> $log_dir/${filename}_rRNA_log.txt
done

# Align to tRNA fasta
for filename in $RPF_filenames
do
bbmap.sh in=$fastq_dir/${filename}_non_rRNA.fastq out=$SAM_dir/${filename}_tRNA.sam ref=$tRNA_fasta outm=$fastq_dir/${filename}_tRNA.fastq outu=$fastq_dir/${filename}_non_rRNA_tRNA.fastq ambiguous=best nodisk threads=24 2> $log_dir/${filename}_tRNA_log.txt
done

# Align to mito fasta
for filename in $RPF_filenames
do
bbmap.sh in=$fastq_dir/${filename}_non_rRNA_tRNA.fastq out=$SAM_dir/${filename}_mito.sam ref=$mito_fasta outm=$fastq_dir/${filename}_mito.fastq outu=$fastq_dir/${filename}_non_rRNA_tRNA_mito.fastq ambiguous=best nodisk threads=24 2> $log_dir/${filename}_mito_log.txt
done

#Align to most abundant transcript
for filename in $RPF_filenames
do
bbmap.sh in=$fastq_dir/${filename}_non_rRNA_tRNA_mito.fastq out=$SAM_dir/${filename}_pc.sam ref=$most_abundant_fasta outm=$fastq_dir/${filename}_pc.fastq outu=$fastq_dir/${filename}_unaligned.fastq ambiguous=best nodisk threads=24 2> $log_dir/${filename}_pc_log.txt
done

#run fastqc on mapped reads
for filename in $RPF_filenames
do
fastqc $fastq_dir/${filename}_rRNA.fastq --outdir=$fastqc_dir &
fastqc $fastq_dir/${filename}_tRNA.fastq --outdir=$fastqc_dir &
fastqc $fastq_dir/${filename}_mito.fastq --outdir=$fastqc_dir &
fastqc $fastq_dir/${filename}_pc.fastq --outdir=$fastqc_dir &
fastqc $fastq_dir/${filename}_unaligned.fastq --outdir=$fastqc_dir &
done
wait

# This section uses SAMtools to convert <.sam> files to <.bam> files, sort the <.bam> files and make an <.bai> index files.
#When using the sorted <.bam> file as input to the counting_script.py script, ensure the corresponding <.bai> (index file) is in the same directory

# do SAM to BAM for the coding regions
#convert sam to sam
for filename in $RPF_filenames
do
samtools view -bS $SAM_dir/${filename}_pc.sam > $BAM_dir/${filename}_pc.bam &
done
wait

#sort bam
for filename in $RPF_filenames
do
samtools sort $BAM_dir/${filename}_pc.bam -o $BAM_dir/${filename}_pc_sorted.bam -@ 6 -m 1G
done

#index bam
for filename in $RPF_filenames
do
samtools index $BAM_dir/${filename}_pc_sorted.bam $BAM_dir/${filename}_pc_sorted.bai &
done
wait

# process protein coding regions
#make an fai (fasta index) file from the fasta using samtools. This is required for the counting script and needs to exist before running counting_script.py
samtools faidx $most_abundant_fasta

#run the counting_script.py with a range of read lengths (adjust below if required, currently set to 25-35)
for filename in $RPF_filenames
do
for length in $(seq 25 35)
do
counting_script.py -bam $BAM_dir/${filename}_pc_sorted.bam -fasta $most_abundant_fasta -len $length -out_file ${filename}_pc_L${length}_Off0.counts -out_dir $counts_dir &
done
done
wait

#set offset
offset=15

#run summing_region_counts.py script
for filename in $RPF_filenames
do
for length in $(seq 25 35)
do
summing_region_counts.py ${filename}_pc_L${length}_Off0.counts $offset $region_lengths -in_dir $counts_dir -out_dir $region_counts_dir &
done
done
wait

#SPLIT COUNT
#set number of nt to splice
n=50

#run summing_spliced_counts.py script
for filename in $RPF_filenames
do
for length in $(seq 25 35)
do
summing_spliced_counts.py ${filename}_pc_L${length}_Off0.counts $n $region_lengths -in_dir $counts_dir -out_dir $spliced_counts_dir &
done
done
wait

#PERIODICITY
#set offset 
offset=15

#run periodicity.py script
for filename in $RPF_filenames
do
for length in $(seq 25 35)
do
periodicity.py ${filename}_pc_L${length}_Off0.counts $region_lengths -offset $offset -in_dir $counts_dir -out_dir $periodicity_dir &
done
done
wait

# at this stage, QC needs to be done with R to check for periodicity, region coverage, and read amounts

