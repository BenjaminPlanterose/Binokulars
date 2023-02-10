#!/bin/bash

# Index reference genome
biscuit index my_reference.fa

# Align fastq paired-end reads against reference genome, sort and index resulting BAM alignment file (removing duplicate reads with samblaster).
biscuit align my_reference.fa read1.fq.gz read2.fq.gz | samblaster | samtools sort --write-index -o my_output.bam -O BAM -

# Create a pileup
biscuit pileup -o my_pileup.vcf my_reference.fa my_output.bam
bgzip my_pileup.vcf
tabix -p vcf my_pileup.vcf.gz

# Variant calling
biscuit vcf2bed -t snp my_pileup.vcf.gz > snps.bed

# Parse pileup into the Single Fragment Epiread format
mkdir tmp
biscuit epiread -o temp.epiread -O -B snps.bed my_reference.fa my_output.bam > temp.epiread 
sort -T ./tmp -k2,2 -k3,3n temp.epiread | awk 'BEGIN{{ qname="" ; rec="" }} 
			qname == $2 {{ print rec"\t"$5"\t"$6"\t"$7"\t"$8 ; qname="" }} 
			qname != $2 {{ qname=$2 ; rec=$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8 ; pair=$3}}' > single_fragment.epiread
rm temp.epiread

# Count Cs and Ts in Fwd and Rv read
awk '{ print $1 "\t" $3 "\t" gsub("C","C",$4) "\t" gsub("T","T",$4) "\t" gsub("C","C",$8) "\t" gsub("T","T",$8)}' single_fragment.epiread > CT_reads.txt

# Sort CT file by chromosome and location
sort -k1,1 -k2n CT_reads.txt > CT_reads_sorted.txt

# Create a CHR directory, copy the sorted CT reads and split the file into different chromosomes. Remove copy of the sorted reads.
mkdir CHR
cp CT_reads_sorted.txt ./CHR
cd CHR
awk -F '\t' '{print>$1}' CT_reads_sorted.txt
rm CT_reads_sorted.txt


