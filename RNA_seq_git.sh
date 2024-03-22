#!/bin/bash

SECONDS=0

# Path to working directory
cd /home/affine/Documents/DRISHTI/RNASEQ/


# Step-1: Run the FASTQC to check the quality of the reads

mkdir fastqc_output
fastqc data/demo.fastq -o fastqc_output

# Step-2: Run Trimmomatic to trim the reads with poor quality

trimmomatic SE -phred33 data/demo.fastqc demo_trimmed.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
	
# Step-3: # Run hisat2 for alignment

mkdir hisat2

# get genome ndexes from hisat2 downloads or index the fasta using hisat2-build comand

#wget -P hisat2/indexed_genome https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz

# Step-3.1: indexing the reference fasta file

hisat2-build reference/hg38.fa reference/hg38_indexed

# Step-3.2: Run alignment

hisat2 -q --rna-strandness R -x ~/Documents/DRISHTI/RNASEQ/reference/hg38_indexed/hg38_indexed -U ~/Documents/DRISHTI/RNASEQ/data_1/demo.fastq -S ~/Documents/DRISHTI/RNASEQ/hisat2/demo_trimmed.sam

# Step-4.1: run samtools to sort and convert sam to bam

samtools sort -o hisat2/demo_trimmed_sorted.bam ~/Documents/DRISHTI/RNASEQ/hisat2/demo_trimmed.sam
 
# Step-4.2: run samtools to index the bam file

samtools index hisat2/demo_trimmed_sorted.bam 
 
# Step-5: Run feature count
mkdir quants
featureCounts -S 2 -a gtf_files/hg38.refGene.gtf -o quants/featurecounts.txt hisat2/demo_trimmed_sorted.bam 


duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
