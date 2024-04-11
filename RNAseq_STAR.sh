#!/bin/bash

SECONDS=0

# Path to working directory
cd /home/affine/Documents/DRISHTI/RNASEQ/

# Step-1: Run the FASTQC to check the quality of the reads

mkdir fastqc_output
fastqc data/demo.fastq -o fastqc_output

# Step-2: Run Trimmomatic to trim the reads with poor quality

mkdir trimmed_reads
java -jar trimmomatic-0.39.jar SE -phred33 data/demo.fastq trimmed_reads/demo_trimmed.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50

# Step-3: Run STAR for alignment

mkdir star_output

# Step-3.1: Indexing the reference genome for STAR

STAR --runMode genomeGenerate --genomeDir reference/genome_index --genomeFastaFiles reference/hg38.fa --sjdbGTFfile gtf_files/hg38.refGene.gtf --sjdbOverhang 100

# Step-3.2: Run alignment

STAR --genomeDir reference/genome_index --readFilesIn trimmed_reads/demo_trimmed.fastq --outFileNamePrefix star_output/demo_ --runThreadN 8 --outSAMtype BAM SortedByCoordinate

# Step-4: Run RSEM for gene expression quantification

mkdir rsem_output

rsem-calculate-expression --bam --paired-end --num-threads 8 star_output/demo_Aligned.toTranscriptome.out.bam reference/rsem_index/hg38 reference/rsem_output/demo_quantification

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
