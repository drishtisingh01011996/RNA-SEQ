params.reads = 'reads/*fastq'
params.reference = 'reference/genome.fa'
params.outputDir = 'output/'

Channel.fromFilePairs(params.reads).set { input_reads }
Channel.fromPath(params.reference).set { input_reference }

process fastqc {
    // container 'biocontainers/fastqc:v0.11.9_cv4' #specify this container only if you want to the docker
    input:
    file(read) from input_reads

    output:
    file("${params.outputDir}/fastqc_reports/${read.baseName}_fastqc.html") into fastqc_reports

    script:
    """
    mkdir -p ${params.outputDir}/fastqc_reports/
    fastqc ${read} --outdir ${params.outputDir}/fastqc_reports/
    
    """
}

process hisat2_index {
    input:
    file(reference) from input_reference

    output:
    file("${params.outputDir}/indexed_genome/genome") into indexed_genome

    script:
    """
    mkdir -p ${params.outputDir}/indexed_genome/
    hisat2-build ${reference} ${params.outputDir}/indexed_genome/genome
    """
}

process hisat2_align {
    input:
    file(read) from input_reads
    file(indexed_genome) from indexed_genome

    output:
    file("${params.outputDir}/aligned_bam/${read.baseName}.bam") into aligned_bam

    script:
    """
    mkdir -p ${params.outputDir}/aligned_bam/
    hisat2 -x ${indexed_genome} -U ${read} | samtools view -bS - > ${params.outputDir}/aligned_bam/${read.baseName}.bam
    """
}

workflow {
    fastqc
    hisat2_index
    hisat2_align
}
