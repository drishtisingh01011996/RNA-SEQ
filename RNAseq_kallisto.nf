// Define input parameters
params.reads = 'data/*_R{1,2}.fastq'
params.reference = 'reference/reference_transcriptome.fasta'
params.outputDir = 'results'

// Set up the working directory
workDir = 'DRISHTI/RNAseq'

// Define channels for input files
Channel.fromFilePairs(params.reads).into { input_reads }

// Define process for FastQC analysis
process fastqc {
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

// Define process for Kallisto alignment
process kallisto_align {
    input:
    file(reference) from params.reference
    file(read) from input_reads

    output:
    file("${params.outputDir}/kallisto_output/${read.baseName}.bam") into kallisto_output

    script:
    """
    mkdir -p ${params.outputDir}/kallisto_output/
    kallisto quant -i ${reference} -o ${params.outputDir}/kallisto_output/${read.baseName} ${read}
    """
}

// Define workflow
workflow {
    // Run FastQC analysis
    fastqc

    // Run Kallisto alignment
    kallisto_align
}

