#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ────────────────────────────────
// PARAMETERS
// ────────────────────────────────
// --- Input/Output Paths ---
params.processed_glob = "/home/cm/Mutanex/NEW_SRR25739478/Automate/02_Processedreads/processed_data/*_R{1,2}.fastp.fastq.gz"  // Input processed data Path
params.deliverables     = "/home/cm/Mutanex/NEW_SRR25739478/Automate"                                                         // Base directory (Deliverables) for results 

// --- Tool Paths (Singularity Images/Container ) ---
params.bwatool          = "/home/cm/Mutanex/tools/bwa_0.7.18--he4a0461_1.sif"                                                 
params.samtools         = "/home/cm/Mutanex/tools/samtools_1.21--h50ea8bc_0.sif"
params.gatktool         = "/home/cm/Mutanex/tools/gatk_4.1.3.0.sif"

// --- Reference Genome and Known Sites Path ---
// IMPORTANT: GATK requires the reference genome to have a .fai index and a .dict sequence dictionary in the same directory.
params.reference        = "/home/cm/Mutanex/Reference/HG38/hg38.fa"                                                            // Human Reference Genome
params.knownsites_dbsnp = "/home/cm/Mutanex/Reference/dbsnp_ucsc_general/Homo_sapiens_assembly38.dbsnp138.vcf"                 // NCBI Known snp
params.knownsites_mills = "/home/cm/Mutanex/Reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"                        // Insertion and deletion Database

// ────────────────────────────────
// CHANNEL: Find Preprocessed FASTQ Pairs
// ────────────────────────────────
Channel
    .fromFilePairs(params.processed_glob)
    .ifEmpty { error "Cannot find any processed FASTQ files with the glob pattern: ${params.processed_glob}" }
    .set { read_pairs_ch }

// ────────────────────────────────
// WORKFLOW
// ────────────────────────────────
workflow {
    // This workflow is sequential, matching your BASH script.
    // The output of each process becomes the input for the next.
    ALIGNMENT(read_pairs_ch)
    MARK_DUPLICATES(ALIGNMENT.out.sorted_bam)
    BUILD_BAM_INDEX(MARK_DUPLICATES.out.dedup_bam)
    BASE_RECALIBRATOR(BUILD_BAM_INDEX.out.dedup_bam_with_index)
    APPLY_BQSR(BASE_RECALIBRATOR.out.recal_data)
}

// ────────────────────────────────
// PROCESS 1: BWA Alignment and Samtools Sort
// ────────────────────────────────
process ALIGNMENT {
    tag "BWA MEM on ${sampleId}"
    publishDir "${params.deliverables}/03_Alignment", mode: 'copy', pattern: '*.sorted.bam'

    cpus 4
    memory '8 GB'

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), path("${sampleId}.sorted.bam"), emit: sorted_bam

    script:
    """
    set -euo pipefail

    SORTED_BAM="${sampleId}.sorted.bam"
    RG_STRING="@RG\\tID:${sampleId}\\tSM:${sampleId}\\tPL:ILLUMINA\\tLB:${sampleId}\\tPU:${sampleId}"

    ( singularity exec ${params.bwatool} bwa mem -t ${task.cpus} -M -R "\${RG_STRING}" ${params.reference} ${reads[0]} ${reads[1]} | singularity exec ${params.samtools} samtools sort -@ 2 -o \${SORTED_BAM} )
    """
}

// ────────────────────────────────
// PROCESS 2: Mark Duplicates with GATK
// ────────────────────────────────
process MARK_DUPLICATES {
    tag "MarkDuplicates on ${sampleId}"
    publishDir "${params.deliverables}/03_Alignment", mode: 'copy', pattern: '*.dedup.bam'

    cpus 4
    memory '8 GB'

    input:
    tuple val(sampleId), path(sorted_bam)

    output:
    tuple val(sampleId), path("${sampleId}.dedup.bam"), emit: dedup_bam

    script:
    """
    set -euo pipefail
    
    DEDUP_BAM="${sampleId}.dedup.bam"
    singularity exec ${params.gatktool} gatk MarkDuplicatesSpark -I ${sorted_bam} -O \${DEDUP_BAM}
    """
}

// ────────────────────────────────
// PROCESS 3: Build BAM Index
// ────────────────────────────────
process BUILD_BAM_INDEX {
    tag "BuildBamIndex on ${sampleId}"
    // This process produces an index file, which is an intermediate file.
    // We only publish the final BQSR BAM and its index.

    cpus 2
    memory '8 GB'

    input:
    tuple val(sampleId), path(dedup_bam)

    output:
    tuple val(sampleId), path(dedup_bam), path("${dedup_bam.baseName}.bai"), emit: dedup_bam_with_index

    script:
    """
    set -euo pipefail

    singularity exec ${params.gatktool} gatk BuildBamIndex -I ${dedup_bam}
    """
}

// ────────────────────────────────
// PROCESS 4: Base Quality Score Recalibration (BQSR) - Step 1
// ────────────────────────────────
process BASE_RECALIBRATOR {
    tag "BaseRecalibrator on ${sampleId}"
    publishDir "${params.deliverables}/03_Alignment", mode: 'copy', pattern: '*.table'

    cpus 4
    memory '8 GB'

    input:
    // FIX 2: Now correctly receives the BAM and its index file
    tuple val(sampleId), path(dedup_bam), path(index)

    output:
    // FIX 2: Now passes ALL necessary files (BAM, index, and table) to the next step
    tuple val(sampleId), path(dedup_bam), path(index), path("${sampleId}.recal.table"), emit: recal_data

    script:
    """
    set -euo pipefail

    RECAL_TABLE="${sampleId}.recal.table"
    singularity exec ${params.gatktool} gatk --java-options "-Xmx10G" BaseRecalibrator \\
        -I ${dedup_bam} \\
        -R ${params.reference} \\
        --known-sites ${params.knownsites_dbsnp} \\
        --known-sites ${params.knownsites_mills} \\
        -O \${RECAL_TABLE}
    """
}

// ────────────────────────────────
// PROCESS 5: Apply BQSR
// ────────────────────────────────
process APPLY_BQSR {
    tag "ApplyBQSR on ${sampleId}"
    publishDir "${params.deliverables}/03_Alignment", mode: 'copy', pattern: '*.{bam,bai}' // Copies final BAM and its index

    cpus 4
    memory '8 GB'

    input:
    // FIX 2: Now correctly receives the BAM, its index, and the recalibration table
    tuple val(sampleId), path(dedup_bam), path(index), path(recal_table)

    output:
    // This is the final output of the pipeline
    tuple val(sampleId), path("${sampleId}.bqsr.bam"), path("${sampleId}.bqsr.bai")

    script:
    """
    set -euo pipefail

    FINAL_BAM="${sampleId}.bqsr.bam"
    singularity exec ${params.gatktool} gatk ApplyBQSR \\
        -I ${dedup_bam} \\
        -R ${params.reference} \\
        --bqsr-recal-file ${recal_table} \\
        -O \${FINAL_BAM}
    """
}
