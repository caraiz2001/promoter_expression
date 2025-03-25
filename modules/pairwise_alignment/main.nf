#!/usr/bin/env nextflow

process PAIRWISE_ALIGNMENT {
    
    publishDir "${params.outdir}/pairwise_alignment", mode: 'copy'

    input:
    val family_id
    path human_fasta
    path mouse_fasta
    val sequence_type


    output:
    path("*_alignment.tsv") , emit: alignment

    script:
    template 'pairwise_alignment.py'
}