#!/usr/bin/env nextflow

process EXP_SIMILARITY {
    
    publishDir "${params.outdir}/expression_similarity", mode: 'copy'

    input:
    val family_id
    path human_tpm
    path mouse_tpm
    val orthologs_mapped
    val transform_method


    output:
    path("${family_id}_${transform_method}_exp_similarity.tsv") , emit: exp_sim

    script:
    template 'exp_similarity.py'
}