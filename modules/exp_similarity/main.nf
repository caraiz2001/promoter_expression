#!/usr/bin/env nextflow

process EXP_SIMILARITY {
    
    publishDir "${params.outdir}/expression_similarity", mode: 'copy'

    input:
    val family_id
    path human_transformed
    path mouse_transformed
    val transform_method
    val orthologs_mapped


    output:
    path("${family_id}_${transform_method}_exp_similarity.tsv") , emit: exp_sim
    path("${family_id}_${transform_method}_report.txt"), emit: report
    path("*.png"), emit: plots

    script:
    template 'exp_similarity.py'
}
