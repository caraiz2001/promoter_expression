#!/usr/bin/env nextflow

process TRANSFORM_EXPRESSION {
    
    publishDir "${params.outdir}/transform_expression", mode: 'copy'

    input:
    val family_id
    path human_tpm
    path mouse_tpm
    val transform_method
    val orthologs_mapped


    output:
    path("*${family_id}_${transform_method}_exp_human.tsv") , emit: human_transformed
    path("*${family_id}_${transform_method}_exp_mouse.tsv") , emit: mouse_transformed
    path("*${family_id}_${transform_method}_report.txt") , emit: report

    script:
    template 'transform_expression.R'
}