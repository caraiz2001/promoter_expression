#!/usr/bin/env nextflow

process COMBINE_SIMILARITY {
    
    publishDir "${params.outdir}/combined_similarity", mode: 'copy'

    input:
    val family_id
    path promoter_similarity
    path protein_similarity
    path expression_similarity
    val transform_method
    path orthologs_mapped


    output:
    path("${family_id}_${transform_method}_combined_similarity.tsv") , emit: combined_similarity_df
    path("*.png"), emit: plots

    script:
    template 'combine_similarity.py'
}