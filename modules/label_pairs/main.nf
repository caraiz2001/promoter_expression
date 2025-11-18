process LABEL_PAIRS {
    
    publishDir "${params.outdir}/pairs_labeled", mode: 'copy'

    input:
    val family_id
    path combined_df
    val transform_method
    val label_method


    output:
    path("${family_id}_${transform_method}_similarity_${label_method}.tsv") , emit: labeled_pairs
    path("*.png"), emit: plots

    script:
    template 'pos_neg_label.py'
}