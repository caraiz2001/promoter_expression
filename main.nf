#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {ALIGN_SEQUENCES} from './workflows/ALIGN_SEQUENCES.nf'
include {TRANSFORM_EXP_DATA} from './workflows/TRANSFORM_EXP_DATA.nf'
include { EXPRESSION_SIMILARITY } from './workflows/EXPRESSION_SIMILARITY.nf'
include {EXPRESSION} from './workflows/EXPRESSION.nf'
include { COMBINED_SIMILARITY } from './workflows/COMBINED_SIMILARITY.nf'
include { LABELED_PAIRS } from './workflows/LABEL_PAIRS.nf'

workflow {
    // human_fasta = Channel.fromPath("${params.human_fasta_path}")
    // mouse_fasta = Channel.fromPath("${params.mouse_fasta_path}")
    // sequence_type = params.sequence_type
    family_id = params.family_id
    transform_method = params.transform_method
    label_method = params.label_method


    human_tpm = Channel.fromPath("${params.human_exp_path}")
    mouse_tpm = Channel.fromPath("${params.mouse_exp_path}")
    
    orthologs = Channel.fromPath("${params.orthologs_mapped}")
    promoter_similarity = Channel.fromPath("${params.promoter_sim_file}")
    protein_similarity = Channel.fromPath("${params.protein_sim_file}")
    // exp_sim_file = Channel.fromPath("${params.exp_sim_file}")




    // ALIGN_SEQUENCES(family_id, human_fasta, mouse_fasta, sequence_type)
    // TRANSFORM_EXP_DATA(family_id, human_tpm, mouse_tpm, transform_method, orthologs)

    // human_transformed = TRANSFORM_EXP_DATA.out.human_transformed
    // human_transformed.view()

    // mouse_transformed = TRANSFORM_EXP_DATA.out.mouse_transformed
    // mouse_transformed.view()

    // EXPRESSION_SIMILARITY(family_id, TRANSFORM_EXP_DATA.out.human_transformed, TRANSFORM_EXP_DATA.out.mouse_transformed, orthologs, transform_method)

    EXPRESSION(family_id, human_tpm, mouse_tpm, transform_method, orthologs)

    COMBINED_SIMILARITY(family_id, promoter_similarity, protein_similarity, EXPRESSION.out.exp_sim, transform_method, orthologs)

    LABELED_PAIRS(family_id, COMBINED_SIMILARITY.out.combined_similarity, transform_method, label_method)




}
