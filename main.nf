#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {ALIGN_SEQUENCES} from './workflows/ALIGN_SEQUENCES.nf'
include {TRANSFORM_EXP_DATA} from './workflows/TRANSFORM_EXP_DATA.nf'
include { EXPRESSION_SIMILARITY } from './workflows/EXPRESSION_SIMILARITY.nf'

workflow {
    // human_fasta = Channel.fromPath("${params.human_fasta_path}")
    // mouse_fasta = Channel.fromPath("${params.mouse_fasta_path}")
    // sequence_type = params.sequence_type
    family_id = params.family_id

    human_exp = Channel.fromPath("${params.human_exp_path}")
    mouse_exp = Channel.fromPath("${params.mouse_exp_path}")
    transform_method = params.transform_method
    orthologs = Channel.fromPath("${params.orthologs_mapped}")


    // ALIGN_SEQUENCES(family_id, human_fasta, mouse_fasta, sequence_type)
    TRANSFORM_EXP_DATA(family_id, human_exp, mouse_exp, transform_method, orthologs)
    EXPRESSION_SIMILARITY(family_id, TRANSFORM_EXP_DATA.out.human_transformed, TRANSFORM_EXP_DATA.out.mouse_transformed, orthologs)


}
