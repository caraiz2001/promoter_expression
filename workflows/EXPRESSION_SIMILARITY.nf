#!/usr/bin/env nextflow
include {EXP_SIMILARITY} from '../modules/exp_similarity/main.nf'


workflow EXPRESSION_SIMILARITY {

    take:
    family_id
    human_tpm
    mouse_tpm
    orthologs
    transform_method

    main:
    EXP_SIMILARITY(family_id, human_tpm, mouse_tpm, orthologs, transform_method)

    exp_sim = EXP_SIMILARITY.out.exp_sim
}