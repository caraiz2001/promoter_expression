#!/usr/bin/env nextflow
include {EXP_SIMILARITY} from '../modules/exp_similarity/main.nf'


workflow EXPRESSION_SIMILARITY {

    take:
    family_id
    human_transformed
    mouse_transformed
    transform_method
    orthologs

    main:
    EXP_SIMILARITY(family_id, human_transformed, mouse_transformed, transform_method, orthologs)

    exp_sim = EXP_SIMILARITY.out.exp_sim
}
