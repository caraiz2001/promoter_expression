#!/usr/bin/env nextflow

// Include the necessary modules
include { TRANSFORM_EXPRESSION } from '../modules/transform_expression/main.nf'
include { EXP_SIMILARITY } from '../modules/exp_similarity/main.nf'

workflow EXPRESSION {

    take:
    family_id
    human_exp
    mouse_exp
    transform_method
    orthologs

    main:
    // First, run the TRANSFORM_EXPRESSION process
    //  TRANSFORM_EXPRESSION(family_id, human_exp, mouse_exp, transform_method, orthologs)
    transformed_data = TRANSFORM_EXPRESSION(family_id, human_exp, mouse_exp, transform_method, orthologs)

    // Then, use the output from TRANSFORM_EXPRESSION as input for EXP_SIMILARITY
    EXP_SIMILARITY(
        family_id, 
        transformed_data.human_transformed, 
        transformed_data.mouse_transformed, 
        transform_method, 
        orthologs
    )

    emit:
    human_transformed = transformed_data.human_transformed
    mouse_transformed = transformed_data.mouse_transformed
    exp_sim = EXP_SIMILARITY.out.exp_sim
}


// workflow EXPRESSION_SIMILARITY {

//     take:
//     family_id
//     human_transformed
//     mouse_transformed
//     orthologs
//     transform_method

//     main:
//     EXP_SIMILARITY(family_id, human_transformed, mouse_transformed, orthologs, transform_method)

//     exp_sim = EXP_SIMILARITY.out.exp_sim
// }

// workflow TRANSFORM_EXP_DATA {

//     take:
//     family_id
//     human_exp
//     mouse_exp
//     transform_method
//     orthologs

//     main:
//     TRANSFORM_EXPRESSION(family_id, human_exp, mouse_exp, transform_method, orthologs)

//     emit:
//     human_transformed = TRANSFORM_EXPRESSION.out.human_transformed
//     mouse_transformed = TRANSFORM_EXPRESSION.out.mouse_transformed
// }