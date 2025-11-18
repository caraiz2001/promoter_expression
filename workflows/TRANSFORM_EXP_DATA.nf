#!/usr/bin/env nextflow
include {TRANSFORM_EXPRESSION} from '../modules/transform_expression/main.nf'


workflow TRANSFORM_EXP_DATA {

    take:
    family_id
    human_exp
    mouse_exp
    transform_method
    orthologs

    main:
    TRANSFORM_EXPRESSION(family_id, human_exp, mouse_exp, transform_method, orthologs)
    human_transformed = TRANSFORM_EXPRESSION.out.human_transformed
    mouse_transformed = TRANSFORM_EXPRESSION.out.mouse_transformed
    
    emit:
    human_transformed
    mouse_transformed
}