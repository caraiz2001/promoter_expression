#!/usr/bin/env nextflow
include {LABEL_PAIRS} from '../modules/label_pairs/main.nf'


workflow LABELED_PAIRS{

    take:
    family_id
    combined_df
    transform_method
    label_method

    main:
    LABEL_PAIRS(family_id, combined_df, transform_method, label_method)
    labeled_pairs = LABEL_PAIRS.out.labeled_pairs
    labeled_pairs.view()
    
    emit:
    labeled_pairs
}