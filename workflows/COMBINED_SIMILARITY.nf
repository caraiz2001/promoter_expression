#!/usr/bin/env nextflow
include {COMBINE_SIMILARITY} from '../modules/combine_similarity/main.nf'


workflow COMBINED_SIMILARITY {

    take:
    family_id
    promoter_sim_file
    domain_sim_file
    exp_sim_file
    transform_method
    orthologs_mapped

    main:
    COMBINE_SIMILARITY(family_id, promoter_sim_file, domain_sim_file, exp_sim_file, transform_method, orthologs_mapped)
    combined_similarity = COMBINE_SIMILARITY.out.combined_similarity_df
    combined_similarity.view()

    emit:
    combined_similarity
}