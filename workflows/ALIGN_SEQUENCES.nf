#!/usr/bin/env nextflow
include {PAIRWISE_ALIGNMENT} from '../modules/pairwise_alignment/main.nf'


workflow ALIGN_SEQUENCES{

    take:
    family_id
    human_fasta
    mouse_fasta
    sequence_type

    main:
    PAIRWISE_ALIGNMENT(family_id, human_fasta,mouse_fasta, sequence_type)

    alignment = PAIRWISE_ALIGNMENT.out.alignment
}