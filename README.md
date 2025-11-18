# Enhancing Non-Coding DNA Alignment Using Gene Expression Data

## Introduction

The alignment of non-coding DNA sequences presents significant challenges due to their variability and lack of conservation compared to protein-coding regions. This project explores an innovative approach to enhancing the alignment of non-coding DNA, particularly promoter regions, by incorporating gene expression data. The ultimate goal is to develop a more informative similarity measure that better reflects evolutionary relationships and functional similarities.

## Motivation

Alignment methods traditionally rely on sequence similarity to determine evolutionary and functional relationships. However, for non-coding DNA sequences such as promoters, conventional alignment techniques often fail to provide meaningful insights. As an example, let us consider orthologs (= same gene in different species) . 
Very closely related, shared origin and presumably same function.

- **Current alignment methods show low accuracy**, as observed in a baseline comparison:
  - **Protein-coding orthologs:** Standard alignment methods correctly identify orthologous relationships with an accuracy of **99.7%**.
  - **Promoter sequences:** The same alignment methods achieve only **43.5% accuracy** in identifying true orthologous relationships.

By incorporating gene expression data, this project aims to achieve a better representation of the promoter sequences create a similarity measure that improves alignment accuracy, making it easier to distinguish between related and unrelated promoter sequences.




3. **Gene Expression Similarity:**
   - Gene expression profiles will be compared across tissues in human and mouse to assess functional relatedness.
   - The hypothesis is that promoters regulating genes with similar expression patterns are more likely to be functionally related.

## Dataset

The starting dataset includes human and mouse one-to-one orthologs within the kinase family. The key components of the dataset are:
- **Mapping Information:** Ortholog relationships between human and mouse genes.
- **Raw Promoter Sequences:** Extracted upstream sequences of orthologous genes.
- **Protein Domain Sequences:** To assess the similarity of downstream coding regions.
- **Gene Expression Data:** Expression profiles across various tissues in both species.

## Repository Structure

The repository follows an organized structure to facilitate reproducibility and scalability:
```
/promoter_expression
│-- data/                  # Raw and processed data (sequences, mapping info, expression data)
|   |-- db_ids/            # Mapped ensemble gene ids (mouse/human orthologs)
|   |-- expression/        # Gene counts
|   |-- sequence/          # Promoter and protein sequences (human/mouse)
│-- modules/               # Nextflow modules for core tasks
│   │-- pairwise_alignment/  # Align promoter or protein sequences
│   │-- transform_expression/ # Transform expression data using several methods 
│   │-- exp_similarity/ # Compute expression-based similarity (cosine and pearson)
|   |-- combine_similarity/ # Combine sequence and expression similarity
|   |-- label_pairs/ # Label positive/negative pairs for training using different criteria
│-- workflows/             # Workflow scripts integrating multiple modules
│-- scripts/               # Additional scripts (data preprocessing, sequence retrieval, etc.)
│-- exploratory_analysis/  # Exploratory data analysis scripts & plots
│-- results/               # Output data (final similarity matrices, figures, reports)
│-- README.rd              # Project documentation
```


## Implementation Plan

### 1. Development of Nextflow Modules
To ensure modularity and reproducibility, core functionalities will be implemented as **Nextflow** modules:
- **Promoter Alignment Module:** Aligns promoter sequences using traditional sequence-based methods.
- **Protein Domain Alignment Module:** Aligns protein domains to provide baseline evolutionary relationships.
- **Expression Similarity Module:** Computes gene expression similarity across tissues and species.

### 2. Workflow Integration
A set of workflows will be developed to automate the full analysis pipeline, allowing for flexibility in selecting different datasets and similarity metrics. These workflows will:
- Process raw data and prepare inputs for alignment.
- Perform multiple alignment strategies and evaluate their performance.
- Generate final similarity matrices integrating sequence-based and expression-based measures.

### 3. Exploratory Analysis
To assess the effectiveness of the proposed approach, exploratory analysis will be conducted. This includes:
- **Visualization of similarity distributions** using various metrics.
- **Comparative analysis** between sequence-based and expression-based similarity.
- **Evaluation of alignment informativeness** in distinguishing true orthologs from unrelated sequences.

## Expected Challenges

- **Developing scalable Nextflow modules** to handle large genomic datasets efficiently.
- **Defining an optimal similarity metric** for expression-based alignment.
- **Ensuring generalizability** beyond the kinase family dataset to other gene families.

## Next Steps

1. **Implement and test Nextflow modules** for promoter alignment, protein domain alignment, and expression similarity.
2. **Optimize similarity metrics** to improve classification accuracy of promoter alignments.
3. **Perform comparative analyses** between sequence-based and expression-based methods.
4. **Expand the dataset** to include additional gene families and species.

## Contributions & Collaboration

This project is a work in progress, and contributions from others in the field are highly welcome. If you have expertise in:
- Genomic sequence alignment
- Gene expression analysis
- Computational pipeline development

please feel free to contribute! If you have suggestions or would like to collaborate, contact [Your Name] at [Your Email].

---
*Author: Cristina Araiz*
*Affiliation: Center for Genomic Regulation*
*Contact: cristina.araiz@crg.eu*

