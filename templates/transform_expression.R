
#!/usr/bin/env Rscript

# ---------------------
# Load data + define paramters + overview
# ---------------------

human_tpm = "${human_tpm}"
mouse_tpm = "${mouse_tpm}"

transform_method = "${transform_method}"
family_id = "${family_id}"›
output_human = paste0(family_id, "_", transform_method, "_exp_human.tsv")
output_mouse = paste0(family_id, "_", transform_method, "_exp_mouse.tsv")

# Load the data
human_exp <- read.csv(human_tpm, sep ="\\t", header = TRUE)
mouse_exp <- read.csv(mouse_tpm, sep ="\\t", header = TRUE)

# Check that column names match between human and mouse
if (!all(colnames(human_exp) == colnames(mouse_exp))) {
  stop("Column names do not match between human and mouse data")
}

original_genes_human <- nrow(human_exp)
original_genes_mouse <- nrow(mouse_exp)

if (length(human_exp\$gene_id) != length(unique(human_exp\$gene_id))){
  # if there are some rows with duplicated gene ids, keep only the first one
    human_exp <- human_exp[!duplicated(human_exp\$gene_id), ]
    cat("There are duplicated gene ids in the human data. Only the first one is kept.\\n")  
}

if (length(mouse_exp\$gene_id) != length(unique(mouse_exp\$gene_id))){
  # if there are some rows with duplicated gene ids, keep only the first one
    mouse_exp <- mouse_exp[!duplicated(mouse_exp\$gene_id), ]
    cat("There are duplicated gene ids in the mouse data. Only the first one is kept.\\n")  
}

unique_genes_human <- nrow(human_exp)
unique_genes_mouse <- nrow(mouse_exp)

# If the first column is not numeric, set it as row names
if (!is.numeric(human_exp[, 1])) {
  rownames(human_exp) <- human_exp[, 1]
  rownames(mouse_exp) <- mouse_exp[, 1]
  human_exp <- human_exp[, -1]
  mouse_exp <- mouse_exp[, -1]
}

# -------------------------
# Apply transformation
# -------------------------

# We will explore different transformation methods
# 1. Log transformation (natural log)
# 2. Centered log-ratio transformation (clr) 
# 3. TPM (no transformation)


# Transform the data based on the selected method on ALL genes

if (transform_method == "log"){
    # Just take the natural log of (tpm + 1)
    human_transformed <- log(human_exp + 1)
    mouse_transformed <- log(mouse_exp + 1)

    print(human_transformed[1:5,1:5])
    
    } else if (transform_method == "clr") {
    # Centered log-ratio transformation
    human_exp_t <- t(human_exp) # Rows are samples and genes are columns
    mouse_exp_t <- t(mouse_exp)

    lgm_human <- apply(log(human_exp_t + 1), 1, mean)
    lgm_mouse <- apply(log(mouse_exp_t + 1), 1, mean)
    
    human_transformed <- log(human_exp_t + 1) - lgm_human
    mouse_transformed <- log(mouse_exp_t + 1) - lgm_mouse

    human_transformed <- t(human_transformed)
    mouse_transformed <- t(mouse_transformed)

    print(human_transformed[1:5,1:5])

    } else if (transform_method == "zscore") {
        stop("Z-score transformation not implemented yet")

    } else if (transform_method == "tpm"){
      human_transformed <- human_exp
      mouse_transformed <- mouse_exp
      print("No transformation applied")

    } else {
    stop("Invalid transformation method")
}

#  --------------------------------
# Filter out non-expressed genes
# --------------------------------

# IMPORTANT: The filtering is done on the TPM values, not on the transformed values

# Define non-expressed genes as those in which all samples have expression below a threshold
expression_threshold = 1

human_non_expressed <- rowSums(human_exp <= expression_threshold) == ncol(human_exp)
mouse_non_expressed <- rowSums(mouse_exp <= expression_threshold) == ncol(mouse_exp)

# human_non_expressed <- FALSE
# mouse_non_expressed <- FALSE

# Filter out non-expressed genes
human_transformed_filtered <- human_transformed[!human_non_expressed, ]
mouse_transformed_filtered <- mouse_transformed[!mouse_non_expressed, ]

# Print the number of non-expressed genes
cat(sprintf("Number of non-expressed genes in human data: %d\\n", sum(human_non_expressed)))
cat(sprintf("Number of non-expressed genes in mouse data: %d\\n", sum(mouse_non_expressed)))

# Print the number of genes after filtering
cat(sprintf("Number of genes after filtering in human data: %d\\n", nrow(human_transformed_filtered)))
cat(sprintf("Number of genes after filtering in mouse data: %d\\n", nrow(mouse_transformed_filtered)))

# Add an extra column (called gene_id) at the beginning of the df to keep the gene ids
human_transformed_filtered <- cbind(gene_id = rownames(human_transformed_filtered), human_transformed_filtered)
mouse_transformed_filtered <- cbind(gene_id = rownames(mouse_transformed_filtered), mouse_transformed_filtered)

print(human_transformed[1:5,1:5])

# Save the transformed data without row names
write.table(human_transformed_filtered, output_human, sep = "\\t", row.names = FALSE)
write.table(mouse_transformed_filtered, output_mouse, sep = "\\t", row.names = FALSE)

# -------------------------------
# Now to evaluate the performance of the transformation
# --------------------------------

evaluate_expression = TRUE

if (evaluate_expression == TRUE){
    orthologs_file = "${orthologs_mapped}"
    orthologs <- read.csv(orthologs_file, sep = "\\t", header = TRUE)

    orthologs_filtered_human <- intersect(orthologs[,1], rownames(human_transformed_filtered))
    orthologs_filtered_mouse <- intersect(orthologs[,2], rownames(mouse_transformed_filtered))

    human_transformed_orthologs <- human_transformed_filtered[orthologs_filtered_human, ]
    mouse_transformed_orthologs <- mouse_transformed_filtered[orthologs_filtered_mouse, ]

    # Print the number of orthologs
    cat(sprintf("Number of human orthologs (filtered): %d\\n", length(orthologs_filtered_human)))
    cat(sprintf("Number of mouse orthologs (filtered): %d\\n", length(orthologs_filtered_mouse)))


    # Convert Gene ID to gene name (using ortholog file)
    rownames(human_transformed_orthologs) <- orthologs[match(rownames(human_transformed_orthologs), orthologs[,1]), 8]
    rownames(mouse_transformed_orthologs) <- orthologs[match(rownames(mouse_transformed_orthologs), orthologs[,2]), 8]

    # Print the number of common orthologs
    common_orthologs <- intersect(rownames(human_transformed_orthologs), rownames(mouse_transformed_orthologs))
    cat(sprintf("Number of common orthologs (after filtering): %d\\n", length(common_orthologs)))

    # Remove rows with NA (genes not mapped to names)
    human_transformed_orthologs <- human_transformed_orthologs[!is.na(rownames(human_transformed_orthologs)), ]
    mouse_transformed_orthologs <- mouse_transformed_orthologs[!is.na(rownames(mouse_transformed_orthologs)), ]

    # Ensure the row names match between human and mouse
    human_transformed_orthologs <- human_transformed_orthologs[common_orthologs, ]
    mouse_transformed_orthologs <- mouse_transformed_orthologs[common_orthologs, ]

    # Print the number of rows in the transformed data
    cat(sprintf("Number of genes in the final human and mouse data: %d\\n", nrow(human_transformed_orthologs)))

    # Compute Pearson correlation
    human_values <- as.numeric(as.matrix(human_transformed_orthologs[, -1]))
    mouse_values <- as.numeric(as.matrix(mouse_transformed_orthologs[, -1]))

    correlation <- cor(human_values, mouse_values, method = "pearson")
    r_squared <- correlation^2

    cat(sprintf("Pearson correlation: %.3f, R² = %.3f  ", correlation, r_squared))
    

}

# ------------------------
# Write report
# ------------------------

report_file = paste0(family_id, "_", transform_method, "_report.txt")
report_lines <- c(
  sprintf("Transformation method: %s", transform_method),
  sprintf("Expression threshold: %.2f", expression_threshold),
  sprintf("Number of genes in the original df (human): %d", original_genes_human),
  sprintf("Number of genes in the original df (mouse): %d", original_genes_mouse),
  sprintf("Number of unique genes (human/mouse): %d, %d", unique_genes_human, unique_genes_mouse),
  sprintf("Number of non-expressed genes in human data: %d", sum(human_non_expressed)),
  sprintf("Number of non-expressed genes in mouse data: %d", sum(mouse_non_expressed)),
  sprintf("Number of genes after filtering in human data: %d", nrow(human_transformed_filtered)),
  sprintf("Number of genes after filtering in mouse data: %d", nrow(mouse_transformed_filtered)),
  sprintf("Number of human orthologs (after filtering): %d", length(orthologs_filtered_human)),
  sprintf("Number of mouse orthologs (after filtering): %d", length(orthologs_filtered_mouse)),
  sprintf("Number of common orthologs (after filtering): %d", length(common_orthologs)),
  sprintf("Number of genes in the final human and mouse data: %d", nrow(human_transformed_orthologs)),
  sprintf("Pearson correlation: %.3f, R² = %.3f", correlation, r_squared)
)
writeLines(report_lines, report_file)
