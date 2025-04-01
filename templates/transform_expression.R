
#!/usr/bin/env Rscript


human_tpm = "${human_tpm}"
mouse_tpm = "${mouse_tpm}"

# Load the data
human_exp <- read.csv(human_tpm, sep ="\\t", header = TRUE)
mouse_exp <- read.csv(mouse_tpm, sep ="\\t", header = TRUE)

# Check that column names match between human and mouse
if (!all(colnames(human_exp) == colnames(mouse_exp))) {
  stop("Column names do not match between human and mouse data")
}

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

# If the first column is not numeric, set it as row names
if (!is.numeric(human_exp[, 1])) {
  rownames(human_exp) <- human_exp[, 1]
  rownames(mouse_exp) <- mouse_exp[, 1]
  human_exp <- human_exp[, -1]
  mouse_exp <- mouse_exp[, -1]
}



# We will explore different transformation methods
# 1. Log transformation (natural log)
# 2. Centered log-ratio transformation (clr) 


transform_method = "${transform_method}"
family_id = "${family_id}"
output_human = paste0(family_id, "_", transform_method, "_exp_human.tsv")
output_mouse = paste0(family_id, "_", transform_method, "_exp_mouse.tsv")


if (transform_method == "log"){
    human_transformed <- log(human_exp + 1)
    mouse_transformed <- log(mouse_exp + 1)
    
    } else if (transform_method == "clr") {
    # Centered log-ratio transformation
    human_exp_t <- t(human_exp) # Rows are samples and genes are columns
    mouse_exp_t <- t(mouse_exp) 
    lgm_human <- apply(log(human_exp_t + 1), 1, mean)
    lgm_mouse <- apply(log(mouse_exp_t + 1), 1, mean)

    print(human_exp_t[1:5,1:5])
    print(lgm_human)

    human_transformed <- log(human_exp_t + 1) - lgm_human
    mouse_transformed <- log(mouse_exp_t + 1) - lgm_mouse

    human_transformed <- t(human_transformed)
    mouse_transformed <- t(mouse_transformed)

    print(human_transformed[1:5,1:5])
    print(mouse_transformed[1:5,1:5])

    } else if (transform_method == "zscore") {
        stop("Z-score transformation not implemented yet")
    } else if (transform_method == "tpm"){
      human_transformed <- human_exp
      mouse_transformed <- mouse_exp
      print("No transformation applied")

    } else {
    stop("Invalid transformation method")
}

# Add an extra column (called gene_id) at the beginning of the df to keep the gene ids
human_transformed <- cbind(gene_id = rownames(human_transformed), human_transformed)
mouse_transformed <- cbind(gene_id = rownames(mouse_transformed), mouse_transformed)

print(human_transformed[1:5,1:5])

# Save the transformed data without row names
write.table(human_transformed, output_human, sep = "\\t", row.names = FALSE)
write.table(mouse_transformed, output_mouse, sep = "\\t", row.names = FALSE)

evaluate_expression = TRUE

if (evaluate_expression == TRUE){
    orthologs_file = "${orthologs_mapped}"
    orthologs <- read.csv(orthologs_file, sep = "\\t", header = TRUE)

    # Filter transformed data to keep only the genes in the orthologs file
    orthologs_human <- intersect(orthologs[,1], rownames(human_transformed))
    orthologs_mouse <- intersect(orthologs[,2], rownames(mouse_transformed))

    human_transformed <- human_transformed[orthologs_human, ]
    mouse_transformed <- mouse_transformed[orthologs_mouse, ]

    print(human_transformed[1:5,1:5])

    # Convert Gene ID to gene name (using ortholog file)
    # Map gene names using the ortholog file
    rownames(human_transformed) <- orthologs[match(rownames(human_transformed), orthologs[,1]), 8]
    rownames(mouse_transformed) <- orthologs[match(rownames(mouse_transformed), orthologs[,2]), 8]

    # Remove rows with NA (genes not mapped to names)
    human_transformed <- human_transformed[!is.na(rownames(human_transformed)), ]
    mouse_transformed <- mouse_transformed[!is.na(rownames(mouse_transformed)), ]

    print(human_transformed[1:5,1:5])

    # Ensure the row names match between human and mouse
    common_genes <- intersect(rownames(human_transformed), rownames(mouse_transformed))
    human_transformed <- human_transformed[common_genes, ]
    mouse_transformed <- mouse_transformed[common_genes, ]

    print(dim(human_transformed))
    print(dim(mouse_transformed))
    print(human_transformed[1:5,1:5])

    # Compute Pearson correlation
    human_values <- as.numeric(as.matrix(human_transformed[, -1]))
    mouse_values <- as.numeric(as.matrix(mouse_transformed[, -1]))

    print(human_values[1:5])
    print(mouse_values[1:5])

    correlation <- cor(human_values, mouse_values, method = "pearson")
    r_squared <- correlation^2

    print(correlation)
    print(r_squared)

    cat(sprintf("Pearson correlation: %.3f, R² = %.3f  ", correlation, r_squared))
}
# Write a report with the results
report_file = paste0(family_id, "_", transform_method, "_report.txt")
writeLines(sprintf("Pearson correlation: %.3f, R² = %.3f", correlation, r_squared), report_file)
