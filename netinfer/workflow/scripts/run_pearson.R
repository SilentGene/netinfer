#!/usr/bin/env Rscript
# Pearson correlation w/ CLR transformation network inference (Manual thresholding)

# Print R version
message(sprintf("Using: %s", R.version.string))

# Load required libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(reshape2)
    library(Hmisc) # For correlation and p-values
})

# Read command line arguments from snakemake
log_file <- snakemake@log[[1]]
abundance_file <- snakemake@input[["abundance"]]
network_file <- snakemake@output[["network"]]
correlation_file <- snakemake@output[["correlation"]]
fdr_threshold <- snakemake@params[["fdr_threshold"]]
rho_threshold <- snakemake@params[["rho_threshold"]]

# Set up logging
log_con <- file(log_file, open="wt")
sink(log_con, type="message")
sink(log_con, type="output")

# Main function
main <- function() {
    message("Starting Pearson correlation network inference")
    
    # Read abundance data
    abundance_data <- read_tsv(abundance_file, show_col_types = FALSE)
    message(sprintf("Loaded abundance table with %d (OTUs) x %d (samples)", 
                   nrow(abundance_data), ncol(abundance_data) - 1))
    
    # Ensure output directory exists
    out_dir <- dirname(network_file)
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # Transpose data: rows as samples, columns as OTUs
    abundance_data_t <- abundance_data |>
        pivot_longer(cols = -`#OTU ID`, names_to = "Sample", values_to = "Abundance") |>
        pivot_wider(names_from = `#OTU ID`, values_from = Abundance)
    
    # Prepare matrix
    data_mat <- as.matrix(abundance_data_t[,-1])
    rownames(data_mat) <- abundance_data_t$Sample
    
    # Apply CLR transformation
    message("Applying CLR transformation...")
    # Add pseudocount to handle zeros
    pseudocount <- 1e-6
    data_mat_input <- data_mat + pseudocount
    
    # Calculate CLR: log(x) - mean(log(x)) per sample (row)
    # transposing because apply returns (cols x rows) when margin=1
    data_mat_clr <- t(apply(data_mat_input, 1, function(x) {
        lx <- log(x)
        lx - mean(lx)
    }))
    
    # Compute Pearson correlations and p-values using Hmisc::rcorr
    message("Computing Pearson correlations and p-values (on CLR data)...")
    cor_results <- rcorr(data_mat_clr, type = "pearson")
    cor_matrix <- cor_results$r
    p_matrix <- cor_results$P
    
    # Convert matrices to long format data frames
    cor_df <- melt(cor_matrix, varnames = c("source", "target"), value.name = "Pearson")
    pval_df <- melt(p_matrix, varnames = c("source", "target"), value.name = "P_value")
    
    # Merge correlation and p-value data
    merged_df <- merge(cor_df, pval_df, by = c("source", "target"))
    
    # Remove self-correlations and duplicated pairs
    merged_df <- merged_df[merged_df$source != merged_df$target, ]
    merged_df <- merged_df[!duplicated(t(apply(merged_df[, 1:2], 1, sort))), ]
    
    # Adjust p-values using False Discovery Rate (FDR) correction
    merged_df$FDR <- p.adjust(merged_df$P_value, method = "fdr")
    
    # Filter for statistically significant (FDR < fdr_threshold) and strong (abs(Pearson) >= rho_threshold) correlations
    filtered_df <- merged_df[which(merged_df$FDR < fdr_threshold & abs(merged_df$Pearson) >= rho_threshold), ]
    
    # Sort by absolute correlation strength
    filtered_df_sorted <- filtered_df[order(abs(filtered_df$Pearson), decreasing = TRUE), ]
    
    # Diagnostic info
    max_cor <- max(abs(cor_df$Pearson), na.rm = TRUE)
    message(sprintf("Maximum absolute correlation in data: %f", max_cor))
    message(sprintf("FDR threshold: %f, Correlation threshold: %f", fdr_threshold, rho_threshold))
    
    # Save network (edge list)
    write_tsv(filtered_df_sorted, network_file)
    message(sprintf("Completed. Network has %d edges.", nrow(filtered_df_sorted)))
    
    # Save full correlation matrix
    cor_matrix_df <- as.data.frame(cor_matrix)
    cor_matrix_df <- cbind(Taxon=rownames(cor_matrix_df), cor_matrix_df)
    write.table(cor_matrix_df, file=correlation_file, 
                sep="\t", row.names=FALSE, quote=FALSE)
}

# Execute main function with error handling
tryCatch({
    main()
}, error = function(e) {
    message(sprintf("Error in Pearson analysis: %s", e$message))
    quit(status=1)
})

# Close log connection
sink(type="message")
sink(type="output")
close(log_con)
