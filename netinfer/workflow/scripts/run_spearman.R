#!/usr/bin/env Rscript
# Spearman correlation network inference

# Print R version
message(sprintf("Using: %s", R.version.string))
# Print script arguments
message("Script arguments:")
message(paste(commandArgs(trailingOnly = FALSE), collapse = " "))

# Load required libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(Hmisc)
    library(reshape2)
    library(jsonlite)
})

# Read command line arguments from snakemake
log_file <- snakemake@log[[1]]
abundance_file <- snakemake@input[["abundance"]]
network_file <- snakemake@output[["network"]]
correlation_file <- snakemake@output[["correlation"]]
pvalues_file <- snakemake@output[["pvalues"]]
fdr_threshold <- snakemake@params[["fdr_threshold"]]
rho_threshold <- snakemake@params[["rho_threshold"]]
include_negative <- snakemake@params[["include_negative"]]
threads <- snakemake@threads

# Set up logging
log_con <- file(log_file, open="wt")
sink(log_con, type="message")
sink(log_con, type="output")

# Main function
main <- function() {
    message("Starting Spearman correlation network inference")
    
    # Read abundance data
    abundance_data <- read_tsv(abundance_file)
    message(sprintf("Loaded abundance table with %d (OTUs) x %d (samples)", 
                   nrow(abundance_data), ncol(abundance_data)))
    
    # Ensure output directory exists
    out_dir <- dirname(network_file)
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # transpose data
    abundance_data_t <- abundance_data |>
        pivot_longer(cols = -`#OTU ID`, names_to = "Sample", values_to = "Abundance") |>
        pivot_wider(names_from = `#OTU ID`, values_from = Abundance)
    abundance_data_t
    
    data_mat <- as.matrix(abundance_data_t[,-1])
    rownames(data_mat) <- abundance_data_t$Sample

    # Drop samples without any abundance (cannot be normalised)
    empty_samples <- rowSums(data_mat) <= 0
    if (any(empty_samples)) {
        message(sprintf("Removing %d sample(s) with zero total abundance", sum(empty_samples)))
        data_mat <- data_mat[!empty_samples, , drop = FALSE]
    }

    # Total-sum scaling: convert each sample to relative abundance so that
    # sequencing depth differences do not create spurious rank correlations.
    # (A pseudocount is unnecessary for rank-based correlation.)
    message("Normalising each sample to relative abundance (total-sum scaling)")
    data_spearman <- data_mat / rowSums(data_mat)

    # Compute Spearman correlations and p-values
    cor_results <- rcorr(data_spearman, type = "spearman")
    cor_matrix <- cor_results$r
    p_matrix <- cor_results$P

    # Convert matrices to long format data frames
    cor_df <- melt(cor_matrix, varnames = c("source", "target"), value.name = "Spearman")
    pval_df <- melt(p_matrix, varnames = c("source", "target"), value.name = "P_value")

    # Merge correlation and p-value data
    merged_df <- merge(cor_df, pval_df, by = c("source", "target"))

    # Remove self-correlations and duplicated pairs
    merged_df <- merged_df[merged_df$`source` != merged_df$`target`, ]
    merged_df <- merged_df[!duplicated(t(apply(merged_df[, 1:2], 1, sort))), ]

    # Adjust p-values using False Discovery Rate (FDR) correction
    merged_df$FDR <- p.adjust(merged_df$P_value, method = "fdr")

    # Keep positive correlations by default. In include-negative mode,
    # retain sufficiently strong negative correlations as well.
    threshold_values <- if (include_negative) abs(merged_df$Spearman) else merged_df$Spearman
    filtered_df <- merged_df[which(merged_df$FDR < fdr_threshold & threshold_values >= rho_threshold), ]

    # Sort the results by the absolute strength of the correlation
    filtered_df_sorted <- filtered_df[order(abs(filtered_df$Spearman), decreasing = TRUE), ]
    
    # Save network
    write_tsv(filtered_df_sorted, network_file)

    message(sprintf("Completed. Network has %d edges.", nrow(filtered_df_sorted)))
    message(sprintf("FDR threshold: %f, Correlation threshold: %f, Include negative: %s",
                   fdr_threshold, rho_threshold, include_negative))
    
    # Save full correlation and p-value matrices
    cor_matrix_df <- as.data.frame(cor_matrix)
    cor_matrix_df <- cbind(Taxon=rownames(cor_matrix_df), cor_matrix_df)
    write.table(cor_matrix_df, file=correlation_file, 
                sep="\t", row.names=FALSE, quote=FALSE)
    p_matrix_df <- as.data.frame(p_matrix)
    p_matrix_df <- cbind(Taxon=rownames(p_matrix_df), p_matrix_df)
    write.table(p_matrix_df, file=pvalues_file, 
                sep="\t", row.names=FALSE, quote=FALSE)
}

# Execute main function with error handling
tryCatch({
    main()
}, error = function(e) {
    message(sprintf("Error in Spearman analysis: %s", e$message))
    quit(status=1)
})

# Close log connection
sink(type="message")
sink(type="output")
close(log_con)
