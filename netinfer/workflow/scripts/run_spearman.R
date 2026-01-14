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
    
    # Add a small pseudocount to avoid issues with zeros
    pseudocount <- 1e-6
    data_spearman <- abundance_data_t[,-1] + pseudocount
    data_spearman

    # Compute Spearman correlations and p-values
    cor_results <- rcorr(as.matrix(data_spearman), type = "spearman")
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

    # Filter for statistically significant (FDR < fdr_threshold) and strong (abs(Spearman) > rho_threshold) correlations
    filtered_df <- merged_df[which(merged_df$FDR < fdr_threshold & abs(merged_df$Spearman) > rho_threshold), ]

    # Sort the results by the absolute strength of the correlation
    filtered_df_sorted <- filtered_df[order(abs(filtered_df$Spearman), decreasing = TRUE), ]
    
    # Save network
    write_tsv(filtered_df_sorted, network_file)

    message(sprintf("Completed. Network has %d edges.", nrow(filtered_df_sorted)))
    message(sprintf("FDR threshold: %f, Correlation threshold: %f", 
                   fdr_threshold, rho_threshold))
    
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
