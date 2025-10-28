#!/usr/bin/env Rscript
"""
Spearman correlation network inference
"""

# Load required libraries
suppressPackageStartupMessages({
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
    abundance_data <- read.table(abundance_file, sep="\t", header=TRUE, row.names=1)
    message(sprintf("Loaded abundance data with dimensions: %d x %d", 
                   nrow(abundance_data), ncol(abundance_data)))
    
    # Add small pseudocount to avoid rank ties
    pseudocount <- 1e-6
    abundance_data <- abundance_data + pseudocount
    
    # Ensure all numeric
    stopifnot(all(sapply(abundance_data, is.numeric)))
    
    # Compute Spearman correlations
    message("Computing Spearman correlations...")
    cor_results <- rcorr(as.matrix(t(abundance_data)), type="spearman")
    cor_matrix <- cor_results$r
    p_matrix <- cor_results$P
    
    # Save correlation and p-value matrices
    write.table(cor_matrix, file=correlation_file, 
                sep="\t", quote=FALSE)
    write.table(p_matrix, file=pvalues_file,
                sep="\t", quote=FALSE)
    
    # Convert to long format
    cor_df <- melt(cor_matrix, varnames=c("source", "target"), value.name="correlation")
    pval_df <- melt(p_matrix, varnames=c("source", "target"), value.name="pvalue")
    
    # Merge correlation and p-values
    network_df <- merge(cor_df, pval_df, by=c("source", "target"))
    
    # Remove self-correlations and duplicates
    network_df <- network_df[network_df$source != network_df$target, ]
    network_df <- network_df[!duplicated(t(apply(network_df[,1:2], 1, sort))), ]
    
    # FDR correction
    network_df$fdr <- p.adjust(network_df$pvalue, method="fdr")
    
    # Filter by FDR and correlation thresholds
    network_df <- network_df[network_df$fdr < fdr_threshold & 
                           abs(network_df$correlation) > rho_threshold, ]
    
    # Sort by absolute correlation
    network_df <- network_df[order(abs(network_df$correlation), decreasing=TRUE), ]
    
    # Save network
    write.table(network_df, file=network_file, 
                sep="\t", row.names=FALSE, quote=FALSE)
    
    message(sprintf("Completed. Network has %d edges.", nrow(network_df)))
    message(sprintf("FDR threshold: %f, Correlation threshold: %f", 
                   fdr_threshold, rho_threshold))
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
