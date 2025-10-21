#!/usr/bin/env Rscript
"""
propR network inference
"""

# Load required libraries
suppressPackageStartupMessages({
    library(propr)
    library(jsonlite)
    library(tidyverse)
})

# Read command line arguments from snakemake
log_file <- snakemake@log[[1]]
abundance_file <- snakemake@input[["abundance"]]
network_file <- snakemake@output[["network"]]
stats_file <- snakemake@output[["stats"]]
rho_matrix_file <- snakemake@output[["rho_matrix"]]
rho_threshold <- snakemake@params[["rho_threshold"]]
threads <- snakemake@threads

# Set up logging
log_con <- file(log_file, open="wt")
sink(log_con, type="message")
sink(log_con, type="output")

# Main function
main <- function() {
    message("Starting propR network inference")
    
    # Read abundance data
    abundance_data <- read.table(abundance_file, sep="\t", header=TRUE, row.names=1)
    message(sprintf("Loaded abundance data with dimensions: %d x %d", 
                   nrow(abundance_data), ncol(abundance_data)))
    
    # Calculate proportionality
    prop_obj <- propr(t(abundance_data), 
                     metric="rho", 
                     cutoff=rho_threshold, 
                     p=100)
    
    # Get network edges
    network_df <- data.frame(
        source = prop_obj@pairs$Partner,
        target = prop_obj@pairs$Pair,
        rho = prop_obj@pairs$prop
    )
    
    # Save results
    write.table(network_df, file=network_file, 
                sep="\t", row.names=FALSE, quote=FALSE)
    
    # Save rho matrix
    rho_matrix <- updateCutoff(prop_obj, cutoff=0)@matrix
    colnames(rho_matrix) <- rownames(rho_matrix) <- colnames(abundance_data)
    write.table(rho_matrix, file=rho_matrix_file,
                sep="\t", quote=FALSE)
    
    # Save statistics
    stats <- list(
        edges_count = nrow(network_df),
        rho_threshold = rho_threshold,
        min_rho = min(network_df$rho),
        max_rho = max(network_df$rho),
        mean_rho = mean(network_df$rho)
    )
    write_json(stats, stats_file, pretty=TRUE)
    
    message(sprintf("Completed. Network has %d edges.", nrow(network_df)))
}

# Execute main function with error handling
tryCatch({
    main()
}, error = function(e) {
    message(sprintf("Error in propR analysis: %s", e$message))
    quit(status=1)
})

# Close log connection
sink(type="message")
sink(type="output")
close(log_con)