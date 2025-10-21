#!/usr/bin/env Rscript
"""
SPIEC-EASI network inference
"""

# Load required libraries
suppressPackageStartupMessages({
    library(SpiecEasi)
    library(jsonlite)
    library(tidyverse)
})

# Read command line arguments from snakemake
log_file <- snakemake@log[[1]]
abundance_file <- snakemake@input[["abundance"]]
network_file <- snakemake@output[["network"]]
stats_file <- snakemake@output[["stats"]]
method <- snakemake@params[["method"]]
weight_threshold <- snakemake@params[["weight_threshold"]]
threads <- snakemake@threads

# Set up logging
log_con <- file(log_file, open="wt")
sink(log_con, type="message")
sink(log_con, type="output")

# Main function
main <- function() {
    message("Starting SPIEC-EASI network inference")
    
    # Read abundance data
    abundance_data <- read.table(abundance_file, sep="\t", header=TRUE, row.names=1)
    message(sprintf("Loaded abundance data with dimensions: %d x %d", 
                   nrow(abundance_data), ncol(abundance_data)))
    
    # Run SPIEC-EASI
    se.out <- spiec.easi(t(abundance_data), 
                         method=method,
                         lambda.min.ratio=1e-3,
                         nlambda=50,
                         pulsar.params=list(thresh=0.05,
                                          subsample.ratio=0.8,
                                          rep.num=50,
                                          ncores=threads))
    
    # Get stability metrics
    stab_metrics <- getStability(se.out)
    
    # Get network matrix
    net_matrix <- getOptMerge(se.out)
    colnames(net_matrix) <- rownames(net_matrix) <- colnames(abundance_data)
    
    # Create edge list
    edges <- which(abs(net_matrix) >= weight_threshold, arr.ind=TRUE)
    network_df <- data.frame(
        source = rownames(net_matrix)[edges[,1]],
        target = colnames(net_matrix)[edges[,2]],
        weight = net_matrix[edges]
    )
    network_df <- network_df[network_df$source < network_df$target,]
    
    # Save results
    write.table(network_df, file=network_file, 
                sep="\t", row.names=FALSE, quote=FALSE)
    
    # Save statistics
    stats <- list(
        optimal_lambda = se.out$opt.lambda,
        stability_metrics = stab_metrics,
        edges_count = nrow(network_df),
        method = method
    )
    write_json(stats, stats_file, pretty=TRUE)
    
    message(sprintf("Completed. Network has %d edges.", nrow(network_df)))
}

# Execute main function with error handling
tryCatch({
    main()
}, error = function(e) {
    message(sprintf("Error in SPIEC-EASI analysis: %s", e$message))
    quit(status=1)
})

# Close log connection
sink(type="message")
sink(type="output")
close(log_con)