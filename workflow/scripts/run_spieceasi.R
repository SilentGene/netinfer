#!/usr/bin/env Rscript
# SPIEC-EASI network inference

# Print R version
message(sprintf("Using: %s", R.version.string))
# Print script arguments
message("Script arguments:")
message(paste(commandArgs(trailingOnly = FALSE), collapse = " "))

# Load required libraries
suppressPackageStartupMessages({
    library(SpiecEasi)
    library(jsonlite)
    library(tidyverse)
    library(reshape2)
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

    # Convert to matrix
    otu_mat <- as.matrix(abundance_data_t[,-1])

    # Data integrity checks for SPIEC-EASI
    otu_mat[is.na(otu_mat)] <- 0 # Replace any NAs with 0
    otu_mat <- otu_mat[rowSums(otu_mat) > 0, ] # Remove samples with zero counts
    otu_mat <- otu_mat[, colSums(otu_mat) > 0] # Remove features with zero counts

    message(sprintf(
        "Dimensions after cleaning for SPIEC-EASI: %d samples, %d OTUs",
        nrow(otu_mat), ncol(otu_mat)
    ))


    # Run SPIEC-EASI
    message("Running SPIEC-EASI...")
    # Use 'mb' for many features (n > 1000), 'glasso' for fewer (n < 1000)
    se.out <- spiec.easi(
        otu_mat,
        method = method,
        nlambda = 20,
        lambda.min.ratio = 1e-2,
        sel.criterion = 'stars',
        pulsar.params = list(rep.num = 50, ncores = threads),
        verbose = TRUE
    )
    
    # Function to extract and format the edgelist
    extract_spieceasi_edges <- function(se_model, genomes) {
        adj_matrix <- getRefit(se_model)
        weight_matrix <- symBeta(getOptBeta(se_model), mode = "maxabs")
        
        adj_dense <- as.matrix(adj_matrix)
        weight_dense <- as.matrix(weight_matrix)
        
        adj_long <- melt(adj_dense, varnames = c("source", "target"), value.name = "Connected")
        weights_long <- melt(weight_dense, varnames = c("source", "target"), value.name = "Weight")
        
        adj_long$`source` <- genomes[adj_long$`source`]
        adj_long$`target` <- genomes[adj_long$`target`]
        weights_long$`source` <- genomes[weights_long$`source`]
        weights_long$`target` <- genomes[weights_long$`target`]
        
        merged_edges <- merge(adj_long, weights_long, by = c("source", "target"))
        edges_only <- merged_edges[merged_edges$Connected == 1 & merged_edges$`source` != merged_edges$`target`, ]
        edges_unique <- edges_only[!duplicated(t(apply(edges_only[, 1:2], 1, sort))), ]
        edges_unique <- edges_unique[order(edges_unique$Weight, decreasing = TRUE), ]
        edges_final <- edges_unique[, c("source", "target", "Weight")]
        
        # Filter for strong correlations
        edges_filtered <- edges_final[edges_final$Weight > 0.5, ]
        print(paste("Found", nrow(edges_filtered), "edges with weight > 0.5"))
        return(edges_filtered)
    }

    # Extract edges and save to file
    genomes_used <- colnames(otu_mat)
    edges_spieceasi <- extract_spieceasi_edges(se.out, genomes = genomes_used)
    
    # Save results
    write.table(edges_spieceasi, file=network_file, 
                sep="\t", row.names=FALSE, quote=FALSE)
    
    # Save statistics
    stab_metrics <- getStability(se.out)
    stats <- list(
        optimal_lambda = se.out$opt.lambda,
        stability_metrics = stab_metrics,
        edges_count = nrow(edges_spieceasi),
        method = method
    )
    write_json(stats, stats_file, pretty=TRUE)
    
    message(sprintf("Completed. Network has %d edges.", nrow(edges_spieceasi)))
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