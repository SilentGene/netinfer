#!/usr/bin/env Rscript
# propR network inference

# Load required libraries
suppressPackageStartupMessages({
    library(propr)
    library(jsonlite)
    library(tidyverse)
    library(reshape2)
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
    data_propr <- abundance_data_t[,-1] + pseudocount
    
    # Run propr with Centered Log-Ratio (CLR) transformation
    # The 'propr' function computes proportionality matrix 'rho'
    rho_results <- propr(data_propr, metric = "rho", symmetrize = TRUE)  # use clr transformation by default if no `ivar` is provided
    rho_matrix <- slot(rho_results, "matrix")

    # Convert matrix to a long edgelist format
    rho_long <- melt(rho_matrix, varnames = c("source", "target"), value.name = "Rho")

    # Remove self-correlations and duplicated pairs
    rho_long <- rho_long[rho_long$`source` != rho_long$`target`, ]
    rho_long <- rho_long[!duplicated(t(apply(rho_long[, 1:2], 1, sort))), ]

    # Filter for strong proportionality (abs(rho) > rho_threshold)
    rho_filtered <- rho_long[abs(rho_long$Rho) > rho_threshold, ]
    message(sprintf("Found %d edges with abs(Rho) > %f", nrow(rho_filtered), rho_threshold))

    # Sort by absolute proportionality value
    network_df <- rho_filtered[order(abs(rho_filtered$Rho), decreasing = TRUE), ]

    # Save results
    write.table(network_df, file=network_file, 
                sep="\t", row.names=FALSE, quote=FALSE)

    # Save rho matrix
    rho_matrix_df <- as.data.frame(rho_matrix)
    rho_matrix_df <- cbind(Taxon=rownames(rho_matrix_df), rho_matrix_df)
    write.table(rho_matrix_df, file=rho_matrix_file, 
                sep="\t", row.names=FALSE, quote=FALSE)
    
    # Save statistics
    stats <- list(
        edges_count = nrow(network_df),
        rho_threshold = rho_threshold,
        min_rho = ifelse(nrow(network_df) > 0, min(network_df$rho), NA),
        max_rho = ifelse(nrow(network_df) > 0, max(network_df$rho), NA),
        mean_rho = ifelse(nrow(network_df) > 0, mean(network_df$rho), NA)
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