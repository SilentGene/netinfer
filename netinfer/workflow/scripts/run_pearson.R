#!/usr/bin/env Rscript
# Pearson (RMT-based) network inference using RMThreshold

# Print R version
message(sprintf("Using: %s", R.version.string))

# Load required libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(RMThreshold)
    library(reshape2)
})

# Read command line arguments from snakemake
log_file <- snakemake@log[[1]]
abundance_file <- snakemake@input[["abundance"]]
network_file <- snakemake@output[["network"]]
correlation_file <- snakemake@output[["correlation"]]
# Optional parameters from config
# RMThreshold parameters could be added here if needed

# Set up logging
log_con <- file(log_file, open="wt")
sink(log_con, type="message")
sink(log_con, type="output")

# Main function
main <- function() {
    message("Starting Pearson (RMT-based) network inference")
    
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
    
    # Compute Pearson correlations
    message("Computing Pearson correlations...")
    cor_matrix <- cor(data_mat, method = "pearson")
    diag(cor_matrix) <- 0  # Set diagonal to 0 for RMT analysis
    
    # Use RMThreshold to find the optimal threshold
    message("Finding optimal threshold using Random Matrix Theory (RMThreshold)...")
    
    # rm.get.threshold with robust non-interactive parameters
    rmt_res <- tryCatch({
        rm.get.threshold(
            cor_matrix, 
            interval = c(0.3, 0.8), 
            discard.zeros = TRUE, 
            unfold.method = "spline",
            plot.comp = FALSE,
            plot.spacing = FALSE,
            interactive = FALSE,
            save.fit = FALSE
        )
    }, error = function(e) {
        message(sprintf("RMThreshold (rm.get.threshold) error: %s", e$message))
        return(NULL)
    })
    
    # Robustly handle the result
    optimal_threshold <- if (!is.null(rmt_res) && !is.null(rmt_res$threshold)) {
        rmt_res$threshold
    } else {
        message("RMThreshold failed to find an optimal threshold. Falling back to default: 0.7")
        0.7
    }
    
    # Diagnostic message
    # Sometimes rmt_res$threshold might be a vector, take the first one
    optimal_threshold <- optimal_threshold[1]
    message(sprintf("Optimal RMT threshold determined: %f", optimal_threshold))
    
    # Convert matrices to long format data frames
    cor_df <- melt(cor_matrix, varnames = c("source", "target"), value.name = "Pearson")
    
    # Add diagnostic: show max correlation in the data
    max_cor <- max(abs(cor_df$Pearson), na.rm = TRUE)
    message(sprintf("Maximum absolute correlation in data: %f", max_cor))
    
    # Remove self-correlations and duplicated pairs
    cor_df <- cor_df[cor_df$source != cor_df$target, ]
    cor_df <- cor_df[!duplicated(t(apply(cor_df[, 1:2], 1, sort))), ]
    
    # Filter by the RMT threshold
    filtered_df <- cor_df[abs(cor_df$Pearson) >= optimal_threshold, ]
    
    # Sort by absolute correlation strength
    filtered_df <- filtered_df[order(abs(filtered_df$Pearson), decreasing = TRUE), ]
    
    # Save network (edge list)
    write_tsv(filtered_df, network_file)
    message(sprintf("Completed. Network has %d edges.", nrow(filtered_df)))
    
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
    message(sprintf("Error in Pearson RMT analysis: %s", e$message))
    quit(status=1)
})

# Close log connection
sink(type="message")
sink(type="output")
close(log_con)
