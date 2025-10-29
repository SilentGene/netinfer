# Spearman with FDR correction network inference script

# Load necessary packages
if (!require("Hmisc")) install.packages("Hmisc")
if (!require("reshape2")) install.packages("reshape2")
library(Hmisc)
library(reshape2)

# Load OTU/ASV table
# Rows = features (e.g., species), columns = samples
input_file <- ".//combined3proj//combined3proj_MQMAGs-coverM.tsv-abundance.tsv"
otu <- read.table(input_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
# remove the row with "unmapped" as rowname
otu_filtered <- otu[rownames(otu) != "unmapped", ]
# transpose the data
otu_t <- t(otu_filtered)

# Remove rows with blank sample names
data <- otu_t[otu_t[[1]] != "", ]

# Add a small pseudocount to avoid rank ties caused by zero inflation
pseudocount <- 1e-6
data_pseudo <- data + pseudocount

# Ensure all numeric
stopifnot(all(sapply(data_pseudo, is.numeric)))

# Compute Spearman correlations
cor_results <- rcorr(as.matrix(data_pseudo), type = "spearman")
cor_matrix <- cor_results$r
p_matrix <- cor_results$P

# Save full matrices
#write.csv(cor_matrix, "spearman_correlation_matrix.csv")
#write.csv(p_matrix, "spearman_pvalue_matrix.csv")

# Melt to long format
cor_df <- melt(cor_matrix, varnames = c("Taxon A", "Taxon B"), value.name = "Spearman")
pval_df <- melt(p_matrix, varnames = c("Taxon A", "Taxon B"), value.name = "P_value")

# Merge and clean
merged_df <- merge(cor_df, pval_df, by = c("Taxon A", "Taxon B"))
merged_df <- merged_df[merged_df$Taxon A != merged_df$Taxon B, ]
merged_df <- merged_df[!duplicated(t(apply(merged_df[,1:2], 1, sort))), ]

# FDR correction
merged_df$FDR <- p.adjust(merged_df$P_value, method = "fdr")

# Save full result
#write.csv(merged_df, "spearman_full_long_format.csv", row.names = FALSE)

# Filter to FDR < 0.05, spearman > 0.7
# FDR is the adjusted p-value
filtered_df <- merged_df[merged_df$FDR < 0.05 & merged_df$Spearman > 0.7, ]

# Sort by absolute Spearman and save top 1 million
filtered_df_sorted <- filtered_df[order(filtered_df$Spearman, decreasing = TRUE), ]
write.table(filtered_df_sorted, paste0(input_file, "_Spearman0.7_p0.05.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
