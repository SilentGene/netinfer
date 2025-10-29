"""
Spearman correlation network inference rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

rule spearman_network:
    input:
        abundance = f"{outdir}/preprocessed/filtered_abundance.tsv"
    output:
        network = f"{outdir}/networks/spearman/network.tsv",      # Edge list with correlations and FDR
        correlation = f"{outdir}/networks/spearman/correlation.tsv", # Full correlation matrix
        pvalues = f"{outdir}/networks/spearman/pvalues.tsv"      # Full p-value matrix
    params:
        fdr_threshold = config["spearman"]["fdr_threshold"],   # FDR significance threshold
        rho_threshold = config["spearman"]["rho_threshold"]    # Minimum absolute correlation
    threads: 1
    log:
        f"{outdir}/logs/spearman.log"
    script:
        "../scripts/run_spearman.R"