"""
Pearson correlation (RMT-based) network inference rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

rule pearson_network:
    input:
        abundance = f"{outdir}/preprocessed_data/filtered_abundance.tsv"
    output:
        network = f"{outdir}/subtool_outputs/pearson/network.tsv",
        correlation = f"{outdir}/subtool_outputs/pearson/correlation.tsv"
    params:
        fdr_threshold = config["pearson"]["fdr_threshold"],
        rho_threshold = config["pearson"]["rho_threshold"]
    threads: 1
    log:
        f"{outdir}/logs/pearson.log"
    script:
        "../scripts/run_pearson.R"
