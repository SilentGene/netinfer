"""
propR network inference rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

rule propr_network:
    input:
        abundance = f"{outdir}/preprocessed/filtered_abundance.tsv"
    output:
        network = f"{outdir}/networks/propr/network.tsv",
        stats = f"{outdir}/networks/propr/stats.json",
        rho_matrix = f"{outdir}/networks/propr/rho_matrix.tsv"
    params:
        rho_threshold = config["propr"]["rho_threshold"]
    threads: config["propr"]["threads"]
    log:
        f"{outdir}/logs/propr.log"
    script:
        "../scripts/run_propr.R"

rule propr_plots:
    input:
        network = f"{outdir}/networks/propr/network.tsv",
        rho_matrix = f"{outdir}/networks/propr/rho_matrix.tsv"
    output:
        heatmap = f"{outdir}/networks/propr/heatmap.pdf",
        pca = f"{outdir}/networks/propr/pca_plot.pdf"
    threads: 1
    log:
        f"{outdir}/logs/propr_plots.log"
    script:
        "../scripts/plot_propr.R"