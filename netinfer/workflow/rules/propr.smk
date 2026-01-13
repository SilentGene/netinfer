"""
propR network inference rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

rule propr_network:
    input:
        abundance = f"{outdir}/preprocessed_data/filtered_abundance.tsv"
    output:
        network = f"{outdir}/subtool_outputs/propr/network.tsv",
        stats = f"{outdir}/subtool_outputs/propr/stats.json",
        rho_matrix = f"{outdir}/subtool_outputs/propr/rho_matrix.tsv"
    params:
        rho_threshold = config["propr"]["rho_threshold"]
    threads: 1
    log:
        f"{outdir}/logs/propr.log"
    script:
        "../scripts/run_propr.R"

rule propr_plots:
    input:
        network = f"{outdir}/subtool_outputs/propr/network.tsv",
        rho_matrix = f"{outdir}/subtool_outputs/propr/rho_matrix.tsv"
    output:
        heatmap = f"{outdir}/subtool_outputs/propr/heatmap.pdf",
        pca = f"{outdir}/subtool_outputs/propr/pca_plot.pdf"
    threads: 1
    log:
        f"{outdir}/logs/propr_plots.log"
    script:
        "../scripts/plot_propr.R"