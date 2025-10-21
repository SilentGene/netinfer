"""
propR network inference rules
"""

rule propr_network:
    input:
        abundance = "results/preprocessed/filtered_abundance.tsv"
    output:
        network = "results/networks/propr/network.tsv",
        stats = "results/networks/propr/stats.json",
        rho_matrix = "results/networks/propr/rho_matrix.tsv"
    params:
        rho_threshold = config["propr"]["rho_threshold"]
    conda:
        "../envs/r.yaml"
    threads: config["propr"]["threads"]
    log:
        "results/logs/propr.log"
    script:
        "../scripts/run_propr.R"

rule propr_plots:
    input:
        network = "results/networks/propr/network.tsv",
        rho_matrix = "results/networks/propr/rho_matrix.tsv"
    output:
        heatmap = "results/networks/propr/heatmap.pdf",
        pca = "results/networks/propr/pca_plot.pdf"
    conda:
        "../envs/r.yaml"
    threads: 1
    log:
        "results/logs/propr_plots.log"
    script:
        "../scripts/plot_propr.R"