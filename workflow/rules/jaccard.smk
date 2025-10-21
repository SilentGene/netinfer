"""
Jaccard index network inference rules
"""

rule jaccard_network:
    input:
        abundance = "results/preprocessed/filtered_abundance.tsv"
    output:
        network = "results/networks/jaccard/network.tsv",
        similarity_matrix = "results/networks/jaccard/similarity_matrix.tsv"
    params:
        weight_threshold = config["jaccard"]["weight_threshold"]
    conda:
        "../envs/python.yaml"
    threads: config["jaccard"]["threads"]
    log:
        "results/logs/jaccard.log"
    script:
        "../scripts/run_jaccard.py"

rule jaccard_plots:
    input:
        network = "results/networks/jaccard/network.tsv",
        similarity_matrix = "results/networks/jaccard/similarity_matrix.tsv"
    output:
        heatmap = "results/networks/jaccard/heatmap.pdf",
        distribution = "results/networks/jaccard/similarity_distribution.pdf"
    conda:
        "../envs/python.yaml"
    threads: 1
    log:
        "results/logs/jaccard_plots.log"
    script:
        "../scripts/plot_jaccard.py"