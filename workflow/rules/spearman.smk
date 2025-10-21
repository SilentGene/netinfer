"""
Spearman correlation network inference rules
"""

rule spearman_network:
    input:
        abundance = "results/preprocessed/filtered_abundance.tsv"
    output:
        network = "results/networks/spearman/network.tsv",
        correlation = "results/networks/spearman/correlation.tsv",
        pvalues = "results/networks/spearman/pvalues.tsv"
    params:
        fdr_threshold = config["spearman"]["fdr_threshold"],
        rho_threshold = config["spearman"]["rho_threshold"]
    conda:
        "../envs/python.yaml"
    threads: config["spearman"]["threads"]
    log:
        "results/logs/spearman.log"
    script:
        "../scripts/run_spearman.py"