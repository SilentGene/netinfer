"""
FlashWeave network inference rules
"""

rule flashweave_network:
    input:
        abundance = "results/preprocessed/filtered_abundance.tsv"
    output:
        network = "results/networks/flashweave/network.tsv",
        stats = "results/networks/flashweave/stats.json"
    params:
        pvalue = config["flashweave"]["pvalue_threshold"],
        weight = config["flashweave"]["weight_threshold"],
        he_mode = config["flashweave"]["heterogeneous"]
    conda:
        "../envs/julia.yaml"
    threads: config["flashweave"]["threads"]
    log:
        "results/logs/flashweave.log"
    script:
        "../scripts/run_flashweave.jl"