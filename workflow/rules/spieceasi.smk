"""
SPIEC-EASI network inference rules
"""

rule spieceasi_network:
    input:
        abundance = "results/preprocessed/filtered_abundance.tsv"
    output:
        network = "results/networks/spieceasi/network.tsv",
        stats = "results/networks/spieceasi/stats.json"
    params:
        method = config["spieceasi"]["method"],
        weight_threshold = config["spieceasi"]["weight_threshold"]
    conda:
        "../envs/r.yaml"
    threads: config["spieceasi"]["threads"]
    log:
        "results/logs/spieceasi.log"
    script:
        "../scripts/run_spieceasi.R"

rule spieceasi_plots:
    input:
        network = "results/networks/spieceasi/network.tsv",
        stats = "results/networks/spieceasi/stats.json"
    output:
        stability = "results/networks/spieceasi/stability_plot.pdf",
        network_viz = "results/networks/spieceasi/network_plot.pdf"
    conda:
        "../envs/r.yaml"
    threads: 1
    log:
        "results/logs/spieceasi_plots.log"
    script:
        "../scripts/plot_spieceasi.R"