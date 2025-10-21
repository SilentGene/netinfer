"""
FlashWeave network inference rules
"""

def get_flashweave_output():
    """Determine FlashWeave output files based on mode"""
    outputs = []
    if config["flashweave"]["enabled"]:
        if config["flashweave"]["heterogeneous"]:
            outputs.extend([
                "results/networks/flashweave/network_HE.tsv",
                "results/networks/flashweave/stats_HE.json"
            ])
        else:
            outputs.extend([
                "results/networks/flashweave/network.tsv",
                "results/networks/flashweave/stats.json"
            ])
    return outputs

rule flashweave_network:
    input:
        abundance = "results/preprocessed/filtered_abundance.tsv"
    output:
        network = "results/networks/flashweave/network{mode}.tsv",
        stats = "results/networks/flashweave/stats{mode}.json"
    params:
        pvalue = config["flashweave"]["pvalue_threshold"],
        weight = config["flashweave"]["weight_threshold"],
        he_mode = lambda wildcards: wildcards.mode == "_HE"
    conda:
        "../envs/julia.yaml"
    threads: config["flashweave"]["threads"]
    log:
        "results/logs/flashweave{mode}.log"
    script:
        "../scripts/run_flashweave.jl"