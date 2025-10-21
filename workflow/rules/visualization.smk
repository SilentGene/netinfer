"""
Network visualization rules
"""

# Only include these rules if visualization is enabled
if config.get("visualization", {}).get("enabled", True):

    rule aggregate_networks:
    input:
        flashweave = "results/networks/flashweave/network.tsv" if config["flashweave"]["enabled"] else [],
        fastspar = "results/networks/fastspar/network.tsv" if config["fastspar"]["enabled"] else [],
        spearman = "results/networks/spearman/network.tsv" if config["spearman"]["enabled"] else [],
        spieceasi = "results/networks/spieceasi/network.tsv" if config["spieceasi"]["enabled"] else [],
        propr = "results/networks/propr/network.tsv" if config["propr"]["enabled"] else [],
        jaccard = "results/networks/jaccard/network.tsv" if config["jaccard"]["enabled"] else [],
        taxonomy = "results/preprocessed/processed_taxonomy.tsv",
        abundance = "results/preprocessed/filtered_abundance.tsv"
    output:
        network = "results/networks/aggregated_network.tsv",
        stats = "results/networks/network_statistics.json"
    conda:
        "../envs/python.yaml"
    threads: 1
    log:
        "results/logs/aggregate_networks.log"
    script:
        "../scripts/aggregate_networks.py"

    rule generate_html_viewer:
        input:
            network = "results/networks/aggregated_network.tsv",
        stats = "results/networks/network_statistics.json",
        abundance = "results/preprocessed/filtered_abundance.tsv",
        taxonomy = "results/preprocessed/processed_taxonomy.tsv"
    output:
        html = "results/visualization/network_viewer.html",
        data = "results/visualization/network_data.json"
    params:
        max_edges = config["visualization"]["max_edges"],
        node_size_by = config["visualization"]["node_size_by"],
        edge_width_by = config["visualization"]["edge_width_by"],
        color_by = config["visualization"]["color_by"]
    conda:
        "../envs/python.yaml"
    threads: 1
    log:
        "results/logs/generate_html.log"
    script:
        "../scripts/generate_html.py"