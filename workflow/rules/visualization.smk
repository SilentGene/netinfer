"""
Network visualization rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified
taxonomy_in = str(config.get("input", {}).get("taxonomy_table", "")).strip()

# Only include these rules if visualization is enabled
if config.get("visualization", {}).get("enabled", True):

    rule aggregate_networks:
        input:
            flashweave = f"{outdir}/networks/flashweave/he/network.tsv" if (config["flashweave"]["enabled"] and config["flashweave"]["heterogeneous"]) else 
                          f"{outdir}/networks/flashweave/normal/network.tsv" if config["flashweave"]["enabled"] else [],
            fastspar = f"{outdir}/networks/fastspar/network.tsv" if config["fastspar"]["enabled"] else [],
            spearman = f"{outdir}/networks/spearman/network.tsv" if config["spearman"]["enabled"] else [],
            spieceasi = f"{outdir}/networks/spieceasi/network.tsv" if config["spieceasi"]["enabled"] else [],
            propr = f"{outdir}/networks/propr/network.tsv" if config["propr"]["enabled"] else [],
            jaccard = f"{outdir}/networks/jaccard/network.tsv" if config["jaccard"]["enabled"] else [],
            taxonomy = f"{outdir}/preprocessed/processed_taxonomy.tsv" if taxonomy_in else [],
            abundance = f"{outdir}/preprocessed/filtered_abundance.tsv"
        output:
            network = f"{outdir}/networks/aggregated_network.tsv",
            stats = f"{outdir}/networks/network_statistics.json"
        threads: 1
        log:
            f"{outdir}/logs/aggregate_networks.log"
        script:
            "../scripts/aggregate_networks.py"

    rule generate_html_viewer:
        input:
            network = f"{outdir}/networks/aggregated_network.tsv",
            stats = f"{outdir}/networks/network_statistics.json",
            abundance = f"{outdir}/preprocessed/filtered_abundance.tsv",
            taxonomy = f"{outdir}/preprocessed/processed_taxonomy.tsv" if taxonomy_in else []
        output:
            html = f"{outdir}/visualization/network_viewer.html",
            data = f"{outdir}/visualization/network_data.json"
        params:
            max_edges = config["visualization"]["max_edges"],
            node_size_by = config["visualization"]["node_size_by"],
            edge_width_by = config["visualization"]["edge_width_by"],
            color_by = config["visualization"]["color_by"]
        threads: 1
        log:
            f"{outdir}/logs/generate_html.log"
        script:
            "../scripts/generate_html.py"