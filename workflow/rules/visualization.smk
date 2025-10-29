"""
Network visualization rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified
taxonomy_in = str(config.get("input", {}).get("taxonomy_table", "")).strip()

# Only include these rules if visualization is enabled
if config.get("visualization", {}).get("enabled", True):

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