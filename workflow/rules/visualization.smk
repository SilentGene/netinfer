"""
Network visualization rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified
suffix = str(config.get("suffix", "")).strip()
suffix_tag = f"-{suffix}" if suffix else ""
taxonomy_in = str(config.get("input", {}).get("taxonomy_table", "")).strip()
metadata_in = str(config.get("input", {}).get("metadata_table", "")).strip()

# Only include these rules if visualization is enabled
if config.get("visualization", {}).get("enabled", True):

    rule generate_html_viewer:
        input:
            network = f"{outdir}/networks/merged_edges{suffix_tag}.tsv",
            stats = f"{outdir}/networks/network_info{suffix_tag}.json",
            abundance = f"{outdir}/preprocessed/filtered_abundance.tsv",
            taxonomy = f"{outdir}/preprocessed/processed_taxonomy.tsv" if taxonomy_in else [],
            metadata = f"{outdir}/preprocessed/processed_metadata.tsv" if metadata_in else []
        output:
            html = f"{outdir}/visualization/network_viewer.html"
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