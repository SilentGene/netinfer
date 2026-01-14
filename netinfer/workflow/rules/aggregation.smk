"""
Network aggregation rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified
# Optional suffix for outputs
suffix = str(config.get("suffix", "")).strip()
suffix_tag = f"-{suffix}" if suffix else ""
taxonomy_in = str(config.get("input", {}).get("taxonomy_table", "")).strip()
infer_tax = config.get("infer_taxonomy", False)
has_taxonomy = bool(taxonomy_in) or infer_tax

rule aggregate_networks:
    input:
        flashweave = f"{outdir}/subtool_outputs/flashweave/normal/network.tsv" if config["flashweave"]["enabled"] else [],
        flashweaveHE = f"{outdir}/subtool_outputs/flashweave/HE/network.tsv" if config["flashweaveHE"]["enabled"] else [],
        fastspar = f"{outdir}/subtool_outputs/fastspar/network.tsv" if config["fastspar"]["enabled"] else [],
        spearman = f"{outdir}/subtool_outputs/spearman/network.tsv" if config["spearman"]["enabled"] else [],
        spieceasi = f"{outdir}/subtool_outputs/spieceasi/network.tsv" if config["spieceasi"]["enabled"] else [],
        propr = f"{outdir}/subtool_outputs/propr/network.tsv" if config["propr"]["enabled"] else [],
        jaccard = f"{outdir}/subtool_outputs/jaccard/network.tsv" if config["jaccard"]["enabled"] else [],
        pearson = f"{outdir}/subtool_outputs/pearson/network.tsv" if config["pearson"]["enabled"] else [],
        taxonomy = f"{outdir}/preprocessed_data/processed_taxonomy.tsv" if has_taxonomy else [],
        abundance = f"{outdir}/preprocessed_data/filtered_abundance.tsv"
    output:
        combined_table = f"{outdir}/final_results/merged_edges{suffix_tag}.tsv",
        combined_graph = f"{outdir}/final_results/merged_network{suffix_tag}.gml",
        combined_gexf = f"{outdir}/final_results/merged_network{suffix_tag}.gexf",
        stats = f"{outdir}/final_results/network_info{suffix_tag}.json"
    threads: 1
    params:
        trusted_methods = config.get("trusted_methods")
    log:
        f"{outdir}/logs/aggregate_networks.log"
    script:
        "../scripts/aggregate_networks.py"
