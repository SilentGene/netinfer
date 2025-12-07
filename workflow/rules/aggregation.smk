"""
Network aggregation rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified
taxonomy_in = str(config.get("input", {}).get("taxonomy_table", "")).strip()
infer_tax = config.get("infer_taxonomy", False)
has_taxonomy = bool(taxonomy_in) or infer_tax

rule aggregate_networks:
    input:
        flashweave = f"{outdir}/networks/flashweave/normal/network.tsv" if config["flashweave"]["enabled"] else [],
        flashweaveHE = f"{outdir}/networks/flashweave/HE/network.tsv" if config["flashweaveHE"]["enabled"] else [],
        fastspar = f"{outdir}/networks/fastspar/network.tsv" if config["fastspar"]["enabled"] else [],
        spearman = f"{outdir}/networks/spearman/network.tsv" if config["spearman"]["enabled"] else [],
        spieceasi = f"{outdir}/networks/spieceasi/network.tsv" if config["spieceasi"]["enabled"] else [],
        propr = f"{outdir}/networks/propr/network.tsv" if config["propr"]["enabled"] else [],
        jaccard = f"{outdir}/networks/jaccard/network.tsv" if config["jaccard"]["enabled"] else [],
        taxonomy = f"{outdir}/preprocessed/processed_taxonomy.tsv" if has_taxonomy else [],
        abundance = f"{outdir}/preprocessed/filtered_abundance.tsv"
    output:
        combined_table = f"{outdir}/networks/aggregated_network.tsv",
        combined_graph = f"{outdir}/networks/aggregated_network.gml",
        stats = f"{outdir}/networks/network_statistics.json",
        diff_phyla_table = f"{outdir}/networks/different_phyla_associations.tsv" if has_taxonomy else []
    threads: 1
    params:
        trusted_methods = config.get("trusted_methods")
    log:
        f"{outdir}/logs/aggregate_networks.log"
    script:
        "../scripts/aggregate_networks.py"
