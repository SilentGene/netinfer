"""
Jaccard index network inference rules

Note: Jaccard index is a measure of similarity between two sets, and it is defined as the size of the intersection (presence/absence) divided by the size of the union of the two sets. It does not account for the relative abundance of the taxa.
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

rule jaccard_network:
    input:
        abundance = f"{outdir}/preprocessed_data/filtered_abundance.tsv"
    output:
        network = f"{outdir}/subtool_outputs/jaccard/network.tsv",
        similarity_matrix = f"{outdir}/subtool_outputs/jaccard/similarity_matrix.tsv"
    params:
        weight_threshold = config["jaccard"]["weight_threshold"]
    threads: 1
    log:
        f"{outdir}/logs/jaccard.log"
    script:
        "../scripts/run_jaccard.py"

rule jaccard_plots:
    input:
        network = f"{outdir}/subtool_outputs/jaccard/network.tsv",
        similarity_matrix = f"{outdir}/subtool_outputs/jaccard/similarity_matrix.tsv"
    output:
        heatmap = f"{outdir}/subtool_outputs/jaccard/heatmap.pdf",
        distribution = f"{outdir}/subtool_outputs/jaccard/similarity_distribution.pdf"
    threads: 1
    log:
        f"{outdir}/logs/jaccard_plots.log"
    script:
        "../scripts/plot_jaccard.py"