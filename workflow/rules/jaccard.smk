"""
Jaccard index network inference rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

rule jaccard_network:
    input:
        abundance = f"{outdir}/preprocessed/filtered_abundance.tsv"
    output:
        network = f"{outdir}/networks/jaccard/network.tsv",
        similarity_matrix = f"{outdir}/networks/jaccard/similarity_matrix.tsv"
    params:
        weight_threshold = config["jaccard"]["weight_threshold"]
    threads: 1
    log:
        f"{outdir}/logs/jaccard.log"
    script:
        "../scripts/run_jaccard.py"

rule jaccard_plots:
    input:
        network = f"{outdir}/networks/jaccard/network.tsv",
        similarity_matrix = f"{outdir}/networks/jaccard/similarity_matrix.tsv"
    output:
        heatmap = f"{outdir}/networks/jaccard/heatmap.pdf",
        distribution = f"{outdir}/networks/jaccard/similarity_distribution.pdf"
    threads: 1
    log:
        f"{outdir}/logs/jaccard_plots.log"
    script:
        "../scripts/plot_jaccard.py"