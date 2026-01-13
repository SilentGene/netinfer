"""
SPIEC-EASI network inference rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

rule spieceasi_network:
    input:
        abundance = f"{outdir}/preprocessed_data/filtered_abundance.tsv"
    output:
        network = f"{outdir}/subtool_outputs/spieceasi/network.tsv",
        stats = f"{outdir}/subtool_outputs/spieceasi/stats.json"
    params:
        method = config["spieceasi"]["method"],
        weight_threshold = config["spieceasi"]["weight_threshold"]
    threads: config["spieceasi"]["threads"]
    log:
        f"{outdir}/logs/spieceasi.log"
    script:
        "../scripts/run_spieceasi.R"

rule spieceasi_plots:
    input:
        network = f"{outdir}/subtool_outputs/spieceasi/network.tsv",
        stats = f"{outdir}/subtool_outputs/spieceasi/stats.json"
    output:
        stability = f"{outdir}/subtool_outputs/spieceasi/stability_plot.pdf",
        network_viz = f"{outdir}/subtool_outputs/spieceasi/network_plot.pdf"
    threads: 1
    log:
        f"{outdir}/logs/spieceasi_plots.log"
    script:
        "../scripts/plot_spieceasi.R"