"""
FlashWeave network inference rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

rule flashweave_network:
    input:
        abundance = f"{outdir}/preprocessed/filtered_abundance.tsv"
    output:
        network = f"{outdir}/networks/flashweave/normal/network.tsv",
        graph = f"{outdir}/networks/flashweave/normal/network.gml"
    params:
        pvalue = config["flashweave"]["pvalue_threshold"],
        weight = config["flashweave"]["weight_threshold"],
        script = "scripts/run_flashweave.jl"
    threads: 1
    log:
        f"{outdir}/logs/flashweave_normal.log"
    shell:
        """
        julia {params.script} \
            --data {input.abundance} \
            --transposed \
            --alpha {params.pvalue} \
            --outtable {output.network} \
            --outgraph {output.graph} \
            --no_heterogeneous \
            2> {log}
        """