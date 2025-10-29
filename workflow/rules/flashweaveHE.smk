"""
FlashWeave network inference rules (Heterogeneous mode)
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

rule flashweave_he_network:
    input:
        abundance = f"{outdir}/preprocessed/filtered_abundance.tsv"
    output:
        network = f"{outdir}/networks/flashweave/HE/network.tsv",
        graph = f"{outdir}/networks/flashweave/HE/network.gml"
    params:
        pvalue = config["flashweaveHE"]["pvalue_threshold"],
        weight = config["flashweaveHE"]["weight_threshold"],
        script = workflow.source_path("../scripts/run_flashweave.jl")
    threads: 1
    log:
        f"{outdir}/logs/flashweave_HE.log"
    shell:
        """
        julia {params.script} \
            --data {input.abundance} \
            --transposed \
            --alpha {params.pvalue} \
            --outtable {output.network} \
            --outgraph {output.graph} \
            2> {log}
        """