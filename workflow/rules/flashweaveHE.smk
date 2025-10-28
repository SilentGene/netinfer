"""
FlashWeave network inference rules (Heterogeneous mode)
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

rule flashweave_he_network:
    input:
        abundance = f"{outdir}/preprocessed/filtered_abundance.tsv"
    output:
        network = f"{outdir}/networks/flashweave/he/network.tsv",
        stats = f"{outdir}/networks/flashweave/he/stats.json"
    params:
        pvalue = config["flashweaveHE"]["pvalue_threshold"],
        weight = config["flashweaveHE"]["weight_threshold"]
        # No --no_heterogeneous flag, so it will use default heterogeneous mode
    threads: config["flashweaveHE"]["threads"]
    log:
        f"{outdir}/logs/flashweave_he.log"
    shell:
        """
        julia {workflow.basedir}/scripts/run_flashweave.jl \
            --data {input.abundance} \
            --alpha {params.pvalue} \
            --output network \
            2> {log}
        
        # Move output files to correct location
        mkdir -p $(dirname {output.network})
        mv *_edges_info.tsv {output.network}
        mv *.gml {output.stats}
        """