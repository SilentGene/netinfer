"""
FlashWeave network inference rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

"""
FlashWeave network inference rules (Normal mode)
"""

rule flashweave_network:
    input:
        abundance = f"{outdir}/preprocessed/filtered_abundance.tsv"
    output:
        network = f"{outdir}/networks/flashweave/normal/network.tsv",
        graph = f"{outdir}/networks/flashweave/normal/network.gml"
    params:
        pvalue = config["flashweave"]["pvalue_threshold"],
        weight = config["flashweave"]["weight_threshold"],
        extra_args = "--no_heterogeneous"  # Force normal mode
    threads: config["flashweave"]["threads"]
    log:
        f"{outdir}/logs/flashweave_normal.log"
    shell:
        """
        julia {workflow.basedir}/scripts/run_flashweave.jl \
            --data {input.abundance} \
            --alpha {params.pvalue} \
            --output fwNetwork \
            {params.extra_args} \
            2> {log}
        
        # Move output files to correct location
        mkdir -p $(dirname {output.network})
        mv *_edges_info.tsv {output.network}
        mv *.gml {output.graph}
        """