"""
FlashWeave network inference rules
"""

"""
FlashWeave network inference rules (Normal mode)
"""

rule flashweave_network:
    input:
        abundance = "results/preprocessed/filtered_abundance.tsv"
    output:
        network = "results/networks/flashweave/normal/network.tsv",
        graph = "results/networks/flashweave/normal/network.gml"
    params:
        pvalue = config["flashweave"]["pvalue_threshold"],
        weight = config["flashweave"]["weight_threshold"],
        extra_args = "--no_heterogeneous"  # Force normal mode
    conda:
        "../envs/julia.yaml"
    threads: config["flashweave"]["threads"]
    log:
        "results/logs/flashweave_normal.log"
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