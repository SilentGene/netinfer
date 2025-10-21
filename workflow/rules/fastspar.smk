"""
FastSpar network inference rules
"""

rule fastspar_correlation:
    input:
        abundance = "results/preprocessed/filtered_abundance.tsv"
    output:
        correlation = "results/networks/fastspar/correlation.tsv",
        covariance = "results/networks/fastspar/covariance.tsv"
    params:
        iterations = config["fastspar"]["iterations"]
    conda:
        "../envs/fastspar.yaml"
    threads: config["fastspar"]["threads"]
    log:
        "results/logs/fastspar_correlation.log"
    shell:
        """
        fastspar --otu_table {input.abundance} \
                 --correlation {output.correlation} \
                 --covariance {output.covariance} \
                 --iterations {params.iterations} \
                 --threads {threads} \
                 > {log} 2>&1
        """

rule fastspar_bootstrap:
    input:
        abundance = "results/preprocessed/filtered_abundance.tsv"
    output:
        bootstraps = directory("results/networks/fastspar/bootstrap")
    params:
        iterations = config["fastspar"]["iterations"],
        num_bootstraps = 1000
    conda:
        "../envs/fastspar.yaml"
    threads: 1
    log:
        "results/logs/fastspar_bootstrap.log"
    shell:
        """
        mkdir -p {output.bootstraps}
        fastspar_bootstrap --otu_table {input.abundance} \
                          --number {params.num_bootstraps} \
                          --prefix {output.bootstraps}/boot \
                          > {log} 2>&1
        """

rule fastspar_pvalues:
    input:
        abundance = "results/preprocessed/filtered_abundance.tsv",
        correlation = "results/networks/fastspar/correlation.tsv",
        bootstraps = "results/networks/fastspar/bootstrap"
    output:
        pvalues = "results/networks/fastspar/pvalues.tsv"
    params:
        iterations = config["fastspar"]["iterations"]
    conda:
        "../envs/fastspar.yaml"
    threads: config["fastspar"]["threads"]
    log:
        "results/logs/fastspar_pvalues.log"
    shell:
        """
        fastspar_pvalues --otu_table {input.abundance} \
                         --correlation {input.correlation} \
                         --prefix {input.bootstraps}/boot \
                         --permutations {params.iterations} \
                         --threads {threads} \
                         --outfile {output.pvalues} \
                         > {log} 2>&1
        """

rule fastspar_network:
    input:
        correlation = "results/networks/fastspar/correlation.tsv",
        pvalues = "results/networks/fastspar/pvalues.tsv"
    output:
        network = "results/networks/fastspar/network.tsv"
    params:
        pvalue_threshold = config["fastspar"]["pvalue_threshold"],
        weight_threshold = config["fastspar"]["weight_threshold"]
    conda:
        "../envs/python.yaml"
    threads: 1
    log:
        "results/logs/fastspar_network.log"
    script:
        "../scripts/process_fastspar.py"