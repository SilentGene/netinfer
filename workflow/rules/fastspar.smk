"""
FastSpar network inference rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

rule fastspar_correlation:
    input:
        abundance = f"{outdir}/preprocessed/filtered_abundance.tsv"
    output:
        correlation = f"{outdir}/networks/fastspar/correlation.tsv",
        covariance = f"{outdir}/networks/fastspar/covariance.tsv"
    params:
        iterations = config["fastspar"]["iterations"]
    threads: config["fastspar"]["threads"]
    log:
        f"{outdir}/logs/fastspar_correlation.log"
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
        abundance = f"{outdir}/preprocessed/filtered_abundance.tsv"
    output:
        bootstraps = directory(f"{outdir}/networks/fastspar/bootstrap")
    params:
        iterations = config["fastspar"]["iterations"],
        num_bootstraps = 1000
    threads: 1
    log:
        f"{outdir}/logs/fastspar_bootstrap.log"
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
        abundance = f"{outdir}/preprocessed/filtered_abundance.tsv",
        correlation = f"{outdir}/networks/fastspar/correlation.tsv",
        bootstraps = f"{outdir}/networks/fastspar/bootstrap"
    output:
        pvalues = f"{outdir}/networks/fastspar/pvalues.tsv"
    params:
        iterations = config["fastspar"]["iterations"]
    threads: config["fastspar"]["threads"]
    log:
        f"{outdir}/logs/fastspar_pvalues.log"
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
        correlation = f"{outdir}/networks/fastspar/correlation.tsv",
        pvalues = f"{outdir}/networks/fastspar/pvalues.tsv"
    output:
        network = f"{outdir}/networks/fastspar/network.tsv"
    params:
        pvalue_threshold = config["fastspar"]["pvalue_threshold"],
        weight_threshold = config["fastspar"]["weight_threshold"]
    threads: 1
    log:
        f"{outdir}/logs/fastspar_network.log"
    script:
        "../scripts/process_fastspar.py"