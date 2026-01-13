"""
FastSpar network inference rules (aligned with the reference script steps)
Steps:
1) Run FastSpar to compute median correlation/covariance (input already preprocessed upstream)
2) Generate bootstrap tables
3) For each bootstrap, run FastSpar to compute correlations
4) Compute p-values from bootstrap correlations
5) Build final network edges (thresholds applied)
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

# Paths
FS_DIR = f"{outdir}/subtool_outputs/fastspar"
BOOT_DIR = f"{FS_DIR}/bootstrap"
COR_DIR = f"{FS_DIR}/correlations"
fastspar_log_dir = f"{outdir}/logs/fastspar_corr"

# FastSpar parameters (read from config; no internal defaults)
BOOT_N = int(config["fastspar"]["bootstraps"])  # number of bootstrap tables
ITER_N = int(config["fastspar"]["iterations"])   # iterations for fastspar
THREADS_N = int(config["fastspar"]["threads"])   # threads for fastspar
COR_THR = float(config["fastspar"]["correlation_threshold"])  # fastspar --threshold
P_THR = float(config["fastspar"]["pvalue_threshold"])         # final filtering
W_THR = float(config["fastspar"]["weight_threshold_filter"])         # final |R| threshold


rule fastspar_correlation:
    input:
        abundance=f"{outdir}/preprocessed_data/filtered_abundance.tsv"
    output:
        correlation=f"{FS_DIR}/median_correlation.tsv",
        covariance=f"{FS_DIR}/median_covariance.tsv"
    params:
        iterations=ITER_N,
        cor_thr=COR_THR
    threads: THREADS_N
    log:
        f"{outdir}/logs/fastspar_correlation.log"
    shell:
        """
        mkdir -p $(dirname {output.correlation})
        fastspar --otu_table {input.abundance} \
                 --correlation {output.correlation} \
                 --covariance {output.covariance} \
                 --iterations {params.iterations} \
                 --threshold {params.cor_thr} \
                 --threads {threads} \
                 --yes \
                 > {log} 2>&1
        """


rule fastspar_bootstrap:
    input:
        abundance=f"{outdir}/preprocessed_data/filtered_abundance.tsv"
    output:
        boots=expand(f"{BOOT_DIR}/boot_{{i}}.tsv", i=range(BOOT_N))
    params:
        boot=BOOT_N,
        prefix=BOOT_DIR
    threads: 1
    log:
        f"{outdir}/logs/fastspar_bootstrap.log"
    shell:
        """
        mkdir -p {params.prefix}
        mkdir -p {COR_DIR}  # create this for the next rule
        fastspar_bootstrap --otu_table {input.abundance} \
                           --number {params.boot} \
                           --prefix {params.prefix}/boot \
                           > {log} 2>&1
        """


rule fastspar_bootstrap_correlation:
    input:
        boot=f"{BOOT_DIR}/boot_{{i}}.tsv"
    output:
        corr=f"{COR_DIR}/cor_boot_{{i}}.tsv",
        cov=f"{COR_DIR}/cov_boot_{{i}}.tsv"
    params:
        iterations=ITER_N,
        prefix=COR_DIR
    threads: 1
    log:
        f"{fastspar_log_dir}/fastspar_bootcorr_{{i}}.log"
    shell:
        """
        fastspar --otu_table {input.boot} \
                 --correlation {output.corr} \
                 --covariance {output.cov} \
                 --iterations {params.iterations} \
                 --threads 1 \
                 --yes \
                 > {log} 2>&1
        """



rule fastspar_pvalues:
    """Compute p-values from precomputed bootstrap correlations."""
    input:
        abundance=f"{outdir}/preprocessed_data/filtered_abundance.tsv",
        median_corr=f"{FS_DIR}/median_correlation.tsv",
        corr=expand(f"{COR_DIR}/cor_boot_{{i}}.tsv", i=range(BOOT_N))
    output:
        pvalues=f"{FS_DIR}/pvalues.tsv"
    params:
        permutations=BOOT_N,
        threads=THREADS_N
    log:
        f"{outdir}/logs/fastspar_pvalues.log"
    shell:
        """
        fastspar_pvalues --otu_table {input.abundance} \
                         --correlation {input.median_corr} \
                         --prefix {COR_DIR}/cor_boot_ \
                         --permutations {params.permutations} \
                         --threads {params.threads} \
                         --outfile {output.pvalues} \
                         > {log} 2>&1
        """


rule fastspar_network:
    input:
        correlation=f"{FS_DIR}/median_correlation.tsv",
        pvalues=f"{FS_DIR}/pvalues.tsv"
    output:
        network=f"{FS_DIR}/network.tsv"
    params:
        pvalue_threshold=P_THR,
        weight_threshold=W_THR
    threads: 1
    log:
        f"{outdir}/logs/fastspar_network.log"
    script:
        "../scripts/process_fastspar.py"