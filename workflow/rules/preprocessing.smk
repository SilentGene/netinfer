"""
Preprocessing rules for NetInfer pipeline
"""

rule preprocess_abundance:
    input:
        abundance = config["input"]["abundance_table"]
    output:
        filtered = "results/preprocessed/filtered_abundance.tsv",
        stats = "results/preprocessed/abundance_stats.tsv"
    conda:
        "../envs/python.yaml"
    threads: 1
    log:
        "results/logs/preprocess_abundance.log"
    script:
        "../scripts/preprocess_abundance.py"

rule preprocess_taxonomy:
    input:
        taxonomy = config["input"]["taxonomy_table"]
    output:
        processed = "results/preprocessed/processed_taxonomy.tsv"
    conda:
        "../envs/python.yaml"
    threads: 1
    log:
        "results/logs/preprocess_taxonomy.log"
    script:
        "../scripts/preprocess_taxonomy.py"

rule preprocess_metadata:
    input:
        metadata = config["input"]["metadata_table"]
    output:
        processed = "results/preprocessed/processed_metadata.tsv"
    conda:
        "../envs/python.yaml"
    threads: 1
    log:
        "results/logs/preprocess_metadata.log"
    script:
        "../scripts/preprocess_metadata.py"