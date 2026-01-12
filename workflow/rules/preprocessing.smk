"""
Preprocessing rules for NetInfer pipeline
"""

# Get output directory from config and ensure no whitespace
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

# Resolve inputs from config
from snakemake.exceptions import WorkflowError

abundance_in = str(config.get("input", {}).get("abundance_table", "")).strip()
taxonomy_in = str(config.get("input", {}).get("taxonomy_table", "")).strip()
metadata_in = str(config.get("input", {}).get("metadata_table", "")).strip()

if not abundance_in:
    raise WorkflowError(
        "Required config 'input.abundance_table' is empty. Set it in config/config.yaml or run via the 'netinfer' CLI with --input."
    )

infer_tax = config.get("infer_taxonomy", False)

rule preprocess_abundance:
    input:
        abundance = abundance_in
    output:
        filtered = f"{outdir}/preprocessed_data/filtered_abundance.tsv",
        stats = f"{outdir}/preprocessed_data/abundance_stats.tsv"
    threads: 1
    log:
        f"{outdir}/logs/preprocess_abundance.log"
    script:
        "../scripts/preprocess_abundance.py"

if infer_tax:
    rule infer_taxonomy:
        input:
            abundance = abundance_in
        output:
            processed = f"{outdir}/preprocessed_data/processed_taxonomy.tsv"
        threads: 1
        log:
            f"{outdir}/logs/infer_taxonomy.log"
        script:
            "../scripts/infer_taxonomy.py"
elif taxonomy_in:
    rule preprocess_taxonomy:
        input:
            taxonomy = taxonomy_in
        output:
            processed = f"{outdir}/preprocessed_data/processed_taxonomy.tsv"
        threads: 1
        log:
            f"{outdir}/logs/preprocess_taxonomy.log"
        script:
            "../scripts/preprocess_taxonomy.py"

if metadata_in:
    rule preprocess_metadata:
        input:
            metadata = metadata_in
        output:
            processed = f"{outdir}/preprocessed_data/processed_metadata.tsv"
        threads: 1
        log:
            f"{outdir}/logs/preprocess_metadata.log"
        script:
            "../scripts/preprocess_metadata.py"