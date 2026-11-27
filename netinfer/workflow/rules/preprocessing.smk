"""
Preprocessing rules for NetInfer pipeline
"""

# Get output directory from config and ensure no whitespace
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

# Resolve inputs from config
from snakemake.exceptions import WorkflowError

count_in = str(config.get("input", {}).get("count_table", "") or "").strip()
relab_in = str(config.get("input", {}).get("relative_abundance_table", "") or "").strip()
other_in = str(config.get("input", {}).get("other_table", "") or "").strip()
taxonomy_in = str(config.get("input", {}).get("taxonomy_table", "")).strip()
metadata_in = str(config.get("input", {}).get("metadata_table", "")).strip()

if str(config.get("input", {}).get("abundance_table", "") or "").strip():
    raise WorkflowError(
        "Config key 'input.abundance_table' is no longer supported. Use 'input.count_table' "
        "(raw read counts) and/or 'input.relative_abundance_table' instead "
        "(CLI: --input_count / --input_relative_abundance)."
    )
if not count_in and not relab_in and not other_in:
    raise WorkflowError(
        "No abundance table given. Set 'input.count_table', 'input.relative_abundance_table' and/or "
        "'input.other_table' in the config, or run via the 'netinfer' CLI with --input_count, "
        "--input_relative_abundance and/or --input_others."
    )
if count_in and other_in:
    raise WorkflowError(
        "'input.count_table' and 'input.other_table' cannot be combined: read counts are always "
        "preferred for network inference. Drop the other table (--input_others)."
    )

# Which values the outputs display: relative abundance if available, otherwise
# read counts, otherwise the other values
display_value_type = "relative_abundance" if relab_in else ("count" if count_in else "other")
# Any input table works for reading feature IDs (taxonomy inference)
abundance_in = count_in or other_in or relab_in

infer_tax = config.get("infer_taxonomy", False)

rule preprocess_abundance:
    input:
        counts = count_in if count_in else [],
        relab = relab_in if relab_in else [],
        others = other_in if other_in else []
    output:
        filtered = f"{outdir}/preprocessed_data/filtered_abundance.tsv",
        display = f"{outdir}/preprocessed_data/display_values.tsv",
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
            metadata = metadata_in,
            abundance = f"{outdir}/preprocessed_data/filtered_abundance.tsv",
            stats = f"{outdir}/preprocessed_data/abundance_stats.tsv"
        output:
            processed = f"{outdir}/preprocessed_data/processed_metadata.tsv"
        threads: 1
        log:
            f"{outdir}/logs/preprocess_metadata.log"
        script:
            "../scripts/preprocess_metadata.py"