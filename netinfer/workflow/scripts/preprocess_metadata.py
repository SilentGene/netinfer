#!/usr/bin/env python3
"""Validate and align sample metadata for FlashWeave, omitting sample IDs."""

import csv
import logging
from collections import Counter
from pathlib import Path


def validate_sample_ids(sample_ids, source):
    """Require nonempty, unique IDs without changing their spelling."""
    if not sample_ids:
        raise ValueError(f"{source} contains no samples")
    if any(not sample_id.strip() for sample_id in sample_ids):
        raise ValueError(f"{source} contains empty sample names")
    duplicates = [name for name, count in Counter(sample_ids).items() if count > 1]
    if duplicates:
        raise ValueError(f"{source} contains duplicate sample names: {duplicates}")


def preprocess_metadata(metadata_path, abundance_path, output_path, logger=None, removed_samples=()):
    """Match metadata rows exactly to the filtered abundance table's columns.

    Read all fields as text to preserve identifiers (e.g. 001 or NA) and
    metadata values. Metadata uses commas for .csv files and tabs otherwise.
    The filtered abundance table is TSV with feature IDs in its first column.
    Metadata rows for samples removed during abundance preprocessing
    (``removed_samples``) are allowed and dropped.
    """
    with open(abundance_path, newline="", encoding="utf-8-sig") as handle:
        abundance_header = next(csv.reader(handle, delimiter="\t"), [])
    samples = abundance_header[1:]
    validate_sample_ids(samples, "Abundance table")

    delimiter = "," if Path(metadata_path).suffix.lower() == ".csv" else "\t"
    with open(metadata_path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader, [])
        if len(header) < 2:
            raise ValueError("Metadata must be a headered TSV or CSV with sample IDs and at least one metadata column")
        if any(not column.strip() for column in header[1:]):
            raise ValueError("Metadata contains empty metadata column names")
        if len(set(header[1:])) != len(header[1:]):
            raise ValueError("Metadata contains duplicate metadata column names")
        rows = []
        for row in reader:
            if len(row) != len(header):
                raise ValueError(
                    f"Metadata line {reader.line_num} has {len(row)} fields; expected {len(header)}"
                )
            rows.append(row)

    metadata_samples = [row[0] for row in rows]
    validate_sample_ids(metadata_samples, "Metadata table")
    abundance_set = set(samples)
    metadata_set = set(metadata_samples)
    missing = [sample for sample in samples if sample not in metadata_set]
    removed_set = set(removed_samples)
    dropped = [sample for sample in metadata_samples if sample in removed_set]
    extra = [sample for sample in metadata_samples
             if sample not in abundance_set and sample not in removed_set]
    if missing or extra:
        raise ValueError(
            f"Metadata sample names do not match abundance table. "
            f"Missing from metadata: {missing}; extra in metadata: {extra}"
        )

    by_sample = {row[0]: row[1:] for row in rows}
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header[1:])
        writer.writerows(by_sample[sample] for sample in samples)
    if logger and dropped:
        logger.info("Dropped %d metadata row(s) for samples removed during preprocessing: %s", len(dropped), dropped)
    if logger:
        logger.info("Validated %d samples; reordered=%s", len(samples), metadata_samples != samples)
        logger.info("Saved metadata without sample IDs to %s", output_path)


def main(snakemake):
    """Run from the preprocess_metadata Snakemake rule."""
    Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("preprocess_metadata")
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(snakemake.log[0])
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(handler)
    try:
        removed_samples = []
        if hasattr(snakemake.input, "stats"):
            import yaml
            with open(snakemake.input.stats) as handle:
                stats = yaml.safe_load(handle) or {}
            removed_samples = [str(s) for s in stats.get("removed_samples") or []]
        preprocess_metadata(
            snakemake.input.metadata,
            snakemake.input.abundance,
            snakemake.output.processed,
            logger,
            removed_samples,
        )
    except Exception:
        logger.exception("Metadata preprocessing failed")
        raise
    finally:
        logger.removeHandler(handler)
        handler.close()


if __name__ == "__main__":
    main(snakemake)
