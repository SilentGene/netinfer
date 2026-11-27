#!/usr/bin/env python3
"""
Preprocess abundance tables for network inference.

Inputs (at least one; features x samples, TSV/CSV/BIOM):
- a read count table (recommended),
- a relative abundance table (fraction 0-1 or percentage 0-100, detected
  automatically),
- an "other" table with non-negative values that are neither counts nor
  relative abundance (e.g. CoverM trimmed_mean, sandpiper full_coverage).
  It cannot be combined with a count table.

- Inference table: read counts if given, otherwise the other table, otherwise
  relative abundance.
- Relative abundance used for the abundance filter: the given relative
  abundance table if any, otherwise computed per sample from the inference table.
- Display table (dashboard and abundance columns in the final tables):
  relative abundance in percent if given, otherwise the inference table values.

Steps (all thresholds come from the config):
1. Drop 'unmapped'/'undefined' rows.
2. Sample filtering: drop samples with zero total and, if counts are given,
   samples with fewer than `min_sample_depth` reads.
3. Feature filtering:
   - prevalence: fraction of samples with abundance > 0 >= `min_prevalence`
     and number of such samples >= `min_prevalent_samples`
   - peak relative abundance (max across samples) >= `min_peak_rel_abundance`
   - optional legacy raw-value threshold `min_peak_abundance` (deprecated)
4. Drop features with constant values across samples (zero variance).

The inference table keeps its ORIGINAL values.
"""

import logging
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import yaml

DEFAULTS = {
    "min_sample_depth": 1000,
    "min_prevalence": 0.1,
    "min_prevalent_samples": 0,
    "min_peak_rel_abundance": 0.001,
}


def setup_logger(log_file: str) -> logging.Logger:
    """Set up logging to both file and console."""
    logger = logging.getLogger('preprocess_abundance')
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    for handler in (logging.FileHandler(log_file), logging.StreamHandler()):
        handler.setLevel(logging.INFO)
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def load_abundance_table(file_path: str, logger: Optional[logging.Logger] = None) -> pd.DataFrame:
    """Load a Feature x Sample table and return it as Sample x Feature."""
    if file_path.endswith('.biom'):
        try:
            import biom
        except ImportError:
            raise ImportError("biom-format package required for BIOM files")
        table = biom.load_table(file_path)
        # BIOM matrix is observation (feature) x sample, same orientation as TSV/CSV input
        df = pd.DataFrame(
            table.matrix_data.toarray(),
            index=table.ids('observation'),
            columns=table.ids('sample')
        )
    else:
        sep = ',' if file_path.endswith('.csv') else '\t'
        df = pd.read_csv(file_path, sep=sep, index_col=0)

    if logger:
        logger.info(f"Loaded {file_path}: {df.shape[0]} features x {df.shape[1]} samples")

    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    df = df.fillna(0)
    if (df.to_numpy(dtype=float) < 0).any():
        raise ValueError(f"{file_path} contains negative values; provide untransformed data")
    # Transpose to Sample x Feature format for processing
    return df.T


def drop_unmapped(df: pd.DataFrame, logger: Optional[logging.Logger] = None) -> pd.DataFrame:
    """Drop 'unmapped'/'undefined' features from a Sample x Feature table."""
    mask = df.columns.str.lower().isin(['unmapped', 'undefined'])
    if mask.any() and logger:
        logger.info(f"Removed {int(mask.sum())} unmapped/undefined row(s)")
    return df.loc[:, ~mask]


def relab_to_fraction(df: pd.DataFrame, logger: Optional[logging.Logger] = None) -> Tuple[pd.DataFrame, str]:
    """Detect whether relative abundance is a fraction or percentage; return fractions.

    Based on the largest per-sample total: <= 1 means fractions (0-1),
    > 1 means percentages (0-100). Totals > 100 are not relative abundance.
    A tiny tolerance absorbs floating-point rounding.
    """
    max_total = float(df.sum(axis=1).max())
    if max_total <= 1.0 + 1e-6:
        unit = 'fraction'
        frac = df
    elif max_total <= 100.0 * (1 + 1e-3):
        unit = 'percent'
        frac = df / 100.0
    else:
        raise ValueError(
            f"Relative abundance table has sample totals up to {max_total:g} (> 100). "
            "This does not look like relative abundance; did you mean --input_count?"
        )
    if logger:
        logger.info(f"Relative abundance unit detected: {unit} (largest sample total {max_total:g})")
    return frac, unit


def is_integer_table(df: pd.DataFrame) -> bool:
    """True if every value is a whole number (i.e. looks like read counts)."""
    values = df.to_numpy(dtype=float)
    return bool(np.all(values == np.round(values)))


def filter_samples(df: pd.DataFrame, min_sample_depth: float, is_counts: bool,
                   logger: Optional[logging.Logger] = None) -> Tuple[pd.DataFrame, list]:
    """Drop empty samples and, for read counts, low-depth samples."""
    depth = df.sum(axis=1)
    to_drop = depth <= 0
    if is_counts and min_sample_depth > 0:
        to_drop |= depth < min_sample_depth
    elif not is_counts and min_sample_depth > 0 and logger:
        logger.info("No read count table given: skipping min_sample_depth filter")
    removed = [str(s) for s in df.index[to_drop]]
    if removed and logger:
        logger.info(f"Removed {len(removed)} sample(s) with total < {min_sample_depth} "
                    f"(or zero total): {removed}")
    return df.loc[~to_drop], removed


def filter_features(df: pd.DataFrame,
                    min_prevalence: float,
                    min_prevalent_samples: int,
                    min_peak_rel_abundance: float,
                    min_peak_abundance: Optional[float] = None,
                    rel: Optional[pd.DataFrame] = None) -> Tuple[pd.DataFrame, dict]:
    """Filter features of a Sample x Feature table; returns filtered table and step counts.

    `rel` is the relative abundance (fractions) aligned to `df`; if None it is
    computed per sample from `df` using the full feature set.
    """
    num_samples = df.shape[0]
    if rel is None:
        rel = df.div(df.sum(axis=1), axis=0)

    present = (df > 0).sum(axis=0)
    keep_prev = (present / num_samples >= min_prevalence) & (present >= min_prevalent_samples)
    keep_abund = rel.max(axis=0) >= min_peak_rel_abundance
    keep = keep_prev & keep_abund
    if min_peak_abundance is not None:
        keep &= df.max(axis=0) >= min_peak_abundance

    filtered = df.loc[:, keep]
    constant = filtered.nunique(axis=0) <= 1
    filtered = filtered.loc[:, ~constant]

    counts = {
        'failed_prevalence': int((~keep_prev).sum()),
        'failed_peak_rel_abundance': int((~keep_abund).sum()),
        'removed_constant': int(constant.sum()),
    }
    return filtered, counts


def get_param(config: dict, key: str):
    value = config.get(key)
    return DEFAULTS[key] if value is None else value


def input_path(snakemake, name: str) -> Optional[str]:
    """Return the path of an optional named Snakemake input, or None."""
    value = getattr(snakemake.input, name, None)
    if value is None or (not isinstance(value, str) and len(value) == 0):
        return None
    return str(value) if isinstance(value, str) else str(value[0])


def align_relab(relab: pd.DataFrame, ref: pd.DataFrame) -> pd.DataFrame:
    """Align relative abundance (Sample x Feature) to the samples/features of `ref`."""
    missing_samples = [s for s in ref.index if s not in relab.index]
    missing_features = [f for f in ref.columns if f not in relab.columns]
    if missing_samples or missing_features:
        raise ValueError(
            "Read count and relative abundance tables do not match. "
            f"Samples missing from the relative abundance table: {missing_samples[:10]} "
            f"({len(missing_samples)} total); features missing: {missing_features[:10]} "
            f"({len(missing_features)} total)"
        )
    return relab.loc[ref.index, ref.columns]


def main(snakemake):
    """Main processing function."""
    logger = setup_logger(snakemake.log[0])
    logger.info("Starting abundance table preprocessing")
    config = snakemake.config

    min_sample_depth = float(get_param(config, 'min_sample_depth'))
    min_prevalence = float(get_param(config, 'min_prevalence'))
    min_prevalent_samples = int(get_param(config, 'min_prevalent_samples'))
    min_peak_rel_abundance = float(get_param(config, 'min_peak_rel_abundance'))
    legacy_peak = config.get('min_peak_abundance')
    if legacy_peak is not None:
        legacy_peak = float(legacy_peak)
        logger.warning("'min_peak_abundance' is deprecated (it is compared with raw values); "
                       "use 'min_peak_rel_abundance' instead. Still applying it: "
                       f"peak raw value >= {legacy_peak}")

    try:
        count_path = input_path(snakemake, 'counts')
        relab_path = input_path(snakemake, 'relab')
        other_path = input_path(snakemake, 'others')
        if not count_path and not relab_path and not other_path:
            raise ValueError("Provide a read count, relative abundance and/or other abundance table")
        if count_path and other_path:
            raise ValueError("A read count table and an other abundance table cannot be combined")

        counts = drop_unmapped(load_abundance_table(count_path, logger), logger) if count_path else None
        relab_frac, relab_unit = None, None
        if relab_path:
            relab_raw = load_abundance_table(relab_path, logger)
            relab_frac, relab_unit = relab_to_fraction(relab_raw, logger)  # totals include unmapped
            relab_frac = drop_unmapped(relab_frac, logger)

        if counts is not None:
            df = counts
            inference_type = 'read_count'
            logger.info("Network inference will use the read count table")
            if not is_integer_table(counts):
                logger.warning("The read count table contains non-integer values. "
                               "FastSpar and SPIEC-EASI assume integer read counts.")
        elif other_path:
            df = drop_unmapped(load_abundance_table(other_path, logger), logger)
            inference_type = 'other'
            logger.warning("No read count table given: network inference will use the other abundance "
                           "table. FastSpar results are unreliable and SPIEC-EASI results are biased on "
                           "non-count data; raw read counts are recommended.")
        else:
            inference_type = 'relative_abundance'
            df = relab_frac * (100.0 if relab_unit == 'percent' else 1.0)  # original values
            logger.warning("No read count table given: network inference will use relative abundance. "
                           "FastSpar results are unreliable and SPIEC-EASI results are biased on "
                           "non-count data; raw read counts are recommended.")
        input_samples, input_features = df.shape

        df, removed_samples = filter_samples(df, min_sample_depth, counts is not None, logger)
        if df.shape[0] == 0:
            raise ValueError("No samples left after sample filtering")

        rel = align_relab(relab_frac, df) if relab_frac is not None else None
        filtered, step_counts = filter_features(
            df, min_prevalence, min_prevalent_samples, min_peak_rel_abundance, legacy_peak, rel
        )
        n_samples, n_features = filtered.shape
        logger.info(f"Kept {n_features}/{input_features} features and {n_samples}/{input_samples} samples")
        if n_features < 2:
            raise ValueError(f"Only {n_features} feature(s) left after filtering; "
                             "relax the filtering thresholds")
        if n_features > n_samples:
            logger.warning(f"More features ({n_features}) than samples ({n_samples}); "
                           "network inference has limited power. Consider stricter filtering.")

        def save(table: pd.DataFrame, path: str) -> None:
            out = table.T  # back to Feature x Sample
            out.index.name = "#OTU ID"
            out.to_csv(path, sep='\t')

        save(filtered, snakemake.output.filtered)
        logger.info(f"Saved filtered inference table to {snakemake.output.filtered}")

        if relab_frac is not None:
            display = relab_frac.loc[filtered.index, filtered.columns] * 100.0
            display_type = 'relative_abundance_percent'
        else:
            display = filtered
            display_type = 'read_count' if counts is not None else 'other'
        save(display, snakemake.output.display)
        logger.info(f"Saved display table ({display_type}) to {snakemake.output.display}")

        stats = {
            'count_table': count_path,
            'relative_abundance_table': relab_path,
            'relative_abundance_unit_detected': relab_unit,
            'other_table': other_path,
            'inference_value_type': inference_type,
            'display_value_type': display_type,
            'input_samples': input_samples,
            'samples_after_filter': n_samples,
            'removed_samples': removed_samples,
            'input_features': input_features,
            'features_after_filter': n_features,
            'removed_features': input_features - n_features,
            **step_counts,
            'min_sample_depth': min_sample_depth,
            'min_prevalence': min_prevalence,
            'min_prevalent_samples': min_prevalent_samples,
            'min_peak_rel_abundance': min_peak_rel_abundance,
            'min_peak_abundance_deprecated': legacy_peak,
        }
        with open(snakemake.output.stats, 'w') as f:
            yaml.dump(stats, f, sort_keys=False)
        logger.info("Saved preprocessing statistics")

    except Exception as e:
        logger.error(f"Error during preprocessing: {str(e)}")
        raise


if __name__ == "__main__":
    main(snakemake)
