#!/usr/bin/env python3
"""
Preprocess abundance table for network inference
This script filters features in an abundance table based on prevalence and abundance thresholds.
It supports input files in TSV, CSV, or BIOM format and outputs a filtered table along with statistics.
"""

import pandas as pd
import logging
from typing import Tuple
import yaml

def setup_logger(log_file: str) -> logging.Logger:
    """Set up logging to both file and console."""
    logger = logging.getLogger('preprocess_abundance')
    logger.setLevel(logging.INFO)
    
    # File handler
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    # Formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

def load_abundance_table(file_path: str, logger: logging.Logger = None) -> pd.DataFrame:
    """Load abundance table from TSV/CSV/BIOM file and remove unmapped rows."""
    if file_path.endswith('.biom'):
        try:
            import biom
            table = biom.load_table(file_path)
            df = pd.DataFrame(
                table.matrix_data.toarray().T,
                index=table.ids('sample'),
                columns=table.ids('observation')
            )
        except ImportError:
            raise ImportError("biom-format package required for BIOM files")
    else:
        # Detect separator based on file extension
        sep = ',' if file_path.endswith('.csv') else '\t'
        df = pd.read_csv(file_path, sep=sep, index_col=0)
    
    if logger:
        logger.info(f"Loaded table with initial shape {df.shape}")
    
    # Remove rows starting with 'unmapped' or 'Unmapped'
    unmapped_mask = df.index.str.lower().str.startswith('unmapped')
    unmapped_count = unmapped_mask.sum()
    df = df[~unmapped_mask]
    
    if logger and unmapped_count > 0:
        logger.info(f"Removed {unmapped_count} unmapped row(s)")
        logger.info(f"Table shape after removing unmapped rows: {df.shape}")
    
    # Transpose to Sample x Feature format for processing
    return df.T

def filter_features(df: pd.DataFrame, 
                   min_prevalence: float,
                   min_abundance: float) -> Tuple[pd.DataFrame, dict]:
    """Filter features by prevalence and abundance."""
    # Calculate statistics
    total_reads = df.sum(axis=1)
    rel_abundance = df.div(total_reads, axis=0)
    
    prevalence = (df > 0).mean()
    max_abundance = rel_abundance.max()
    
    # Apply filters
    mask = (prevalence >= min_prevalence) & (max_abundance >= min_abundance)
    filtered_df = df.loc[:, mask]
    
    # Gather statistics
    stats = {
        'input_features': df.shape[1],
        'features_after_filter': filtered_df.shape[1],
        'removed_features': df.shape[1] - filtered_df.shape[1],
        'specified_min_prevalence': min_prevalence,
        'specified_min_abundance': min_abundance
    }
    
    return filtered_df, stats

def main(snakemake):
    """Main processing function."""
    # Set up logging
    logger = setup_logger(snakemake.log[0])
    logger.info("Starting abundance table preprocessing")
    
    try:
        # Load abundance table
        logger.info(f"Loading abundance table from {snakemake.input.abundance}")
        abundance_df = load_abundance_table(snakemake.input.abundance, logger=logger)
        
        # Filter features
        filtered_df, stats = filter_features(
            abundance_df,
            min_prevalence=snakemake.config['min_prevalence'],
            min_abundance=snakemake.config['min_abundance']
        )
        
        # Save filtered table with "#OTU ID" as the index name
        # Transpose back to Feature x Sample format for output
        filtered_df = filtered_df.T
        filtered_df.index.name = "#OTU ID"
        filtered_df.to_csv(snakemake.output.filtered, sep='\t')
        logger.info(f"Saved filtered table with {filtered_df.shape[0]} features")
        
        # Save statistics
        with open(snakemake.output.stats, 'w') as f:
            yaml.dump(stats, f)
        logger.info("Saved preprocessing statistics")
        
    except Exception as e:
        logger.error(f"Error during preprocessing: {str(e)}")
        raise e

if __name__ == "__main__":
    main(snakemake)