#!/usr/bin/env python3
"""
Preprocess taxonomy table for network inference
This script processes the taxonomy mapping file, ensuring proper formatting
and consistency with the abundance table features.
"""

import pandas as pd
import logging
import yaml
from typing import Dict, Any

def setup_logger(log_file: str) -> logging.Logger:
    """Set up logging to both file and console."""
    logger = logging.getLogger('preprocess_taxonomy')
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

def load_taxonomy_table(file_path: str, logger: logging.Logger = None) -> pd.DataFrame:
    """Load taxonomy table from TSV/CSV file."""
    # Detect separator based on file extension
    sep = ',' if file_path.endswith('.csv') else '\t'
    df = pd.read_csv(file_path, sep=sep, index_col=0)
    
    if logger:
        logger.info(f"Loaded taxonomy table with shape {df.shape}")
        
    # Remove rows starting with 'unmapped' or 'Unmapped'
    unmapped_mask = df.index.str.lower().str.startswith('unmapped')
    unmapped_count = unmapped_mask.sum()
    df = df[~unmapped_mask]
    
    if logger and unmapped_count > 0:
        logger.info(f"Removed {unmapped_count} unmapped row(s)")
        logger.info(f"Table shape after removing unmapped rows: {df.shape}")
    
    return df

def standardize_taxonomy(df: pd.DataFrame, logger: logging.Logger = None) -> pd.DataFrame:
    """Standardize taxonomy format and fill missing values."""
    # Expected taxonomy levels
    tax_levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    
    # If columns are unnamed (numeric), assume standard taxonomy levels
    if all(str(col).isdigit() for col in df.columns):
        df.columns = tax_levels[:len(df.columns)]
        if logger:
            logger.info("Renamed numeric columns to standard taxonomy levels")
    
    # Fill missing values with "unclassified"
    df = df.fillna('unclassified')
    
    # Ensure all taxonomy strings are stripped of whitespace
    df = df.apply(lambda x: x.str.strip())
    
    if logger:
        logger.info(f"Standardized taxonomy table with columns: {', '.join(df.columns)}")
    
    return df

def main(snakemake):
    """Main processing function."""
    # Set up logging
    logger = setup_logger(snakemake.log[0])
    logger.info("Starting taxonomy table preprocessing")
    
    try:
        # Load taxonomy table
        logger.info(f"Loading taxonomy table from {snakemake.input.taxonomy}")
        taxonomy_df = load_taxonomy_table(snakemake.input.taxonomy, logger=logger)
        
        # Standardize taxonomy format
        taxonomy_df = standardize_taxonomy(taxonomy_df, logger=logger)
        
        # Save processed table
        taxonomy_df.to_csv(snakemake.output.processed, sep='\t')
        logger.info(f"Saved processed taxonomy table to {snakemake.output.processed}")
        
    except Exception as e:
        logger.error(f"Error during preprocessing: {str(e)}")
        raise e

if __name__ == "__main__":
    main(snakemake)