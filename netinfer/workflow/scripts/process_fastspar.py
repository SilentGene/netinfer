#!/usr/bin/env python3
"""
Process FastSpar results and create network file
"""

import pandas as pd
import numpy as np
import logging
import sys
from typing import Tuple

def setup_logger(log_file: str) -> logging.Logger:
    """Set up logging to both file and console."""
    logger = logging.getLogger('process_fastspar')
    logger.setLevel(logging.INFO)
    
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

def process_fastspar_results(correlation_file: str,
                           pvalue_file: str,
                           pvalue_threshold: float,
                           weight_threshold: float) -> pd.DataFrame:
    """Process FastSpar correlation and p-value matrices into network edges."""
    # Load correlation and p-value matrices
    corr_df = pd.read_csv(correlation_file, sep='\t', index_col=0)
    pval_df = pd.read_csv(pvalue_file, sep='\t', index_col=0)
    
    # Create empty list for edges
    edges = []
    
    # Process upper triangle of the matrices
    for i in range(len(corr_df.index)):
        for j in range(i + 1, len(corr_df.columns)):
            weight = corr_df.iloc[i, j]
            pvalue = pval_df.iloc[i, j]
            
            # Apply thresholds
            if abs(weight) >= weight_threshold and pvalue <= pvalue_threshold:
                edges.append({
                    'source': corr_df.index[i],
                    'target': corr_df.columns[j],
                    'weight': weight,
                    'pvalue': pvalue
                })
    
    # Convert to DataFrame
    network_df = pd.DataFrame(edges)
    return network_df

def main(snakemake):
    """Main processing function."""
    logger = setup_logger(snakemake.log[0])
    logger.info("Starting FastSpar network processing")
    
    try:
        # Process results
        network_df = process_fastspar_results(
            correlation_file=snakemake.input.correlation,
            pvalue_file=snakemake.input.pvalues,
            pvalue_threshold=snakemake.params.pvalue_threshold,
            weight_threshold=snakemake.params.weight_threshold
        )
        
        # Save network
        network_df.to_csv(snakemake.output.network, sep='\t', index=False)
        logger.info(f"Saved network with {len(network_df)} edges")
        
    except Exception as e:
        logger.error(f"Error during FastSpar processing: {str(e)}")
        raise e

if __name__ == "__main__":
    main(snakemake)