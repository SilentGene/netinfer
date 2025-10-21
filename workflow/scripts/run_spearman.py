#!/usr/bin/env python3
"""
Spearman correlation network inference
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import logging
from typing import Tuple
import time

def setup_logger(log_file: str) -> logging.Logger:
    """Set up logging to both file and console."""
    logger = logging.getLogger('spearman_network')
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

def calculate_spearman_parallel(abundance_df: pd.DataFrame,
                              chunk_size: int = 100) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate Spearman correlation and p-values in chunks for memory efficiency."""
    features = abundance_df.columns
    n_features = len(features)
    
    # Initialize result matrices
    correlation = pd.DataFrame(np.zeros((n_features, n_features)),
                             index=features, columns=features)
    pvalues = pd.DataFrame(np.ones((n_features, n_features)),
                          index=features, columns=features)
    
    # Process upper triangle in chunks
    for i in range(0, n_features, chunk_size):
        for j in range(i, n_features, chunk_size):
            # Get chunk ranges
            i_end = min(i + chunk_size, n_features)
            j_end = min(j + chunk_size, n_features)
            
            # Calculate correlations for chunk
            for feat1 in features[i:i_end]:
                for feat2 in features[j:j_end]:
                    if feat1 != feat2:
                        rho, pval = stats.spearmanr(abundance_df[feat1], 
                                                  abundance_df[feat2])
                        correlation.loc[feat1, feat2] = rho
                        correlation.loc[feat2, feat1] = rho
                        pvalues.loc[feat1, feat2] = pval
                        pvalues.loc[feat2, feat1] = pval
    
    return correlation, pvalues

def create_network(correlation: pd.DataFrame,
                  pvalues: pd.DataFrame,
                  fdr_threshold: float,
                  rho_threshold: float) -> pd.DataFrame:
    """Create network from correlation and p-value matrices."""
    # Apply FDR correction
    flat_pvals = pvalues.values[np.triu_indices_from(pvalues.values, k=1)]
    _, pvals_adj, _, _ = multipletests(flat_pvals, method='fdr_bh')
    
    # Create pvalue matrix with FDR corrected values
    pvalues_adj = pd.DataFrame(np.zeros_like(pvalues),
                              index=pvalues.index,
                              columns=pvalues.columns)
    pvalues_adj.values[np.triu_indices_from(pvalues_adj.values, k=1)] = pvals_adj
    pvalues_adj = pvalues_adj + pvalues_adj.T
    
    # Create edge list
    edges = []
    for i in range(len(correlation.index)):
        for j in range(i + 1, len(correlation.columns)):
            rho = correlation.iloc[i, j]
            pval = pvalues_adj.iloc[i, j]
            
            if abs(rho) >= rho_threshold and pval <= fdr_threshold:
                edges.append({
                    'source': correlation.index[i],
                    'target': correlation.columns[j],
                    'correlation': rho,
                    'pvalue': pval
                })
    
    return pd.DataFrame(edges)

def main(snakemake):
    """Main processing function."""
    logger = setup_logger(snakemake.log[0])
    logger.info("Starting Spearman correlation network inference")
    
    try:
        # Load abundance data
        abundance_df = pd.read_csv(snakemake.input.abundance, sep='\t', index_col=0)
        logger.info(f"Loaded abundance data with shape {abundance_df.shape}")
        
        # Calculate correlations
        start_time = time.time()
        correlation, pvalues = calculate_spearman_parallel(abundance_df.T)
        logger.info(f"Calculated correlations in {time.time() - start_time:.2f} seconds")
        
        # Create network
        network_df = create_network(
            correlation=correlation,
            pvalues=pvalues,
            fdr_threshold=snakemake.params.fdr_threshold,
            rho_threshold=snakemake.params.rho_threshold
        )
        logger.info(f"Created network with {len(network_df)} edges")
        
        # Save results
        correlation.to_csv(snakemake.output.correlation, sep='\t')
        pvalues.to_csv(snakemake.output.pvalues, sep='\t')
        network_df.to_csv(snakemake.output.network, sep='\t', index=False)
        
    except Exception as e:
        logger.error(f"Error during Spearman network inference: {str(e)}")
        raise e

if __name__ == "__main__":
    main(snakemake)