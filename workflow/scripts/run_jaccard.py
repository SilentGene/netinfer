#!/usr/bin/env python3
"""
Jaccard index network inference
"""

import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import logging
from typing import Tuple
import time

def setup_logger(log_file: str) -> logging.Logger:
    """Set up logging to both file and console."""
    logger = logging.getLogger('jaccard_network')
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

def calculate_jaccard_similarity(abundance_df: pd.DataFrame) -> pd.DataFrame:
    """Calculate Jaccard similarity matrix."""
    # Convert abundance to presence/absence
    presence_absence = (abundance_df > 0).astype(int)
    
    # Calculate Jaccard distances
    jaccard_dist = pdist(presence_absence.T, metric='jaccard')
    
    # Convert distances to similarities
    jaccard_sim = 1 - jaccard_dist
    
    # Convert to square matrix
    similarity_matrix = pd.DataFrame(squareform(jaccard_sim),
                                   index=abundance_df.columns,
                                   columns=abundance_df.columns)
    
    return similarity_matrix

def create_network(similarity_matrix: pd.DataFrame,
                  weight_threshold: float) -> pd.DataFrame:
    """Create network from similarity matrix."""
    edges = []
    for i in range(len(similarity_matrix.index)):
        for j in range(i + 1, len(similarity_matrix.columns)):
            similarity = similarity_matrix.iloc[i, j]
            
            if similarity >= weight_threshold:
                edges.append({
                    'source': similarity_matrix.index[i],
                    'target': similarity_matrix.columns[j],
                    'similarity': similarity
                })
    
    return pd.DataFrame(edges)

def main(snakemake):
    """Main processing function."""
    logger = setup_logger(snakemake.log[0])
    logger.info("Starting Jaccard similarity network inference")
    
    try:
        # Load abundance data
        abundance_df = pd.read_csv(snakemake.input.abundance, sep='\t', index_col=0)
        logger.info(f"Loaded abundance data with shape {abundance_df.shape}")
        
        # Calculate Jaccard similarities
        start_time = time.time()
        similarity_matrix = calculate_jaccard_similarity(abundance_df)
        logger.info(f"Calculated similarities in {time.time() - start_time:.2f} seconds")
        
        # Create network
        network_df = create_network(
            similarity_matrix=similarity_matrix,
            weight_threshold=snakemake.params.weight_threshold
        )
        logger.info(f"Created network with {len(network_df)} edges")

        # order network_df by weighted similarity descending
        network_df = network_df.sort_values(by='similarity', ascending=False)
        
        # Save results
        similarity_matrix.to_csv(snakemake.output.similarity_matrix, sep='\t')
        network_df.to_csv(snakemake.output.network, sep='\t', index=False)
        
    except Exception as e:
        logger.error(f"Error during Jaccard network inference: {str(e)}")
        raise e

if __name__ == "__main__":
    main(snakemake)