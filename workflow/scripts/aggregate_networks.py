#!/usr/bin/env python3
"""
Aggregate networks from different methods and calculate consensus scores
"""

import pandas as pd
import numpy as np
import networkx as nx
import json
import logging
from typing import Dict, List, Tuple
from pathlib import Path

def setup_logger(log_file: str) -> logging.Logger:
    """Set up logging to both file and console."""
    logger = logging.getLogger('aggregate_networks')
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

def load_network_file(file_path: str, method: str) -> pd.DataFrame:
    """Load network file and standardize format."""
    if not Path(file_path).exists():
        return pd.DataFrame()
        
    df = pd.read_csv(file_path, sep='\t')
    
    # Standardize column names
    if 'source' not in df.columns or 'target' not in df.columns:
        if all(col in df.columns for col in ['from', 'to']):
            df = df.rename(columns={'from': 'source', 'to': 'target'})
        else:
            raise ValueError(f"Network file {file_path} must have source/target or from/to columns")
    
    # Add method column
    df['method'] = method
    
    # Ensure source < target for consistency
    df[['source', 'target']] = np.sort(df[['source', 'target']], axis=1)
    
    return df

def calculate_edge_statistics(edge_df: pd.DataFrame,
                            abundance_df: pd.DataFrame) -> pd.DataFrame:
    """Calculate edge statistics from abundance data."""
    stats = []
    
    for _, row in edge_df.iterrows():
        taxon_a = row['source']
        taxon_b = row['target']
        
        # Get abundance vectors
        abund_a = abundance_df[taxon_a]
        abund_b = abundance_df[taxon_b]
        
        # Calculate statistics
        stats.append({
            'source': taxon_a,
            'target': taxon_b,
            'prevalence_a': (abund_a > 0).mean(),
            'prevalence_b': (abund_b > 0).mean(),
            'mean_abundance_a': abund_a.mean(),
            'mean_abundance_b': abund_b.mean(),
            'max_abundance_a': abund_a.max(),
            'max_abundance_b': abund_b.max(),
            'abundance_b_at_max_a': abund_b[abund_a.idxmax()],
            'abundance_a_at_max_b': abund_a[abund_b.idxmax()],
            'sample_at_max_a': abund_a.idxmax(),
            'sample_at_max_b': abund_b.idxmax()
        })
    
    return pd.DataFrame(stats)

def calculate_network_metrics(network_df: pd.DataFrame) -> Dict:
    """Calculate network-level metrics."""
    G = nx.from_pandas_edgelist(network_df, 'source', 'target')
    
    metrics = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'density': nx.density(G),
        'components': nx.number_connected_components(G),
        'largest_component_size': len(max(nx.connected_components(G), key=len)),
        'avg_degree': sum(dict(G.degree()).values()) / G.number_of_nodes(),
        'avg_clustering': nx.average_clustering(G),
        'methods_used': network_df['method'].nunique()
    }
    
    return metrics

def main(snakemake):
    """Main processing function."""
    logger = setup_logger(snakemake.log[0])
    logger.info("Starting network aggregation")
    
    try:
        # Load networks from each method
        networks = {
            'FlashWeave': load_network_file(snakemake.input.flashweave, 'FlashWeave'),
            'FastSpar': load_network_file(snakemake.input.fastspar, 'FastSpar'),
            'Spearman': load_network_file(snakemake.input.spearman, 'Spearman'),
            'SPIEC-EASI': load_network_file(snakemake.input.spieceasi, 'SPIEC-EASI'),
            'propR': load_network_file(snakemake.input.propr, 'propR'),
            'Jaccard': load_network_file(snakemake.input.jaccard, 'Jaccard')
        }
        
        # Combine networks
        combined_df = pd.concat([df for df in networks.values() if not df.empty])
        logger.info(f"Combined {len(combined_df)} edges from all methods")
        
        # Calculate edge consensus
        edge_counts = combined_df.groupby(['source', 'target'])['method'].agg(['unique', 'count'])
        edge_counts.columns = ['methods', 'methods_count']
        edge_counts['methods'] = edge_counts['methods'].apply(lambda x: ','.join(sorted(x)))
        
        # Load abundance data and calculate edge statistics
        abundance_df = pd.read_csv(snakemake.input.abundance, sep='\t', index_col=0)
        edge_stats = calculate_edge_statistics(edge_counts.reset_index(), abundance_df)
        
        # Load taxonomy data if available
        if Path(snakemake.input.taxonomy).exists():
            taxonomy_df = pd.read_csv(snakemake.input.taxonomy, sep='\t')
            taxonomy_dict = taxonomy_df.set_index('Feature')['Taxonomy'].to_dict()
            edge_stats['taxonomy_a'] = edge_stats['source'].map(taxonomy_dict)
            edge_stats['taxonomy_b'] = edge_stats['target'].map(taxonomy_dict)
        
        # Combine all information
        final_df = pd.merge(
            edge_stats,
            edge_counts.reset_index(),
            on=['source', 'target']
        )
        
        # Calculate network-level metrics
        network_stats = calculate_network_metrics(final_df)
        
        # Save results
        final_df.to_csv(snakemake.output.network, sep='\t', index=False)
        with open(snakemake.output.stats, 'w') as f:
            json.dump(network_stats, f, indent=2)
        
        logger.info(f"Completed network aggregation. Final network has {len(final_df)} edges.")
        
    except Exception as e:
        logger.error(f"Error during network aggregation: {str(e)}")
        raise e

if __name__ == "__main__":
    main(snakemake)