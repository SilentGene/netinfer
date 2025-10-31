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
from functools import reduce

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

def load_network_file(file_path: str, method: str, logger: logging.Logger) -> pd.DataFrame:
    """Load network file and standardize format."""
    if not Path(file_path).exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(file_path, sep='\t')
    except pd.errors.EmptyDataError:
        df = df = pd.DataFrame(columns=['source', 'target', method])
    # only keep first three columns
    df = df.iloc[:, :3]
    # Standardize column names with 'source', 'target', method
    df.columns = ['source', 'target', method]
    # Ensure source < target for consistency
    df[['source', 'target']] = np.sort(df[['source', 'target']], axis=1)
    # print log 
    logger.info(f"Loaded {len(df)} edges from {file_path}")
    return df

def combine_methods(networks: Dict[str, pd.DataFrame], trusted_methods: List[str]) -> pd.DataFrame:
    """Combine networks from different methods into a single DataFrame."""
    # First, combine trusted methods
    trusted_dfs = [df for method, df in networks.items() if method in trusted_methods and not df.empty]
    if not trusted_dfs:
        raise ValueError("No trusted methods available for aggregation. Please check the config.yaml file.")
    combined_trusted = reduce(
        lambda left, right: pd.merge(left, right, on=['source', 'target'], how='outer'),
        trusted_dfs
    )
    
    # Then, append non-trusted methods
    non_trusted_dfs = [df for method, df in networks.items() if method not in trusted_methods and not df.empty]
    if non_trusted_dfs:
        combined_df = pd.concat([combined_trusted] + non_trusted_dfs, ignore_index=True)
    else:
        combined_df = combined_trusted

    return combined_df


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

def calculate_network_metrics(G: nx.Graph) -> Dict:
    """Calculate network-level metrics."""
    
    metrics = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'density': nx.density(G),
        'components': nx.number_connected_components(G),
        'largest_component_size': len(max(nx.connected_components(G), key=len)),
        'avg_degree': sum(dict(G.degree()).values()) / G.number_of_nodes(),
        'avg_clustering': nx.average_clustering(G),
        'methods_used': len(set(nx.get_edge_attributes(G, 'method').values()))
    }
    
    return metrics

def main(snakemake):
    """Main processing function."""
    logger = setup_logger(snakemake.log[0])
    logger.info("Starting network aggregation")
    
    try:
        # Load networks from each method
        networks = {
            'FlashWeave': load_network_file(snakemake.input.flashweave, 'FlashWeave', logger),
            'FlashWeaveHE': load_network_file(snakemake.input.flashweaveHE, 'FlashWeaveHE', logger),
            'FastSpar': load_network_file(snakemake.input.fastspar, 'FastSpar', logger),
            'Spearman': load_network_file(snakemake.input.spearman, 'Spearman', logger),
            'SPIEC-EASI': load_network_file(snakemake.input.spieceasi, 'SPIEC-EASI', logger),
            'propR': load_network_file(snakemake.input.propr, 'propR', logger),
            'Jaccard': load_network_file(snakemake.input.jaccard, 'Jaccard', logger)
        }
        
        combined_df = combine_methods(networks, snakemake.params.trusted_methods)

        logger.info(f"Combined {len(combined_df)} edges from trusted methods: {snakemake.params.trusted_methods} and appended weights from other methods")
        # export combined_df
        combined_df.to_csv("../test/testout/test.out.tsv", sep='\t', index=False)
        
        # add a column "n_methods" indicating how many methods detected each edge. That is the count of non-null values across method columns.
        method_cols = [col for col in combined_df.columns if col not in ['source', 'target']]
        combined_df['n_methods'] = combined_df[method_cols].notnull().sum(axis=1)

        # Load abundance data and calculate edge statistics
        abundance_df = pd.read_csv(snakemake.input.abundance, sep='\t', index_col=0)
        edge_stats = calculate_edge_statistics(combined_df.reset_index(), abundance_df)
        
        # Load taxonomy data if available
        if Path(snakemake.input.taxonomy).exists():
            taxonomy_df = pd.read_csv(snakemake.input.taxonomy, sep='\t')
            taxonomy_dict = taxonomy_df.set_index('Feature')['Taxonomy'].to_dict()
            edge_stats['taxonomy_a'] = edge_stats['source'].map(taxonomy_dict)
            edge_stats['taxonomy_b'] = edge_stats['target'].map(taxonomy_dict)
        
        # Combine all information
        final_df = pd.merge(
            combined_df,
            edge_stats,
            on=['source', 'target']
        )

        # order final_df by n_methods descending and they by max weight across trusted methods
        final_df = final_df.sort_values(by=['n_methods'] + snakemake.params.trusted_methods, ascending=False)
        
        # Save combined table
        final_df.to_csv(snakemake.output.combined_table, sep='\t', index=False)
        
        # Save gml graph
        G = nx.from_pandas_edgelist(final_df, 'source', 'target', edge_attr=True)
        nx.write_gml(G, snakemake.output.combined_graph)

        # Calculate network-level metrics
        network_stats = calculate_network_metrics(G)

        # Save network statistics
        with open(snakemake.output.stats, 'w') as f:
            json.dump(network_stats, f, indent=2)


        logger.info(f"Completed network aggregation. Final network has {len(final_df)} edges.")
        
    except Exception as e:
        logger.exception(f"Error during network aggregation: {str(e)}")
        raise

if __name__ == "__main__":
    main(snakemake)