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
    # Preserve the intended method order from the input dict
    method_order = list(networks.keys())
    # First, combine trusted methods
    trusted_dfs = [df for method, df in networks.items() if method in trusted_methods and not df.empty]
    if not trusted_dfs:
        raise ValueError("No trusted methods available for aggregation. Please check the config.yaml file.")
    
    # Outer join for trusted methods to keep union of all trusted edges
    combined_df = reduce(
        lambda left, right: pd.merge(left, right, on=['source', 'target'], how='outer'),
        trusted_dfs
    )
    
    # Then, left join non-trusted methods
    # This ensures we only keep edges present in trusted methods, but annotate them with other methods' info
    non_trusted_dfs = [df for method, df in networks.items() if method not in trusted_methods and not df.empty]
    for df in non_trusted_dfs:
        combined_df = pd.merge(combined_df, df, on=['source', 'target'], how='left')

    # Ensure columns exist for all methods whose files were present, even if empty
    # If a method produced no edges, keep its column with NaNs to indicate no associations
    for method in networks.keys():
        if method not in combined_df.columns:
            combined_df[method] = np.nan

    # Reorder columns to keep a stable order: source, target, then methods in method_order
    ordered_cols = ['source', 'target'] + method_order
    combined_df = combined_df[[col for col in ordered_cols if col in combined_df.columns]]

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
            'prevalence_a': (abund_a > 0).mean() * 100,
            'prevalence_b': (abund_b > 0).mean() * 100,
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

def get_diff_phyla_edges(final_df: pd.DataFrame) -> pd.DataFrame:
    """Get edges connecting different phyla based on taxonomy columns."""
    if 'Taxonomy A' not in final_df.columns or 'Taxonomy B' not in final_df.columns:
        return pd.DataFrame()  # Return empty DataFrame if taxonomy columns are missing

    def extract_phylum(taxonomy: str) -> str:
        """ Extract the string after p__ in the taxonomy string. """
        if pd.isna(taxonomy) or taxonomy == "":
            return "Unknown"
        for part in taxonomy.split(';'):
            part = part.strip()
            if part.startswith('p__'):
                return part[3:].strip()
        return "Unknown"

    final_df['Phylum A'] = final_df['Taxonomy A'].apply(extract_phylum)
    final_df['Phylum B'] = final_df['Taxonomy B'].apply(extract_phylum)

    diff_phyla_df = final_df[final_df['Phylum A'] != final_df['Phylum B']].copy()
    
    # Drop the temporary phylum columns
    diff_phyla_df.drop(columns=['Phylum A', 'Phylum B'], inplace=True)
    
    return diff_phyla_df

def main(snakemake):
    """Main processing function."""
    logger = setup_logger(snakemake.log[0])
    logger.info("Starting network aggregation")
    
    try:
        # Load networks from each method
        networks = {
            'flashweave': load_network_file(snakemake.input.flashweave, 'flashweave', logger),
            'flashweaveHE': load_network_file(snakemake.input.flashweaveHE, 'flashweaveHE', logger),
            'fastspar': load_network_file(snakemake.input.fastspar, 'fastspar', logger),
            'spearman': load_network_file(snakemake.input.spearman, 'spearman', logger),
            'spieceasi': load_network_file(snakemake.input.spieceasi, 'spieceasi', logger),
            'propr': load_network_file(snakemake.input.propr, 'propr', logger),
            'jaccard': load_network_file(snakemake.input.jaccard, 'jaccard', logger)
        }
        
        combined_df = combine_methods(networks, snakemake.params.trusted_methods)

        # Ensure all trusted methods are present as columns (fill with NaN if missing)
        # This prevents KeyError during sorting if a trusted method produced no edges
        for method in snakemake.params.trusted_methods:
            if method not in combined_df.columns:
                combined_df[method] = np.nan

        logger.info(f"Combined {len(combined_df)} edges from trusted methods: {snakemake.params.trusted_methods} and appended weights from other methods")
        # export combined_df
        #combined_df.to_csv("../test/testout/test.out.tsv", sep='\t', index=False)

        
        # add a column "n_methods" indicating how many methods detected each edge. That is the count of non-null values across method columns.
        method_cols = [col for col in combined_df.columns if col not in ['source', 'target']]
        combined_df['n_methods'] = combined_df[method_cols].notnull().sum(axis=1)

        # Load abundance data and calculate edge statistics
        abundance_df = pd.read_csv(snakemake.input.abundance, sep='\t', index_col=0)
        abundance_df = abundance_df.T
        edge_stats = calculate_edge_statistics(combined_df.reset_index(), abundance_df)
        
        # Load taxonomy data if available
        taxonomy_input = snakemake.input.taxonomy
        # Handle Namedlist or list input from snakemake
        if isinstance(taxonomy_input, (list, tuple)) or 'Namedlist' in str(type(taxonomy_input)):
            taxonomy_path = str(taxonomy_input[0]) if len(taxonomy_input) > 0 else ""
        else:
            taxonomy_path = str(taxonomy_input)

        if taxonomy_path and Path(taxonomy_path).exists():
            taxonomy_df = pd.read_csv(taxonomy_path, sep='\t')
            # Use first column as Feature and second as Taxonomy, ignoring column names
            taxonomy_dict = dict(zip(taxonomy_df.iloc[:, 0], taxonomy_df.iloc[:, 1]))
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
        
        # Save gml graph (Save before renaming to avoid spaces in GML keys)
        # Create a copy for GML export to handle NaNs without affecting the CSV output
        gml_df = final_df.copy()
        
        # Fill NaNs in numeric columns with 0.0 for GML compatibility (Gephi doesn't like NaNs)
        numeric_cols = gml_df.select_dtypes(include=[np.number]).columns
        gml_df[numeric_cols] = gml_df[numeric_cols].fillna(0.0)
        
        # Fill NaNs in other columns with empty string
        gml_df = gml_df.fillna("")
        
        G = nx.from_pandas_edgelist(gml_df, 'source', 'target', edge_attr=True)
        nx.write_gml(G, snakemake.output.combined_graph)

        # Calculate network-level metrics
        network_stats = calculate_network_metrics(G)

        # Rename columns to formal names for output
        column_mapping = {
            'source': 'Taxon A',
            'target': 'Taxon B',
            'flashweave': 'FlashWeave',
            'flashweaveHE': 'FlashWeaveHE',
            'fastspar': 'FastSpar',
            'spearman': 'Spearman',
            'spieceasi': 'SpiecEasi',
            'propr': 'PropR',
            'jaccard': 'Jaccard',
            'prevalence_a': 'Prevalence A (%)',
            'prevalence_b': 'Prevalence B (%)',
            'mean_abundance_a': 'Mean Abundance A',
            'mean_abundance_b': 'Mean Abundance B',
            'max_abundance_a': 'Max Abundance A',
            'max_abundance_b': 'Max Abundance B',
            'abundance_b_at_max_a': 'Abundance B at Max A',
            'abundance_a_at_max_b': 'Abundance A at Max B',
            'sample_at_max_a': 'Sample at Max A',
            'sample_at_max_b': 'Sample at Max B',
            'taxonomy_a': 'Taxonomy A',
            'taxonomy_b': 'Taxonomy B',
            'n_methods': 'Num supporting methods'
        }
        final_df = final_df.rename(columns=column_mapping)

        # Save combined table
        final_df.to_csv(snakemake.output.combined_table, sep='\t', index=False)

        # Association between different phyla
        if snakemake.output.diff_phyla_table:
            diff_phyla_df = get_diff_phyla_edges(final_df)
            diff_phyla_df.to_csv(snakemake.output.diff_phyla_table, sep='\t', index=False)
        
        # Save network statistics
        with open(snakemake.output.stats, 'w') as f:
            json.dump(network_stats, f, indent=2)


        logger.info(f"Completed network aggregation. Final network has {len(final_df)} edges.")
        
    except Exception as e:
        logger.exception(f"Error during network aggregation: {str(e)}")
        raise

if __name__ == "__main__":
    main(snakemake)