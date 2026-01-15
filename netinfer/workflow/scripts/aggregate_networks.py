#!/usr/bin/env python3
"""
Aggregate networks from different methods and calculate consensus scores
"""

import pandas as pd
import numpy as np
import networkx as nx
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

def calculate_network_metrics(G: nx.Graph, trusted_methods_list: List[str] = None, other_methods_list: List[str] = None) -> Dict:
    """Calculate network-level metrics."""
    
    metrics = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'density': nx.density(G),
        'num_components': nx.number_connected_components(G),
        'largest_component_size': len(max(nx.connected_components(G), key=len)),
        'avg_degree': sum(dict(G.degree()).values()) / G.number_of_nodes(),
        'avg_clustering': nx.average_clustering(G),
        'modularity_label_propagation': nx.algorithms.community.modularity(G, nx.algorithms.community.label_propagation_communities(G)),
        'modularity_greedy': nx.algorithms.community.modularity(G, nx.algorithms.community.greedy_modularity_communities(G)),
        'modularity_louvain': nx.algorithms.community.modularity(G, nx.algorithms.community.louvain_communities(G)),
        'num_louvain_communities': len(nx.algorithms.community.louvain_communities(G)),
        'trusted_methods_used': ", ".join(trusted_methods_list) if trusted_methods_list else "None",
        'other_methods_used': ", ".join(other_methods_list) if other_methods_list else "None"
    }
    
    return metrics

def add_phylum_info(final_df: pd.DataFrame) -> pd.DataFrame:
    """Add boolean column indicating if edge connects different phyla."""
    if 'Taxonomy A' not in final_df.columns or 'Taxonomy B' not in final_df.columns:
        return final_df

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

    # Add Boolean column
    final_df['Inter phylum'] = (final_df['Phylum A'] != final_df['Phylum B']) & \
                                     (final_df['Phylum A'] != "Unknown") & \
                                     (final_df['Phylum B'] != "Unknown")
    
    # Drop the temporary phylum columns
    final_df.drop(columns=['Phylum A', 'Phylum B'], inplace=True)
    
    return final_df

def calculate_zi_pi(G: nx.Graph, node_to_community: Dict) -> Dict[str, Dict]:
    """
    Calculate Within-module connectivity (Zi) and Among-module connectivity (Pi)
    for each node based on the Zi-Pi topology framework.
    
    Roles:
    - Peripheral nodes (Specialists): Zi < 2.5 and Pi < 0.62
    - Connectors (Generalists): Zi < 2.5 and Pi > 0.62
    - Module hubs (Core nodes): Zi > 2.5 and Pi < 0.62
    - Network hubs (Super generalists): Zi > 2.5 and Pi > 0.62
    """
    node_metrics = {}
    
    # Pre-calculate degrees for all nodes
    degrees = dict(G.degree())
    
    # Group nodes by community
    community_nodes = {}
    for node, comm_id in node_to_community.items():
        if comm_id not in community_nodes:
            community_nodes[comm_id] = []
        community_nodes[comm_id].append(node)
        
    # Calculate stats for each module (community)
    module_stats = {}
    for comm_id, nodes in community_nodes.items():
        # Get degrees of nodes within this module (Kis)
        kis_values = []
        for node in nodes:
            # Degree of node within its own module
            kis = 0
            for neighbor in G.neighbors(node):
                if node_to_community.get(neighbor) == comm_id:
                    kis += 1
            kis_values.append(kis)
        
        module_stats[comm_id] = {
            'mean_kis': np.mean(kis_values) if kis_values else 0,
            'std_kis': np.std(kis_values) if kis_values else 0
        }
        
    # Calculate Zi and Pi for each node
    for node in G.nodes():
        if node not in node_to_community:
            continue
            
        comm_id = node_to_community[node]
        ki = degrees.get(node, 0) # Total degree
        
        # Calculate Kis (degree within own module)
        kis = 0
        kis_by_module = {} # Count links to each module for Pi calculation
        
        for neighbor in G.neighbors(node):
            neighbor_comm = node_to_community.get(neighbor)
            if neighbor_comm is not None:
                if neighbor_comm not in kis_by_module:
                    kis_by_module[neighbor_comm] = 0
                kis_by_module[neighbor_comm] += 1
                
                if neighbor_comm == comm_id:
                    kis += 1
        
        # Calculate Zi
        # Zi = (kis - mean_ks) / std_ks
        mean_ks = module_stats[comm_id]['mean_kis']
        std_ks = module_stats[comm_id]['std_kis']
        
        if std_ks == 0:
            zi = 0.0
        else:
            zi = (kis - mean_ks) / std_ks
            
        # Calculate Pi
        # Pi = 1 - sum((kis_m / ki)^2) for all modules m
        if ki == 0:
            pi = 0.0
        else:
            sum_sq_ratio = sum([(k_m / ki) ** 2 for k_m in kis_by_module.values()])
            pi = 1 - sum_sq_ratio
            
        # Determine Role
        if zi < 2.5:
            if pi < 0.62:
                role = "Peripheral node" # Specialist
            else:
                role = "Connector" # Generalist
        else:
            if pi < 0.62:
                role = "Module hub" # Core
            else:
                role = "Network hub" # Super generalist
                
        node_metrics[node] = {
            'Zi': round(zi, 4),
            'Pi': round(pi, 4),
            'Role': role
        }
        
    return node_metrics

def main(snakemake):
    """Main processing function."""
    logger = setup_logger(snakemake.log[0])
    logger.info("Starting network aggregation")
    
    try:
        # This is the order that methods will appear in the final table
        networks = {
            'flashweave': load_network_file(snakemake.input.flashweave, 'flashweave', logger),
            'flashweaveHE': load_network_file(snakemake.input.flashweaveHE, 'flashweaveHE', logger),
            'spieceasi': load_network_file(snakemake.input.spieceasi, 'spieceasi', logger),
            'fastspar': load_network_file(snakemake.input.fastspar, 'fastspar', logger),
            'propr': load_network_file(snakemake.input.propr, 'propr', logger),
            'spearman': load_network_file(snakemake.input.spearman, 'spearman', logger),
            'pearson': load_network_file(snakemake.input.pearson, 'pearson', logger),
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
        
        # Identify node-specific columns (those ending in _a or _b)
        node_cols = [c for c in gml_df.columns if c.endswith(('_a', '_b'))]
        
        node_attr_mapping = {} # node -> {attr_name: value}
        for _, row in gml_df.iterrows():
            src, tgt = row['source'], row['target']
            
            # Map attributes to source node (Node A)
            if src not in node_attr_mapping:
                node_attr_mapping[src] = {}
                for c in node_cols:
                    # Standard source attributes end in _a and don't describe B
                    if c.endswith('_a') and 'abundance_b' not in c:
                        node_attr_mapping[src][c] = row[c]
                    # Cross-attribute: A's abundance when B is max
                    if c == 'abundance_a_at_max_b':
                        node_attr_mapping[src][c] = row[c]
            
            # Map attributes to target node (Node B)
            if tgt not in node_attr_mapping:
                node_attr_mapping[tgt] = {}
                for c in node_cols:
                    # Standard target attributes end in _b and don't describe A
                    if c.endswith('_b') and 'abundance_a' not in c:
                        node_attr_mapping[tgt][c] = row[c]
                    # Cross-attribute: B's abundance when A is max
                    if c == 'abundance_b_at_max_a':
                        node_attr_mapping[tgt][c] = row[c]
                
        # Drop all node-related columns from edge dataframe
        gml_edge_df = gml_df.drop(columns=node_cols)
        
        # Create Graph: only remaining columns will be edge attributes
        G = nx.from_pandas_edgelist(gml_edge_df, 'source', 'target', edge_attr=True)
        
        # Add attributes to nodes in the network
        for node, attrs in node_attr_mapping.items():
            if node in G:
                for attr_name, val in attrs.items():
                    G.nodes[node][attr_name] = val

        # Calculate Louvain community detection
        louvain_communities = nx.algorithms.community.louvain_communities(G)
        
        # Create a mapping from node to community ID
        node_to_community = {}
        for community_id, community_nodes in enumerate(louvain_communities):
            for node in community_nodes:
                node_to_community[node] = community_id
        
        # Add community information as node attribute
        for node in G.nodes():
            G.nodes[node]['louvain_community'] = node_to_community.get(node, -1)
            
        # Calculate Zi-Pi metrics and add to node attributes
        zi_pi_metrics = calculate_zi_pi(G, node_to_community)
        for node, metrics in zi_pi_metrics.items():
            for key, value in metrics.items():
                G.nodes[node][key] = value
        
        # Add community information to the final_df for TSV export
        # Map community IDs for source and target nodes
        final_df['community_a'] = final_df['source'].map(node_to_community)
        final_df['community_b'] = final_df['target'].map(node_to_community)
        
        # Add a boolean column indicating if edge connects different communities
        final_df['inter_community'] = final_df['community_a'] != final_df['community_b']

        nx.write_gml(G, snakemake.output.combined_graph)
        
        # Also export as GEXF for Gephi Lite compatibility
        nx.write_gexf(G, snakemake.output.combined_gexf)


        # Calculate network-level metrics
        # Separate active methods into trusted and other
        trusted_active = [m for m, df in networks.items() 
                         if not df.empty and m in snakemake.params.trusted_methods]
        other_active = [m for m, df in networks.items() 
                       if not df.empty and m not in snakemake.params.trusted_methods]
        
        network_stats = calculate_network_metrics(G, trusted_active, other_active)

        # Rename columns to formal names for output
        column_mapping = {
            'source': 'Taxon A',
            'target': 'Taxon B',
            'flashweave': 'FlashWeave',
            'flashweaveHE': 'FlashWeaveHE',
            'fastspar': 'FastSpar',
            'spearman': 'Spearman',
            'pearson': 'Pearson',
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
            'community_a': 'Community A',
            'community_b': 'Community B',
            'inter_community': 'Inter community',
            'n_methods': 'Num supporting methods'
        }
        final_df = final_df.rename(columns=column_mapping)

        # Add between phyla info and reorder columns
        if 'Taxonomy A' in final_df.columns and 'Taxonomy B' in final_df.columns:
            final_df = add_phylum_info(final_df)
            
            # Reorder columns: "Inter phylum" between "Sample at Max B" and "Taxonomy A"
            cols = list(final_df.columns)
            if 'Inter phylum' in cols:
                cols.remove('Inter phylum')
                idx = cols.index('Taxonomy A')
                cols.insert(idx, 'Inter phylum')
                final_df = final_df[cols]

        # Save combined table
        final_df.to_csv(snakemake.output.combined_table, sep='\t', index=False)
        
        # Export Node Table (merged_nodes.tsv)
        # Collect all attributes for each node
        node_data = []
        for node, attrs in G.nodes(data=True):
            data_row = {'Node': node}
            data_row.update(attrs)
            node_data.append(data_row)
            
        nodes_df = pd.DataFrame(node_data)
        
        # Renaissance columns if needed (e.g. taxonomy_a -> Taxonomy)
        # We know from node_attr_mapping logic that attributes were mapped from source/target columns.
        # Common attributes: 'taxonomy_a'/'taxonomy_b' -> we should standardize this.
        # However, in G.nodes, we just assigned whatever was in node_attr_mapping.
        # Let's check what keys are usually there. 
        # They come from columns ending in _a or _b.
        # For a node, if it was a source, it got _a attrs. If target, _b attrs.
        # We need to clean up these keys to remove _a and _b suffixes for the node table.
        
        # Create a clean DataFrame
        cleaned_data = []
        
        # Helper to get first non-null value from possible keys
        def get_val(row_dict, keys):
            for k in keys:
                if k in row_dict and pd.notnull(row_dict[k]) and row_dict[k] != "":
                    return row_dict[k]
            return None

        # Iterate through node data and clean it up
        for item in node_data:
            clean_item = {'Node': item['Node']}
            
            # Transfer standard metrics directly
            for k in ['Zi', 'Pi', 'Role', 'louvain_community']:
                if k in item:
                    clean_item[k] = item[k]
            
            # identifying keys for taxonomy and abundance
            # They usually end with _a or _b or have no suffix if they came from node_attr_mapping directly
            # Let's look at all keys in the item
            all_keys = list(item.keys())
            
            # Find unique base names (e.g. 'taxonomy', 'mean_abundance')
            base_names = set()
            for k in all_keys:
                if k in ['Node', 'Zi', 'Pi', 'Role', 'louvain_community']:
                    continue
                if k.endswith('_a'):
                    base_names.add(k[:-2])
                elif k.endswith('_b'):
                    base_names.add(k[:-2])
                else:
                    base_names.add(k)
            
            # For each base name, try to find a value from _a, _b, or bare versions
            for base in base_names:
                # Prioritize keys that might be relevant
                possible_keys = [base, f"{base}_a", f"{base}_b"]
                
                # Special handling for cross-attributes which shouldn't be merged lightly?
                # e.g. abundance_b_at_max_a -> This is specific to a relationship, maybe not a node property?
                # Actually, in calculate_edge_statistics:
                # 'abundance_b_at_max_a': 'Abundance B at Max A'
                # These were edge attributes. When we mapped them to nodes in line 246 (original code):
                # if c == 'abundance_a_at_max_b': node_attr_mapping[src][c] = row[c]
                # So 'Abundance A at Max B' was assigned to Node A. 
                # This attribute effectively describes Node A.
                # So it is a valid node attribute. 
                # Its name in the DF was 'abundance_a_at_max_b'.
                # Base name would be 'abundance_a_at_max_b' (no _a/_b suffix).
                # Wait, my previous rename logic strips _b ! 
                # abundance_a_at_max_b -> ends with _b -> becomes abundance_a_at_max.
                # This might collide if there is 'abundance_a_at_max'. 
                # Let's check calculate_edge_statistics names: 
                # 'abundance_b_at_max_a'
                # 'abundance_a_at_max_b'
                # neither of these end in just _a or _b in a way that implies simple node attribute duplication.
                # BUT, wait.
                # 'mean_abundance_a' -> ends in _a -> base 'mean_abundance'.
                # 'mean_abundance_b' -> ends in _b -> base 'mean_abundance'.
                # This is correct.
                
                # However, 'abundance_a_at_max_b' -> ends in _b -> base 'abundance_a_at_max' ??
                # This logic is dangerous for keys that naturally end in a/b but aren't just copies.
                # But 'abundance_a_at_max_b' IS the value for Node A.
                # 'abundance_b_at_max_a' IS the value for Node B.
                # If we have Node A, we want 'abundance_a_at_max_b'. 
                # If we convert it to 'abundance_a_at_max', is that meaningful? 
                # Maybe "Abundance at Max Partner"?
                
                # The user output shows columns: "abundance_a_at_max", "abundance_b_at_max".
                # And duplication of "prevalence", "mean_abundance".
                
                # Let's stick to merging "mean_abundance_a" and "mean_abundance_b" into "mean_abundance".
                # And keeping unique complex metrics if possible, or simplifying them.
                
                val = get_val(item, possible_keys)
                if val is not None:
                    clean_item[base] = val
                    
            cleaned_data.append(clean_item)
            
        nodes_df = pd.DataFrame(cleaned_data)
        
        # Rename specific columns to match user preference and formal names
        column_mapping_nodes = {
            'louvain_community': 'Louvain Community',
            'Zi': 'Within-module connectivity (Zi)',
            'Pi': 'Among-module connectivity (Pi)',
            'Role': 'Zi-Pi Role',
            'taxonomy': 'Taxonomy',
            'prevalence': 'Prevalence (%)',
            'mean_abundance': 'Mean Abundance',
            'max_abundance': 'Max Abundance',
            'sample_at_max': 'Sample at Max Abundance'
        }
        
        nodes_df = nodes_df.rename(columns=column_mapping_nodes)
        
        # Remove columns that might be confusing or redundant if they were stripped incorrectly
        # e.g. 'abundance_a_at_max' which came from 'abundance_a_at_max_b'
        # Let's keep them but maybe rename if they are clear?
        # Actually 'abundance_a_at_max_b' means "Abundance of A when B is max".
        # This is an edge property (dependent on B). Mapping it to Node A is only valid in the context of specific edges.
        # But here we are aggregating to a Node table.
        # If Node A has many edges, which 'abundance_a_at_max_b' do we keep?
        # The code at line 246 said:
        # if src not in node_attr_mapping: ...
        # So it takes the value from the FIRST edge encountered for that node.
        # This is somewhat arbitrary for edge-dependent properties.
        # Standard node properties (abundance, taxonomy) should be constant across edges.
        # So merging _a and _b for abundance/taxonomy is correct.
        # But 'abundance_a_at_max_b' varies by edge. Exporting it for a node is misleading.
        # I should probably DROP edge-specific attributes from the node table to avoid confusion.
        
        cols_to_drop = [c for c in nodes_df.columns if 'at_max' in c]
        nodes_df = nodes_df.drop(columns=cols_to_drop, errors='ignore')
        
        # Reorder columns
        preferred_start = ['Node', 'Taxonomy', 'Zi-Pi Role', 'Within-module connectivity (Zi)', 'Among-module connectivity (Pi)', 'Louvain Community']
        preferred_start = [c for c in preferred_start if c in nodes_df.columns]
        
        other_cols = [c for c in nodes_df.columns if c not in preferred_start]
        nodes_df = nodes_df[preferred_start + other_cols]
        
        # Save node table
        nodes_df.to_csv(snakemake.output.nodes_table, sep='\t', index=False)

        # Save network statistics as tab-separated text file with descriptions
        metric_descriptions = {
            'nodes': 'Total number of nodes (taxa) in the network',
            'edges': 'Total number of edges (associations) in the network',
            'density': 'Network density, ratio of actual edges to possible edges (0-1)',
            'num_components': 'Number of connected components in the network',
            'largest_component_size': 'Number of nodes in the largest connected component',
            'avg_degree': 'Average degree (number of connections) per node',
            'avg_clustering': 'Average clustering coefficient, measure of local connectivity (0-1)',
            'modularity_label_propagation': 'Modularity score using label propagation algorithm (-1 to 1)',
            'modularity_greedy': 'Modularity score using greedy optimization algorithm (-1 to 1)',
            'modularity_louvain': 'Modularity score using Louvain algorithm (-1 to 1)',
            'num_louvain_communities': 'Number of communities detected by Louvain algorithm',
            'trusted_methods_used': 'List of trusted methods that contributed to the network',
            'other_methods_used': 'List of other methods that were run and included in the output'
        }
        
        with open(snakemake.output.stats, 'w') as f:
            # Write header
            f.write('Item\tValue\tDescription\n')
            # Write each metric
            for key, value in network_stats.items():
                description = metric_descriptions.get(key, 'No description available')
                f.write(f'{key}\t{value}\t{description}\n')


        logger.info(f"Completed network aggregation. Final network has {len(final_df)} edges.")
        
    except Exception as e:
        logger.exception(f"Error during network aggregation: {str(e)}")
        raise

if __name__ == "__main__":
    main(snakemake)