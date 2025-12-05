#!/usr/bin/env python3
"""
Generate interactive HTML network visualization using Pyvis
"""

import pandas as pd
import json
import networkx as nx
from pathlib import Path
import logging
import os
from typing import Dict, Tuple
from pyvis.network import Network

def setup_logger(log_file: str) -> logging.Logger:
    """Set up logging to both file and console."""
    logger = logging.getLogger('generate_html')
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

def get_color_from_string(s: str) -> str:
    """Generate a consistent color from a string."""
    if not s or pd.isna(s):
        return '#97C2FC' # Default pyvis blue
    
    idx = abs(hash(s)) % len(COLORS)
    return COLORS[idx]

def prepare_network_data(network_df: pd.DataFrame,
                        max_edges: int = 10000) -> Tuple[Dict, nx.Graph]:
    """Prepare network data for visualization."""
    
    # Identify available edge attributes
    potential_methods = ['FlashWeave', 'FlashWeaveHE', 'FastSpar', 'Spearman', 'SpiecEasi', 'PropR', 'Jaccard']
    available_methods = [col for col in potential_methods if col in network_df.columns]
    
    # Handle column name variations for method count
    if 'Num supporting methods' in network_df.columns:
        network_df['methods_count'] = network_df['Num supporting methods']
    elif 'n_methods' in network_df.columns:
        network_df['methods_count'] = network_df['n_methods']
    elif 'methods_count' not in network_df.columns:
        # Fallback: calculate it if missing
        network_df['methods_count'] = network_df[available_methods].notnull().sum(axis=1)
    
    edge_attrs = ['methods_count'] + available_methods
    if 'methods' in network_df.columns:
        edge_attrs.append('methods')

    # Create graph
    G = nx.from_pandas_edgelist(
        network_df,
        'Taxon A',
        'Taxon B',
        edge_attr=edge_attrs
    )
    
    # Calculate node metrics
    degrees = dict(G.degree())
    betweenness = nx.betweenness_centrality(G)
    closeness = nx.closeness_centrality(G)
    
    nx.set_node_attributes(G, degrees, 'degree')
    nx.set_node_attributes(G, betweenness, 'betweenness')
    nx.set_node_attributes(G, closeness, 'closeness')
    
    # Reconstruct methods list if missing in edge attributes
    for u, v, data in G.edges(data=True):
        if 'methods' not in data:
            found_methods = [m for m in available_methods if pd.notnull(data.get(m))]
            data['methods'] = "; ".join(found_methods)

    # Prepare JSON data (legacy requirement)
    nodes_data = []
    for node in G.nodes():
        node_data = {
            'id': node,
            'degree': degrees[node],
            'betweenness': betweenness[node],
            'closeness': closeness[node]
        }
        nodes_data.append(node_data)
        
    edges_data = []
    for u, v, data in G.edges(data=True):
        edge_data = {
            'source': u,
            'target': v,
            'methods_count': data.get('methods_count', 0),
            'methods': data.get('methods', '')
        }
        edges_data.append(edge_data)

    # Filter edges if too many
    if G.number_of_edges() > max_edges:
        sorted_edges = sorted(G.edges(data=True), key=lambda x: x[2].get('methods_count', 0), reverse=True)
        edges_to_keep = sorted_edges[:max_edges]
        
        # Create new graph with filtered edges
        G_filtered = nx.Graph()
        G_filtered.add_edges_from(edges_to_keep)
        
        # Copy node attributes
        for node in G_filtered.nodes():
            if node in G.nodes:
                G_filtered.nodes[node].update(G.nodes[node])
        
        G = G_filtered
        
        # Update JSON data to match filtered graph
        nodes_set = set(G.nodes())
        nodes_data = [n for n in nodes_data if n['id'] in nodes_set]
        edges_data = [e for e in edges_data if G.has_edge(e['source'], e['target'])]

    return {'nodes': nodes_data, 'edges': edges_data}, G

def style_network(G: nx.Graph, params: Dict):
    """Apply Pyvis styles to the graph."""
    
    node_size_by = params.get('node_size_by', 'degree')
    edge_width_by = params.get('edge_width_by', 'weight')
    
    # Node styling
    for node in G.nodes():
        attrs = G.nodes[node]
        
        # Size
        base_size = 10
        if node_size_by in attrs:
            val = attrs[node_size_by]
            # Normalize or scale? For now just simple scaling
            attrs['size'] = base_size + (val * 2) 
        else:
            attrs['size'] = base_size
            
        # Color and Group
        # Default color since taxonomy is removed
        attrs['color'] = '#97C2FC'
        attrs['title'] = f"ID: {node}<br>Degree: {attrs.get('degree',0)}"
            
        attrs['label'] = str(node)

    # Edge styling
    for u, v, data in G.edges(data=True):
        # Width
        base_width = 1
        if edge_width_by == 'weight':
            # Use methods_count as proxy
            val = data.get('methods_count', 1)
            data['width'] = val
        elif edge_width_by in data:
            data['width'] = data[edge_width_by]
        else:
            data['width'] = base_width
            
        # Tooltip
        methods = data.get('methods', '')
        count = data.get('methods_count', 0)
        data['title'] = f"Methods ({count}): {methods}"

def main(snakemake):
    logger = setup_logger(snakemake.log[0])
    logger.info("Starting Pyvis visualization generation")
    
    try:
        # Load inputs
        network_df = pd.read_csv(snakemake.input.network, sep='\t')
        
        # Prepare data
        network_data, G = prepare_network_data(
            network_df, 
            max_edges=snakemake.params.max_edges
        )
        
        # Style graph for Pyvis
        style_network(G, snakemake.params)
        
        # Generate Pyvis network
        # use_DOT=False to avoid graphviz dependency if possible, though pyvis doesn't use it by default
        net = Network(height="800px", width="100%", bgcolor="#ffffff", font_color="black", select_menu=True, filter_menu=True)
        net.from_nx(G)
        
        # Add physics controls
        net.show_buttons(filter_=['physics'])
        
        # Use Barnes-Hut and disable stabilization via valid JSON options (prevents 0% progress hang)
        net.set_options("""
        {
            "physics": {
                "enabled": true,
                "solver": "barnesHut",
                "stabilization": { "enabled": false },
                "barnesHut": {
                    "gravitationalConstant": -2000,
                    "centralGravity": 0.3,
                    "springLength": 95,
                    "springConstant": 0.04,
                    "damping": 0.09,
                    "avoidOverlap": 0
                }
            }
        }
        """)
        
        # Save HTML
        # Change directory to output folder so lib/ is generated there
        output_path = Path(snakemake.output.html)
        output_dir = output_path.parent
        output_file = output_path.name
        
        current_dir = os.getcwd()
        try:
            if not output_dir.exists():
                output_dir.mkdir(parents=True, exist_ok=True)
            os.chdir(output_dir)
            net.save_graph(output_file)
        finally:
            os.chdir(current_dir)
        
        # Save JSON data
        with open(snakemake.output.data, 'w') as f:
            json.dump(network_data, f, indent=2)
            
        logger.info(f"Saved visualization to {snakemake.output.html}")
        
    except Exception as e:
        logger.error(f"Error during visualization: {str(e)}")
        raise e

if __name__ == "__main__":
    main(snakemake)