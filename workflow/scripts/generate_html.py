#!/usr/bin/env python3
"""
Generate interactive HTML network visualization
"""

import pandas as pd
import numpy as np
import json
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import networkx as nx
from pathlib import Path
import logging
from typing import Dict, List, Tuple
import jinja2

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

def prepare_network_data(network_df: pd.DataFrame,
                        taxonomy_df: pd.DataFrame = None,
                        max_edges: int = 10000) -> Tuple[Dict, nx.Graph]:
    """Prepare network data for visualization."""
    # Create networkx graph
    # Identify available edge attributes
    potential_methods = ['FlashWeave', 'FlashWeaveHE', 'FastSpar', 'Spearman', 'SpiecEasi', 'PropR', 'Jaccard']
    available_methods = [col for col in potential_methods if col in network_df.columns]
    
    edge_attrs = ['Num supporting methods'] + available_methods
    # Also include 'methods' and 'methods_count' if they exist (legacy support)
    if 'methods' in network_df.columns:
        edge_attrs.append('methods')
    if 'methods_count' in network_df.columns:
        edge_attrs.append('methods_count')
    
    G = nx.from_pandas_edgelist(
        network_df,
        'Taxon A',
        'Taxon B',
        edge_attr=edge_attrs
    )
    
    # Calculate node metrics
    node_metrics = {
        'degree': dict(G.degree()),
        'betweenness': nx.betweenness_centrality(G),
        'closeness': nx.closeness_centrality(G)
    }
    
    # Add taxonomy information if available
    if taxonomy_df is not None:
        taxonomy_dict = taxonomy_df.set_index('Feature')['Taxonomy'].to_dict()
        nx.set_node_attributes(G, taxonomy_dict, 'taxonomy')
    
    # Prepare JSON-serializable data
    nodes_data = []
    for node in G.nodes():
        node_data = {
            'id': node,
            'degree': node_metrics['degree'][node],
            'betweenness': node_metrics['betweenness'][node],
            'closeness': node_metrics['closeness'][node]
        }
        if taxonomy_df is not None and node in taxonomy_dict:
            node_data['taxonomy'] = taxonomy_dict[node]
        nodes_data.append(node_data)
    
    edges_data = []
    for edge in G.edges(data=True):
        data = edge[2]
        
        # Get method count
        count = data.get('methods_count', data.get('Num supporting methods', 0))
        
        # Get list of methods
        if 'methods' in data:
            methods_list = data['methods']
            if isinstance(methods_list, str):
                methods_list = methods_list.split(';')
        else:
            # Reconstruct from available method columns
            methods_list = [m for m in available_methods if pd.notnull(data.get(m))]
            
        edges_data.append({
            'source': edge[0],
            'target': edge[1],
            'methods_count': int(count),
            'methods': methods_list
        })
    
    # Limit edges if needed
    if len(edges_data) > max_edges:
        edges_data.sort(key=lambda x: x['methods_count'], reverse=True)
        edges_data = edges_data[:max_edges]
    
    return {'nodes': nodes_data, 'edges': edges_data}, G

def create_abundance_plot(source: str,
                        target: str,
                        abundance_df: pd.DataFrame) -> go.Figure:
    """Create abundance correlation plot for a pair of taxa."""
    fig = make_subplots(rows=2, cols=2,
                       subplot_titles=['Correlation Plot',
                                     'Abundance Distribution',
                                     'Time Series'])
    
    # Correlation plot
    fig.add_trace(
        go.Scatter(x=abundance_df[source],
                  y=abundance_df[target],
                  mode='markers',
                  name='Samples'),
        row=1, col=1
    )
    
    # Abundance distributions
    fig.add_trace(
        go.Box(y=abundance_df[source], name=source),
        row=1, col=2
    )
    fig.add_trace(
        go.Box(y=abundance_df[target], name=target),
        row=1, col=2
    )
    
    # Time series
    fig.add_trace(
        go.Scatter(y=abundance_df[source],
                  name=source,
                  mode='lines'),
        row=2, col=1
    )
    fig.add_trace(
        go.Scatter(y=abundance_df[target],
                  name=target,
                  mode='lines'),
        row=2, col=1
    )
    
    return fig

def generate_html(network_data: Dict,
                 stats: Dict,
                 abundance_df: pd.DataFrame) -> str:
    """Generate HTML visualization."""
    template = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>NetInfer Network Visualization</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/vis-network/9.1.2/vis-network.min.js"></script>
        <style>
            body { margin: 0; padding: 20px; font-family: Arial, sans-serif; }
            #network { height: 600px; border: 1px solid lightgray; }
            #details { margin-top: 20px; }
            .stats { margin-bottom: 20px; }
            .plot { height: 400px; }
        </style>
    </head>
    <body>
        <h1>Network Visualization</h1>
        <div class="stats">
            <h2>Network Statistics</h2>
            <pre>{{ stats | tojson(indent=2) }}</pre>
        </div>
        <div id="network"></div>
        <div id="details">
            <h2>Edge Details</h2>
            <div id="abundance-plot" class="plot"></div>
        </div>
        <script>
            const network_data = {{ network_data | tojson }};
            // Network visualization code here
            // (Will be implemented in the full version)
        </script>
    </body>
    </html>
    """
    
    return jinja2.Template(template).render(
        network_data=network_data,
        stats=stats
    )

def main(snakemake):
    """Main processing function."""
    logger = setup_logger(snakemake.log[0])
    logger.info("Starting HTML visualization generation")
    
    try:
        # Load network data
        network_df = pd.read_csv(snakemake.input.network, sep='\t')
        
        # Load taxonomy data if available
        taxonomy_df = None
        taxonomy_input = snakemake.input.taxonomy
        # Handle Namedlist or list input from snakemake
        if isinstance(taxonomy_input, (list, tuple)) or 'Namedlist' in str(type(taxonomy_input)):
            taxonomy_path = str(taxonomy_input[0]) if len(taxonomy_input) > 0 else ""
        else:
            taxonomy_path = str(taxonomy_input)

        if taxonomy_path and Path(taxonomy_path).exists():
            taxonomy_df = pd.read_csv(taxonomy_path, sep='\t')
        
        # Prepare network data
        network_data, G = prepare_network_data(
            network_df,
            taxonomy_df,
            max_edges=snakemake.params.max_edges
        )
        
        # Load network statistics
        with open(snakemake.input.stats) as f:
            network_stats = json.load(f)
        
        # Load abundance data
        abundance_df = pd.read_csv(snakemake.input.abundance, sep='\t', index_col=0)
        
        # Generate HTML
        html_content = generate_html(network_data, network_stats, abundance_df)
        
        # Save results
        with open(snakemake.output.html, 'w') as f:
            f.write(html_content)
        
        with open(snakemake.output.data, 'w') as f:
            json.dump(network_data, f, indent=2)
        
        logger.info("Completed HTML visualization generation")
        
    except Exception as e:
        logger.error(f"Error during HTML generation: {str(e)}")
        raise e

if __name__ == "__main__":
    main(snakemake)