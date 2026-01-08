#!/usr/bin/env python3
"""
Generate interactive HTML network visualization dashboard
"""

import pandas as pd
import json
import logging
from pathlib import Path
import numpy as np
from scipy import stats

def setup_logger(log_file: str) -> logging.Logger:
    """Set up logging."""
    logger = logging.getLogger('generate_html')
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(log_file)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)
    return logger

def get_method_columns(df: pd.DataFrame) -> list:
    """Identify method columns in the network file."""
    potential_methods = ['FlashWeave', 'FlashWeaveHE', 'FastSpar', 'Spearman', 'SpiecEasi', 'PropR', 'Jaccard']
    return [col for col in potential_methods if col in df.columns]

def prepare_data(network_df: pd.DataFrame, abundance_df: pd.DataFrame, logger: logging.Logger) -> dict:
    """Prepare data structures for the dashboard."""
    
    # 1. Process Network Data
    method_cols = get_method_columns(network_df)
    logger.info(f"Found method columns: {method_cols}")
    
    # Ensure required columns exist
    required_cols = ['Taxon A', 'Taxon B']
    if not all(col in network_df.columns for col in required_cols):
        raise ValueError(f"Network file missing required columns: {required_cols}")

    # Create a simplified list of edges for the table
    edges_data = []
    
    # Collect all unique taxons involved in the network
    involved_taxa = set(network_df['Taxon A'].unique()) | set(network_df['Taxon B'].unique())
    logger.info(f"Unique taxa in network: {len(involved_taxa)}")

    for idx, row in network_df.iterrows():
        edge = {
            'id': idx,
            'TaxonA': row['Taxon A'],
            'TaxonB': row['Taxon B'],
            'TaxonomyA': row.get('Taxonomy A', 'N/A'),
            'TaxonomyB': row.get('Taxonomy B', 'N/A'),
            'methods': {m: row[m] for m in method_cols if pd.notna(row[m])}
        }
        edges_data.append(edge)

    # 2. Process Abundance Data
    # Filter abundance matrix to include only taxa present in the network
    # Assumes abundance_df index is Taxon ID
    
    # Check if index name matches or if first column is the ID
    if abundance_df.index.name != 'Taxon A' and abundance_df.index.name != 'featureid':
        # heuristic to find the taxon column if it's not the index
        if 'featureid' in abundance_df.columns:
            abundance_df.set_index('featureid', inplace=True)
        elif '#OTU ID' in abundance_df.columns:
             abundance_df.set_index('#OTU ID', inplace=True)

    # Intersection of taxa
    valid_taxa = [t for t in involved_taxa if t in abundance_df.index]
    missing_taxa = involved_taxa - set(valid_taxa)
    
    if missing_taxa:
        logger.warning(f"{len(missing_taxa)} taxa from network not found in abundance file. Example: {list(missing_taxa)[:3]}")

    filtered_abundance = abundance_df.loc[valid_taxa].copy()
    
    # Convert to dictionary conforming to: { 'TaxonID': [v1, v2, v3...] }
    # And keep sample names
    samples = list(filtered_abundance.columns)
    abundance_map = filtered_abundance.to_dict(orient='index')
    
    # Clean up numpy types for JSON serialization
    for taxon, data in abundance_map.items():
        # data is {sample: value} because orient='index'
        # we want a sorted list of values corresponding to 'samples' list
        # Actually to preserve order, let's just grab the values in order of 'samples'
        # But wait, to_dict(orient='index') gives {col: val}, so we can just look them up
        values = [data.get(s, 0) for s in samples]
        # Replace NaN with 0 or null
        abundance_map[taxon] = [float(v) if pd.notnull(v) else 0 for v in values]
        
    return {
        'edges': edges_data,
        'samples': samples,
        'abundance': abundance_map
    }

def generate_html(data_json: str) -> str:
    """Generate the full HTML content."""
    
    html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>NetInfer Result</title>
    
    <!-- Dependencies -->
    <link href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.3.0/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdn.datatables.net/1.13.4/css/dataTables.bootstrap5.min.css" rel="stylesheet">
    
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.4/jquery.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.3.0/js/bootstrap.bundle.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.4/js/dataTables.bootstrap5.min.js"></script>
    <script src="https://cdn.plot.ly/plotly-2.20.0.min.js"></script>
    <!-- Simple Linear Regression Library (Inline for simplicity) -->
    <script>
    function linearRegression(x, y) {
        var n = y.length;
        var sum_x = 0;
        var sum_y = 0;
        var sum_xy = 0;
        var sum_xx = 0;
        var sum_yy = 0;

        for (var i = 0; i < n; i++) {
            sum_x += x[i];
            sum_y += y[i];
            sum_xy += (x[i] * y[i]);
            sum_xx += (x[i] * x[i]);
            sum_yy += (y[i] * y[i]);
        }

        var slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
        var intercept = (sum_y - slope * sum_x) / n;
        var r2 = Math.pow((n * sum_xy - sum_x * sum_y) / Math.sqrt((n * sum_xx - sum_x * sum_x) * (n * sum_yy - sum_y * sum_y)), 2);

        return { slope: slope, intercept: intercept, r2: r2 };
    }
    </script>

    <style>
        body { height: 100vh; overflow: hidden; display: flex; flex-direction: column; }
        .header { background-color: #f8f9fa; padding: 10px 20px; border-bottom: 1px solid #dee2e6; }
        .main-content { flex: 1; display: flex; overflow: hidden; }
        .left-panel { width: 45%; padding: 10px; overflow-y: auto; border-right: 1px solid #dee2e6; }
        .right-panel { width: 55%; padding: 10px; overflow-y: auto; display: flex; flex-direction: column; gap: 15px; }
        
        .plot-container { border: 1px solid #eee; border-radius: 5px; padding: 5px; background: white; min-height: 250px; }
        .info-box { background: #e9ecef; padding: 10px; border-radius: 5px; margin-bottom: 10px; }
        
        /* Table styles */
        #networkTable tbody tr { cursor: pointer; }
        #networkTable tbody tr.selected { background-color: #0d6efd !important; color: white; }
        /* Wrap long text in table */
        #networkTable { table-layout: fixed; width: 100%; }
        #networkTable td { 
            white-space: normal !important; 
            word-wrap: break-word !important; 
            word-break: break-word !important;
            overflow-wrap: break-word !important;
        }
    </style>
</head>
<body>

<div class="header">
    <h3>NetInfer Result</h3>
</div>

<div class="main-content">
    <div class="left-panel">
        <table id="networkTable" class="table table-striped table-hover table-sm" style="width:100%">
            <thead>
                <tr>
                    <th>Taxon A</th>
                    <th>Taxon B</th>
                    <th>Taxonomy A</th>
                    <th>Taxonomy B</th>
                    <!-- Hidden columns for searching/data storage if needed -->
                    <th style="display:none;">ID</th> 
                </tr>
            </thead>
            <tbody>
            </tbody>
        </table>
    </div>
    
    <div class="right-panel">
        <div id="infoBox" class="info-box">
            <strong>Select a row to view details</strong>
        </div>
        
        <div id="lineChart" class="plot-container"></div>
        <div id="scatterChart" class="plot-container"></div>
        <div id="barChart" class="plot-container"></div>
    </div>
</div>

<script>
    // Embed data
    const DATA = ##DATA_PLACEHOLDER##;
    
    $(document).ready(function() {
        // Initialize DataTable
        const table = $('#networkTable').DataTable({
            data: DATA.edges,
            columns: [
                { data: 'TaxonA' },
                { data: 'TaxonB' },
                { data: 'TaxonomyA' },
                { data: 'TaxonomyB' },
                { data: 'id', visible: false }
            ],
            scrollY: 'calc(100vh - 150px)',
            scrollCollapse: true,
            paging: true,
            lengthMenu: [[20, 50, 100, -1], [20, 50, 100, "All"]],
            select: true
        });

        // Row Selection Handler
        $('#networkTable tbody').on('click', 'tr', function () {
            if ($(this).hasClass('selected')) {
                $(this).removeClass('selected');
            } else {
                table.$('tr.selected').removeClass('selected');
                $(this).addClass('selected');
                
                const trData = table.row(this).data();
                if (trData) {
                    updateDashboard(trData);
                }
            }
        });
        
        // Auto-select first row
        // if (DATA.edges.length > 0) {
        //     $('#networkTable tbody tr:first').click();
        // }
    });

    function updateDashboard(edge) {
        // 1. Update Info Box
        const html = `
            <h5>Selected Interaction</h5>
            <div class="row">
                <div class="col-md-6">
                    <strong>Taxon A:</strong> ${edge.TaxonA}<br>
                    <small class="text-muted text-break">${edge.TaxonomyA}</small>
                </div>
                <div class="col-md-6">
                    <strong>Taxon B:</strong> ${edge.TaxonB}<br>
                    <small class="text-muted text-break">${edge.TaxonomyB}</small>
                </div>
            </div>
        `;
        $('#infoBox').html(html);

        const abA = DATA.abundance[edge.TaxonA];
        const abB = DATA.abundance[edge.TaxonB];
        const samples = DATA.samples;

        if (!abA || !abB) {
            console.error("Abundance data missing for one of the taxa");
            return;
        }

        // 2. Line Chart (Abundance Profile)
        const traceA = {
            x: samples,
            y: abA,
            mode: 'lines+markers',
            name: 'Taxon A',
            line: {color: '#1f77b4'}
        };
        const traceB = {
            x: samples,
            y: abB,
            mode: 'lines+markers',
            name: 'Taxon B',
            line: {color: '#ff7f0e'}
        };
        
        const layoutLine = {
            title: 'Abundance Profile',
            margin: { t: 30, b: 40, l: 50, r: 20 },
            showlegend: true,
            legend: { x: 0, y: 1.1, orientation: 'h' }
        };
        
        Plotly.newPlot('lineChart', [traceA, traceB], layoutLine, {responsive: true});

        // 3. Scatter Plot + Regression
        // Calculate regression
        const reg = linearRegression(abA, abB);
        
        // Generate regression line points
        const minX = Math.min(...abA);
        const maxX = Math.max(...abA);
        const lineX = [minX, maxX];
        const lineY = [reg.intercept + reg.slope * minX, reg.intercept + reg.slope * maxX];
        
        const scatterTrace = {
            x: abA,
            y: abB,
            text: samples,
            mode: 'markers',
            type: 'scatter',
            name: 'Samples',
            marker: { color: 'rgba(0,0,0,0.5)' }
        };
        
        const regLineTrace = {
            x: lineX,
            y: lineY,
            mode: 'lines',
            type: 'scatter',
            name: 'Regression',
            line: { color: 'red', dash: 'dash' }
        };
        
        const layoutScatter = {
            title: `Correlation (R² = ${reg.r2.toFixed(3)}, Slope = ${reg.slope.toFixed(3)})`,
            xaxis: { title: 'Abundance Taxon A' },
            yaxis: { title: 'Abundance Taxon B' },
            margin: { t: 30, b: 40, l: 50, r: 20 },
            showlegend: false
        };
        
        Plotly.newPlot('scatterChart', [scatterTrace, regLineTrace], layoutScatter, {responsive: true});

        // 4. Bar Chart (Method Weights)
        const methods = Object.keys(edge.methods);
        const weights = Object.values(edge.methods);
        
        const barTrace = {
            x: methods,
            y: weights,
            text: weights.map(w => typeof w === 'number' ? w.toFixed(4) : w),
            textposition: 'auto',
            type: 'bar',
            marker: {
                color: '#66c2a5'
            }
        };
        
        const layoutBar = {
            title: 'Inference Methods Weights',
            xaxis: { title: 'Method' },
            yaxis: { title: 'Weight/Correlation' },
            margin: { t: 30, b: 40, l: 50, r: 20 }
        };
        
        Plotly.newPlot('barChart', [barTrace], layoutBar, {responsive: true});
    }
</script>

</body>
</html>
    """
    
    return html_template.replace('##DATA_PLACEHOLDER##', data_json)

def main(snakemake):
    logger = setup_logger(snakemake.log[0])
    logger.info("Starting dashboard generation")
    
    try:
        # Load Inputs
        logger.info(f"Reading network file: {snakemake.input.network}")
        network_df = pd.read_csv(snakemake.input.network, sep='\t')
        
        logger.info(f"Reading abundance file: {snakemake.input.abundance}")
        abundance_df = pd.read_csv(snakemake.input.abundance, sep='\t', index_col=0)
        
        # Prepare Data
        data = prepare_data(network_df, abundance_df, logger)
        
        # Generate HTML
        data_json = json.dumps(data)
        html_content = generate_html(data_json)
        
        # Output
        out_path = Path(snakemake.output.html)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(out_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
            
        logger.info("Successfully generated dashboard")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        raise e

if __name__ == "__main__":
    main(snakemake)