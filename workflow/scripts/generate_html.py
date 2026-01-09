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
            'PrevalenceA': row.get('Prevalence A (%)', 0.0),
            'PrevalenceB': row.get('Prevalence B (%)', 0.0),
            'methods': {m: (float(row[m]) if pd.notna(row[m]) else 0.0) for m in method_cols}
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

def generate_html(data_json: str, logo_svg: str = "") -> str:
    """Generate the full HTML content."""
    
    html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <link rel="icon" href="https://raw.githubusercontent.com/SilentGene/netinfer/refs/heads/main/docs/logo.svg" sizes="any" type="image/svg+xml">
    <title>NetInfer</title>
    
    <!-- Dependencies -->
    <link href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.3.0/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdn.datatables.net/1.13.4/css/dataTables.bootstrap5.min.css" rel="stylesheet">
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600&display=swap" rel="stylesheet">
    
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.4/jquery.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.3.0/js/bootstrap.bundle.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.4/js/dataTables.bootstrap5.min.js"></script>
    <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
    <!-- Simple Linear Regression Library -->
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
        :root {
            --primary-color: #2c3e50;
            --secondary-color: #546e7a;
            --accent-color: #009688; /* Teal */
            --bg-color: #f4f7f6;
            --card-bg: #ffffff;
            --border-color: #e0e0e0;
            --text-color: #37474f;
            --text-light: #78909c;
            --shadow: 0 2px 8px rgba(0,0,0,0.04);
            --shadow-hover: 0 4px 12px rgba(0,0,0,0.08);
            --font-family: 'Inter', system-ui, -apple-system, sans-serif;
        }

        body { 
            height: 100vh; 
            overflow: hidden; 
            display: flex; 
            flex-direction: column; 
            font-family: var(--font-family);
            background-color: var(--bg-color);
            color: var(--text-color);
        }

        .header { 
            padding: 15px 25px; 
            z-index: 10;
            display: flex;
            align-items: center;
        }

        .header-logo {
            width: 60px;
            margin-right: 15px;
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .header-logo svg {
            width: 100%;
            height: 100%;
        }

        .header h3 {
            margin: 0;
            font-weight: 600;
            font-size: 1.2rem;
            color: var(--primary-color);
            letter-spacing: -0.02em;
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            display: inline-block;
        }

        .header-link {
            text-decoration: none;
            display: flex;
            align-items: center;
        }

        .header-link:hover h3 {
            color: var(--accent-color);
            transform: translateY(-1px) scale(1.02);
            text-shadow: 0 4px 8px rgba(0, 150, 136, 0.15);
        }

        .badge-science {
            background-color: rgba(0, 150, 136, 0.1);
            color: var(--accent-color);
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 0.75rem;
            font-weight: 600;
            margin-left: 12px;
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }

        .main-content { 
            flex: 1; 
            display: flex; 
            flex-direction: column;
            overflow-y: auto; 
            padding: 20px;
            gap: 20px;
        }

        .top-section {
            display: flex;
            gap: 20px;
        }

        .bottom-section {
            height: auto;
            flex-shrink: 0;
            display: flex;
            flex-direction: column;
        }

        .left-panel { 
            width: 50%; 
            min-width: 450px;
            min-height: 600px;
            background: var(--card-bg);
            border-radius: 8px;
            box-shadow: var(--shadow);
            border: 1px solid var(--border-color);
            display: flex;
            flex-direction: column;
        }
        
        .left-panel-content {
            padding: 15px;
            height: 100%;
        }

        .right-panel { 
            width: 50%; 
            display: flex; 
            flex-direction: column; 
            gap: 20px;
            padding-right: 5px; /* Space for scrollbar */
        }
        
        /* Cards */
        .chart-wrapper, .info-box, .scroll-container { 
            border-radius: 8px; 
            padding: 20px; 
            background: var(--card-bg); 
            min-height: 300px;
            box-shadow: var(--shadow);
            border: 1px solid var(--border-color);
        }

        .info-box { 
            min-height: auto;
            background: var(--card-bg);
            border: none;
            box-shadow: none;
            padding: 20px;
        }
        
        .section-label {
            font-size: 0.75rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-light);
            margin-bottom: 5px;
            font-weight: 600;
        }

        .scroll-container {
            overflow-x: auto;
            width: 100%;
        }

        /* Table Aesthetics */
        #networkTable {
            width: 100% !important;
            border-collapse: separate;
            border-spacing: 0;
            font-size: 0.65rem;
        }
        
        #networkTable thead th {
            font-weight: 700;
            color: var(--accent-color);
            background-color: #f8f8f8;
            border-bottom: none;
            padding: 12px 10px;
            text-transform: uppercase;
            font-size: 0.7rem;
            letter-spacing: 0.05em;
        }

        #networkTable td { 
            white-space: normal !important; 
            word-wrap: break-word !important; 
            word-break: break-word !important;
            overflow-wrap: break-word !important;
            border-bottom: 1px solid #f0f0f0;
            color: var(--text-color);
        }
        
        #networkTable tbody tr { transition: background-color 0.15s; }
        #networkTable tbody tr:hover { background-color: #f1f8e9; }
        #networkTable tbody tr.selected { 
            background-color: #b2dfdb !important; 
            color: #004d40 !important;
        }
        /* Fix for DataTables Bootstrap 5 which uses box-shadow on cells */
        #networkTable tbody tr.selected > * {
            box-shadow: inset 0 0 0 9999px #b2dfdb !important;
            color: #004d40 !important;
        }

        /* Pagination Styling */
        .dataTables_paginate {
            font-size: 0.75rem;
        }
        .page-link {
            padding: 5px 10px !important;
            color: var(--accent-color) !important;
            font-size: 1rem;
            border: none;
        }
        .page-item.active .page-link {
            background-color: var(--accent-color) !important;
            border-color: var(--accent-color) !important;
            color: #fff !important;
        }
        .page-link:hover {
            color: var(--secondary-color) !important;
            background-color: #e0f2f1 !important;
        }
        .page-link:focus {
            box-shadow: 0 0 0 0.25rem rgba(0, 150, 136, 0.25) !important;
        }
        .pagination {
            --bs-pagination-disabled-bg: none;
            --bs-pagination-disabled-color: var(--text-color);
        }

        /* Scrollbar Styling */
        ::-webkit-scrollbar { width: 8px; height: 8px; }
        ::-webkit-scrollbar-track { background: transparent; }
        ::-webkit-scrollbar-thumb { background: #cfd8dc; border-radius: 4px; }
        ::-webkit-scrollbar-thumb:hover { background: #b0bec5; }
        
        .dataTables_wrapper .dataTables_info {
            font-size: 0.8rem !important;
            color: var(--text-light) !important;
            padding-top: 10px;
        }

        .dataTables_wrapper .dataTables_filter input {
            border: 1px solid var(--border-color);
            border-radius: 4px;
            padding: 6px 12px;
            margin-left: 8px;
            font-size: 0.9rem;
            outline: none;
        }
        .dataTables_wrapper .dataTables_filter input:focus {
            border-color: var(--accent-color);
            box-shadow: 0 0 0 2px rgba(0, 150, 136, 0.1);
        }

        .table-controls {
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 10px 0;
            gap: 15px;
        }

        .table-controls-left {
            display: flex;
            align-items: center;
            gap: 15px;
        }

        .table-controls-right {
            display: flex;
            align-items: center;
            gap: 15px;
        }

        .dataTables_wrapper .dataTables_length, 
        .dataTables_wrapper .dataTables_filter {
            margin: 0 !important;
            float: none !important;
            display: flex;
            align-items: center;
        }

        .dataTables_wrapper .dataTables_filter input {
            width: 200px;
            height: 36px;
            border: 1px solid transparent;
            border-radius: 8px;
            padding: 0 14px;
            font-size: 0.85rem;
            transition: all 0.25s cubic-bezier(0.4, 0, 0.2, 1);
            background-color: #f0f2f5;
            color: #1a1a1a;
        }

        .dataTables_wrapper .dataTables_filter input:hover {
            background-color: #e4e6e9;
        }

        .dataTables_wrapper .dataTables_filter input:focus {
            border-color: var(--accent-color);
            background-color: #fff;
            box-shadow: 0 0 0 3px rgba(0, 150, 136, 0.15);
            width: 240px;
        }

        .method-selector {
            height: 32px;
            padding: 0 30px 0 12px;
            border-radius: 6px;
            border: 1px solid #ddd;
            background: #fff url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='8' height='8' viewBox='0 0 8 8'%3E%3Cpath fill='%23666' d='M0 2l4 4 4-4H0z'/%3E%3C/svg%3E") no-repeat right 10px center;
            background-size: 8px 8px;
            appearance: none;
            font-size: 0.85rem;
            font-weight: bold;
            color: var(--accent-color);
            cursor: pointer;
            transition: border-color 0.2s;
        }

        .method-selector:focus {
            border-color: var(--accent-color);
            outline: none;
            box-shadow: 0 0 0 3px rgba(0, 150, 136, 0.15);
        }

        .form-control-sm {
            background-color: #f8f8f8;
        }

        .form-control-sm:focus {
            border-color: var(--accent-color);
            outline: none;
            box-shadow: 0 0 0 .25rem rgba(49, 165, 154, 0.25);
        }

        .form-select-sm:focus {
            border-color: var(--accent-color);
            outline: none;
            box-shadow: 0 0 0 .25rem rgba(49, 165, 154, 0.25);
        }

        .dataTables_wrapper .dataTables_length select {
            height: 36px;
            border-radius: 6px;
            border: 1px solid #ddd;
            padding: 0 8px;
            font-size: 0.85rem;
            transition: border-color 0.2s;
        }

        .dataTables_wrapper .dataTables_length select:focus {
            border-color: var(--accent-color);
            outline: none;
            box-shadow: 0 0 0 3px rgba(0, 150, 136, 0.15);
        }

        .control-label {
            font-size: 0.7rem;
            font-weight: 700;
            color: #999;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            margin-right: 8px;
        }
    </style>
</head>
<body>

<div class="header">
    <a href="https://github.com/SilentGene/netinfer" target="_blank" class="header-link">
        <div class="header-logo">##LOGO_PLACEHOLDER##</div>
        <h3>NetInfer Result</h3>
    </a>
    <span class="badge-science">Analysis Dashboard</span>
</div>

<div class="main-content">
    <div class="top-section">
        <div class="left-panel">
            <div class="left-panel-content">
                <div class="table-controls">
                    <div class="table-controls-left">
                        <b style="color: #2c3e50; font-size: 1.1rem; margin-right: 10px;">Network Interactions</b>
                    </div>
                    <div class="table-controls-right" id="customTableTools" style="display: flex; align-items: center; gap: 15px;">
                        <div class="d-flex align-items-center">
                            <select id="methodSelect" class="method-selector"></select>
                        </div>
                        <!-- DataTables Search and Length flow into here -->
                    </div>
                </div>
                <table id="networkTable" class="table table-hover" style="width:100%">
                    <thead>
                        <tr>
                            <th style="width: 20%">Taxon A</th>
                            <th style="width: 20%">Taxon B</th>
                            <th style="width: 20%">Taxonomy A</th>
                            <th style="width: 20%">Taxonomy B</th>
                            <th id="methodHeader" style="width: 20%">Score</th>
                            <th style="display:none;">ID</th> 
                        </tr>
                    </thead>
                    <tbody>
                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="right-panel">
            <div id="infoBox" class="info-box">
                <div class="text-center text-muted" style="padding: 20px;">
                    <em>Select an interaction from the table to start visualizing.</em>
                </div>
            </div>
            
            <div class="row">
                <div class="col-md-6">
                    <!-- Wrapper for styling (padding/bg) -->
                    <div class="chart-wrapper">
                        <div id="scatterChart"></div>
                    </div>
                </div>
                <div class="col-md-6">
                    <!-- Wrapper for styling (padding/bg) -->
                    <div class="chart-wrapper">
                        <div id="barChart"></div>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <div class="bottom-section">
        <div class="scroll-container">
             <div id="lineChart"></div>
        </div>
    </div>
</div>

<script>
    // Embed data
    const DATA = ##DATA_PLACEHOLDER##;
    
    // Plotly common config
    const PLOT_CONFIG = {
        responsive: true,
        displayModeBar: 'hover'
    };
    
    const PLOT_TEMPLATE = 'plotly_white';
    const COLOR_A = '#26a69a'; // Teal Light
    const COLOR_B = '#f5a233'; // Yellow
    const COLOR_REG = '#37474f'; // Dark Grey

    $(document).ready(function() {
        // Initialize Methods Dropdown
        const methodNames = Object.keys(DATA.edges[0].methods);
        let currentMethod = methodNames[0];
        
        const $methodSelect = $('#methodSelect');
        methodNames.forEach(m => {
            $methodSelect.append(`<option value="${m}">${m}</option>`);
        });
        
        $('#methodHeader').text(currentMethod);

        // Initialize DataTable
        const table = $('#networkTable').DataTable({
            data: DATA.edges,
            columns: [
                { data: 'TaxonA' },
                { data: 'TaxonB' },
                { data: 'TaxonomyA' },
                { data: 'TaxonomyB' },
                { 
                    data: null, 
                    render: function(data, type, row) {
                        const val = row.methods[currentMethod];
                        return (typeof val === 'number') ? val.toFixed(4) : val;
                    }
                },
                { data: 'id', visible: false }
            ],
            paging: true,
            pageLength: 10,
            lengthMenu: [[10, 20, 50, 100, -1], [10, 20, 50, 100, "All"]],
            select: {
                style: 'single'
            },
            language: {
                search: "",
                searchPlaceholder: "Search anything...",
                lengthMenu: "_MENU_",
                info: "Showing _START_ to _END_ of _TOTAL_ edges",
                'paginate': {
                    'previous': '<',
                    'next': '>'
                }
            },
            dom: '<"table-header-custom"lf>rtip'
        });

        // Move the DataTables controls into our custom container
        $('.dataTables_length').prepend('<span class="control-label">Show</span>');
        $('.table-header-custom').appendTo('#customTableTools');
        $('.table-header-custom').css({'display': 'flex', 'gap': '15px', 'align-items': 'center'});

        // Handle Method Change
        $methodSelect.on('change', function() {
            currentMethod = $(this).val();
            $('#methodHeader').text(currentMethod);
            table.rows().invalidate().draw();
        });

        // Row Selection Handler
        $('#networkTable tbody').on('click', 'tr', function () {
            if ($(this).hasClass('selected')) {
                // Do nothing if already selected to avoid clearing
            } else {
                table.$('tr.selected').removeClass('selected');
                $(this).addClass('selected');
                
                const trData = table.row(this).data();
                if (trData) {
                    updateDashboard(trData);
                }
            }
        });
    });

    function updateDashboard(edge) {
        // 1. Update Info Box
        const html = `
            <div style="padding-bottom: 15px; font-size: 16px;">
                <b style="color: #455a64;">Selected Interaction Overview</b>
            </div>
            <div class="row g-3">
                <div class="col-md-6">
                    <div style="background: ${COLOR_A}; color: white; padding: 25px; border-radius: 12px; height: 100%; box-shadow: 0 4px 15px rgba(38, 166, 154, 0.15); display: flex; flex-direction: column;">
                        <div class="text-uppercase mb-2" style="font-size: 0.7rem; font-weight: 800; letter-spacing: 0.1em; opacity: 0.8;">Taxon A</div>
                        <h4 class="mb-1" style="font-weight: 700; word-break: break-all;">${edge.TaxonA}</h4>
                        <div class="text-break mb-4" style="font-size: 0.85rem; opacity: 0.9; line-height: 1.4;">${edge.TaxonomyA}</div>
                        <div class="mt-auto d-flex flex-column align-items-center">
                            <div id="prevChartA" style="width: 100px; height: 100px;"></div>
                            <div style="font-size: 0.7rem; font-weight: 700; text-transform: uppercase; margin-top: 5px; opacity: 0.8;">Prevalence</div>
                        </div>
                    </div>
                </div>
                <div class="col-md-6">
                    <div style="background: ${COLOR_B}; color: white; padding: 25px; border-radius: 12px; height: 100%; box-shadow: 0 4px 15px rgba(245, 162, 51, 0.15); display: flex; flex-direction: column;">
                        <div class="text-uppercase mb-2" style="font-size: 0.7rem; font-weight: 800; letter-spacing: 0.1em; opacity: 0.8;">Taxon B</div>
                        <h4 class="mb-1" style="font-weight: 700; word-break: break-all;">${edge.TaxonB}</h4>
                        <div class="text-break mb-4" style="font-size: 0.85rem; opacity: 0.9; line-height: 1.4;">${edge.TaxonomyB}</div>
                        <div class="mt-auto d-flex flex-column align-items-center">
                            <div id="prevChartB" style="width: 100px; height: 100px;"></div>
                            <div style="font-size: 0.7rem; font-weight: 700; text-transform: uppercase; margin-top: 5px; opacity: 0.8;">Prevalence</div>
                        </div>
                    </div>
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
            line: {color: COLOR_A, width: 2},
            marker: {size: 6}
        };
        const traceB = {
            x: samples,
            y: abB,
            mode: 'lines+markers',
            name: 'Taxon B',
            line: {color: COLOR_B, width: 2},
            marker: {size: 6}
        };
        
        // Dynamic width calculation
        const pxPerSample = 30; 
        const containerWidth = $('.scroll-container').width() - 40; // Adjust for padding
        const calcWidth = samples.length * pxPerSample;
        const finalWidth = Math.max(containerWidth, calcWidth);

        const layoutLine = {
            title: { text: '<b style="color: #455a64;">Abundance Profile</b>', font: { size: 16 }, x: 0, xanchor: 'left'},
            width: finalWidth,
            height: 450,
            template: PLOT_TEMPLATE,
            margin: { t: 40, r: 20, l: 50, b: 80 },
            xaxis: {
                automargin: true,
                tickangle: 45,
                gridcolor: '#f0f0f0'
            },
            yaxis: {
                gridcolor: '#f0f0f0',
                zerolinecolor: '#e0e0e0'
            },
            showlegend: true,
            legend: { 
                x: 0.02, 
                y: 0.98, 
                xanchor: 'left', 
                yanchor: 'top', 
                bgcolor: 'rgba(255, 255, 255, 0.8)',
                bordercolor: '#000000',
                borderwidth: 1
            },
            font: { family: 'Inter' }
        };
        
        Plotly.newPlot('lineChart', [traceA, traceB], layoutLine, {responsive: false, displayModeBar: 'hover'});

        // 2.5 Prevalence Donut Charts
        function renderDonut(divId, value, color, name) {
            const data = [{
                values: [value, 100 - value],
                labels: ['Present', 'Absent'],
                type: 'pie',
                hole: 0.75,
                marker: { colors: ['#ffffff', 'rgba(255, 255, 255, 0.2)'] },
                textinfo: 'none',
                hoverinfo: 'label+percent',
                showlegend: false
            }];
            const layout = {
                height: 100,
                width: 100,
                margin: { t: 0, b: 0, l: 0, r: 0 },
                annotations: [{
                    text: value.toFixed(1) + '%',
                    showarrow: false,
                    font: { size: 12, weight: 'bold', color: '#ffffff' }
                }],
                paper_bgcolor: 'rgba(0,0,0,0)',
                plot_bgcolor: 'rgba(0,0,0,0)'
            };
            Plotly.newPlot(divId, data, layout, {displayModeBar: false, staticPlot: false});
        }
        
        renderDonut('prevChartA', edge.PrevalenceA, COLOR_A, 'Taxon A');
        renderDonut('prevChartB', edge.PrevalenceB, COLOR_B, 'Taxon B');

        // 3. Scatter Plot + Regression
        const reg = linearRegression(abA, abB);
        
        const minX = Math.min(...abA);
        const maxX = Math.max(...abA);
        // Extend line slightly for visuals
        const span = maxX - minX;
        const lineX = [minX, maxX];
        const lineY = [reg.intercept + reg.slope * minX, reg.intercept + reg.slope * maxX];
        
        const scatterTrace = {
            x: abA,
            y: abB,
            text: samples,
            mode: 'markers',
            type: 'scatter',
            name: 'Samples',
            marker: { 
                color: '#546e7a',
                size: 8,
                opacity: 0.6,
                line: { color: 'white', width: 1 }
            }
        };
        
        const regLineTrace = {
            x: lineX,
            y: lineY,
            mode: 'lines',
            type: 'scatter',
            name: 'Regression',
            line: { color: COLOR_REG, dash: 'dash', width: 2 }
        };
        
        const layoutScatter = {

            title: { 
                text: `<b style="color: #455a64;">Correlation</b><br><span style="font-size: 13px; color: #546e7a;">R² = ${reg.r2.toFixed(3)}, Slope = ${reg.slope.toFixed(3)}</span>`, 
                font: { size: 16 }, 
                x: 0, 
                xanchor: 'left',
                y: 0.96,
                yanchor: 'top'
            },
            height: 300,
            template: PLOT_TEMPLATE,
            xaxis: { title: 'Abundance Taxon A', gridcolor: '#f0f0f0', zerolinecolor: '#e0e0e0' },
            yaxis: { title: 'Abundance Taxon B', gridcolor: '#f0f0f0', zerolinecolor: '#e0e0e0' },
            margin: { t: 50, b: 40, l: 50, r: 20 },
            showlegend: false,
            font: { family: 'Inter' }
        };
        
        Plotly.newPlot('scatterChart', [scatterTrace, regLineTrace], layoutScatter, PLOT_CONFIG);

        // 4. Bar Chart
        const methods = Object.keys(edge.methods);
        const weights = Object.values(edge.methods);
        
        const barTrace = {
            x: methods,
            y: weights,
            text: weights.map(w => typeof w === 'number' ? w.toFixed(2) : w),
            textposition: 'auto',
            type: 'bar',
            marker: {
                color: COLOR_A,
                opacity: 0.8
            }
        };
        
        const layoutBar = {
            title: { 
                text: '<b style="color: #455a64;">Inference Methods Weights</b>', 
                font: { size: 16 }, 
                x: 0, 
                xanchor: 'left',
                y: 0.96,
                yanchor: 'top'
            },
            height: 300,
            template: PLOT_TEMPLATE,
            xaxis: { gridcolor: '#f0f0f0', automargin: true },
            yaxis: { title: 'Weight/Correlation', gridcolor: '#f0f0f0' },
            margin: { t: 50, b: 40, l: 50, r: 20 },
            font: { family: 'Inter' }
        };
        
        Plotly.newPlot('barChart', [barTrace], layoutBar, PLOT_CONFIG);
    }
</script>

</body>
</html>
    """
    
    return html_template.replace('##DATA_PLACEHOLDER##', data_json).replace('##LOGO_PLACEHOLDER##', logo_svg)

def main(snakemake):
    logger = setup_logger(snakemake.log[0])
    logger.info("Starting dashboard generation")
    
    try:
        # Load logo if exists
        logo_svg = ""
        try:
            # Path relative to script: ../../docs/logo.svg
            script_path = Path(__file__).resolve()
            logo_path = script_path.parent.parent.parent / 'docs' / 'logo.svg'
            if logo_path.exists():
                with open(logo_path, 'r', encoding='utf-8') as f:
                    logo_svg = f.read()
                    # Remove XML declaration if present
                    if '<?xml' in logo_svg:
                        logo_svg = logo_svg[logo_svg.find('?>')+2:].strip()
                logger.info(f"Loaded logo from {logo_path}")
            else:
                logger.warning(f"Logo not found at {logo_path}")
        except Exception as e:
            logger.warning(f"Could not load logo: {e}")

        # Load Inputs
        logger.info(f"Reading network file: {snakemake.input.network}")
        network_df = pd.read_csv(snakemake.input.network, sep='\t')
        
        logger.info(f"Reading abundance file: {snakemake.input.abundance}")
        abundance_df = pd.read_csv(snakemake.input.abundance, sep='\t', index_col=0)
        
        # Prepare Data
        data = prepare_data(network_df, abundance_df, logger)
        
        # Generate HTML
        data_json = json.dumps(data)
        html_content = generate_html(data_json, logo_svg)
        
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