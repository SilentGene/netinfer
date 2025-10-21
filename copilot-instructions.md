# NetInfer - Copilot Instructions

## Project Overview

**NetInfer** is a Snakemake-based bioinformatics pipeline for inferring co-abundance networks from microbiome data (OTU tables, ASV tables, or any feature abundance tables).

## Core Functionality

### Input
- Microbiome abundance table (OTU/ASV/feature × samples)
  - Formats: TSV, CSV, BIOM
  - Rows: taxonomic features/OTUs
  - Columns: samples
  - Values: abundance/counts
- Optional taxonomy table for features
  - Formats: TSV, CSV
  - Rows: taxonomic features/OTUs
  - Columns: 1: feature ID, 2: taxonomy string
- Optional metadata table for samples
  - Formats: TSV, CSV
  - Rows: samples
  - Columns: metadata variables

### Processing
The pipeline integrates multiple co-abundance network inference methods:
- **FlashWeave (HE mode/normal mode)** (Julia-based)
- **FastSpar** (Python/C++)
- **Jaccard Index** (Python-based)
- **Spearman correlation w/ FDR correction** (R-based)
- **SPIEC-EASI** (R-based)
- **propR** (R-based)
- Additional methods can be added modularly

### Prefiltering parameters for different methods:
- FlashWeave: p<=0.001, weight>=4
- FastSpar: p<=0.05, weight>=0.3
- Spearman: FDR<=0.05, rho>=0.7
- SPIEC-EASI: weight>=0.5
- propR: rho>=0.5
- Jaccard: weight>=0.3

### Output
1. **Summary table**: Aggregated network edges from all methods with consensus scores
   - Taxon A
   - Taxonomy A
   - Taxon B
   - Taxonomy B
   - Prevalence A
   - Prevalence B
   - Mean Abundance A (%)
   - Mean Abundance B (%)
   - Max Abundance A (%)
   - Abundance B at Max A (%)
   - Sample at Max A
   - Max Abundance B (%)
   - Abundance A at Max B (%)
   - Sample at Max B
   - TaxA->TaxB: merge two taxa names into one string for easy searching
   - FlashWeave-hetF Weight
   - FlashWeaveHE Weight
   - SpiecEasi Weight
   - Spearman Cor
   - propR Rho
   - FastSpar Weight
   - Jaccard Index
   - (Any other methods added...)
   - Supporting methods Count

2. **Interactive HTML visualization**: 
   - Final aggregated table embedded
   - Click on a row to view abundance correlation plots across samples
   - Network graph visualization (using vis.js or similar)
   - Filter edges by weight, taxonomy, or any attribute
   - Highlight nodes by degree or other metrics
   - Export filtered networks

## Technical Stack

### Core Technologies
- **Workflow**: Snakemake (Python-based)
- **Languages**: Python (primary), R, Julia
- **Configuration**: YAML files for parameters
- **Environment**: Conda/Mamba for dependency management
- **GitHub**: version control and collaboration

### Key Python Libraries
- `pandas`: data manipulation
- `numpy`: numerical operations
- `networkx`: network analysis
- `plotly`: interactive plotting
- `jinja2`: HTML template rendering
- `pyyaml`: configuration parsing

### Visualization
- **Network**: vis.js, cytoscape.js, or pyvis
- **Plots**: Plotly.js for interactive abundance curves
- **Output**: Single standalone HTML file (no server required)

## Project Structure

```
NetInfer/
├── workflow/
│   ├── Snakefile              # Main workflow
│   ├── rules/                 # Modular rule files
│   │   ├── preprocessing.smk
│   │   ├── flashweave.smk
│   │   ├── fastspar.smk
│   │   ├── spearman.smk
│   │   ├── spieceasi.smk
│   │   ├── (anyothermethods.smk)
│   │   └── visualization.smk
│   ├── scripts/               # Analysis scripts
│   │   ├── aggregate_networks.py
│   │   ├── generate_html.py
│   │   └── utils.py
│   └── envs/                  # Conda environment files
│       ├── julia.yaml
│       ├── python.yaml
│       └── r.yaml
├── config/
│   └── config.yaml            # Pipeline configuration
├── resources/                 # Reference data (if needed)
├── results/                   # Output directory
└── docs/                      # Documentation
```

## Coding Guidelines

### Snakemake Rules
- Each method should be in a separate rule file
- Use `conda:` directive for environment management
- Parameterize all thresholds and options via config.yaml
- Log all outputs with timestamps
- Use `temp()` for intermediate files when appropriate

### Python Scripts
- Use type hints for function signatures
- Follow PEP 8 style guidelines
- Add docstrings for all functions
- Handle errors gracefully with try-except blocks
- Log progress and warnings

### Configuration
- All user-adjustable parameters in `config.yaml`
- Include sensible defaults
- Document each parameter with comments

### Testing
- Include test data (small subset)
- Validate input formats
- Check for required columns/formats

## Key Design Principles

1. **Modularity**: Each inference method is independent
2. **Flexibility**: Easy to add new methods
3. **Reproducibility**: Fixed random seeds, versioned environments
4. **User-friendly**: Clear error messages, progress tracking
5. **Performance**: Parallel execution where possible
6. **Portability**: Single HTML output, no external dependencies for viewing

## Common Tasks for AI Assistant

### When writing Snakemake rules:
- Define clear input/output with proper wildcards
- Use appropriate resources (threads, memory)
- Include shell commands or script directives
- Add log files for debugging

### When writing Python scripts:
- Parse command-line arguments with argparse
- Read/write data with pandas
- Create modular, reusable functions
- Generate publication-quality plots

### When generating HTML:
- Create responsive layouts
- Include interactive controls (sliders, dropdowns, checkboxes)
- Embed all data in the HTML (no external files)
- Optimize for performance (large networks)

### When handling network data:
- Store edges as DataFrames with columns: node1, node2, weight, method, p-value
- Standardize edge directions (e.g., alphabetically sorted nodes)
- Filter by statistical significance
- Calculate network metrics (degree, betweenness, etc.)

## Example Workflow Commands

```bash
# Simple run
netinfer --input abundance_table.tsv --output results_dir --threads 4

# Specify methods
netinfer --input abundance_table.tsv --output results_dir --threads 4 --methods flashweave,fastspar,spearman

# Skip visualization
netinfer --input abundance_table.tsv --output results_dir --threads 4 --no-visual
```

## Notes for AI
- When suggesting code, prefer pandas over loops for data manipulation
- For network visualization, prioritize vis.js for simplicity
- Always validate input data formats before processing
- Consider memory usage for large networks (>10,000 nodes)
- Use multiprocessing for embarrassingly parallel tasks
- Keep the final HTML file size reasonable (<50MB if possible)
