# NetInfer

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io)
[![Conda](https://img.shields.io/badge/conda-compatible-brightgreen.svg)](https://docs.conda.io/en/latest/)

A Snakemake-based bioinformatics pipeline for inferring co-abundance networks from microbiome data.

> [!NOTE]
> NetInfer is still under development. Please use with caution and report any issues on the GitHub repository.

## Features

- Multiple network inference methods:
  - Direct Association Methods:
      - [FlashWeave (HE mode/normal mode)](https://github.com/meringlab/FlashWeave.jl)
      - [SPIEC-EASI](https://github.com/zdk123/SpiecEasi)
  - Correlation-based Methods:
      - [propr](https://www.nature.com/articles/s41598-017-16520-0)
      - [Spearman correlation with FDR correction](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient)
      - [Jaccard Index](https://en.wikipedia.org/wiki/Jaccard_index)
      - [FastSpar](https://github.com/scwatts/fastspar)


- Interactive HTML visualization
- Flexible input formats (TSV, CSV, BIOM)
- Comprehensive output with consensus scores
- Parallel execution support
- Conda environment management

## Quick Start

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html) (or [Mamba](https://github.com/mamba-org/mamba))

### Installation

```bash
# Clone the repository
git clone https://github.com/SilentGene/NetInfer.git
cd NetInfer

# Create and activate conda environment
conda env create -f environment.yaml  # it takes a while
conda activate netinfer

# Install FlashWeave (Julia package)
julia workflow/scripts/install_flashweave.jl
```

### Basic Usage

1. Prepare your data:
   - Required: abundance table `abundance_table.tsv`
   - (Optional) taxonomy table `taxonomy.tsv` for additional annotations
   - (Optional) sample metadata table `metadata.tsv` for `FlashWeave` and additional annotations

2. Run the pipeline:
```bash
# Simple run with all methods enabled and using default settings
netinfer --input abundance_table.tsv --output results_dir --threads 6

# Specify methods
netinfer --input abundance_table.tsv --output results_dir --threads 6 --methods flashweave,fastspar,spearman

# Skip visualization
netinfer --input abundance_table.tsv --output results_dir --threads 6 --no-visual

# Specify every detail via my own config file
netinfer --input abundance_table.tsv --output results_dir --threads 6 --config my_config.yaml
```

### Output

The pipeline generates:
1. Filtered and processed input data (`results/preprocessed/`)
2. Individual network files for each method (`results/networks/`)
3. Aggregated network with consensus scores (`results/networks/aggregated_network.tsv`)
4. Interactive HTML visualization (`results/visualization/network_viewer.html`)

## Input File Formats

### Abundance Table
- Format: TSV/CSV/BIOM
- Rows: Features (OTUs/ASVs)
- Columns: Samples
- Values: Raw counts or relative abundances

Example:
```
Feature         Sample1  Sample2  Sample3
Otu1           100      150      80
Otu2           50       60       40
...
```

### Taxonomy Table (Optional)
- Format: TSV/CSV
- Required columns: 
  - Feature ID (matching abundance table)
  - Taxonomy string

Example:
```
Feature  Taxonomy
Otu1     k__Bacteria;p__Firmicutes;c__Clostridia
Otu2     k__Bacteria;p__Bacteroidetes;c__Bacteroidia
...
```

### Metadata Table (Optional)
- Format: TSV/CSV
- Rows: Samples (matching abundance table)
- Columns: Metadata variables

## Configuration

See `config/config.yaml` for all available parameters and their descriptions.

Key parameters:
- Method-specific thresholds
- Filtering criteria
- Visualization options
- Resource allocation

## Method Details

### FlashWeave
- P-value threshold: ≤0.001
- Weight threshold: ≥4.0
- Supports heterogeneous mode

### FastSpar
- P-value threshold: ≤0.05
- Iterations: 1000
- Correlation threshold: 0.2
- Absolute value threshold for correlation coefficient: 0.3

### Spearman
- FDR threshold: ≤0.05
- Correlation threshold: ≥0.7

### SPIEC-EASI
- Weight threshold: ≥0.5
- Methods: MB (Meinshausen-Buhlmann) or GLasso

### propR
- Correlation threshold: ≥0.5

### Jaccard
- Similarity threshold: ≥0.3

## Visualization

The interactive HTML viewer provides:
- Network graph visualization
- Edge filtering by weight/method
- Node highlighting by degree
- Abundance correlation plots
- Taxonomy-based coloring
- Export capabilities


## Troubleshooting

If you see an error like:

```
LockException:
Error: Directory cannot be locked. Please make sure that no other Snakemake process is trying to create the same files...
```

it usually means a previous Snakemake run was interrupted and left a lock behind. You can safely unlock and continue.

```bash
netinfer <original_args> --snake_args="--unlock"
```


## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Citation

If you use NetInfer in your research, please cite:

```bibtex
@software{netinfer2025,
  author = {Heyu Lin},
  title = {NetInfer: A Snakemake Pipeline for Microbiome Network Inference},
  year = {2025},
  url = {https://github.com/SilentGene/NetInfer}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

Heyu Lin - heyu.lin@qut.edu.au

Project Link: [https://github.com/SilentGene/NetInfer](https://github.com/SilentGene/NetInfer)