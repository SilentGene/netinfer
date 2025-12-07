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
conda env create -f environment.yaml  # creates 'netinfer' environment and installs dependencies
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
```powershell
# Simple run with all methods enabled and using default settings
netinfer --input abundance_table.tsv --output results_dir --threads 6

# Include taxonomy
netinfer --input abundance_table.tsv --output results_dir --threads 6 --taxonomy taxonomy.tsv

# I don't have taxonomy, but want to infer from feature IDs
netinfer --input abundance_table.tsv --output results_dir --threads 6 --infer-taxonomy

# Only use my favorite methods
netinfer --input abundance_table.tsv --output results_dir --threads 6 --methods flashweave,fastspar,spearman

# Skip visualization
netinfer --input abundance_table.tsv --output results_dir --threads 6 --no-visual

# Expert: Specify every detail via my own config file
netinfer --input abundance_table.tsv --output results_dir --threads 6 --config my_config.yaml
```

### Output

The pipeline generates:
1. Filtered and processed input data (`results/preprocessed/`)
2. Individual network files for each method (`results/networks/`)
3. Aggregated network with consensus scores (`results/networks/aggregated_network.tsv`)
4. Associations between different phyla if `--taxonomy` or `--infer-taxonomy` is specified (`results/networks/phyla_associations.tsv`)
5. Interactive HTML visualization (`results/visualization/network_viewer.html`)

## Input File Formats

### Abundance Table
> No matter what column names are, the first column will always be treated as Feature IDs, and the rest columns as Sample IDs.
- Format: TSV/CSV/BIOM
- Rows: Features (OTUs/ASVs)
- Columns: Samples
- Values: Raw counts or relative abundances

Example:
```
Feature        Sample1  Sample2  Sample3
OTU1           100      150      80
OTU2           50       60       40
...
```

### Taxonomy Table (Optional)
> No matter what column names are, the first column will always be treated as Feature IDs, and the second column as Taxonomy strings.
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
- Supports normal and heterogeneous modes

### FastSpar
- P-value threshold: ≤0.05
- Iterations: 1000
- Correlation threshold: 0.2
- Absolute value threshold for correlation coefficient: 0.3
- This method is computationally intensive; consider using more threads

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

## ALL flags
``` 
$ netinfer --help
usage: netinfer [-h] [--input INPUT] [--output OUTPUT] [--taxonomy TAXONOMY] 
                [--infer-taxonomy] [--metadata METADATA] [--methods METHODS] 
                [--config CONFIG] [--threads THREADS] [--no-visual] 
                [--snake-args SNAKE_ARGS]

NetInfer: Microbiome Network Inference Pipeline

options:
  -h, --help            show this help message and exit
  --input INPUT         Input abundance table file (TSV/CSV/BIOM format). Overrides config file if specified.
  --output OUTPUT       Output directory for results. Overrides config file if specified.
  --taxonomy TAXONOMY   Taxonomy mapping file (optional)
  --infer-taxonomy      Infer taxonomy from feature IDs in the abundance table (looks for 'd__' or 'p__')
  --metadata METADATA   Sample metadata file (optional)
  --methods METHODS     Comma-separated list of methods to use (default: all). Available methods: flashweave, flashweaveHE, fastspar,
                        spearman, spieceasi, propr, jaccard
  --config CONFIG       Path to a base config YAML file. CLI arguments will override settings in this file. The final merged config will be
                        saved to the output directory.
  --threads THREADS     Number of threads to use (default: 1)
  --no-visual           Skip visualization generation
  --snake-args, --snake_args SNAKE_ARGS
                        Additional Snakemake command-line arguments as a single string, e.g. 
                        --snake-args "--unlock --rerun-incomplete --dry-run"
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

Heyu Lin - heyu.lin(SHIFT+2)qut.edu.au

Project Link: [https://github.com/SilentGene/NetInfer](https://github.com/SilentGene/NetInfer)