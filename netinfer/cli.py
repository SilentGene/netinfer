#!/usr/bin/env python3
"""
NetInfer: A command-line tool for microbiome network inference
"""

import argparse
import os
import sys
import yaml
import logging
import subprocess
from pathlib import Path
from typing import List, Optional
import shlex
from importlib import resources

from netinfer import __version__

def setup_logger() -> logging.Logger:
    """Set up logging configuration."""
    logger = logging.getLogger('netinfer')
    logger.setLevel(logging.INFO)
    
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    return logger

def print_snakemake_version(logger: logging.Logger) -> None:
    """Log the installed Snakemake version early so users see environment info."""
    try:
        import snakemake  # type: ignore
        version = getattr(snakemake, "__version__", None) or getattr(snakemake, "version", None)
        if version is None:
            # Fallback: some rare builds expose version via module attribute
            version = str(snakemake)
        logger.info(f"Snakemake version: {version}")
    except Exception as e:
        logger.warning(f"Could not determine Snakemake version: {e}")

def print_python_version(logger: logging.Logger) -> None:
    """Log the Python runtime version early for reproducibility."""
    try:
        vinfo = sys.version_info
        version = f"{vinfo.major}.{vinfo.minor}.{vinfo.micro}"
        # Optional: include implementation (CPython, PyPy, etc.)
        import platform
        impl = platform.python_implementation()
        logger.info(f"Python version: {version} ({impl})")
    except Exception as e:
        logger.warning(f"Could not determine Python version: {e}")

def find_config_file() -> Path:
    """Find the config file using modern importlib.resources."""
    # First try current directory
    cwd_config = Path.cwd() / "config.yaml"
    if cwd_config.is_file():
        return cwd_config
        
    # Then try using importlib.resources (Standard package layout)
    try:
        config_path = resources.files("netinfer").joinpath("config/config.yaml")
        if config_path.is_file():
            return Path(config_path)
    except (ImportError, TypeError):
        pass
        
    # If not found, provide detailed error message
    searched_paths = [
        str(cwd_config),
        "netinfer/config/config.yaml (package resource)"
    ]
    raise FileNotFoundError(
        "Could not find config.yaml in any of the expected locations:\n" +
        "\n".join(f"- {p}" for p in searched_paths)
    )

def find_workflow_dir() -> Path:
    """Find the workflow directory using modern importlib.resources."""
    # First try current directory
    cwd_workflow = Path.cwd() / "workflow"
    if (cwd_workflow / "Snakefile").exists():
        return cwd_workflow
        
    # Then try using importlib.resources (Standard package layout)
    try:
        workflow_dir = resources.files("netinfer").joinpath("workflow")
        if (workflow_dir / "Snakefile").exists():
            return Path(workflow_dir)
    except (ImportError, TypeError):
        pass
        
    # If not found, provide detailed error message
    searched_paths = [
        str(cwd_workflow),
        "netinfer/workflow (package resource)"
    ]
    raise FileNotFoundError(
        "Could not find workflow directory with Snakefile in any of the expected locations:\n" +
        "\n".join(f"- {p}" for p in searched_paths)
    )

def create_config(count_file: Optional[str] = None,
                 relab_file: Optional[str] = None,
                 other_file: Optional[str] = None,
                 output_dir: Optional[str] = None,
                 taxonomy_file: Optional[str] = None,
                 metadata_file: Optional[str] = None,
                 methods: Optional[List[str]] = None,
                 no_visual: bool = False,
                 infer_taxonomy: bool = False,
                 include_negative: Optional[bool] = None,
                 base_config_path: Optional[str] = None,
                 suffix: Optional[str] = None,
                 filter_params: Optional[dict] = None) -> str:
    """Create a config file for the pipeline run.
    
    Args:
        count_file: Read count table path (overrides config if provided)
        relab_file: Relative abundance table path (overrides config if provided)
        other_file: Other abundance table path, e.g. coverage (overrides config if provided)
        output_dir: Output directory path (overrides config if provided)
        taxonomy_file: Taxonomy table path (overrides config if provided)
        metadata_file: Metadata table path (overrides config if provided)
        methods: List of methods to enable (overrides config if provided)
        no_visual: Whether to disable visualization (overrides config if provided)
        infer_taxonomy: Whether to infer taxonomy from abundance table
        include_negative: Whether signed association methods should filter
            on absolute values, retaining sufficiently strong negative edges.
            If None, preserve the base configuration value.
        base_config_path: Path to base config file (uses default if None)
        filter_params: Preprocessing thresholds to override (keys as in config;
            None values are ignored)
    
    Returns:
        Path to the final config file written to output directory
    """
    # Load base config (user-specified or default)
    if base_config_path:
        config_file = Path(base_config_path)
        if not config_file.is_file():
            raise FileNotFoundError(f"Specified config file not found: {base_config_path}")
    else:
        config_file = find_config_file()
    
    logger = logging.getLogger('netinfer')
    logger.info(f"Using base configuration from: {config_file}")
    
    with open(config_file) as f:
        config = yaml.safe_load(f)
    
    # Update input/output paths (CLI overrides config)
    config.setdefault("input", {})
    if str(config["input"].get("abundance_table", "") or "").strip():
        raise ValueError(
            "Config key 'input.abundance_table' is no longer supported. Use 'input.count_table' "
            "(raw read counts) and/or 'input.relative_abundance_table' instead."
        )
    config["input"].pop("abundance_table", None)
    if count_file:
        config["input"]["count_table"] = os.path.abspath(count_file)
    if relab_file:
        config["input"]["relative_abundance_table"] = os.path.abspath(relab_file)
    if other_file:
        config["input"]["other_table"] = os.path.abspath(other_file)
    given = {k: str(config["input"].get(k, "") or "").strip()
             for k in ("count_table", "relative_abundance_table", "other_table")}
    if not any(given.values()):
        raise ValueError(
            "No abundance table given. Use --input_count (raw read counts, recommended), "
            "--input_relative_abundance and/or --input_others, or set them in the config file."
        )
    if given["count_table"] and given["other_table"]:
        raise ValueError(
            "--input_count and --input_others cannot be combined: read counts are always "
            "preferred for network inference. Drop --input_others."
        )
    if taxonomy_file:
        config["input"]["taxonomy_table"] = os.path.abspath(taxonomy_file)
    if metadata_file:
        config["input"]["metadata_table"] = os.path.abspath(metadata_file)
    
    # Set infer_taxonomy flag
    config["infer_taxonomy"] = infer_taxonomy

    # By default, retain only positive associations. When requested, apply
    # method thresholds to absolute values so strong negative edges are kept.
    if include_negative is not None:
        config["include_negative"] = include_negative

    # Preprocessing (filtering) thresholds: CLI overrides config
    for key, value in (filter_params or {}).items():
        if value is not None:
            config[key] = value

    # Sanitize and set optional output suffix
    if suffix:
        try:
            import re
            # Replace any char not alphanumeric, dash or underscore with underscore
            clean = re.sub(r"[^A-Za-z0-9_-]", "_", suffix.strip())
            # Collapse multiple underscores
            clean = re.sub(r"_+", "_", clean)
            # Trim underscores from ends
            clean = clean.strip("_")
        except Exception:
            clean = suffix.strip()
        config["suffix"] = clean
    else:
        # Ensure key exists for downstream logic
        config["suffix"] = ""
    
    # Determine output directory (required either from CLI or config)
    if output_dir:
        output_dir_abs = os.path.abspath(output_dir).strip()
        config["output_dir"] = output_dir_abs
    else:
        # Use output_dir from config file
        if "output_dir" not in config or not config["output_dir"]:
            raise ValueError("output_dir must be specified either via --output or in config file")
        output_dir_abs = os.path.abspath(config["output_dir"]).strip()
        config["output_dir"] = output_dir_abs
    
    # Update method selection
    if methods:
        # Define available methods and their config mappings
        method_mappings = {
            "flashweave": {"config_key": "flashweave", "he_mode": False},
            "flashweaveHE": {"config_key": "flashweave", "he_mode": True},
            "fastspar": {"config_key": "fastspar"},
            "pearson": {"config_key": "pearson"},
            "spearman": {"config_key": "spearman"},
            "spieceasi": {"config_key": "spieceasi"},
            "propr": {"config_key": "propr"},
            "jaccard": {"config_key": "jaccard"}
        }
        
        # Validate methods
        invalid_methods = [m for m in methods if m not in method_mappings]
        if invalid_methods:
            raise ValueError(f"Invalid method(s): {', '.join(invalid_methods)}. "
                           f"Available methods: {', '.join(method_mappings.keys())}")
        
        # Disable all methods by default
        for method in method_mappings.values():
            config_key = method["config_key"]
            if config_key in config:
                config[config_key]["enabled"] = False
        
        # Enable only selected methods
        for method in methods:
            mapping = method_mappings[method]
            config_key = mapping["config_key"]
            config[config_key]["enabled"] = True
            
            # Handle special cases
            if config_key == "flashweave":
                config[config_key]["heterogeneous"] = mapping.get("he_mode", False)
    
    # Update visualization
    if "visualization" not in config:
        config["visualization"] = {}
    config["visualization"]["enabled"] = not no_visual
    
    # Record which NetInfer version produced this run
    config["netinfer_version"] = __version__

    # Create temporary config file with clean YAML output (use absolute paths)
    os.makedirs(output_dir_abs, exist_ok=True)
    config_path = os.path.abspath(os.path.join(output_dir_abs, "config.yaml"))
    with open(config_path, "w") as f:
        # Use default_flow_style=False to ensure clean formatting
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    
    return config_path

def run_pipeline(config_path: str, threads: int, extra_snake_args: Optional[List[str]] = None, working_dir: Optional[str] = None) -> int:
    """Run the Snakemake pipeline with given configuration.

    Uses `python -m snakemake` for compatibility with Snakemake ≥9 where the
    old programmatic API (snakemake.snakemake) is no longer available.
    """
    workflow_dir = find_workflow_dir()
    snakefile = workflow_dir / "Snakefile"

    # Determine Snakemake working directory
    # Prefer the explicit working_dir argument (user-specified output directory)
    snakemake_directory: Optional[str] = None
    if working_dir:
        snakemake_directory = os.path.abspath(working_dir)

    # Build base command
    cmd = [
        sys.executable, "-m", "snakemake",
        "--snakefile", str(snakefile),
        "--cores", str(threads),
        "--configfiles", config_path,
        "--printshellcmds",
        "--show-failed-logs",
        "--keep-incomplete"
    ]

    # Force Snakemake working directory to output_dir
    if snakemake_directory:
        cmd.extend(["--directory", snakemake_directory])

    if extra_snake_args:
        cmd.extend(extra_snake_args)

    logger = logging.getLogger('netinfer')
    logger.info("Invoking Snakemake via module runner…")
    logger.debug("Command: %s", " ".join(cmd))

    # First, attempt to clear any stale lock in the workflow directory.
    # Safe when no lock exists; required if a previous run was interrupted.
    unlock_cmd = [
        sys.executable, "-m", "snakemake",
        "--snakefile", str(snakefile),
        "--unlock",
    ]
    logger.debug("Pre-run unlock (best-effort): %s", " ".join(unlock_cmd))
    try:
        subprocess.run(unlock_cmd, cwd=str(workflow_dir), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as e:
        logger.debug("Unlock attempt ignored due to error: %s", e)

    proc = subprocess.run(cmd, cwd=str(workflow_dir))
    return proc.returncode

def main():
    """Main entry point for the command-line interface."""
    parser = argparse.ArgumentParser(
        description="NetInfer: Microbiome Network Inference Pipeline"
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"NetInfer {__version__}"
    )
    
    parser.add_argument(
        "--input_count", "--input-count",
        dest="input_count",
        help=(
            "Read count table (TSV/CSV/BIOM; features x samples), e.g. from `coverm genome -m count`. "
            "Recommended; used for network inference. Overrides config file if specified."
        )
    )
    parser.add_argument(
        "--input_relative_abundance", "--input-relative-abundance",
        dest="input_relative_abundance",
        help=(
            "Relative abundance table (TSV/CSV/BIOM; features x samples; fraction or percentage), "
            "e.g. from `coverm genome -m relative_abundance`. Used for abundance filtering and display; "
            "used for network inference only if --input_count is not given. Overrides config file if specified."
        )
    )
    parser.add_argument(
        "--input_others", "--input-others",
        dest="input_others",
        help=(
            "Other abundance table (TSV/CSV/BIOM; features x samples) with non-negative values that are "
            "neither read counts nor relative abundance, e.g. CoverM `trimmed_mean` or sandpiper "
            "`full_coverage`. Not optimal: use only if read counts are unavailable. Used for network "
            "inference; cannot be combined with --input_count."
        )
    )
    # Removed option: kept only to give a clear error instead of passing it to Snakemake
    parser.add_argument("--input", dest="legacy_input", help=argparse.SUPPRESS)
    
    parser.add_argument(
        "--output",
        help="Output directory for results. Overrides config file if specified."
    )
    
    parser.add_argument(
        "--taxonomy",
        help="Taxonomy mapping file (optional)"
    )
    
    parser.add_argument(
        "--infer-taxonomy",
        action="store_true",
        help="Infer taxonomy from feature IDs in the abundance table (looks for 'd__' or 'p__')"
    )
    
    parser.add_argument(
        "--metadata",
        help="Sample metadata file (optional)"
    )
    
    parser.add_argument(
        "--methods",
        help="Comma-separated list of methods to use (default: all). "
             "Available methods: flashweave, flashweaveHE, fastspar, "
             "pearson, spearman, spieceasi, propr, jaccard"
    )
    
    parser.add_argument(
        "--config",
        help=(
            "Path to a base config YAML file. CLI arguments will override settings in this file. "
            "The final merged config will be saved to the output directory."
        )
    )
    
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use (default: 1)"
    )
    
    parser.add_argument(
        "--no-visual",
        action="store_true",
        help="Skip visualization generation"
    )
    parser.add_argument(
        "--include-negative",
        action="store_true",
        default=None,
        help=(
            "Apply signed association thresholds to absolute values, retaining "
            "strong negative as well as positive associations"
        )
    )
    filt = parser.add_argument_group("abundance filtering (overrides config; see README 'Preprocessing')")
    filt.add_argument("--min-sample-depth", type=float, default=None,
                      help="Drop samples with fewer total reads (integer counts only; 0 disables). Default: 1000")
    filt.add_argument("--min-prevalence", type=float, default=None,
                      help="Minimum fraction of samples in which a feature is present (> 0). Default: 0.1")
    filt.add_argument("--min-prevalent-samples", type=int, default=None,
                      help="Minimum number of samples in which a feature is present (0 disables). Default: 0")
    filt.add_argument("--min-peak-rel-abundance", type=float, default=None,
                      help="Minimum peak (max per-sample) relative abundance of a feature, e.g. 0.001 = 0.1%%. Default: 0.001")
    parser.add_argument(
        "--snake-args", "--snake_args",
        dest="snake_args",
        help=(
            "Additional Snakemake command-line arguments as a single string, "
            "e.g. --snake-args \"--unlock --rerun-incomplete --dry-run\""
        )
    )
    parser.add_argument(
        "--suffix",
        help=(
            "Optional suffix to append to final aggregated outputs (inserted before file extension). "
            "Invalid filename characters will be replaced with '_'"
        )
    )
    
    # Parse known args but allow passthrough of unknown ones (to Snakemake)
    args, unknown_snake_args = parser.parse_known_args()
    logger = setup_logger()
    # Print environment versions before any heavy work
    logger.info(f"NetInfer version: {__version__}")
    print_python_version(logger)
    print_snakemake_version(logger)
    
    try:
        # Parse methods if specified
        methods = None
        if args.methods:
            methods = [m.strip().lower() for m in args.methods.split(",")]
        
        if args.legacy_input:
            raise ValueError(
                "--input is no longer supported. Use --input_count for raw read counts (recommended) "
                "and/or --input_relative_abundance for relative abundance (or --input_others for e.g. coverage)."
            )
        # Validate that either an input table and --output, or --config is provided
        if not args.config and (not (args.input_count or args.input_relative_abundance or args.input_others)
                                or not args.output):
            raise ValueError(
                "Either --config must be specified, or --output and at least one of "
                "--input_count / --input_relative_abundance / --input_others must be provided"
            )
        
        # Create configuration (with CLI overrides applied)
        logger.info("Creating pipeline configuration...")
        config_path = create_config(
            count_file=args.input_count,
            relab_file=args.input_relative_abundance,
            other_file=args.input_others,
            output_dir=args.output,
            taxonomy_file=args.taxonomy,
            metadata_file=args.metadata,
            methods=methods,
            no_visual=args.no_visual,
            infer_taxonomy=args.infer_taxonomy,
            include_negative=args.include_negative,
            base_config_path=args.config,
            suffix=args.suffix,
            filter_params={
                "min_sample_depth": args.min_sample_depth,
                "min_prevalence": args.min_prevalence,
                "min_prevalent_samples": args.min_prevalent_samples,
                "min_peak_rel_abundance": args.min_peak_rel_abundance,
            }
        )
        logger.info(f"Configuration file created at: {config_path}")
        
        # Prepare additional Snakemake args (passthrough)
        extra_snake_args: Optional[List[str]] = None
        if args.snake_args:
            try:
                extra_snake_args = shlex.split(args.snake_args)
            except Exception:
                # Fallback: naive split by whitespace
                extra_snake_args = args.snake_args.split()
        # Also include any unknown args from our parser to forward to Snakemake
        if unknown_snake_args:
            if extra_snake_args is None:
                extra_snake_args = []
            extra_snake_args.extend(unknown_snake_args)

        # Run pipeline
        logger.info("Starting pipeline execution...")
        # Force Snakemake working directory to the user-specified output directory when provided
        result = run_pipeline(config_path, args.threads, extra_snake_args, working_dir=args.output)
        
        if result == 0:
            logger.info("Pipeline completed successfully!")
        else:
            logger.error("Pipeline execution failed.")
        
        return result
        
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
