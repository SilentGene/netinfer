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

def create_config(input_file: Optional[str] = None,
                 output_dir: Optional[str] = None,
                 taxonomy_file: Optional[str] = None,
                 metadata_file: Optional[str] = None,
                 methods: Optional[List[str]] = None,
                 no_visual: bool = False,
                 infer_taxonomy: bool = False,
                 base_config_path: Optional[str] = None,
                 suffix: Optional[str] = None) -> str:
    """Create a config file for the pipeline run.
    
    Args:
        input_file: Input abundance table path (overrides config if provided)
        output_dir: Output directory path (overrides config if provided)
        taxonomy_file: Taxonomy table path (overrides config if provided)
        metadata_file: Metadata table path (overrides config if provided)
        methods: List of methods to enable (overrides config if provided)
        no_visual: Whether to disable visualization (overrides config if provided)
        infer_taxonomy: Whether to infer taxonomy from abundance table
        base_config_path: Path to base config file (uses default if None)
    
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
    if input_file:
        config["input"]["abundance_table"] = os.path.abspath(input_file)
    if taxonomy_file:
        config["input"]["taxonomy_table"] = os.path.abspath(taxonomy_file)
    if metadata_file:
        config["input"]["metadata_table"] = os.path.abspath(metadata_file)
    
    # Set infer_taxonomy flag
    config["infer_taxonomy"] = infer_taxonomy

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
        "--use-conda",
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
        "--input",
        help="Input abundance table file (TSV/CSV/BIOM format). Overrides config file if specified."
    )
    
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
    print_python_version(logger)
    print_snakemake_version(logger)
    
    try:
        # Parse methods if specified
        methods = None
        if args.methods:
            methods = [m.strip().lower() for m in args.methods.split(",")]
        
        # Validate that either --input/--output or --config is provided
        if not args.config and (not args.input or not args.output):
            raise ValueError(
                "Either --config must be specified, or both --input and --output must be provided"
            )
        
        # Create configuration (with CLI overrides applied)
        logger.info("Creating pipeline configuration...")
        config_path = create_config(
            input_file=args.input,
            output_dir=args.output,
            taxonomy_file=args.taxonomy,
            metadata_file=args.metadata,
            methods=methods,
            no_visual=args.no_visual,
            infer_taxonomy=args.infer_taxonomy,
            base_config_path=args.config,
            suffix=args.suffix
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