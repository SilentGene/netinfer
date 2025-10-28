#!/usr/bin/env python3
"""
NetInfer: A command-line tool for microbiome network inference
"""

import argparse
import os
import sys
import yaml
import shutil
import logging
from pathlib import Path
from typing import List, Optional
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

def find_config_file() -> Path:
    """Find the config file using modern importlib.resources."""
    # First try current directory
    cwd_config = Path.cwd() / "config" / "config.yaml"
    if cwd_config.exists():
        return cwd_config
        
    # Then try using importlib.resources
    try:
        with resources.files("netinfer").joinpath("config/config.yaml") as config_path:
            if config_path.exists():
                return config_path
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
        
    # Then try using importlib.resources
    try:
        with resources.files("netinfer").joinpath("workflow") as workflow_dir:
            if (workflow_dir / "Snakefile").exists():
                return workflow_dir
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

def create_config(input_file: str,
                 output_dir: str,
                 taxonomy_file: Optional[str] = None,
                 metadata_file: Optional[str] = None,
                 methods: Optional[List[str]] = None,
                 no_visual: bool = False) -> str:
    """Create a temporary config file for the pipeline run."""
    # Load default config
    config_file = find_config_file()
    with open(config_file) as f:
        config = yaml.safe_load(f)
    
    # Update input/output paths
    config["input"]["abundance_table"] = os.path.abspath(input_file)
    if taxonomy_file:
        config["input"]["taxonomy_table"] = os.path.abspath(taxonomy_file)
    if metadata_file:
        config["input"]["metadata_table"] = os.path.abspath(metadata_file)
    
    config["output_dir"] = os.path.abspath(output_dir)
    
    # Update method selection
    if methods:
        # Define available methods and their config mappings
        method_mappings = {
            "flashweave": {"config_key": "flashweave", "he_mode": False},
            "flashweaveHE": {"config_key": "flashweave", "he_mode": True},
            "fastspar": {"config_key": "fastspar"},
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
    
    # Create temporary config file
    os.makedirs(output_dir, exist_ok=True)
    config_path = os.path.join(output_dir, "config.yaml")
    with open(config_path, "w") as f:
        yaml.dump(config, f)
    
    return config_path

def run_pipeline(config_path: str, threads: int) -> int:
    """Run the Snakemake pipeline with given configuration."""
    import snakemake
    
    workflow_dir = find_workflow_dir()
    snakefile = workflow_dir / "Snakefile"
    
    success = snakemake.snakemake(
        str(snakefile),
        cores=threads,
        configfiles=[config_path],
        use_conda=True,
        conda_frontend="mamba",
        printshellcmds=True,
        show_failed_logs=True,
        keep_incomplete=True
    )
    
    return 0 if success else 1

def main():
    """Main entry point for the command-line interface."""
    parser = argparse.ArgumentParser(
        description="NetInfer: Microbiome Network Inference Pipeline"
    )
    
    parser.add_argument(
        "--input",
        required=True,
        help="Input abundance table file (TSV/CSV/BIOM format)"
    )
    
    parser.add_argument(
        "--output",
        required=True,
        help="Output directory for results"
    )
    
    parser.add_argument(
        "--taxonomy",
        help="Taxonomy mapping file (optional)"
    )
    
    parser.add_argument(
        "--metadata",
        help="Sample metadata file (optional)"
    )
    
    parser.add_argument(
        "--methods",
        help="Comma-separated list of methods to use (default: all). "
             "Available methods: flashweave, flashweaveHE, fastspar, "
             "spearman, spieceasi, propr, jaccard"
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
    
    args = parser.parse_args()
    logger = setup_logger()
    
    try:
        # Parse methods if specified
        methods = None
        if args.methods:
            methods = [m.strip().lower() for m in args.methods.split(",")]
        
        # Create configuration
        logger.info("Creating pipeline configuration...")
        config_path = create_config(
            input_file=args.input,
            output_dir=args.output,
            taxonomy_file=args.taxonomy,
            metadata_file=args.metadata,
            methods=methods,
            no_visual=args.no_visual
        )
        
        # Run pipeline
        logger.info("Starting pipeline execution...")
        result = run_pipeline(config_path, args.threads)
        
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