"""
NetInfer: A Snakemake pipeline for microbiome network inference
"""

from importlib.metadata import PackageNotFoundError, version

try:
    # Single source of truth: the version field in pyproject.toml
    __version__ = version("netinfer")
except PackageNotFoundError:  # running from source without installing
    __version__ = "unknown"
