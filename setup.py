from setuptools import setup, find_packages

setup(
    name="netinfer",
    version="0.1.0",
    description="A Snakemake pipeline for microbiome network inference",
    author="Heyu Lin",
    author_email="heyu.lin@qut.edu.au",
    url="https://github.com/SilentGene/NetInfer",
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "netinfer=netinfer.cli:main",
        ],
    },
    install_requires=[
        "snakemake>=6.0.0",
        "pyyaml",
        "pandas",
        "numpy",
    ],
    python_requires=">=3.8",
    package_data={
        "netinfer": [
            "workflow/*",
            "workflow/rules/*",
            "workflow/scripts/*",
            "workflow/envs/*",
            "config/*",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)