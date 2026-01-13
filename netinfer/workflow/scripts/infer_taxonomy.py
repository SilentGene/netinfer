#!/usr/bin/env python3
"""
Infer taxonomy from abundance table feature IDs.
"""

import pandas as pd
import logging
import sys

def setup_logger(log_file: str) -> logging.Logger:
    """Set up logging."""
    logger = logging.getLogger('infer_taxonomy')
    logger.setLevel(logging.INFO)
    
    fh = logging.FileHandler(log_file)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)
    
    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(ch)
    
    return logger

def main(snakemake):
    logger = setup_logger(snakemake.log[0])
    abundance_path = snakemake.input.abundance
    output_path = snakemake.output.processed
    
    logger.info(f"Inferring taxonomy from abundance table: {abundance_path}")
    
    try:
        # Load abundance table
        if str(abundance_path).endswith('.biom'):
            try:
                import biom
                table = biom.load_table(abundance_path)
                features = table.ids('observation')
            except ImportError:
                logger.error("biom-format package required for BIOM files")
                raise
        else:
            # Assume TSV/CSV
            sep = ',' if str(abundance_path).endswith('.csv') else '\t'
            try:
                # Read only the index (first column)
                df = pd.read_csv(abundance_path, sep=sep, index_col=0, usecols=[0])
                features = df.index.astype(str)
            except Exception:
                # Fallback
                df = pd.read_csv(abundance_path, sep=sep, index_col=0)
                features = df.index.astype(str)

        taxonomy_data = []
        for feature in features:
            # Logic: find 'd__' or 'k__' as domain info, otherwise find 'p__' as phylum info
            start_idx = feature.find('d__')
            if start_idx == -1:
                start_idx = feature.find('k__')
            if start_idx == -1:
                start_idx = feature.find('p__')
                
            if start_idx != -1:
                tax_str = feature[start_idx:]
            else:
                tax_str = "Unclassified"
                
            taxonomy_data.append({'Feature': feature, 'Taxonomy': tax_str})
            
        tax_df = pd.DataFrame(taxonomy_data)
        
        # Save directly to the processed path
        tax_df.to_csv(output_path, sep='\t', index=False)
        logger.info(f"Inferred taxonomy saved to: {output_path}")
        
    except Exception as e:
        logger.exception(f"Error during taxonomy inference: {e}")
        raise

if __name__ == "__main__":
    main(snakemake)
