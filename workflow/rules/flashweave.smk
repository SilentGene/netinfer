"""
FlashWeave network inference rules
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

rule flashweave_network:
    input:
        abundance = f"{outdir}/preprocessed/filtered_abundance.tsv"
    output:
        network = f"{outdir}/networks/flashweave/normal/network.tsv",
        graph = f"{outdir}/networks/flashweave/normal/network.gml"
    params:
        pvalue = config["flashweave"]["pvalue_threshold"],
        weight = config["flashweave"]["weight_threshold"],
        script = "scripts/run_flashweave.jl"
    threads: 1
    log:
        f"{outdir}/logs/flashweave_normal.log"
    shell:
        """
        julia {params.script} \
            --data {input.abundance} \
            --transposed \
            --alpha {params.pvalue} \
            --outtable {output.network} \
            --outgraph {output.graph} \
            --no_heterogeneous \
            2> {log}

        # Post-filter by weight threshold to keep associations with sufficient weight
        python - << 'PY'
import pandas as pd
import numpy as np
path = r"{output.network}"
thr = float({params.weight})
df = pd.read_csv(path, sep='\t')
# Guard: if file unexpectedly has fewer than 3 columns, keep as-is
if df.shape[1] >= 3:
    weight_col = df.columns[2]
    # Coerce to numeric, drop rows below threshold
    df[weight_col] = pd.to_numeric(df[weight_col], errors='coerce')
    df = df[df[weight_col] >= thr]
df.to_csv(path, sep='\t', index=False)
PY
        """