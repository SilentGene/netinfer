"""
FlashWeave network inference rules (Heterogeneous mode)
"""

# Get output directory from config
outdir = config.get("output_dir", "results").strip()  # Default to "results" if not specified

rule flashweave_he_network:
    input:
        abundance = f"{outdir}/preprocessed_data/filtered_abundance.tsv"
    output:
        network = f"{outdir}/subtool_outputs/flashweave/HE/network.tsv",
        graph = f"{outdir}/subtool_outputs/flashweave/HE/network.gml"
    params:
        pvalue = config["flashweaveHE"]["pvalue_threshold"],
        weight = config["flashweaveHE"]["weight_threshold"],
        script = f"{workflow.basedir}/scripts/run_flashweave.jl"
    threads: 1
    log:
        f"{outdir}/logs/flashweave_HE.log"
    shell:
        """
        julia {params.script} \
            --data {input.abundance} \
            --transposed \
            --alpha {params.pvalue} \
            --outtable {output.network} \
            --outgraph {output.graph} \
            2> {log}

        # Post-filter by weight threshold (HE mode)
        python - << 'PY'
import pandas as pd
import numpy as np
path = r"{output.network}"
thr = float({params.weight})
df = pd.read_csv(path, sep='\t')
if df.shape[1] >= 3:
    weight_col = df.columns[2]
    df[weight_col] = pd.to_numeric(df[weight_col], errors='coerce')
    df = df[df[weight_col] >= thr]
df.to_csv(path, sep='\t', index=False)
PY
        """