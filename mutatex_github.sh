#!/bin/bash

# === Get the script's directory for absolute paths ===
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# === User settings (hard-coded) ===
PDBFILE="$SCRIPT_DIR/fold_recruitment_of_smelgif_to_swi_model_0.pdb"
OUTDIR="$SCRIPT_DIR/mutatex_results"
THREADS=32
NRUNS=5

FOLDX_BIN="$SCRIPT_DIR/modules/foldx_20270131_5.1"
ROTABASE="$SCRIPT_DIR/modules/rotabase.txt"
INTERFACE_TEMPLATE="$SCRIPT_DIR/modules/interface_runfile_template.txt"
REPAIR_TEMPLATE="$SCRIPT_DIR/modules/repair_runfile_template.txt"
MUTATE_TEMPLATE="$SCRIPT_DIR/modules/mutate_runfile_template.txt"

# === Create output folder ===
mkdir -p "$OUTDIR"
cd "$OUTDIR" || exit

# === Run MutateX ===
mutatex \
  -p "$THREADS" \
  -n "$NRUNS" \
  -B \
  -R "$REPAIR_TEMPLATE" \
  -M "$MUTATE_TEMPLATE" \
  -I "$INTERFACE_TEMPLATE" \
  -x "$FOLDX_BIN" \
  -b "$ROTABASE" \
  "$PDBFILE"

# === Return to original folder ===
cd - > /dev/null

# === Generate heatmaps using Python ===
python3 << EOF
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

outdir = "$OUTDIR"
pdbname = os.path.basename("$PDBFILE").replace(".pdb","")
results_dir = os.path.join(outdir, pdbname, "results")

# Folding ΔΔG heatmap
folding_csv = os.path.join(results_dir, "ddg_global.csv")
if os.path.exists(folding_csv):
    df = pd.read_csv(folding_csv)
    pivot = df.pivot(index="mutation", columns="position", values="ddg")
    plt.figure(figsize=(18,6))
    sns.heatmap(pivot, cmap="coolwarm", center=0)
    plt.title(f"{pdbname} Folding ΔΔG Heatmap")
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, f"{pdbname}_folding_heatmap.png"))
    plt.close()

# Binding ΔΔG heatmap
binding_csv = os.path.join(results_dir, "ddg_binding.csv")
if os.path.exists(binding_csv):
    df = pd.read_csv(binding_csv)
    pivot = df.pivot(index="mutation", columns="position", values="ddg")
    plt.figure(figsize=(18,6))
    sns.heatmap(pivot, cmap="coolwarm", center=0)
    plt.title(f"{pdbname} Binding ΔΔG Heatmap")
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, f"{pdbname}_binding_heatmap.png"))
    plt.close()
EOF

echo "MutateX run and heatmaps completed. Results in $OUTDIR/$pdbname/results/"
