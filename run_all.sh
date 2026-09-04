#!/usr/bin/env bash
# Regenerate every table in results/ and every figure in figures/ from data/.
# Usage:  ./run_all.sh            (uses `python` on PATH; override with PYTHON=/path/to/python)
# Runtime: about 1-2 minutes on a laptop CPU. No GPU, no trajectories, no MD engine needed.
set -euo pipefail
cd "$(dirname "$0")"
PYTHON="${PYTHON:-python}"
t0=$(date +%s)
step() { printf '\n==> [%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

step "environment: $($PYTHON --version 2>&1) at $(command -v "$PYTHON")"
$PYTHON -c "import numpy, scipy, pandas, matplotlib, sklearn, yaml, typer" || { echo "missing dependency -- see requirements.txt"; exit 1; }
mkdir -p results figures

step "0/7 sanity tests (dihedral sign convention, basin assignment, data shape)"
$PYTHON -m pytest -q tests

step "1/7 chi2 populations per temperature + block-bootstrap CIs -> results/populations.csv, figures/populations_vs_T.png"
$PYTHON scripts/01_populations.py

step "2/7 Van't Hoff on PBC-clean populations (several basin definitions) -> results/vant_hoff.csv, figures/vant_hoff.png"
$PYTHON scripts/02_vant_hoff.py

step "3/7 chi2 vs indole SASA per T -> results/sasa_by_state.csv, figures/chi2_vs_sasa.png"
$PYTHON scripts/03_chi2_sasa.py

step "4/7 GMM of the indole micro-environment -> results/gmm_*.csv, figures/gmm_pockets.png"
$PYTHON scripts/04_gmm_pockets.py

step "5/7 WTMetaD FES (6 T) + block convergence -> results/metad_dG.csv, figures/metad_fes_6T.png"
$PYTHON scripts/05_metad_fes.py

step "6/7 PBC / NoJump demonstration -> results/pbc_basin_changes.csv, figures/pbc_demo.png"
$PYTHON scripts/06_pbc_demo.py

step "7/7 composite README figure -> figures/hero.png"
$PYTHON scripts/07_hero_figure.py

# --- pending / guarded steps (not run here) ------------------------------------------------
# * Re-extraction of data/md_clean/clean_perframe.npz from the trajectories requires the
#   6 x 1 us xtc files (not shipped; see scripts/fetch_trajectories.md) and MDAnalysis:
#     python scripts/prepare_data.py --source-root <thesis tree> --pdb <8ong.pdb>
# * ORCA TD-DFT and excited-state MD branches of the thesis are NOT part of this repository.
# * Multiple-walker WTMetaD / OPES re-sampling (the fix for the unconverged FES) is pending.

step "done in $(( $(date +%s) - t0 )) s; data/ size: $(du -sh data | cut -f1); results: $(ls results | wc -l) files; figures: $(ls figures | wc -l) files"
