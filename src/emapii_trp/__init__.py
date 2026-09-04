"""emapii_trp -- re-analysis package for the Trp128 (Trp125) chi2 switch in EMAP II.

All numbers in the README are produced by the scripts in ``scripts/`` that import this
package. Residue numbering in code follows PDB 8ONG (Trp128); publications use the
mature-protein numbering (Trp125 = Trp128 - 3).
"""
from __future__ import annotations

from pathlib import Path

__version__ = "0.1.0"

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = REPO_ROOT / "figures"

R_KJ = 8.314462618e-3  # gas constant, kJ/(mol K)


def celsius_to_kelvin(t_c: float) -> float:
    """Convert Celsius to Kelvin (IUPAC offset 273.15)."""
    return t_c + 273.15
