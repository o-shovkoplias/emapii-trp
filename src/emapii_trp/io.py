"""Readers for the small data files shipped in ``data/``."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import DATA_DIR

CLEAN_KEYS = ("chi", "dgln", "dlys", "dne", "dmin", "sasa", "KKK", "T")


def load_clean_perframe(path: Path | None = None) -> pd.DataFrame:
    """Load the PBC-clean per-frame table (6006 frames = 6 T x 1001 frames, 1 ns apart).

    Columns
    -------
    T      : temperature, C
    frame  : frame index within its trajectory (= time in ns)
    chi    : Trp128 chi2 (CA-CB-CG-CD1), degrees, GROMACS/IUPAC sign, min-image corrected
    dgln   : indole centroid -> Gln132 side-chain centroid, A
    dlys   : indole centroid -> Lys124-126 CA centroid, A
    dne    : Trp128 NE1 -> Gln132 OE1, A
    dmin   : minimum heavy-atom distance indole -> Gln132 side chain, A
    sasa   : indole SASA (NoJump-unwrapped protein, freesasa), A^2
    KKK    : summed side-chain SASA of Lys124-126, A^2
    """
    npz = np.load(path or DATA_DIR / "md_clean/clean_perframe.npz")
    df = pd.DataFrame({k: npz[k] for k in CLEAN_KEYS})
    df["T"] = df["T"].astype(int)
    df["frame"] = df.groupby("T").cumcount()
    return df[["T", "frame", *[k for k in CLEAN_KEYS if k != "T"]]]


def load_xvg(path: Path, usecols: tuple[int, ...] | None = None) -> np.ndarray:
    """Load a GROMACS/xmgrace .xvg (lines starting with # or @ are comments)."""
    return np.loadtxt(path, comments=["#", "@"], usecols=usecols)


def load_raw_chi(T: int) -> pd.DataFrame:
    """RAW (no min-image) chi1/chi2 at 1 ns, for the PBC demonstration only."""
    arr = load_xvg(DATA_DIR / f"md_raw_pbc_broken/chi_mda_{T}_pbc_broken.xvg")
    return pd.DataFrame({"time_ns": arr[:, 0], "chi1_raw": arr[:, 1], "chi2_raw": arr[:, 2]})


def load_gmx_chi2_100ps(T: int) -> pd.DataFrame | None:
    """gmx angle chi2 at 100 ps (5 temperatures; none exists for 25 C)."""
    p = DATA_DIR / f"md_clean/chi2_gmx_100ps_{T}.xvg"
    if not p.exists():
        return None
    arr = load_xvg(p)
    return pd.DataFrame({"time_ps": arr[:, 0], "chi2": arr[:, 1]})


def load_kordysh(path: Path | None = None) -> pd.DataFrame:
    """Digitised Kordysh 2005 lambda_max(T) table."""
    return pd.read_csv(path or DATA_DIR / "experiment/kordysh2005_lambda_max.csv", comment="#")


def load_fes_grid(T: int) -> dict:
    """FES grid reduced from the PLUMED HILLS file for temperature ``T`` (C)."""
    npz = np.load(DATA_DIR / f"metad/fes_grid_{T}.npz")
    return {k: npz[k] for k in npz.files}


def load_fes_basins(name: str = "fes_basins_6T.json") -> dict:
    """Thesis-era basin summary JSON (kept for cross-checking)."""
    with (DATA_DIR / "metad" / name).open() as fh:
        return json.load(fh)
