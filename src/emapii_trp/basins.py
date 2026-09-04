"""chi2 basin assignment and periodic-angle helpers.

Sign convention: GROMACS ``gmx angle -type dihedral`` == MDAnalysis ``calc_dihedrals`` ==
IUPAC right-handed. S1 (flip-in) ~ +100 deg, S2 ~ -60 deg, S3 (flip-out) ~ -120 deg.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

Window = Sequence[float]


def in_window(chi: np.ndarray, window: Window) -> np.ndarray:
    """Boolean mask for lo <= chi <= hi (degrees, no wrapping across +/-180)."""
    lo, hi = window
    return (chi >= lo) & (chi <= hi)


def assign_states(chi: np.ndarray, windows: Mapping[str, Window]) -> np.ndarray:
    """Label each frame with the first matching window name, else 'other'."""
    out = np.full(len(chi), "other", dtype=object)
    for name, win in windows.items():
        m = in_window(chi, win) & (out == "other")
        out[m] = name
    return out.astype(str)


def wrap_deg(x: np.ndarray) -> np.ndarray:
    """Wrap an angle difference into (-180, 180]."""
    return ((np.asarray(x) + 180.0) % 360.0) - 180.0


def dihedral_deg(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """Signed dihedral (degrees, IUPAC right-handed) for four points.

    Implementation note: ``m = cross(b2_hat, n1)`` (not ``cross(n1, b2_hat)``) --
    the other order silently inverts the sign. Checked against ``gmx angle``.
    """
    b0, b1, b2 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
    b1_hat = b1 / np.linalg.norm(b1)
    m = np.cross(b1_hat, n1)
    return float(np.degrees(np.arctan2(np.dot(m, n2), np.dot(n1, n2))))
