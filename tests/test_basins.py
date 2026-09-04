"""Minimal sanity tests: dihedral sign convention, basin assignment, data integrity."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from emapii_trp.basins import assign_states, dihedral_deg, wrap_deg  # noqa: E402
from emapii_trp.io import load_clean_perframe  # noqa: E402


def test_dihedral_is_right_handed() -> None:
    p0 = np.array([1.0, 0.0, 0.0]); p1 = np.zeros(3)
    p2 = np.array([0.0, 0.0, 1.0]); p3 = np.array([0.0, 1.0, 1.0])
    assert abs(dihedral_deg(p0, p1, p2, p3) - 90.0) < 1e-6  # wrong cross order gives -90


def test_assign_states() -> None:
    chi = np.array([108.0, -61.0, -122.0, 0.0])
    lab = assign_states(chi, {"S1": (50, 150), "S2": (-90, -30), "S3": (-160, -90)})
    assert list(lab) == ["S1", "S2", "S3", "other"]


def test_wrap() -> None:
    assert wrap_deg(np.array([190.0]))[0] == -170.0


def test_clean_table_shape() -> None:
    df = load_clean_perframe()
    assert len(df) == 6006
    assert sorted(df["T"].unique().tolist()) == [25, 30, 35, 40, 45, 50]
    assert (df.groupby("T").size() == 1001).all()
