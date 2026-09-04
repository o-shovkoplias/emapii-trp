"""Shared matplotlib style: one palette, one font scale, readable in print."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COLORS = {
    "S1": "#2c7fb8",     # flip-in, blue
    "S2": "#7a7a7a",     # intermediate, grey
    "S3": "#d95f0e",     # flip-out, orange
    "other": "#bdbdbd",
    "exp": "#222222",    # experiment (Kordysh 2005)
    "accent": "#6a51a3",
}
T_COLORS = {25: "#4575b4", 30: "#91bfdb", 35: "#a6d96a", 40: "#fdae61", 45: "#f46d43", 50: "#d73027"}


def set_style() -> None:
    """Apply the repository plotting style."""
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 200, "savefig.bbox": "tight",
        "font.size": 9.5, "axes.titlesize": 10, "axes.labelsize": 9.5,
        "legend.fontsize": 8, "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False,
    })


def shade_basins(ax, axis: str = "x", windows: dict | None = None) -> None:
    """Shade the S1/S2/S3 chi2 windows along ``axis`` ('x' or 'y')."""
    windows = windows or {"S1": (50, 150), "S2": (-90, -30), "S3": (-160, -90)}
    span = ax.axvspan if axis == "x" else ax.axhspan
    for name, (lo, hi) in windows.items():
        span(lo, hi, color=COLORS[name], alpha=0.10, lw=0)
