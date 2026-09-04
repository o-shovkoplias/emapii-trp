#!/usr/bin/env python3
"""Analysis 6 -- why the minimum-image / NoJump correction matters for a dihedral.

Compares, frame by frame, the RAW chi2 (MDAnalysis dihedral without the box; the file
family that produced the retracted T_c = 41.4 C) with the PBC-clean chi2 and with the
gmx-angle chi2 at the same times. Reports the fraction of frames whose basin label
changes per temperature and shows the 25 C histograms.

Outputs: results/pbc_basin_changes.csv, results/convention_check.json, figures/pbc_demo.png
"""
from __future__ import annotations

import json
from pathlib import Path

import _common  # noqa: F401
import numpy as np
import pandas as pd
import typer

from emapii_trp import FIGURES_DIR, RESULTS_DIR
from emapii_trp.basins import assign_states, wrap_deg
from emapii_trp.config import load_config
from emapii_trp.io import load_clean_perframe, load_gmx_chi2_100ps, load_raw_chi
from emapii_trp.plotting import COLORS, set_style, shade_basins

app = typer.Typer(add_completion=False)


@app.command()
def main(config: Path = typer.Option(Path("config.yaml"))) -> None:
    """Raw vs clean chi2: basin changes, convention check, demonstration figure."""
    cfg = load_config(config); win = cfg["basins"]["published"]
    df = load_clean_perframe()
    rows, conv = [], {}
    raw_all = {}
    for T, g in df.groupby("T"):
        raw = load_raw_chi(T)
        assert len(raw) == len(g)
        c_raw = raw.chi2_raw.to_numpy(); c_clean = g.chi.to_numpy()
        raw_all[T] = c_raw
        lab_raw = assign_states(c_raw, win); lab_clean = assign_states(c_clean, win)
        d = np.abs(wrap_deg(c_raw - c_clean))
        changed = lab_raw != lab_clean
        rows.append({"T_C": int(T), "n": len(g), "frames_differ_gt_5deg": int((d > 5).sum()), "frac_differ_gt_5deg": round(float((d > 5).mean()), 4),
                     "frames_basin_changed": int(changed.sum()), "frac_basin_changed": round(float(changed.mean()), 4),
                     "raw_S3_count": int((lab_raw == "S3").sum()), "clean_S3_count": int((lab_clean == "S3").sum()),
                     "raw_P_S3": round(float((lab_raw == "S3").mean()), 4), "clean_P_S3": round(float((lab_clean == "S3").mean()), 4),
                     "median_abs_shift_on_changed_deg": round(float(np.median(d[changed])), 0) if changed.any() else 0.0})
        gmx = load_gmx_chi2_100ps(T)
        if gmx is not None:
            at_ns = gmx[np.isclose(gmx.time_ps % 1000, 0)].chi2.to_numpy()[: len(c_clean)]
            dd = np.abs(wrap_deg(at_ns - c_clean[: len(at_ns)]))
            conv[int(T)] = {"n_compared": int(len(at_ns)), "frac_within_1deg_clean_vs_gmx": round(float((dd < 1).mean()), 4),
                            "frac_within_1deg_raw_vs_gmx": round(float((np.abs(wrap_deg(at_ns - c_raw[: len(at_ns)])) < 1).mean()), 4)}
    tab = pd.DataFrame(rows)
    RESULTS_DIR.mkdir(exist_ok=True)
    tab.to_csv(RESULTS_DIR / "pbc_basin_changes.csv", index=False)
    pooled = {"frames_basin_changed_total": int(tab.frames_basin_changed.sum()), "frac_basin_changed_pooled": round(float(tab.frames_basin_changed.sum() / tab.n.sum()), 4),
              "raw_P_S3_25C": float(tab.loc[tab.T_C == 25, "raw_P_S3"].iloc[0]), "clean_P_S3_25C": float(tab.loc[tab.T_C == 25, "clean_P_S3"].iloc[0]),
              "raw_P_S3_30C": float(tab.loc[tab.T_C == 30, "raw_P_S3"].iloc[0]), "clean_P_S3_30C": float(tab.loc[tab.T_C == 30, "clean_P_S3"].iloc[0])}
    (RESULTS_DIR / "convention_check.json").write_text(json.dumps({"clean_vs_gmx_angle_at_1ns": conv, "pooled": pooled,
        "sign_convention": "clean chi2 == gmx angle -type dihedral (right-handed IUPAC): S1 ~ +100 deg, S3 ~ -120 deg"}, indent=1) + "\n")
    typer.echo(tab.to_string(index=False)); typer.echo(json.dumps(conv)); typer.echo(json.dumps(pooled))

    set_style()
    import matplotlib.pyplot as plt
    g25 = df[df["T"] == 25]
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.7))
    ax = axes[0]
    bins = np.arange(-180, 181, 6)
    ax.hist(raw_all[25], bins=bins, color="#999999", alpha=0.9, label="RAW dihedral (no min-image)")
    ax.hist(g25.chi, bins=bins, histtype="step", color="k", lw=1.4, label="PBC-clean (box-aware)")
    shade_basins(ax, "x"); ax.set_yscale("log")
    ax.set_xlabel("Trp128 chi2 at 25 C (deg)"); ax.set_ylabel("frames (log)")
    r25 = tab[tab.T_C == 25].iloc[0]
    ax.set_title(f"(a) 25 C: RAW puts {int(r25.raw_S3_count)} frames in S3, clean puts {int(r25.clean_S3_count)}")
    ax.legend(fontsize=7.5, loc="upper center")
    ax = axes[1]
    ax.scatter(raw_all[25], g25.chi, s=6, color="k", alpha=0.5, rasterized=True)
    ax.plot([-180, 180], [-180, 180], color=COLORS["S1"], lw=0.8)
    ax.set_xlabel("RAW chi2 (deg)"); ax.set_ylabel("PBC-clean chi2 (deg)"); ax.set_title("(b) off-diagonal frames = side chain split across the box")
    ax.set_xlim(-180, 180); ax.set_ylim(-180, 180)
    ax = axes[2]
    xpos = np.arange(len(tab))
    ax.bar(xpos, 100 * tab.frac_basin_changed, color=COLORS["S3"], width=0.6)
    for x, v in zip(xpos, 100 * tab.frac_basin_changed.to_numpy()):
        ax.text(float(x), float(v) + 0.2, f"{v:.1f}%", ha="center", fontsize=8)
    ax.set_xticks(xpos); ax.set_xticklabels(tab.T_C.astype(str))
    ax.set_xlabel("temperature (C)"); ax.set_ylabel("% frames whose basin label changes")
    ax.set_title(f"(c) mislabelled frames, pooled {100 * pooled['frac_basin_changed_pooled']:.1f}%")
    fig.tight_layout()
    FIGURES_DIR.mkdir(exist_ok=True)
    fig.savefig(FIGURES_DIR / "pbc_demo.png")
    typer.echo(f"wrote {RESULTS_DIR / 'pbc_basin_changes.csv'}, {FIGURES_DIR / 'pbc_demo.png'}")


if __name__ == "__main__":
    app()
