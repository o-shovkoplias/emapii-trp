#!/usr/bin/env python3
"""Composite README figure: populations + Kordysh | FES at 25 C and 45 C | GMM pockets.

Reads results/ written by scripts 01, 04 and data/metad grids. Run after 01-05.
Output: figures/hero.png
"""
from __future__ import annotations

import json
from pathlib import Path

import _common  # noqa: F401
import numpy as np
import pandas as pd
import typer

from emapii_trp import FIGURES_DIR, RESULTS_DIR
from emapii_trp.config import load_config
from emapii_trp.io import load_fes_grid, load_kordysh
from emapii_trp.plotting import COLORS, set_style

app = typer.Typer(add_completion=False)
CLUSTER_COLORS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a"]


@app.command()
def main(config: Path = typer.Option(Path("config.yaml"))) -> None:
    """Assemble the 1x3 hero figure from existing results."""
    cfg = load_config(config)
    pop = pd.read_csv(RESULTS_DIR / "populations.csv"); pop = pop[pop.window_set == "published"]
    summ = json.loads((RESULTS_DIR / "populations_summary.json").read_text())
    labels = pd.read_csv(RESULTS_DIR / "gmm_labels.csv")
    cent = pd.read_csv(RESULTS_DIR / "gmm_centroids.csv", index_col=0)
    from emapii_trp.io import load_clean_perframe
    df = load_clean_perframe().merge(labels[["T", "frame", "cluster_name"]], on=["T", "frame"])
    kord = load_kordysh()
    set_style()
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(14, 4.3))
    gs = GridSpec(1, 4, figure=fig, width_ratios=[1.25, 0.9, 0.9, 1.15], wspace=0.55)

    ax = fig.add_subplot(gs[0])
    Ts = pop.T_C.to_numpy()
    for s, lab in (("S1", "S1 flip-in"), ("S3", "S3 flip-out")):
        y = pop[f"P_{s}"].to_numpy(); lo = pop[f"P_{s}_lo"].to_numpy(); hi = pop[f"P_{s}_hi"].to_numpy()
        ax.errorbar(Ts, y, yerr=[np.clip(y - lo, 0, None), np.clip(hi - y, 0, None)], fmt="o-", color=COLORS[s], ms=5, lw=1.6, capsize=3, label=lab)
    cross = summ.get("crossover_T_C_linear_interp")
    if cross:
        ax.axvline(cross, color="k", ls=":", lw=0.9); ax.text(cross - 0.4, 0.80, f"P(S3)=P(S1)\n~{cross:.0f} C", fontsize=7.5, ha="right")
    ax.set_xlabel("temperature (C)"); ax.set_ylabel("Trp128 chi2-state population (6 x 1 us MD)"); ax.set_ylim(-0.02, 1.05); ax.set_xticks(Ts)
    ax2 = ax.twinx(); ax2.spines["right"].set_visible(True)
    meas = kord[kord.kind == "measured"]
    ax2.plot(kord.T_C, kord.lambda_max_nm, "--", color=COLORS["exp"], lw=0.9, alpha=0.6); ax2.plot(meas.T_C, meas.lambda_max_nm, "s", color=COLORS["exp"], ms=6, label="Kordysh 2005 lambda_max")
    ax2.set_ylabel("emission max (nm), exp.", labelpad=2); ax2.set_ylim(332, 352)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", bbox_to_anchor=(0.0, -0.17), ncol=3, fontsize=7.5)
    ax.set_title("(a) rotamer populations track the 14 nm red shift")

    S1w = cfg["basins"]["published"]["S1"]; S3w = cfg["basins"]["published"]["S3"]
    for i, T in enumerate((25, 45)):
        ax = fig.add_subplot(gs[1 + i])
        g = load_fes_grid(T); fes = g["fes_kJmol"]; c = g["chi1_deg"]; edges = np.linspace(-180, 180, len(c) + 1)
        pcm = ax.pcolormesh(edges, edges, fes.T, cmap="viridis", vmin=0, vmax=25, shading="flat")
        ax.contour(c, c, fes.T, levels=np.arange(0, 26, 5), colors="white", linewidths=0.35, alpha=0.8)
        ax.axhspan(*S1w, color=COLORS["S1"], alpha=0.15, lw=0); ax.axhspan(*S3w, color=COLORS["S3"], alpha=0.15, lw=0)
        ax.text(150, 100, "S1", color="w", fontweight="bold", ha="center", fontsize=8); ax.text(150, -125, "S3", color="w", fontweight="bold", ha="center", fontsize=8)
        ax.set_xlabel("chi1 (deg)"); ax.set_xticks(range(-180, 181, 90)); ax.set_yticks(range(-180, 181, 90))
        if i == 0:
            ax.set_ylabel("chi2 (deg)")
        ax.set_title(f"({'bc'[i]}) WTMetaD F(chi1,chi2), {T} C\n1 us single walker, not converged", fontsize=9)
    cb = fig.colorbar(pcm, ax=ax, shrink=0.85, pad=0.03); cb.set_label("F (kJ/mol)", fontsize=8)

    ax = fig.add_subplot(gs[3])
    order = list(cent.sort_values("dgln").index)
    for k, col in zip(order, CLUSTER_COLORS):
        name = cent.loc[k, "name"]; g = df[df.cluster_name == name]
        ax.scatter(g.dlys, g.dgln, s=5, alpha=0.35, color=col, label=f"{name} ({100 * len(g) / len(df):.0f}%)", rasterized=True)
    ax.set_xlabel("d(indole -> Lys124-126) (A)"); ax.set_ylabel("d(indole -> Gln132) (A)")
    ax.legend(fontsize=7, markerscale=3, loc="upper left"); ax.set_title("(d) GMM (n=4) of the indole micro-environment")
    fig.savefig(FIGURES_DIR / "hero.png", dpi=170)
    typer.echo(f"wrote {FIGURES_DIR / 'hero.png'}")


if __name__ == "__main__":
    app()
