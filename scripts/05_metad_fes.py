#!/usr/bin/env python3
"""Analysis 5 -- WTMetaD free-energy surfaces F(chi1, chi2) at six temperatures.

Reads the 200x200 grids reduced from the PLUMED HILLS files (data/metad/fes_grid_{T}.npz),
integrates the S1 and S3 basins (chi1 = trans gate, |chi1| > 120 deg), extracts the 1-D
trans-gated profile F(chi2) and the S1->S3 barrier, and re-plots the block-wise dG(t)
stored in data/metad/fes_basins_6T.json. The honest result: the six single-walker 1-us
runs give a NON-monotonic dG(T) -- they are not converged and no T_c is derived from them.

Outputs: results/metad_dG.csv, results/metad_summary.json,
         figures/metad_fes_6T.png, figures/metad_block_convergence.png
"""
from __future__ import annotations

import json
from pathlib import Path

import _common  # noqa: F401
import numpy as np
import pandas as pd
import typer

from emapii_trp import FIGURES_DIR, RESULTS_DIR, R_KJ, celsius_to_kelvin
from emapii_trp.config import load_config
from emapii_trp.io import load_fes_basins, load_fes_grid
from emapii_trp.plotting import COLORS, T_COLORS, set_style
from emapii_trp.stats import vant_hoff_fit

app = typer.Typer(add_completion=False)


def basin_integrals(fes: np.ndarray, c1: np.ndarray, c2: np.ndarray, T_K: float, S1: tuple, S3: tuple, trans_cut: float) -> dict:
    kT = R_KJ * T_K
    X, Y = np.meshgrid(c1, c2, indexing="ij")
    trans = np.abs(X) > trans_cut
    m1 = trans & (Y >= S1[0]) & (Y <= S1[1]); m3 = trans & (Y >= S3[0]) & (Y <= S3[1])
    P = np.exp(-fes / kT)
    prof = -kT * np.log(np.exp(-fes[trans] .reshape(-1, fes.shape[1]) / kT).sum(axis=0))  # trans-gated F(chi2)
    prof -= prof.min()
    i1 = np.argmin(np.where((c2 >= S1[0]) & (c2 <= S1[1]), prof, np.inf)); i3 = np.argmin(np.where((c2 >= S3[0]) & (c2 <= S3[1]), prof, np.inf))
    lo, hi = sorted((i1, i3))
    barrier_inner = prof[lo:hi + 1].max() - prof[i1]           # path through chi2 ~ 0 / -60 (S2 side)
    barrier_wrap = max(prof[hi:].max(), prof[:lo + 1].max()) - prof[i1]  # path through +/-180
    return {"dG_boltzmann_kJmol": float(-kT * np.log(P[m3].sum() / P[m1].sum())),
            "dG_min_kJmol": float(fes[m3].min() - fes[m1].min()),
            "F_S1_min_kJmol": float(fes[m1].min()), "F_S3_min_kJmol": float(fes[m3].min()),
            "chi2_S1_min_deg": float(c2[i1]), "chi2_S3_min_deg": float(c2[i3]),
            "barrier_S1_to_S3_via_S2_kJmol": float(barrier_inner), "barrier_S1_to_S3_via_180_kJmol": float(barrier_wrap),
            "profile_chi2": prof}


@app.command()
def main(config: Path = typer.Option(Path("config.yaml"))) -> None:
    """Integrate basins, extract barriers, plot FES panels and block convergence."""
    cfg = load_config(config); mcfg = cfg["metad"]
    S1 = tuple(cfg["basins"]["published"]["S1"]); S3 = tuple(cfg["basins"]["published"]["S3"])
    rows, grids, profiles = [], {}, {}
    for T in cfg["temperatures_C"]:
        g = load_fes_grid(T)
        fes = g["fes_kJmol"].astype(float); c1 = g["chi1_deg"]; c2 = g["chi2_deg"]
        r = basin_integrals(fes, c1, c2, celsius_to_kelvin(T), S1, S3, mcfg["chi1_trans_cut_deg"])
        profiles[T] = r.pop("profile_chi2"); grids[T] = (fes, c1, c2)
        rows.append({"T_C": T, "n_hills": int(g["n_hills"]), "t_end_ns": round(float(g["t_end_ns"]), 1), **{k: round(v, 2) for k, v in r.items()}})
    tab = pd.DataFrame(rows)
    RESULTS_DIR.mkdir(exist_ok=True)
    tab.to_csv(RESULTS_DIR / "metad_dG.csv", index=False)

    # block-wise dG(t) from the thesis-era JSON (same estimator), for the convergence figure
    basins = load_fes_basins("fes_basins_6T.json")
    blocks = {int(k): v["blocks"] for k, v in basins.items() if k.isdigit()}
    last200 = {}
    for T, bl in blocks.items():
        d = [b["dG_boltzmann"] for b in bl if b["t_end_ns"] >= 600]
        last200[T] = round(float(max(d) - min(d)), 2)
    vh = vant_hoff_fit(celsius_to_kelvin(tab.T_C.to_numpy()), -tab.dG_boltzmann_kJmol.to_numpy() / (R_KJ * celsius_to_kelvin(tab.T_C.to_numpy())))
    sign_changes = int(np.sum(np.diff(np.sign(tab.dG_boltzmann_kJmol.to_numpy())) != 0))
    summary = {
        "dG_boltzmann_kJmol_per_T": dict(zip(tab.T_C.astype(int).tolist(), tab.dG_boltzmann_kJmol.tolist())),
        "dG_sign_changes_across_T": sign_changes,
        "max_abs_step_between_neighbouring_T_kJmol": round(float(np.abs(np.diff(tab.dG_boltzmann_kJmol)).max()), 2),
        "spread_dG_over_last_400ns_blocks_kJmol": last200,
        "barrier_S1_to_S3_via_S2_kJmol_per_T": dict(zip(tab.T_C.astype(int).tolist(), tab.barrier_S1_to_S3_via_S2_kJmol.tolist())),
        "barrier_S1_to_S3_via_180_kJmol_per_T": dict(zip(tab.T_C.astype(int).tolist(), tab.barrier_S1_to_S3_via_180_kJmol.tolist())),
        "vant_hoff_on_metad_dG": {k: (round(v, 2) if isinstance(v, float) else v) for k, v in vh.items()},
        "verdict": "Six single-walker 1-us WTMetaD runs (sigma 0.35 rad, H 1.2 kJ/mol, gamma 10, PACE 500) are not "
                   "converged: dG(S3-S1) is non-monotonic in T and jumps by >10 kJ/mol between neighbouring "
                   "temperatures; each run remains trapped in the S1- or S3-dominated branch it fell into during "
                   "the first ~100-200 ns. No transition temperature is derived from them. Multiple-walker or OPES "
                   "sampling is the pending fix.",
    }
    (RESULTS_DIR / "metad_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    typer.echo(tab[["T_C", "n_hills", "t_end_ns", "dG_boltzmann_kJmol", "dG_min_kJmol", "barrier_S1_to_S3_via_S2_kJmol", "barrier_S1_to_S3_via_180_kJmol"]].to_string(index=False))
    typer.echo(f"Van't Hoff on metad dG: T_c = {vh['Tc_C']:.0f} C, R2 = {vh['R2']:.2f}  (meaningless -- unconverged)")

    set_style()
    import matplotlib.pyplot as plt
    vmax = mcfg["fes_vmax_kJmol"]
    fig, axes = plt.subplots(2, 3, figsize=(11, 7), constrained_layout=True)
    for ax, T in zip(axes.ravel(), cfg["temperatures_C"]):
        fes, c1, c2 = grids[T]
        edges = np.linspace(-180, 180, len(c1) + 1)
        pcm = ax.pcolormesh(edges, edges, fes.T, cmap="viridis", vmin=0, vmax=vmax, shading="flat")
        cs = ax.contour(c1, c2, fes.T, levels=np.arange(0, vmax + 1, 2.5), colors="white", linewidths=0.4, alpha=0.8)
        ax.clabel(cs, inline=True, fontsize=5.5, fmt="%.0f")
        ax.axhspan(*S1, color=COLORS["S1"], alpha=0.12, lw=0); ax.axhspan(*S3, color=COLORS["S3"], alpha=0.12, lw=0)
        for x0 in (-180, 120):
            ax.axvspan(x0, x0 + 60, color="red", alpha=0.05, lw=0)
        ax.text(150, 100, "S1", color="white", fontweight="bold", ha="center"); ax.text(150, -125, "S3", color="white", fontweight="bold", ha="center")
        r = tab[tab.T_C == T].iloc[0]
        ax.set_title(f"T = {T} C    dG(S3-S1) = {r.dG_boltzmann_kJmol:+.1f} kJ/mol")
        ax.set_xlabel("chi1 (deg)"); ax.set_ylabel("chi2 (deg)")
        ax.set_xticks(range(-180, 181, 60)); ax.set_yticks(range(-180, 181, 60))
    fig.colorbar(pcm, ax=axes.ravel().tolist(), shrink=0.8, label="F (kJ/mol)")
    fig.suptitle("WTMetaD F(chi1, chi2) of Trp128, 1 us single walker per T -- dG(T) is non-monotonic: NOT converged\n(red bands: chi1 trans gate |chi1| > 120 deg used for dG; shaded rows: S1 / S3 chi2 windows)", fontsize=10)
    FIGURES_DIR.mkdir(exist_ok=True)
    fig.savefig(FIGURES_DIR / "metad_fes_6T.png")

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9))
    ax = axes[0]
    for T, bl in sorted(blocks.items()):
        ax.plot([b["t_end_ns"] for b in bl], [b["dG_boltzmann"] for b in bl], "o-", color=T_COLORS[T], ms=4, label=f"{T} C")
    ax.axhline(0, color="k", lw=0.6); ax.set_xlabel("hills deposited up to t (ns)"); ax.set_ylabel("dG(S3 - S1) (kJ/mol), chi1 = trans")
    ax.set_title("(a) block-wise dG(t): runs settle into opposite branches"); ax.legend(ncol=2, fontsize=7.5)
    ax = axes[1]
    for T in cfg["temperatures_C"]:
        ax.plot(grids[T][2], profiles[T], color=T_COLORS[T], lw=1.4, label=f"{T} C")
    ax.axvspan(*S1, color=COLORS["S1"], alpha=0.1, lw=0); ax.axvspan(*S3, color=COLORS["S3"], alpha=0.1, lw=0)
    ax.set_xlabel("chi2 (deg)"); ax.set_ylabel("F(chi2) (kJ/mol), chi1 = trans"); ax.set_xlim(-180, 180); ax.set_ylim(0, 30)
    ax.set_title("(b) trans-gated 1-D profiles"); ax.legend(ncol=2, fontsize=7.5)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "metad_block_convergence.png")
    typer.echo(f"wrote {RESULTS_DIR / 'metad_dG.csv'}, figures/metad_fes_6T.png, figures/metad_block_convergence.png")


if __name__ == "__main__":
    app()
