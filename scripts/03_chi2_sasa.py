#!/usr/bin/env python3
"""Analysis 3 -- chi2 vs indole SASA per temperature ("buried flip-out" check).

Compares the indole SASA of the S1 and S3 rotamers (i) at the experimental endpoints
(S1 at 25 C vs S3 at 45 C), (ii) pooled over all six temperatures, and (iii) within each
temperature -- the three comparisons do not agree in sign, which is the point.

Outputs: results/sasa_by_state.csv, results/sasa_summary.json, figures/chi2_vs_sasa.png
"""
from __future__ import annotations

import json
from pathlib import Path

import _common  # noqa: F401
import numpy as np
import pandas as pd
import typer

from emapii_trp import FIGURES_DIR, RESULTS_DIR
from emapii_trp.basins import assign_states
from emapii_trp.config import load_config
from emapii_trp.io import load_clean_perframe
from emapii_trp.plotting import COLORS, set_style, shade_basins

app = typer.Typer(add_completion=False)


def _ci_diff(a: np.ndarray, b: np.ndarray, rng: np.random.Generator, n: int = 5000) -> tuple[float, float]:
    d = np.array([b[rng.integers(0, len(b), len(b))].mean() - a[rng.integers(0, len(a), len(a))].mean() for _ in range(n)])
    return float(np.quantile(d, 0.025)), float(np.quantile(d, 0.975))


@app.command()
def main(config: Path = typer.Option(Path("config.yaml"))) -> None:
    """Per-T x state SASA table, endpoint vs pooled vs fixed-T comparison, 2-D density figure."""
    cfg = load_config(config)
    df = load_clean_perframe()
    df["state"] = assign_states(df["chi"].to_numpy(), cfg["basins"]["published"])
    rng = np.random.default_rng(cfg["bootstrap"]["seed"])

    rows = []
    for (T, s), g in df.groupby(["T", "state"]):
        rows.append({"T_C": int(T), "state": s, "n": len(g), "sasa_mean_A2": round(g.sasa.mean(), 2),
                     "sasa_median_A2": round(g.sasa.median(), 2), "sasa_sd_A2": round(g.sasa.std(ddof=1), 2) if len(g) > 1 else np.nan,
                     "frac_sasa_gt_100": round(float((g.sasa > 100).mean()), 4), "sasa_max_A2": round(g.sasa.max(), 2)})
    for s, g in df.groupby("state"):
        rows.append({"T_C": "pooled", "state": s, "n": len(g), "sasa_mean_A2": round(g.sasa.mean(), 2),
                     "sasa_median_A2": round(g.sasa.median(), 2), "sasa_sd_A2": round(g.sasa.std(ddof=1), 2),
                     "frac_sasa_gt_100": round(float((g.sasa > 100).mean()), 4), "sasa_max_A2": round(g.sasa.max(), 2)})
    tab = pd.DataFrame(rows)
    RESULTS_DIR.mkdir(exist_ok=True)
    tab.to_csv(RESULTS_DIR / "sasa_by_state.csv", index=False)

    S1_25 = df[(df["T"] == 25) & (df.state == "S1")].sasa.to_numpy()
    S3_45 = df[(df["T"] == 45) & (df.state == "S3")].sasa.to_numpy()
    S1_all = df[df.state == "S1"].sasa.to_numpy(); S3_all = df[df.state == "S3"].sasa.to_numpy()
    fixedT = {}
    for T in (35, 40, 50):
        a = df[(df["T"] == T) & (df.state == "S1")].sasa.to_numpy(); b = df[(df["T"] == T) & (df.state == "S3")].sasa.to_numpy()
        if len(a) > 20 and len(b) > 20:
            lo, hi = _ci_diff(a, b, rng)
            fixedT[T] = {"S1_mean": round(a.mean(), 2), "S3_mean": round(b.mean(), 2), "delta_S3_minus_S1": round(b.mean() - a.mean(), 2), "ci95": [round(lo, 2), round(hi, 2)], "n_S1": len(a), "n_S3": len(b)}
    lo_e, hi_e = _ci_diff(S1_25, S3_45, rng); lo_p, hi_p = _ci_diff(S1_all, S3_all, rng)
    s45 = df[df["T"] == 45].sasa
    per_T_mean = df.groupby("T").sasa.mean().round(2).to_dict()
    summary = {
        "endpoint_S1_at_25C_mean_A2": round(S1_25.mean(), 2), "endpoint_S3_at_45C_mean_A2": round(S3_45.mean(), 2),
        "endpoint_delta_S3_minus_S1_A2": round(S3_45.mean() - S1_25.mean(), 2), "endpoint_delta_ci95": [round(lo_e, 2), round(hi_e, 2)],
        "pooled_S1_mean_A2": round(S1_all.mean(), 2), "pooled_S3_mean_A2": round(S3_all.mean(), 2),
        "pooled_delta_S3_minus_S1_A2": round(S3_all.mean() - S1_all.mean(), 2), "pooled_delta_ci95": [round(lo_p, 2), round(hi_p, 2)],
        "pooled_n_S1": int(len(S1_all)), "pooled_n_S3": int(len(S3_all)),
        "fixed_T_comparisons": fixedT,
        "T45_frac_frames_sasa_gt_100": round(float((s45 > 100).mean()), 4), "T45_max_sasa_A2": round(float(s45.max()), 2),
        "mean_indole_sasa_per_T_A2": {int(k): v for k, v in per_T_mean.items()},
        "spearman_chi2_sasa_pooled": round(float(df[["chi", "sasa"]].corr(method="spearman").iloc[0, 1]), 3),
        "reading": "The endpoint comparison (S3 at 45 C more buried than S1 at 25 C) confounds rotamer with "
                   "temperature; pooled by rotamer and at fixed T the sign reverses. Neither rotamer reaches "
                   "bulk-like exposure within the 1-us runs -- a statement about the sampled timescale, not a "
                   "physical prohibition.",
    }
    (RESULTS_DIR / "sasa_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    typer.echo(json.dumps({k: v for k, v in summary.items() if k != "reading"}, indent=1))

    set_style()
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.2), sharex=True, sharey=True)
    for ax, (T, g) in zip(axes.ravel(), df.groupby("T")):
        h = ax.hist2d(g.chi, g.sasa, bins=[np.arange(-180, 181, 10), np.arange(0, 121, 5)], cmap="Greys", norm=LogNorm(vmin=1, vmax=60))
        shade_basins(ax, "x")
        for s in ("S1", "S3"):
            gs = g[g.state == s]
            if len(gs) > 5:
                ax.plot(gs.chi.median(), gs.sasa.mean(), "o", color=COLORS[s], ms=7, mec="k")
                ax.annotate(f"{s}: {gs.sasa.mean():.1f} A2 (n={len(gs)})", (gs.chi.median(), gs.sasa.mean()), xytext=(0, 9),
                            textcoords="offset points", ha="center", fontsize=7.5, color=COLORS[s])
        ax.set_title(f"T = {T} C   <SASA> = {g.sasa.mean():.1f} A2")
        ax.axhline(100, color="r", ls="--", lw=0.7)
    for ax in axes[1]:
        ax.set_xlabel("Trp128 chi2 (deg, GROMACS sign)")
    for ax in axes[:, 0]:
        ax.set_ylabel("indole SASA (A2, NoJump)")
    fig.suptitle("chi2 vs indole SASA per temperature. Blue = S1 window, orange = S3 window; red dashed = 100 A2 (never reached at 45 C)", fontsize=9.5)
    fig.tight_layout()
    FIGURES_DIR.mkdir(exist_ok=True)
    fig.savefig(FIGURES_DIR / "chi2_vs_sasa.png")
    typer.echo(f"wrote {RESULTS_DIR / 'sasa_by_state.csv'}, {FIGURES_DIR / 'chi2_vs_sasa.png'}")


if __name__ == "__main__":
    app()
