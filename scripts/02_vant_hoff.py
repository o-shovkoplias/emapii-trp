#!/usr/bin/env python3
"""Analysis 2 -- Van't Hoff regression of ln[P(S3)/P(S1)] vs 1/T on PBC-clean populations.

The midpoint T_c = dH/dS is reported for several basin definitions and with/without the
non-monotonic 50 C run, together with R^2 and a block-bootstrap CI. Temperatures with
zero frames in either state cannot enter ln K and are listed as excluded.

Outputs: results/vant_hoff.csv, results/vant_hoff_points.csv, figures/vant_hoff.png
"""
from __future__ import annotations

from pathlib import Path

import _common  # noqa: F401
import numpy as np
import pandas as pd
import typer

from emapii_trp import FIGURES_DIR, RESULTS_DIR, R_KJ, celsius_to_kelvin
from emapii_trp.basins import in_window
from emapii_trp.config import load_config
from emapii_trp.io import load_clean_perframe
from emapii_trp.plotting import COLORS, T_COLORS, set_style
from emapii_trp.stats import vant_hoff_fit

app = typer.Typer(add_completion=False)


def _counts(df: pd.DataFrame, S1: tuple, S3: tuple) -> pd.DataFrame:
    rows = []
    for T, g in df.groupby("T"):
        chi = g["chi"].to_numpy()
        rows.append({"T_C": int(T), "T_K": celsius_to_kelvin(T), "n": len(g),
                     "n_S1": int(in_window(chi, S1).sum()), "n_S3": int(in_window(chi, S3).sum())})
    return pd.DataFrame(rows)


def _boot_Tc(df: pd.DataFrame, S1: tuple, S3: tuple, Ts: list[int], *, block: int, n_boot: int,
             rng: np.random.Generator) -> tuple[float, float]:
    """Moving-block bootstrap of the whole fit (resample frames within each T)."""
    series = {T: g["chi"].to_numpy() for T, g in df.groupby("T") if T in Ts}
    tcs = []
    for _ in range(n_boot):
        TK, lnK = [], []
        for T, chi in series.items():
            n = len(chi); nb = int(np.ceil(n / block))
            starts = rng.integers(0, n - block + 1, nb)
            idx = (starts[:, None] + np.arange(block)[None, :]).ravel()[:n]
            c = chi[idx]
            a, b = in_window(c, S1).sum(), in_window(c, S3).sum()
            if a > 0 and b > 0:
                TK.append(celsius_to_kelvin(T)); lnK.append(np.log(b / a))
        if len(TK) >= 3:
            tcs.append(vant_hoff_fit(np.array(TK), np.array(lnK))["Tc_C"])
    tcs = np.array(tcs)
    tcs = tcs[np.isfinite(tcs)]
    return float(np.quantile(tcs, 0.025)), float(np.quantile(tcs, 0.975))


@app.command()
def main(config: Path = typer.Option(Path("config.yaml"))) -> None:
    """Fit ln K vs 1/T for each basin definition x drop-T variant."""
    cfg = load_config(config)
    df = load_clean_perframe()
    vcfg = cfg["vant_hoff"]; bcfg = cfg["bootstrap"]
    rng = np.random.default_rng(bcfg["seed"])
    fits, points = [], []
    for d in vcfg["definitions"]:
        S1, S3 = tuple(d["S1"]), tuple(d["S3"])
        cnt = _counts(df, S1, S3)
        cnt["window_def"] = d["name"]
        cnt["lnK"] = np.where((cnt.n_S1 > 0) & (cnt.n_S3 > 0), np.log(cnt.n_S3.clip(lower=1) / cnt.n_S1.clip(lower=1)), np.nan)
        cnt["dG_S3_S1_kJmol"] = -R_KJ * cnt["T_K"] * cnt["lnK"]
        points.append(cnt)
        for drop in vcfg["drop_T_variants"]:
            use = cnt[(~cnt.T_C.isin(drop)) & cnt.lnK.notna()]
            excluded = cnt[cnt.lnK.isna() & ~cnt.T_C.isin(drop)].T_C.tolist()
            if len(use) < 3:
                continue
            fit = vant_hoff_fit(use.T_K.to_numpy(), use.lnK.to_numpy())
            lo, hi = _boot_Tc(df, S1, S3, use.T_C.tolist(), block=bcfg["block_length_frames"],
                              n_boot=min(bcfg["n_boot"], 1000), rng=rng)
            fits.append({"window_def": d["name"], "S1_window": list(S1), "S3_window": list(S3),
                         "dropped_T": drop, "excluded_T_zero_count": excluded, "n_points": fit["n"],
                         "Tc_C": round(fit["Tc_C"], 1), "Tc_boot95_lo": round(lo, 1), "Tc_boot95_hi": round(hi, 1),
                         "R2": round(fit["R2"], 2), "dH_kJmol": round(fit["dH_kJmol"], 1),
                         "dS_JmolK": round(fit["dS_JmolK"], 1), "T_used": use.T_C.tolist()})
    fits_df = pd.DataFrame(fits)
    pts_df = pd.concat(points, ignore_index=True)
    RESULTS_DIR.mkdir(exist_ok=True)
    fits_df.to_csv(RESULTS_DIR / "vant_hoff.csv", index=False)
    pts_df.to_csv(RESULTS_DIR / "vant_hoff_points.csv", index=False)
    typer.echo(fits_df[["window_def", "dropped_T", "excluded_T_zero_count", "n_points", "Tc_C", "Tc_boot95_lo", "Tc_boot95_hi", "R2", "dH_kJmol", "dS_JmolK"]].to_string(index=False))
    typer.echo(f"T_c range over definitions: {fits_df.Tc_C.min():.1f} .. {fits_df.Tc_C.max():.1f} C; "
               f"R^2 range {fits_df.R2.min():.2f} .. {fits_df.R2.max():.2f}")

    # ---------------- figure ----------------
    set_style()
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))
    ax = axes[0]
    styles = {"published_windows": ("o", COLORS["S1"]), "sign_cut": ("s", COLORS["S3"]), "narrow_windows": ("^", COLORS["accent"])}
    for d in vcfg["definitions"]:
        name = d["name"]; mk, col = styles.get(name, ("o", "k"))
        p = pts_df[(pts_df.window_def == name) & pts_df.lnK.notna()]
        ax.plot(1000 / p.T_K, p.lnK, mk, color=col, ms=6, mfc="white", mew=1.4, label=f"{name}")
        for _, r in p.iterrows():
            ax.annotate(f"{r.T_C:.0f}", (1000 / r.T_K, r.lnK), textcoords="offset points", xytext=(4, 3), fontsize=7)
        f = fits_df[(fits_df.window_def == name) & (fits_df.dropped_T.apply(len) == 0)]
        if len(f):
            use = p; x = np.linspace((1000 / use.T_K).min() - 0.02, (1000 / use.T_K).max() + 0.02, 50)
            fit = vant_hoff_fit(use.T_K.to_numpy(), use.lnK.to_numpy())
            ax.plot(x, fit["slope"] * (x / 1000) + fit["intercept"], "-", color=col, lw=1, alpha=0.7)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("1000 / T  (1/K)"); ax.set_ylabel("ln K = ln[P(S3)/P(S1)]")
    ax.set_title("Van't Hoff plot (all usable T; 25 C and 30 C have P(S3)=0)")
    ax.legend(fontsize=7.5)
    ax = axes[1]
    lab = [f"{r.window_def}, {'all T' if not r.dropped_T else 'drop ' + ','.join(map(str, r.dropped_T)) + ' C'}\nR2 = {r.R2:.2f}, n = {r.n_points}" for r in fits_df.itertuples()]
    y = np.arange(len(fits_df))
    ax.errorbar(fits_df.Tc_C, y, xerr=[np.clip(fits_df.Tc_C - fits_df.Tc_boot95_lo, 0, None), np.clip(fits_df.Tc_boot95_hi - fits_df.Tc_C, 0, None)],
                fmt="o", color="k", capsize=3)
    xlo, xhi = 20.0, 62.0
    for yi, r in zip(y, fits_df.itertuples()):
        if (r.Tc_boot95_lo < xlo) or (r.Tc_boot95_hi > xhi):
            ax.text(xhi - 0.5, yi + 0.22, "bootstrap CI extends beyond the axis", ha="right", fontsize=7, color="#777")
    ax.axvspan(37, 43, color=COLORS["exp"], alpha=0.08, lw=0)
    ax.text(40, -0.9, "experimental midpoint estimates 37-43 C\n(Kordysh 2005, Lozhko 2026 HSQC, Malyna)", ha="center", fontsize=7, color=COLORS["exp"])
    ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=7.5)
    ax.set_xlabel("T_c = dH/dS  (C), block-bootstrap 95% CI")
    ax.set_title("T_c is definition-dependent")
    ax.set_xlim(xlo, xhi); ax.set_ylim(-1.6, len(fits_df) - 0.4)
    fig.tight_layout()
    FIGURES_DIR.mkdir(exist_ok=True)
    fig.savefig(FIGURES_DIR / "vant_hoff.png")
    typer.echo(f"wrote {RESULTS_DIR / 'vant_hoff.csv'}, {FIGURES_DIR / 'vant_hoff.png'}")


if __name__ == "__main__":
    app()
