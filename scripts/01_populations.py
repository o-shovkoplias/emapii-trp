#!/usr/bin/env python3
"""Analysis 1 -- chi2 basin populations per temperature with block-bootstrap CIs.

Inputs : data/md_clean/clean_perframe.npz (PBC-clean chi2, 1 ns, 6 x 1001 frames)
         data/experiment/kordysh2005_lambda_max.csv
Outputs: results/populations.csv, results/populations_summary.json,
         figures/populations_vs_T.png
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
from emapii_trp.io import load_clean_perframe, load_kordysh
from emapii_trp.plotting import COLORS, set_style
from emapii_trp.stats import block_bootstrap, integrated_autocorr_time

app = typer.Typer(add_completion=False)


@app.command()
def main(config: Path = typer.Option(Path("config.yaml"), help="Parameter file.")) -> None:
    """Compute P(S1), P(S2), P(S3), P(other) per T for two window sets and plot vs Kordysh."""
    cfg = load_config(config)
    df = load_clean_perframe()
    rng = np.random.default_rng(cfg["bootstrap"]["seed"])
    bcfg = cfg["bootstrap"]
    rows = []
    for win_name, windows in cfg["basins"].items():
        for T, g in df.groupby("T"):
            lab = assign_states(g["chi"].to_numpy(), windows)
            row = {"window_set": win_name, "T_C": int(T), "n_frames": len(g)}
            for s in ("S1", "S2", "S3", "other"):
                ind = (lab == s).astype(float)
                p, lo, hi = block_bootstrap(ind, np.mean, block=bcfg["block_length_frames"],
                                            n_boot=bcfg["n_boot"], rng=rng, ci=bcfg["ci"])
                row[f"P_{s}"] = round(p, 4)
                row[f"P_{s}_lo"] = round(lo, 4)
                row[f"P_{s}_hi"] = round(hi, 4)
            row["tau_S3_frames"] = round(integrated_autocorr_time((lab == "S3").astype(float)), 1)
            row["tau_S1_frames"] = round(integrated_autocorr_time((lab == "S1").astype(float)), 1)
            row["chi2_mean_S1_deg"] = round(float(g["chi"][lab == "S1"].mean()), 1) if (lab == "S1").any() else np.nan
            row["chi2_mean_S3_deg"] = round(float(g["chi"][lab == "S3"].mean()), 1) if (lab == "S3").any() else np.nan
            rows.append(row)
    out = pd.DataFrame(rows)
    RESULTS_DIR.mkdir(exist_ok=True)
    out.to_csv(RESULTS_DIR / "populations.csv", index=False)

    pub = out[out.window_set == "published"].set_index("T_C")
    # crossover: first T where P_S3 >= P_S1 (linear interpolation between grid points)
    Ts = pub.index.to_numpy(float)
    diff = (pub["P_S3"] - pub["P_S1"]).to_numpy()
    cross = np.nan
    for i in range(len(Ts) - 1):
        if diff[i] < 0 <= diff[i + 1]:
            cross = Ts[i] + (Ts[i + 1] - Ts[i]) * (-diff[i]) / (diff[i + 1] - diff[i])
            break
    pooled = assign_states(df["chi"].to_numpy(), cfg["basins"]["published"])
    summary = {
        "window_set_for_headline": "published",
        "P_S1": pub["P_S1"].to_dict(), "P_S3": pub["P_S3"].to_dict(), "P_S2": pub["P_S2"].to_dict(),
        "P_S1_CI": {int(t): [pub.loc[t, "P_S1_lo"], pub.loc[t, "P_S1_hi"]] for t in pub.index},
        "P_S3_CI": {int(t): [pub.loc[t, "P_S3_lo"], pub.loc[t, "P_S3_hi"]] for t in pub.index},
        "crossover_T_C_linear_interp": None if np.isnan(cross) else round(float(cross), 1),
        "pooled_fraction": {s: round(float(np.mean(pooled == s)), 4) for s in ("S1", "S2", "S3", "other")},
        "chi2_mode_S1_deg": round(float(np.median(df["chi"][pooled == "S1"])), 1),
        "chi2_mode_S3_deg": round(float(np.median(df["chi"][pooled == "S3"])), 1),
        "tau_S3_frames_max": float(pub["tau_S3_frames"].max()),
        "block_length_frames": bcfg["block_length_frames"], "n_boot": bcfg["n_boot"],
        "note": "Single 1-us replica per T (AMBER03+TIP3P). CIs are within-trajectory moving-block "
                "bootstrap and do not include between-replica variance.",
    }
    (RESULTS_DIR / "populations_summary.json").write_text(json.dumps(summary, indent=1) + "\n")

    typer.echo(out[out.window_set == "published"][["T_C", "P_S1", "P_S1_lo", "P_S1_hi", "P_S3", "P_S3_lo", "P_S3_hi", "P_other", "tau_S3_frames"]].to_string(index=False))
    typer.echo(f"S1/S3 crossover (published windows, linear interpolation): {cross:.1f} C")

    # ---------------- figure ----------------
    set_style()
    import matplotlib.pyplot as plt
    kord = load_kordysh()
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    for s in ("S1", "S3", "S2"):
        y = pub[f"P_{s}"].to_numpy(); lo = pub[f"P_{s}_lo"].to_numpy(); hi = pub[f"P_{s}_hi"].to_numpy()
        ax.errorbar(Ts, y, yerr=[np.clip(y - lo, 0, None), np.clip(hi - y, 0, None)], fmt="o-", color=COLORS[s], ms=5, lw=1.6, capsize=3,
                    label={"S1": "S1 flip-in  (chi2 in [50,150])", "S3": "S3 flip-out (chi2 in [-160,-90])", "S2": "S2 (chi2 in [-90,-30])"}[s])
    ax.set_xlabel("temperature (C)")
    ax.set_ylabel("population of Trp128 chi2 state (1 us AMBER03 MD)")
    ax.set_ylim(-0.02, 1.05)
    ax.set_xticks(Ts)
    if not np.isnan(cross):
        ax.axvline(cross, color="k", ls=":", lw=0.9)
        ax.text(cross + 0.4, 0.98, f"P(S3)=P(S1)\n~{cross:.0f} C", fontsize=8, va="top")
    ax2 = ax.twinx()
    meas = kord[kord.kind == "measured"]
    ax2.plot(kord["T_C"], kord["lambda_max_nm"], "--", color=COLORS["exp"], lw=1, alpha=0.6)
    ax2.plot(meas["T_C"], meas["lambda_max_nm"], "s", color=COLORS["exp"], ms=6, label="Kordysh 2005 lambda_max (exp.)")
    ax2.set_ylabel("Trp emission maximum (nm), Kordysh 2005")
    ax2.set_ylim(332, 352)
    ax2.spines["right"].set_visible(True)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize=7.5)
    ax.set_title("Trp128 rotamer populations vs T (MD) and the experimental red shift")
    FIGURES_DIR.mkdir(exist_ok=True)
    fig.savefig(FIGURES_DIR / "populations_vs_T.png")
    typer.echo(f"wrote {RESULTS_DIR / 'populations.csv'}, {FIGURES_DIR / 'populations_vs_T.png'}")


if __name__ == "__main__":
    app()
