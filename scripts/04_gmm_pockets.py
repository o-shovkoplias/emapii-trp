#!/usr/bin/env python3
"""Analysis 4 -- Gaussian-mixture clustering of the indole micro-environment.

Features (per frame, PBC-clean): dgln, dlys, dne, dmin (A) and indole SASA (A2), z-scored.
BIC is scanned over n = 1..8; n = 4 is then fitted (a choice made for interpretability --
BIC keeps decreasing, see results/gmm_bic.csv) and cross-tabulated against the chi2 state
and the temperature. Cluster names are assigned from the fitted centroids by a fixed rule
(printed in results/gmm_centroids.csv) so they are reproducible, not hand-picked.

Outputs: results/gmm_bic.csv, results/gmm_centroids.csv, results/gmm_labels.csv,
         results/gmm_crosstab_cluster_state.csv, results/gmm_crosstab_cluster_T.csv,
         results/gmm_summary.json, figures/gmm_pockets.png
"""
from __future__ import annotations

import json
from pathlib import Path

import _common  # noqa: F401
import numpy as np
import pandas as pd
import typer
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

from emapii_trp import FIGURES_DIR, RESULTS_DIR
from emapii_trp.basins import assign_states
from emapii_trp.config import load_config
from emapii_trp.io import load_clean_perframe
from emapii_trp.plotting import COLORS, set_style

app = typer.Typer(add_completion=False)
CLUSTER_COLORS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"]


def name_clusters(cent: pd.DataFrame) -> dict[int, str]:
    """Rule: closest-to-Gln132 -> 'Gln132 pocket'; then closest-to-Lys124-126 -> 'Lys face';
    then largest SASA -> 'loose / exposed'; then smallest SASA -> 'deep buried'; rest 'other-k'."""
    names: dict[int, str] = {}
    free = set(cent.index)
    k = cent.loc[list(free), "dgln"].idxmin(); names[k] = "Gln132 pocket"; free.remove(k)
    k = cent.loc[list(free), "dlys"].idxmin(); names[k] = "Lys face"; free.remove(k)
    if free:
        k = cent.loc[list(free), "sasa"].idxmax(); names[k] = "loose / exposed"; free.remove(k)
    if free:
        k = cent.loc[list(free), "sasa"].idxmin(); names[k] = "deep buried"; free.remove(k)
    for i, k in enumerate(sorted(free)):
        names[k] = f"other-{i + 1}"
    return names


@app.command()
def main(config: Path = typer.Option(Path("config.yaml"))) -> None:
    """Scan BIC, fit n=4 GMM, cross-tabulate, plot."""
    cfg = load_config(config); gcfg = cfg["gmm"]
    df = load_clean_perframe()
    df["state"] = assign_states(df["chi"].to_numpy(), cfg["basins"]["published"])
    X = StandardScaler().fit_transform(df[gcfg["features"]].to_numpy())

    bic_rows = []
    for n in gcfg["n_components_scan"]:
        gm = GaussianMixture(n, covariance_type="full", random_state=gcfg["seed"], n_init=gcfg["n_init"]).fit(X)
        bic_rows.append({"n_components": n, "BIC": round(gm.bic(X), 1), "AIC": round(gm.aic(X), 1)})
    bic = pd.DataFrame(bic_rows)
    RESULTS_DIR.mkdir(exist_ok=True)
    bic.to_csv(RESULTS_DIR / "gmm_bic.csv", index=False)

    n = gcfg["n_components"]
    gm = GaussianMixture(n, covariance_type="full", random_state=gcfg["seed"], n_init=gcfg["n_init"]).fit(X)
    df["cluster"] = gm.predict(X)
    cent = df.groupby("cluster")[gcfg["features"]].mean().round(2)
    cent["n"] = df.groupby("cluster").size()
    cent["frac_S3"] = df.groupby("cluster").state.apply(lambda s: round((s == "S3").mean(), 3))
    cent["frac_S1"] = df.groupby("cluster").state.apply(lambda s: round((s == "S1").mean(), 3))
    names = name_clusters(cent)
    cent["name"] = pd.Series(names)
    df["cluster_name"] = df.cluster.map(names)
    cent.to_csv(RESULTS_DIR / "gmm_centroids.csv")
    df[["T", "frame", "chi", "state", "cluster", "cluster_name"]].to_csv(RESULTS_DIR / "gmm_labels.csv", index=False)

    ct_state = pd.crosstab(df.cluster_name, df.state)
    ct_state.to_csv(RESULTS_DIR / "gmm_crosstab_cluster_state.csv")
    ct_T = pd.crosstab(df.cluster_name, df["T"], normalize="columns").round(3)
    ct_T.to_csv(RESULTS_DIR / "gmm_crosstab_cluster_T.csv")
    ct_3 = pd.crosstab([df["T"], df.state], df.cluster_name)
    ct_3.to_csv(RESULTS_DIR / "gmm_crosstab_T_state_cluster.csv")

    gln = df.cluster_name == "Gln132 pocket"; lys = df.cluster_name == "Lys face"
    S3 = df.state == "S3"; S1 = df.state == "S1"
    summary = {
        "n_components": n, "features": gcfg["features"], "BIC_argmin_in_scan": int(bic.loc[bic.BIC.idxmin(), "n_components"]),
        "centroids": json.loads(cent.to_json(orient="index")),
        "P(Gln132 pocket | S3)": round(float(gln[S3].mean()), 3), "P(S3 | Gln132 pocket)": round(float(S3[gln].mean()), 3),
        "P(Lys face | S1)": round(float(lys[S1].mean()), 3), "P(S1 | Lys face)": round(float(S1[lys].mean()), 3),
        "Gln132_pocket_fraction_per_T": {int(k): float(v) for k, v in ct_T.loc["Gln132 pocket"].items()},
        "Lys_face_fraction_per_T": {int(k): float(v) for k, v in ct_T.loc["Lys face"].items()} if "Lys face" in ct_T.index else {},
        "note": "n=4 is an interpretability choice; the BIC minimum is at the end of the scan (n=8, not monotonic; see gmm_bic.csv). "
                "Cluster names follow the fixed centroid rule in name_clusters().",
    }
    (RESULTS_DIR / "gmm_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    typer.echo(bic.to_string(index=False)); typer.echo(cent.to_string()); typer.echo(ct_T.to_string())
    typer.echo(f"P(Gln pocket|S3)={summary['P(Gln132 pocket | S3)']}, P(S1|Lys face)={summary['P(S1 | Lys face)']}")

    set_style()
    import matplotlib.pyplot as plt
    order = list(cent.sort_values("dgln").index)
    cmap = {k: CLUSTER_COLORS[i] for i, k in enumerate(order)}
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.4))
    ax = axes[0, 0]
    ax.plot(bic.n_components, bic.BIC, "o-", color="k"); ax.axvline(n, color=COLORS["accent"], ls="--", lw=1)
    ax.set_xlabel("number of Gaussian components"); ax.set_ylabel("BIC (lower is better)"); ax.set_title("(a) BIC scan -- n = 4 chosen for interpretability")
    ax = axes[0, 1]
    for k in order:
        g = df[df.cluster == k]
        ax.scatter(g.dlys, g.dgln, s=5, alpha=0.35, color=cmap[k], label=f"{names[k]} (n={len(g)})", rasterized=True)
    ax.set_xlabel("d(indole -> Lys124-126 face) (A)"); ax.set_ylabel("d(indole -> Gln132 side chain) (A)")
    ax.set_title("(b) clusters in the two pocket distances"); ax.legend(fontsize=7.5, markerscale=3)
    ax = axes[1, 0]
    bottom = np.zeros(len(ct_T.columns))
    for k in order:
        vals = ct_T.loc[names[k]].to_numpy()
        ax.bar(ct_T.columns.astype(str), vals, bottom=bottom, color=cmap[k], label=names[k], width=0.7); bottom += vals
    ax.set_xlabel("temperature (C)"); ax.set_ylabel("fraction of frames"); ax.set_title("(c) cluster occupancy vs temperature"); ax.legend(fontsize=7.5, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax = axes[1, 1]
    comp = pd.crosstab(df.cluster_name, df.state, normalize="index").reindex([names[k] for k in order])
    bottom = np.zeros(len(comp))
    for s in ("S1", "S2", "S3", "other"):
        if s in comp:
            ax.barh(comp.index, comp[s], left=bottom, color=COLORS[s], label=s); bottom += comp[s].to_numpy()
    ax.set_xlabel("fraction of cluster frames"); ax.set_title("(d) chi2 state composition of each cluster")
    ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=4)
    fig.tight_layout()
    FIGURES_DIR.mkdir(exist_ok=True)
    fig.savefig(FIGURES_DIR / "gmm_pockets.png")
    typer.echo(f"wrote results/gmm_*.csv, {FIGURES_DIR / 'gmm_pockets.png'}")


if __name__ == "__main__":
    app()
