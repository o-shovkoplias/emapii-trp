#!/usr/bin/env python3
"""Populate data/ from the original thesis working tree (one-off, documented provenance).

This script is NOT part of run_all.sh: the shipped data/ directory is already the
output of this script. It is kept so that every file in data/ has a machine-readable
origin. It only copies / downsamples / reduces; it never modifies the sources.

Usage (from the repository root):
    python scripts/prepare_data.py --source-root /path/to/diploma \
        --pdb /path/to/8ong.pdb [--skip-hills]

The heavy step is reducing six 162 MB PLUMED HILLS files (1e6 Gaussians each) to
200x200 free-energy grids (data/metad/fes_grid_{T}.npz, ~160 kB each) using the
same FFT-convolution reconstruction as the thesis (well-tempered rescaling
gamma/(gamma-1)). Pass --skip-hills to leave them out.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import typer

app = typer.Typer(add_completion=False)
REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
TS = [25, 30, 35, 40, 45, 50]

HILLS_DIRS = {25: "wtmetad_25", 30: "wtmetad_30", 35: "wtmetad_35",
              40: "wtmetad_40_new", 45: "wtmetad_45", 50: "wtmetad_50"}


def _copy(src: Path, dst: Path, manifest: list[dict], note: str) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    manifest.append({"dest": str(dst.relative_to(REPO)), "source": str(src),
                     "bytes": dst.stat().st_size, "note": note})
    typer.echo(f"  {dst.relative_to(REPO)}  <-  {src}")


def read_hills(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (time_ps, chi1_rad, chi2_rad, height_kJmol) from a PLUMED HILLS file."""
    d = np.loadtxt(path, comments="#", usecols=(0, 1, 2, 5))
    return d[:, 0], d[:, 1], d[:, 2], d[:, 3]


def fes_from_hills(chi1: np.ndarray, chi2: np.ndarray, h: np.ndarray, *, nbins: int,
                   sigma: float, bias_factor: float) -> np.ndarray:
    """Sum periodic Gaussian hills on a grid by FFT convolution; return F (kJ/mol, min = 0).

    F = -(gamma/(gamma-1)) * V_bias, the standard well-tempered estimator
    (Barducci, Bussi, Parrinello 2008). Grid is [-pi, pi)^2, axis 0 = chi1, axis 1 = chi2.
    """
    gmin, gmax = -np.pi, np.pi
    dx = (gmax - gmin) / nbins
    i1 = np.floor((chi1 - gmin) / dx).astype(int) % nbins
    i2 = np.floor((chi2 - gmin) / dx).astype(int) % nbins
    hist = np.zeros((nbins, nbins))
    np.add.at(hist, (i1, i2), h)
    idx = np.arange(nbins)
    dist = np.minimum(idx, nbins - idx) * dx
    k1 = np.exp(-dist**2 / (2 * sigma**2))
    kernel = np.outer(k1, k1)
    bias = np.real(np.fft.ifft2(np.fft.fft2(hist) * np.fft.fft2(kernel)))
    fes = -(bias_factor / (bias_factor - 1.0)) * bias
    return fes - fes.min()


@app.command()
def main(
    source_root: Path = typer.Option(..., help="Root of the thesis working tree (contains ang/, simulations/, papers/)."),
    pdb: Path = typer.Option(..., help="Path to the downloaded 8ONG PDB file (19 NMR models)."),
    skip_hills: bool = typer.Option(False, help="Do not reduce HILLS files to FES grids."),
    chi2_stride: int = typer.Option(10, help="Keep every Nth row of the 10-ps gmx chi2 files (10 -> 100 ps)."),
) -> None:
    """Copy and reduce the source files into data/ and write data/MANIFEST.json."""
    manifest: list[dict] = []
    ang = source_root / "ang"

    typer.echo("[1] PBC-clean per-frame table (6006 frames = 6 T x 1001 frames at 1 ns)")
    _copy(source_root / "papers/buried/figures/clean_perframe.npz",
          DATA / "md_clean/clean_perframe.npz", manifest,
          "Produced by papers/buried/figures/make_figures_clean.py from cloud protein_dt1000.xtc with "
          "box-aware MDAnalysis calc_dihedrals (min-image), NoJump SASA, box-aware distances.")

    typer.echo("[2] NoJump indole SASA + chi2 (1 ns), and 100-ps gmx chi2 (PBC-clean, gmx angle)")
    for T in TS:
        _copy(ang / f"results/{T}/sasa_indole_{T}_fixed.xvg", DATA / f"md_clean/sasa_indole_{T}.xvg",
              manifest, "freesasa on NoJump-unwrapped protein; cols: time_ns SASA_indole SASA_benzene chi2 state")
        g = ang / f"results/{T}/chi2_gromacs_correct.xvg"
        if g.exists():
            raw = np.loadtxt(g, comments=["#", "@"])
            sub = raw[::chi2_stride]
            dst = DATA / f"md_clean/chi2_gmx_{10 * chi2_stride}ps_{T}.xvg"
            hdr = (f"# Trp128 chi2 (CA-CB-CG-CD1), gmx angle -type dihedral on the full 1-us AMBER03 "
                   f"trajectory at T={T} C\n# downsampled every {chi2_stride} rows from a 10-ps file "
                   f"({len(raw)} rows) by scripts/prepare_data.py\n# cols: time_ps  chi2_deg (GROMACS sign)\n")
            dst.write_text(hdr + "\n".join(f"{t:.1f} {c:.3f}" for t, c in sub) + "\n")
            manifest.append({"dest": str(dst.relative_to(REPO)), "source": str(g), "bytes": dst.stat().st_size,
                             "note": f"downsampled x{chi2_stride} from 10-ps gmx angle output"})
            typer.echo(f"  {dst.relative_to(REPO)}  <-  {g} (x{chi2_stride})")
        else:
            typer.echo(f"  (no gmx chi2 file for T={T}; only 1-ns MDAnalysis data exist)")

    typer.echo("[3] RAW (PBC-broken) chi1/chi2 files -- for the 'why NoJump matters' demonstration only")
    for T in TS:
        _copy(ang / f"results/{T}/chi_mda_{T}.xvg", DATA / f"md_raw_pbc_broken/chi_mda_{T}_pbc_broken.xvg",
              manifest, "MDAnalysis dihedral WITHOUT box (min-image) -> ~180 deg flips when the side chain straddles the box; do not use for populations")

    typer.echo("[4] WTMetaD: PLUMED inputs, basin summaries, FES grids")
    for T in TS:
        _copy(source_root / f"simulations/plumed_{T}.dat", DATA / f"metad/plumed_{T}.dat", manifest,
              "PLUMED 2.10 WTMetaD input actually used (phi/psi labels = chi1/chi2 of Trp128)")
    _copy(source_root / "simulations/plumed_opes_25.dat", DATA / "metad/plumed_opes_25.dat", manifest,
          "OPES-Explore input prepared but NOT run (md0021)")
    for name in ("fes_basins_6T.json", "fes_basins_5T.json", "fes_basins_3runs.json"):
        _copy(ang / f"results/metad/{name}", DATA / f"metad/{name}", manifest,
              "basin free energies from ang/scripts/plot_metad_fes_6T.py (FFT hill summation, chi1=trans gate)")
    if not skip_hills:
        for T in TS:
            hills = ang / f"cloud_data/{HILLS_DIRS[T]}/HILLS"
            typer.echo(f"  reducing {hills} ...")
            t_ps, c1, c2, h = read_hills(hills)
            fes = fes_from_hills(c1, c2, h, nbins=200, sigma=0.35, bias_factor=10.0)
            centers = np.degrees(np.linspace(-np.pi, np.pi, 201)[:-1] + np.pi / 200)
            dst = DATA / f"metad/fes_grid_{T}.npz"
            np.savez_compressed(dst, fes_kJmol=fes.astype(np.float32), chi1_deg=centers, chi2_deg=centers,
                                n_hills=len(h), t_end_ns=t_ps[-1] / 1000.0, T_C=T,
                                sigma_rad=0.35, bias_factor=10.0, source=str(hills))
            manifest.append({"dest": str(dst.relative_to(REPO)), "source": str(hills), "bytes": dst.stat().st_size,
                             "note": f"200x200 FES grid from {len(h)} hills (t_end={t_ps[-1]/1000:.0f} ns), FFT sum, WT factor gamma/(gamma-1)"})
            typer.echo(f"    -> {dst.relative_to(REPO)}: {len(h)} hills, {t_ps[-1]/1000:.0f} ns")

    typer.echo("[5] Legacy GMM (n=4) pocket labels + cross-tabs")
    for T in TS:
        _copy(ang / f"results/pocket_classification/pocket_gmm4_{T}.dat", DATA / f"gmm_legacy/pocket_gmm4_{T}.dat",
              manifest, "thesis-era GMM n=4 on WCN/d_S1/d_S3/M/Q_p (ang/scripts/cluster_pocket_n4.py); chi2 column is RAW/PBC-broken")
    for name in ("chi_basin_gmm_crosstab.csv", "chi_6basin_gmm_crosstab.csv", "chi_6basin_3x2_gmm_crosstab.csv"):
        _copy(ang / f"results/pocket_classification/{name}", DATA / f"gmm_legacy/{name}", manifest,
              "thesis-era cross-tab of (chi1,chi2) basins vs GMM clusters")

    typer.echo("[6] Starting structure: PDB 8ONG model 1 only")
    lines = pdb.read_text().splitlines()
    out: list[str] = []
    in_model = False
    for ln in lines:
        if ln.startswith("MODEL") and ln.split()[1] == "1":
            in_model = True
        if in_model:
            out.append(ln)
        if in_model and ln.startswith("ENDMDL"):
            break
    header = [ln for ln in lines if ln.startswith(("HEADER", "TITLE", "COMPND", "EXPDTA", "AUTHOR", "JRNL", "REMARK   2", "SEQRES"))]
    dst = DATA / "structure/8ONG_model1.pdb"
    dst.write_text("\n".join(header + out + ["END"]) + "\n")
    manifest.append({"dest": str(dst.relative_to(REPO)), "source": str(pdb), "bytes": dst.stat().st_size,
                     "note": "model 1 of the 19-model NMR ensemble (wwPDB 8ONG); MD started from this model"})
    typer.echo(f"  {dst.relative_to(REPO)}: {sum(1 for l in out if l.startswith('ATOM'))} ATOM records")

    (DATA / "MANIFEST.json").write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")
    typer.echo(f"wrote data/MANIFEST.json ({len(manifest)} entries)")


if __name__ == "__main__":
    app()
