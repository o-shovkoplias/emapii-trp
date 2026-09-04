# data/ -- provenance

Total size 3.4 MB. Machine-readable list with byte counts and source paths: `MANIFEST.json`
(written by `scripts/prepare_data.py`). Trajectories themselves are not included
(`../scripts/fetch_trajectories.md`).

All MD data come from **six 1 us unbiased simulations, AMBER03 + TIP3P, T = 25/30/35/40/45/50 C**, started from
PDB 8ONG model 1 (rhombic dodecahedron, ~87 A, 14 281 waters, 44 Na+ / 43 Cl-, 45 548 atoms; GROMACS 2024.4/2025.4,
2 fs, LINCS, V-rescale tau_T 0.1 ps, PME 1.2 nm, Parrinello-Rahman tau_P 2 ps), run on the supervisor's cloud
servers. Protein-only frames were saved every 1 ns (1001 frames per T).

| directory / file | content | how it was made | caveat |
|---|---|---|---|
| `structure/8ONG_model1.pdb` | 2618 atoms, model 1 of the 19-model NMR ensemble + header | first MODEL block of the wwPDB file | atom ordering differs from the GROMACS topology; use only as the "starting structure" |
| `md_clean/clean_perframe.npz` | keys `chi dgln dlys dne dmin sasa KKK K T`, 6006 rows | `papers/buried/figures/make_figures_clean.py` (MDAnalysis `calc_dihedrals` with per-frame box; group unwrapping with `minimize_vectors`; SASA/KKK columns joined 1:1 from the NoJump files) | `chi` = Trp128 chi2 in GROMACS sign; agrees with `gmx angle` to < 1 deg (results/convention_check.json) |
| `md_clean/sasa_indole_{T}.xvg` | `time_ns SASA_indole SASA_benzene chi2 state` | freesasa on the NoJump-unwrapped protein | chi2 column is already min-image clean |
| `md_clean/chi2_gmx_100ps_{T}.xvg` | `time_ps chi2_deg`, 10 001 rows, T = 30-50 only | `gmx angle -type dihedral` on the full 10-ps trajectory, every 10th row kept | no such file exists for 25 C |
| `md_raw_pbc_broken/chi_mda_{T}_pbc_broken.xvg` | `time_ns chi1 chi2` **without** min-image correction | thesis-era MDAnalysis script | 6.7 % of frames flipped ~140 deg; used ONLY by `scripts/06_pbc_demo.py` |
| `metad/plumed_{T}.dat` | WTMetaD inputs (`phi`/`psi` labels are chi1/chi2 of Trp128; atoms 1976-1978-1980-1983 and 1978-1980-1983-1984 in the GROMACS topology) | as run | `RESTART` line commented out; TEMP differs per file |
| `metad/plumed_opes_25.dat` | OPES-Explore input | prepared for md0021 | never run |
| `metad/fes_grid_{T}.npz` | `fes_kJmol` (200x200 float32, axis0 = chi1, axis1 = chi2), `chi1_deg`, `chi2_deg`, `n_hills`, `t_end_ns` | `scripts/prepare_data.py::fes_from_hills`: periodic Gaussian hills summed by FFT on [-pi, pi)^2, F = -(gamma/(gamma-1)) V_bias, minimum set to 0 | 1e6 hills per T (1 us); single walker; NOT converged |
| `metad/fes_basins_{6T,5T,3runs}.json` | thesis-era basin dG, block series | `ang/scripts/plot_metad_fes_6T.py` | `scripts/05_metad_fes.py` reproduces the `full` dG values to 3 decimals from the grids |
| `gmm_legacy/pocket_gmm4_{T}.dat`, `*_crosstab.csv` | thesis-era GMM (n = 4) on WCN, d_S1, d_S3, M, Q_p and its cross-tabs | `ang/scripts/cluster_pocket_n4.py` | different feature set from `scripts/04_gmm_pockets.py`; `chi2` column there is RAW (PBC-broken) |
| `experiment/kordysh2005_lambda_max.csv` | lambda_max(T), 9 points, `kind` = measured / interp | digitised from Kordysh et al. 2005 figures | +/- 1 nm |
