# Trajectories

The MD trajectories are not part of this repository (six protein-only `.xtc` files at 1 ns, ~60 MB each;
the full-system trajectories with water are several GB each; six PLUMED `HILLS` files, 162 MB each).

They are available on request from the author (Oleksandr Shovkoplias, o.s.shovkoplias@gmail.com) or via the
supervisor's group (Department of Molecular Physics, Faculty of Physics, Taras Shevchenko National University
of Kyiv). Please state which of the following you need:

| set | force field | length | frames | files |
|---|---|---|---|---|
| unbiased, T = 25/30/35/40/45/50 C | AMBER03 + TIP3P | 1 us each | 1001 (1 ns) protein-only; 100 001 (10 ps) full | `protein_dt1000.xtc`, `md000{1..6}_{T}.tpr`, topology PDB |
| WTMetaD, same six T | AMBER03 + TIP3P, PLUMED 2.10 | 1 us each | 1e6 hills | `HILLS`, `COLVAR`, `plumed_{T}.dat` (inputs are in `data/metad/`) |
| unbiased test, T = 25/37/42 C | CHARMM36m + TIP3P | 200 ns each | 20 001 (10 ps) full system | `md.xtc`, `md.tpr` |

With the six protein-only `.xtc` files and their GROMACS-ordered topology PDB in place, the shipped per-frame
table can be regenerated with MDAnalysis:

```bash
python scripts/prepare_data.py --source-root /path/to/thesis-tree --pdb /path/to/8ong.pdb
```

Mandatory analysis rules learned the hard way (both are enforced in `prepare_data.py` and the thesis scripts):

1. **Minimum-image / NoJump.** 70-97 % of frames have the protein split across the periodic boundary. Dihedrals
   must be computed with the per-frame box (`MDAnalysis.lib.distances.calc_dihedrals(..., box=ts.dimensions)`)
   and SASA on a NoJump-unwrapped protein iterated **sequentially** (NoJump is stateful).
2. **Sign convention.** chi2 = CA-CB-CG-CD1, right-handed IUPAC, identical to `gmx angle -type dihedral`
   (S1 ~ +100 deg, S3 ~ -120 deg). Verify one frame against `gmx angle` before trusting any new code path.
