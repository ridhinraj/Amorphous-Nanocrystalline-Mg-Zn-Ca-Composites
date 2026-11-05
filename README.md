# Amorphous-nanocrystalline Mg-Zn-Ca composites

LAMMPS input files, analysis code and processed data for the molecular dynamics study of peak
strength and local Ca–Zn stress asymmetry in amorphous-nanocrystalline (AN) Mg-Zn-Ca composites.

Forty-eight uniaxial tensile simulations of 389,344-atom cells cover ten compositions, three
pure-glass controls, nanocrystal volume fractions of 5–19 % with packing replicates, two seed radii,
four strain rates, five temperatures and two single-crystal references. All runs use the 2NN
MEAM + ZBL Mg-Zn-Ca potential.

**Paper:** R. Raj and A. Zagabathuni, *Molecular Dynamics of Peak Strength and Local Ca–Zn Stress
Asymmetry in Amorphous-Nanocrystalline Mg-Zn-Ca Composites* (submitted). Section numbers below refer
to that manuscript.

## Layout

| Folder | Contents |
|---|---|
| `md/` | LAMMPS inputs: melt–quench, relax, tension at four rates and five temperatures, GSFE, liquid mixing, lattice relaxation, elastic constants, Voronoi re-run |
| `md/submit/` | SLURM submission and status scripts for the two campaigns |
| `analysis/` | Python analysis: RDF, D²min, bond-orientational and chemical short-range order, interphase shells, per-atom stress regression, SHAP, Gaussian process surrogate, shear-band detection |
| `figures/` | Scripts that draw the manuscript figures, and the rendered files |
| `data/` | Stress–strain curves, the per-case property table of all 48 runs, and the summary tables behind the figures |
| `docs/` | The second-campaign planning document |

## Reproducing a run

1. Put the potential files in `potentials/` (see `potentials/README.md`; they are distributed by NIST,
   not here).
2. Quench a matrix: `lmp -in md/week2_melt_quench.in` in a case directory. `common.mod` and
   `potential.mod` are included by every input and are expected one level up.
3. Build a composite from the quenched matrix: `python3 analysis/build_an_composite.py` (needs
   pymatgen).
4. Relax and load: `lmp -in md/week4_tension.in` (1×10⁹ s⁻¹ baseline; the `_1e8`, `_2e9`, `_5e8` and
   `_310K` variants change rate or temperature).
5. Extract properties: `python3 analysis/extract_props.py`, then `analysis/consolidate_results.py`.

`md/submit/submit_block.sh` and `analysis/setup_new_md.py` generate and submit whole blocks of cases
on SLURM; `md/submit/run_case.sh` runs one case and skips any stage whose output already exists.

## Case names

The paper names composites after their nominal nanocrystal volume fraction; the files here use the
run IDs. The full mapping, with every input of every run, is `data/Table_S1_all48.csv`.

| Paper | Run ID | What it is |
|---|---|---|
| B5, B10, B15 | `vf_1`, `vf_2`, `vf_3` | 5/10/15 % at 12 Å seeds, xCa = 0.05 (Composition 1) |
| B5, B10, B15 (Comp. 5) | `vf_r5_1`–`vf_r5_3` | the same series at xCa = 0.03 |
| B5, B10, B15 (Comp. 6) | `vf_r6_1`–`vf_r6_3` | the same series at xCa = 0.07 |
| B18.9, B18.5 | `vf_20_1`, `vf_25_1` | 20 % and 25 % targets, saturated at 18.91 % and 18.49 % |
| B10-R15, B10-R20 | `vf_rad15_1`, `vf_rad20_1` | 10 % at 15 Å and 20 Å seeds |
| Compositions 1–10 | `run_1`–`run_10` | single-seed composition sweep |
| Pure glasses | `glass0_r5`, `glass0_r1`, `glass0_r6` | the three matrices without seeds |
| Packing replicates | `rep_r5_*`, `rep_r1_*`, `rep_r6_*` | three extra seed placements per matrix at 10 % |
| Temperature series | `temp_250`–`temp_400`, `t310_1` | 250–400 K on the Composition 1 structure |
| Strain-rate series | `rate1e8_1`, `rate5e8_1`, `rate2e9_1`, `rate2e9_run3/4/6` | 1×10⁸ to 2×10⁹ s⁻¹ |
| Crystal references | `hcp_x`, `hcp_z` | defect-free HCP Mg loaded along [21̅1̅0] and [0001] |

Some column headers and `model` values in the data files still say *yield* and *16-run*; the paper
calls these **peak strength** (σmax), **peak strain** (εmax) and the **16-point** model. The quantities
are unchanged.

## Requirements

LAMMPS 22 Jul 2025 (update 4) built with the MEAM and ZBL packages; Python 3.10+ with the packages in
`requirements.txt`. OVITO 3.17 was used for the structure renderings in Fig. 2 only.

## Citing

See `CITATION.cff`. Code is released under the MIT licence (`LICENSE`); the data files in `data/` are
released under CC BY 4.0.
