# Data

Processed results. Raw trajectory dumps (~1.2 GB per tension run) are not included; they are
available from the corresponding author on request.

| File | What it holds |
|---|---|
| `Table_S1_all48.csv` | One row per simulation: composition, volume fraction, seed radius, strain rate, temperature, E, peak strength, peak strain, flow stress at 15 % and softening. This is Supplementary Table S1 of the paper. |
| `stress_strain_properties_all25.csv` | Property table of the first campaign, as written by `analysis/extract_props.py`. Column names are the script's originals (`yield_strength_GPa`, `yield_strain`); these are the peak strength and peak strain of the paper. |
| `new_md_properties_final.csv` | The same for the second campaign (`analysis/collect_new_md.py`). |
| `heldout_predictions.csv` | Gaussian process predictions for every held-out simulation, with σ, the simulated value, the error and z. Supplementary Table S11. The `model` values `13-run`, `16-run` and `25-run` are the paper's 13-, 16- and 25-point models. |
| `tg_summary_all.csv` | Glass transition and fictive temperatures per composition (Supplementary Table S4). |
| `run4_rdf.csv` | Total and partial radial distribution functions of Composition 4 (Fig. 3). |
| `A3_shear_summary_all21.csv` | Shear-band detection: DBSCAN cluster counts per frame for every glass and composite (Supplementary Table S6). |
| `volume_fraction_vs_glass_summary.csv` | Peak strength against volume fraction, relative to the pure glass of the same matrix (Fig. 8). |
| `volume_normalized_regression_summary.csv` | Held-out R² of the per-atom stress regressions and the Ca–Zn log-stress differences (Tables 6 and 7). |
| `stress_strain/` | Engineering stress–strain curves, one file per case. See `MISSING_CURVES.md`. |
| `d2min/` | Per-frame D²min cluster summaries, cumulative and incremental, for the cases analysed in Section 3.5 (Fig. 9). |
| `interphase/` | Per-seed interphase-to-bulk D²min ratios for B5 and B15 at each strain (Fig. 10, Supplementary Table S8). |

## Still to add

- **24 of the 48 stress–strain curves** (the whole second campaign). `stress_strain/MISSING_CURVES.md`
  lists them by case name.
- **Per-atom Voronoi volumes** for the thirteen configurations analysed in Section 3.6. The paper's
  Data availability sentence promises these. They are produced by `md/voronoi_rerun.in` plus
  `analysis/voronoi_py.py`; either add the per-atom files or, if they are too large, change the
  manuscript sentence to say they are available on request.
