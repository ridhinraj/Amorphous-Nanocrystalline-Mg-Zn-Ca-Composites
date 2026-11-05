# Figures

The script names come from an earlier figure numbering and no longer match the figure numbers in the
paper. The mapping is below. All scripts import `style.py` (one palette, colour-blind safe) and
`figdata.py` (loaders), and write a PNG and a PDF at 600 dpi.

| Figure in the paper | Script | Rendered |
|---|---|---|
| 1. Workflow schematic | `schematic/figure01_schematic.js` (node + pptxgenjs) | `rendered/Figure_01.png` |
| 2. Atomistic views (OVITO) | — rendered by hand in OVITO 3.17 | `rendered/Figure_02_PLACEHOLDER.png` |
| 3. Radial distribution functions | — see note below | `rendered/Figure_03.png` |
| 4. Density against zinc content | — see note below | `rendered/Figure_04.png` |
| 5. Thermal history | — see note below | `rendered/Figure_05.png` |
| 6. Stress–strain curves | `fig05.py` | `rendered/Figure_06.pdf` |
| 7. Composition dependence | `fig06.py` | `rendered/Figure_07.pdf` |
| 8. Volume fraction and seed radius | `fig07.py` | `rendered/Figure_08.pdf` |
| 9. D²min, cumulative and incremental | — see note below | `rendered/Figure_09.png` |
| 10. Interphase ratios | `fig_interphase.py` | `rendered/Figure_10.pdf` |
| 11. SHAP and stress asymmetry | `fig10.py` | `rendered/Figure_11.pdf` |
| 12. q6 and Warren–Cowley α | `fig11.py` | `rendered/Figure_12.pdf` |
| 13. Strain rate and temperature | `fig12.py` | `rendered/Figure_13.pdf` |
| 14. Surrogate predictions | `fig13.py` | `rendered/Figure_14.pdf` |

`fig01.py` draws an earlier version of the workflow schematic in matplotlib. It is kept because it is
self-contained; the figure in the paper is the one from `schematic/`.

**Figures 3, 4, 5 and 9** were drawn before the current toolchain and no plotting script survives for
them. The data behind them are in the repository: `data/run4_rdf.csv` (Fig. 3), the density column of
`data/Table_S1_all48.csv` and Supplementary Table S3 (Fig. 4), `data/tg_summary_all.csv` (Fig. 5) and
`data/d2min/` (Fig. 9).

## Paths

The scripts read from `data/` relative to the repository root. If you move them, set `ROOT` at the top
of `figdata.py`.
