# New MD runs for the Section 4.6 limitations (33 cases, 9 blocks)

Everything here reuses the exact templates of `setup_job_array.py` (quench, builder, relax, tension), so the new runs are
protocol-identical to the 25 existing ones. Only seeds, volume fraction, seed radius, temperature, rate or cell size change.

## 1. What to move from cluster 2 to cluster 1

Run `make_transfer_bundle.sh` **on cluster 2, in the project root**. It checks for and packs only these files (≈150–250 MB):

| File (relative to the project root) | Needed by |
|---|---|
| `common.mod`, `potential.mod` | every case |
| `potentials/library.meam`, `potentials/CaZnMg.meam` | every case (`potential.mod` points to `../potentials/...`) |
| `run_5/quenched_matrix_300K.data` | `ca_matched`, `pack_reps`, `glass0`, `indep_q` baseline (xCa 0.03 glass) |
| `run_1/quenched_matrix_300K.data` | `ca_matched`/`pack_reps`/`glass0`, `large` (xCa 0.05 glass) |
| `run_6/quenched_matrix_300K.data` | `pack_reps`, `glass0` (xCa 0.07 glass) |
| `run_1/AN_composite_relaxed_300K.data` | `temp` |
| `run_3/…relaxed…`, `run_4/…relaxed…` | `rate` |
| `pymatgen_offline_bundle.zip` (optional) | only if cluster 1 cannot `pip install` pymatgen |

Not needed: the trajectory dumps, logs, and anything from `vf_*`, `rate*`, `t310_*`. Blocks `clean`, `indep_q` start from a fresh
quench and need only `common.mod`, `potential.mod`, `potentials/`.
The `stress_strain_properties_all25.csv` in this package is the 25 existing yield strengths (used by `collect_new_md.py`).

## 2. Setup on cluster 1 (once)
```bash
tar xzf md_transfer_bundle.tgz             # in the new project root; add the files of THIS package next to them
# edit cluster.env: partition, LMP_EXEC, MODULES (everything cluster-specific is there)
python3 setup_new_md.py                    # writes 33 case directories, cases/<block>.txt, cases_manifest.csv
./verify_inputs.sh                         # atom counts, MEAM+ZBL present, pymatgen import, sbatch
./patch_for_new_lammps.sh                  # ONLY if LAMMPS stops with "Fix halt command before simulation box is defined"
```
`PYMATGEN_SETUP.md` (existing) describes the pymatgen install; the build stage is the only step that needs it.

## 3. Blocks, what each tests, and suggested order

| Order | Block (`./submit_block.sh <block>`) | Cases | Stages | What it addresses | Cost per case (40 ranks) |
|---|---|---|---|---|---|
| 1 | `ca_matched` | 3 | build, relax, tension | Ca series at **matched** geometry (12 Å seeds, 5/10/15 %) for xCa = 0.03 — completes the missing 0.03→0.05 control next to existing `vf_1-3` (0.05) and `vf_r6_1-3` (0.07) | ~1 day |
| 1 | `pack_reps` | 9 | build, relax, tension | 3 extra seed-placement realizations (12 Å, 10 %) for each Ca level → n = 4 per composition with the existing case | ~1 day |
| 1 | `glass0` | 3 | copyglass, relax, tension | 0 % volume fraction (pure glass) for the three Ca levels: real zero point for the volume-fraction trend | ~1 day |
| 2 | `hcp` | 2 | hcp | Defect-free HCP Mg, tension along [2-1-10] and [0001] (339,264 atoms) | ~1 day |
| 2 | `temp` | 4 | copyrelaxed, tension | 250/300/350/400 K at the baseline composition; 300 K is a same-protocol control | ~1 day |
| 2 | `rate` | 2 | copyrelaxed, tension | 2×10⁹ s⁻¹ at xZn = 0.25 and 0.40 (existing check varied Ca only) | ~0.5 day |
| 3 | `indep_q` | 6 | quench, build, relax, tension | **Independent melt-quench realizations** (new composition-assignment and velocity seeds) for xCa 0.03/0.05/0.07 at 12 Å/10 %: the replicates that include glass-structure scatter | **~2 weeks** |
| 3 | `clean` | 2 | quench, build, relax, tension | Compositions 1 and 9 repeated from one uninterrupted quench (same seeds as the originals) | ~2 weeks |
| 4 (optional) | `large` | 2 | replicate, build, relax, tension | 2×2×2 cell (~3.1 M atoms, ~41 nm), 12 Å seeds at 5 and 15 % — shear-band cell-size test; submitted with 160 ranks | ~2–3 days at 160 ranks (estimate; time one first) |

Cost figures are scaled from the comments in your own scripts (“~1 day” for stages 2–4, “~2 weeks” for a quench on 40 ranks);
calibrate against the wall time in `vf_1/log.4_tension`. Set `TIME` in `cluster.env` to cover the longest block you submit
(the quench blocks need > 14 days).
Disk: each tension run writes a per-atom dump (~1.2 GB for 389 k atoms; ~10 GB per `large` case). `hcp` writes none.

```bash
./submit_block.sh ca_matched ; ./submit_block.sh pack_reps ; ./submit_block.sh glass0     # order 1
./submit_block.sh hcp ; ./submit_block.sh temp ; ./submit_block.sh rate                    # order 2
./submit_block.sh <block> --dry                                    # print the sbatch line only
./submit_block.sh <block> --range 1-3 --partition mediumq --time 2-00:00:00   # subset of cases, chosen queue and wall time
```
Every stage is skipped if its output already exists, so a timed-out job can be resubmitted unchanged
(`DRY_RUN=1 ./run_case.sh <block> <line>` shows what would run).

## 4. What the results can and cannot show (read before interpreting)
* **`pack_reps` vs `indep_q`.** `pack_reps` re-places the seeds in the *same* glass, so it measures seed-placement scatter (like the existing ~1 % near-duplicate pair)
  but shares one glass per composition. `indep_q` varies the glass itself. Statistics for the Ca effect should be quoted from the independent set; the placement set is the cheap supplement.
* **`ca_matched` completes the matched series only if the three Ca levels agree on geometry**: all use 12 Å seeds at the same volume fractions.
  Achieved volume fractions come from the builder printout (`log.2_build`); check they match the existing cases (≈ 5/10/15 %).
* **`hcp`** is a *perfect* crystal, so its strength is an upper bound for a 12 Å seed in a glass; use it as a reference ratio, not as the strength of the nanocrystals.
  The axial stress for `hcp_z` is the `Stress_Z_GPa` column (the collector handles this).
* **`temp`** re-equilibrates 50 ps at each temperature (the existing 310 K point used 20 ps), hence the new 300 K control.
* **`large`** replicates a 20 nm glass 2×2×2 and re-thermalises it for 50 ps, so the amorphous topology is eight identical copies. A band that forms is not an artifact of the cell *size*, but the periodic repetition could favour or suppress particular band orientations; state this in any text that uses it. The builder step on 3.1 M atoms is single-core (memory ≈ several GB, tens of minutes).
* Quench blocks use the original-quench `MQ_TEMPLATE` verbatim, including `restart 50000 mq_restart.a/b`. If a quench job times out, resubmit; the existing `submit_resume_stage2.sh` logic from cluster 2 applies.

## 5. Bring back for analysis
Per case: `tension_stress_strain.txt`, `log.4_tension`, `log.3_relax`, `log.2_build` (achieved volume fraction); for D²min/cluster-count analysis also `tension_trajectory.dump` and `AN_composite_relaxed_300K.data`.
Then `python3 collect_new_md.py stress_strain_properties_all25.csv` prints yield strength / modulus / softening per case,
group means and SDs, and the matched-geometry Ca comparison (Welch tests, split into placement and independent sets).

## 6. What was and was not tested here
Tested with a stand-in Lennard-Jones potential on small cells (LAMMPS 22 Jul 2025): the quench, relax, tension, HCP and replicate inputs run, produce the stress–strain file and the data files; the generated default tension input is byte-identical to the template; the shell runner skip logic, verify script and collector were exercised on synthetic data.
**Not tested:** the real Jang 2019 potential (files not available here), the pymatgen builder (unchanged template; only its constants differ), full-size runs, and the SLURM submission itself.
