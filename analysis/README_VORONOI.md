# Per-atom volume correction of the Ca-vs-Zn stress comparison (no new MD, no OVITO)

**Why.** `compute stress/atom` returns stress x volume. The manuscript's "Ca carries ~2x the stress of Zn" used those raw
numbers. Ca atoms are ~2-3x larger than Zn atoms, so the ratio may be mostly volume. Dividing each atom's tensor by its own
Voronoi volume gives stress. LAMMPS does the volumes (`compute voronoi/atom`, VORONOI package) by re-reading stored frames.

## Check the cluster LAMMPS first
    lmp -h | grep -i voronoi          # must show VORONOI in "Installed packages" and voronoi/atom
    head -12 run_5/tension_trajectory.dump | grep "ITEM: ATOMS"     # expect: id type x y z c_peratom_stress[1] ... [6]

## Per case (run_1, run_5, run_6 = the fixed-Zn, varying-Ca series), about 10 minutes each
    cd run_5
    python3 ../voronoi_volume/select_frames.py tension_trajectory.dump --list          # frame 1 should be strain 0.0000
    python3 ../voronoi_volume/select_frames.py tension_trajectory.dump \
            --strains 0 <PEAK_STRAIN> 0.15 --start-step 20000 --erate-ps 1e-3 --dt-ps 0.001 -o frames_selected.dump
    mpirun -np 8 lmp -in ../voronoi_volume/voronoi_rerun.in -var DATA AN_composite_relaxed_300K.data \
            -var SEL frames_selected.dump -var OUT volumes.dump
    python3 ../voronoi_volume/voronoi_merge.py frames_selected.dump volumes.dump --label run_5 \
            --start-step 20000 --erate-ps 1e-3 --dt-ps 0.001 \
            --stress-file tension_stress_strain.txt --save-atoms run_5_atoms.npz

`PEAK_STRAIN` is `yield_strain` in stress_strain_properties_all25.csv: run_1 0.0712, run_5 0.0661, run_6 0.0682.
If the 0% frame is not at step 20000, use the first timestep that `--list` shows as `--start-step`.
`--erate-ps` is the number in `fix deform ... erate`; for 1e9 1/s in metal units it is 1e-3.

## What to look at in the output
1. `sum(V_i)/V_box` must be 1.000 (volumes tile the box). Otherwise do not use the frame.
2. "mean stress from atoms" must match "thermo stress" within a few %. This proves the unit and sign handling.
3. The final block compares, per frame, the Ca/Zn **volume** ratio, the **raw** ratio (what the paper reports) and the
   **stress** ratio. If the stress ratio is near 1 the claim was a volume effect; if it stays well above 1 the claim survives.
`<label>_species_summary.csv` holds the numbers; `<label>_atoms.npz` holds per-atom volume and stress for re-running the
SHAP regression with stress as the target.

Tested on a synthetic 2048-atom Mg/Zn/Ca system with the same dump format (not on your data): volumes tile the box exactly and
the atom-sum stress equals the thermo stress to all printed digits. LAMMPS 22 Jul 2025 with the VORONOI package.
