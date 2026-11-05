"""
multiframe_solute_shap.py

Extends the original final-frame-only solute-clustering SHAP analysis
to run at any specified strain, not just 15%. Built per external
review comment #3: the original analysis was conducted solely on the
final frame, deep into post-yield softening; demonstrating that Ca
atoms carry elevated stress AT OR BEFORE the yield point is necessary
to support the claim that Ca load-partitioning CAUSES the bulk yield
strength enhancement, rather than being a downstream correlate
observed only after yielding has already happened.

For each requested strain milestone, this:
  1. Finds the nearest available dump frame to that strain and reports
     the actual strain of the frame used (frames are recorded every
     5000 steps = 0.005 strain at the baseline rate, so the match is
     not always exact).
  2. Builds the local coordination environment (n_Mg/Zn/Ca neighbors
     within 3.25 Ang, matching the cutoff established from the radial
     distribution function) for every Zn and Ca atom, using the SAME
     nonzero-box-origin-aware periodic wrapping already fixed in
     compute_d2min.py and compute_rdf.py earlier in this project.
  3. Trains a gradient-boosted regression model predicting each
     solute atom's log von Mises stress from its coordination and
     elemental identity (own_type_Ca, n_Mg_neighbors, n_Zn_neighbors,
     n_Ca_neighbors -- the same four features as the original
     analysis), computes SHAP values, and reports the Ca-Zn mean
     log-stress difference.

UNIT NOTE (confirmed against this project's dump header and the
existing detect_shear_bands.py, rather than assumed): the archived
dumps carry LAMMPS's raw compute stress/atom output --
c_peratom_stress[1..6] = sxx, syy, szz, sxy, sxz, syz -- which is
stress*volume (virial units), NOT calibrated per-atom stress in GPa.
Per-atom Voronoi volumes were never computed anywhere in this
pipeline, so there is no rigorous way to convert to absolute stress.
This script therefore computes von Mises directly from these six raw
components using the IDENTICAL formula already used throughout this
project's own detect_shear_bands.py, and works in these same raw
stress-volume units throughout. This is valid for the Ca-vs-Zn
RELATIVE comparison this analysis actually needs (a log-difference
between two groups is unaffected by a shared, uncalibrated overall
scale), exactly as the original single-frame version of this analysis
already did -- it is not valid for quoting an absolute per-atom GPa
value.

Usage:
    python3 multiframe_solute_shap.py <dump_path> --strains 0.065 0.10 0.15
"""

import sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
import shap


def von_mises_stress(sxx, syy, szz, sxy, sxz, syz):
    """Identical formula to detect_shear_bands.py's von_mises_stress,
    reused here rather than reimplemented, so the two analyses stay
    consistent with each other."""
    return np.sqrt(
        0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
        + 3.0 * (sxy ** 2 + sxz ** 2 + syz ** 2)
    )

MG_TYPE, ZN_TYPE, CA_TYPE = 1, 2, 3
T0_STEPS = 20000
TIMESTEP_PS = 0.001


def frame_strain(timestep, erate_per_ps):
    return (timestep - T0_STEPS) * erate_per_ps * TIMESTEP_PS


STRESS_COLS = ["c_peratom_stress[1]", "c_peratom_stress[2]", "c_peratom_stress[3]",
               "c_peratom_stress[4]", "c_peratom_stress[5]", "c_peratom_stress[6]"]


def parse_dump_frames(dump_path, stress_cols=STRESS_COLS):
    """Reads id, type, x, y, z, and the six raw per-atom stress tensor
    components for every frame, computing von Mises from them via
    von_mises_stress (see UNIT NOTE at the top of this file). Column
    names must match the dump's own ITEM: ATOMS header exactly --
    defaults to this project's actual dump format, confirmed via
    `head run_1/tension_trajectory.dump | grep "ITEM: ATOMS"` rather
    than assumed."""
    frames = []
    with open(dump_path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            timestep = int(f.readline().strip())
            f.readline()
            n_atoms = int(f.readline().strip())
            f.readline()
            box = np.array([list(map(float, f.readline().split())) for _ in range(3)])
            header = f.readline().strip().replace("ITEM: ATOMS", "").split()
            idx = {c: header.index(c) for c in ("id", "type", "x", "y", "z", *stress_cols)}
            ids = np.empty(n_atoms, dtype=np.int64)
            types = np.empty(n_atoms, dtype=np.int64)
            pos = np.empty((n_atoms, 3))
            stress = np.empty((n_atoms, 6))
            for i in range(n_atoms):
                parts = f.readline().split()
                ids[i] = int(parts[idx["id"]])
                types[i] = int(parts[idx["type"]])
                pos[i] = [float(parts[idx["x"]]), float(parts[idx["y"]]), float(parts[idx["z"]])]
                stress[i] = [float(parts[idx[c]]) for c in stress_cols]
            order = np.argsort(ids)
            vm = von_mises_stress(*stress[order].T)
            frames.append(dict(timestep=timestep, box=box, ids=ids[order],
                                types=types[order], pos=pos[order], vm=vm))
    return frames


def build_coordination(types, pos, box, cutoff=3.25):
    """Per-atom neighbor counts by element, using nonzero-origin-aware
    periodic wrapping (the same fix applied to compute_d2min.py and
    compute_rdf.py earlier in this project)."""
    L = box[:, 1] - box[:, 0]
    wrapped = (pos - box[:, 0]) % L
    tree = cKDTree(wrapped, boxsize=L)
    pairs = tree.query_pairs(cutoff, output_type="ndarray")

    n = len(types)
    n_mg = np.zeros(n, dtype=int)
    n_zn = np.zeros(n, dtype=int)
    n_ca = np.zeros(n, dtype=int)
    counts = {MG_TYPE: n_mg, ZN_TYPE: n_zn, CA_TYPE: n_ca}
    for i, j in pairs:
        counts[types[j]][i] += 1
        counts[types[i]][j] += 1
    return n_mg, n_zn, n_ca


def analyze_frame(frame, cutoff=3.25, test_size=0.25, random_state=0):
    n_mg, n_zn, n_ca = build_coordination(frame["types"], frame["pos"], frame["box"], cutoff)
    is_solute = np.isin(frame["types"], [ZN_TYPE, CA_TYPE])
    vm = frame["vm"][is_solute]
    valid = vm > 0  # log requires strictly positive stress
    df = pd.DataFrame({
        "own_type_Ca": (frame["types"][is_solute] == CA_TYPE).astype(int)[valid],
        "n_Mg_neighbors": n_mg[is_solute][valid],
        "n_Zn_neighbors": n_zn[is_solute][valid],
        "n_Ca_neighbors": n_ca[is_solute][valid],
        "log_vm": np.log(vm[valid]),
    })

    X, y = df.drop(columns="log_vm"), df["log_vm"]
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=test_size, random_state=random_state)
    model = GradientBoostingRegressor(random_state=random_state)
    model.fit(X_tr, y_tr)
    r2 = r2_score(y_te, model.predict(X_te))

    explainer = shap.TreeExplainer(model)
    shap_values = explainer(X_te)
    mean_abs_shap = pd.Series(np.abs(shap_values.values).mean(axis=0), index=X.columns).sort_values(ascending=False)

    ca_mean = df.loc[df.own_type_Ca == 1, "log_vm"].mean()
    zn_mean = df.loc[df.own_type_Ca == 0, "log_vm"].mean()
    ca_se = df.loc[df.own_type_Ca == 1, "log_vm"].sem()
    zn_se = df.loc[df.own_type_Ca == 0, "log_vm"].sem()

    return dict(r2=r2, mean_abs_shap=mean_abs_shap, ca_zn_diff=ca_mean - zn_mean,
                ca_zn_diff_se=np.sqrt(ca_se ** 2 + zn_se ** 2),
                n_ca=int((df.own_type_Ca == 1).sum()), n_zn=int((df.own_type_Ca == 0).sum()))


def run(dump_path, target_strains, erate_per_ps=0.001, cutoff=3.25):
    frames = parse_dump_frames(dump_path)
    strains = np.array([frame_strain(fr["timestep"], erate_per_ps) for fr in frames])

    results = []
    for target in target_strains:
        i = int(np.argmin(np.abs(strains - target)))
        actual = strains[i]
        print(f"\n=== Target strain {target:.3f} -> nearest available frame at strain {actual:.3f} "
              f"(timestep {frames[i]['timestep']}) ===")
        res = analyze_frame(frames[i], cutoff=cutoff)
        res["target_strain"] = target
        res["actual_strain"] = actual
        results.append(res)
        print(f"  n_Ca={res['n_ca']}, n_Zn={res['n_zn']}, held-out R^2={res['r2']:.3f}")
        print(f"  Ca-Zn log(vm) difference = {res['ca_zn_diff']:+.3f} +/- {res['ca_zn_diff_se']:.3f}")
        print("  Mean |SHAP| by feature:")
        for feat, val in res["mean_abs_shap"].items():
            print(f"    {feat:<18} {val:.4f}")
    return results


if __name__ == "__main__":
    import argparse
    # NOTE: an earlier manual sys.argv-slicing version of this block had a
    # real bug, caught on first synthetic test -- it swept the VALUE of
    # --erate_per_ps (e.g. "0.001") into the strains list as a spurious
    # extra target, since it only excluded tokens starting with "--", not
    # the value immediately following one. argparse handles this correctly
    # by construction.
    parser = argparse.ArgumentParser()
    parser.add_argument("dump_path")
    parser.add_argument("strains", type=float, nargs="*", default=[0.065, 0.10, 0.15])
    parser.add_argument("--erate_per_ps", type=float, default=0.001)
    args = parser.parse_args()
    run(args.dump_path, args.strains, erate_per_ps=args.erate_per_ps)
