#!/usr/bin/env python3
"""
combined_regression.py -- does adding local-structure descriptors improve the per-atom stress regression, and how much of
the Ca-vs-Zn stress premium do they account for?   (guide's comment #6; manuscript Sections 3.6 / 4.1 / 4.6)

Reproduces the original model exactly (model A: own-type Ca + Mg/Zn/Ca neighbor counts at 3.25 A, gradient boosting,
75/25 split, random_state 0) and then adds, on the SAME rows and the SAME split:
    B  + q4, q6 (bond-orientational order, 4.0 A)
    C  + Warren-Cowley own-type alpha (4.0 A)
    D  + size / packing proxies: mean neighbor distance and coordination number (4.0 A)
    E  all of the above + 4.0 A neighbor counts by species
and reports held-out R^2.  Two further analyses address the mechanism question directly:
  * STRUCTURE-ONLY model (no identity, no alpha): out-of-fold predictions give the part of the Ca-Zn log-stress difference that
    local structure alone predicts, and the residual difference that remains;
  * within-species models: how much stress variation structure explains inside Ca atoms and inside Zn atoms.
Without --volumes the target is the raw compute stress/atom von Mises (stress x volume), as in the original analysis.
With --volumes volumes.dump (from voronoi_py.py) every atom's stress tensor is divided by its own Voronoi volume before the
von Mises stress is taken (units GPa), the Voronoi volume is added as a descriptor (models F, G and a second structure-only
model), and the "explained fraction" is re-computed for stress proper rather than stress x volume.

Usage:  python3 combined_regression.py run_1/tension_trajectory.dump 0.070 --erate_per_ps 0.001 --out cr_run1_0070.csv
        python3 combined_regression.py --selftest
Earlier model-A reference values (held-out R^2): run_1 0.212 @0.070, run_5 0.182 @0.065, run_6 0.234 @0.070 (use --expect_r2A).
"""
import argparse, sys, time, os, tempfile
import numpy as np, pandas as pd
import fastlib as fl
from multiframe_solute_shap import build_coordination, von_mises_stress, MG_TYPE, ZN_TYPE, CA_TYPE
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split, cross_val_predict
from sklearn.metrics import r2_score

STRESS = [f"c_peratom_stress[{k}]" for k in range(1, 7)]
BASE = ["own_type_Ca", "n_Mg_neighbors", "n_Zn_neighbors", "n_Ca_neighbors"]          # identical to the original analysis
EXTRA_E = ["q4", "q6", "alpha_own", "r_mean", "N4", "n_Mg4", "n_Zn4", "n_Ca4"]
MODELS = [("A  base (identity + 3.25 A counts)  [original]", BASE),
          ("B  A + q4, q6", BASE + ["q4", "q6"]),
          ("C  A + Warren-Cowley alpha", BASE + ["alpha_own"]),
          ("D  A + mean bond length, coordination", BASE + ["r_mean", "N4"]),
          ("E  all descriptors", BASE + EXTRA_E)]
STRUCT_ONLY = ["n_Mg_neighbors", "n_Zn_neighbors", "n_Ca_neighbors", "n_Mg4", "n_Zn4", "n_Ca4", "N4", "r_mean", "q4", "q6"]

def load_volumes(path, timestep, ids):
    """Voronoi volume (A^3) of every atom, aligned to the id-sorted atoms of the frame, from a voronoi_py.py / LAMMPS dump."""
    info = [f for f in fl.index_dump(path) if f["timestep"] == timestep]
    if not info: sys.exit(f"timestep {timestep} not found in {path}")
    info = info[0]
    with open(path, "rb") as f:
        f.seek(info["offset"]); rows = [f.readline().split() for _ in range(info["n"])]
    arr = np.array(rows, dtype=float); arr = arr[np.argsort(arr[:, info["cols"].index("id")])]
    if not np.array_equal(arr[:, info["cols"].index("id")].astype(int), ids.astype(int)): sys.exit("atom ids differ between dump and volumes file")
    return arr[:, info["cols"].index("c_vor[1]")]

def structural_features(types, pos, box, c_base=3.25, c_big=4.0):
    """features for every atom (rows for non-solutes are computed too but ignored by the caller)"""
    n = len(types); solute = np.isin(types, [ZN_TYPE, CA_TYPE])
    n_mg, n_zn, n_ca = build_coordination(types, pos, box, c_base)          # validated function, exactly the original counts
    i, j, d, r = fl.neighbor_pairs(pos, box, c_big)
    m = solute[i]; i, j, d, r = i[m], j[m], d[m], r[m]
    N4 = np.bincount(i, minlength=n).astype(float)
    cnt = {t: np.bincount(i, weights=(types[j] == t).astype(float), minlength=n) for t in (MG_TYPE, ZN_TYPE, CA_TYPE)}
    q = fl.steinhardt(n, i, d, r)
    x = {t: float(np.mean(types == t)) for t in (MG_TYPE, ZN_TYPE, CA_TYPE)}
    is_ca = types == CA_TYPE
    own = np.where(is_ca, cnt[CA_TYPE], cnt[ZN_TYPE]); xown = np.where(is_ca, x[CA_TYPE], x[ZN_TYPE])
    return pd.DataFrame({"own_type_Ca": is_ca.astype(int), "n_Mg_neighbors": n_mg, "n_Zn_neighbors": n_zn, "n_Ca_neighbors": n_ca,
                         "n_Mg4": cnt[MG_TYPE], "n_Zn4": cnt[ZN_TYPE], "n_Ca4": cnt[CA_TYPE], "N4": N4,
                         "r_mean": np.bincount(i, weights=r, minlength=n) / np.maximum(N4, 1), "q4": q[4], "q6": q[6],
                         "alpha_own": 1.0 - (own / np.maximum(N4, 1)) / xown}), solute

def fit_eval(df, cols, tr, te, seed=0):
    gb = GradientBoostingRegressor(random_state=seed).fit(df.loc[tr, cols], df.loc[tr, "y"])
    return gb, r2_score(df.loc[te, "y"], gb.predict(df.loc[te, cols]))

def analyze(fr, expect_r2A=None, shap_rows=5000, out=None, quiet=False, vol=None):
    t0 = time.time(); types, pos, box = fr["type"], fr["pos"], fr["box"]
    if vol is None:
        vm = von_mises_stress(*[fr[c] for c in STRESS])
    else:
        vm = von_mises_stress(*[fr[c] / vol * 1.0e-4 for c in STRESS])      # bar*A^3 / A^3 -> bar -> GPa
        print(f"TARGET = von Mises of (stress/atom / Voronoi volume), GPa;  mean volume Ca {vol[types == CA_TYPE].mean():.2f}, Zn {vol[types == ZN_TYPE].mean():.2f} A^3")
    F, solute = structural_features(types, pos, box)
    if vol is not None: F["V_vor"] = vol
    keep = solute & (vm > 0)                                  # same filter as the original analysis
    df = F[keep].copy().reset_index(drop=True); df["y"] = np.log(vm[keep])
    nan_rows = df[BASE + EXTRA_E].isna().any(axis=1)
    if nan_rows.any(): print(f"NOTE: dropping {int(nan_rows.sum())} solute atoms with < 4 neighbors at 4.0 A (same rows for every model)"); df = df[~nan_rows].reset_index(drop=True)
    idx = np.arange(len(df)); tr, te = train_test_split(idx, test_size=0.25, random_state=0)      # identical to the original split
    ca, zn = df.own_type_Ca.values == 1, df.own_type_Ca.values == 0
    print(f"solute atoms analyzed: {len(df)} (Ca {ca.sum()}, Zn {zn.sum()}); train {len(tr)} / test {len(te)}   [features {time.time()-t0:.0f} s]")
    rows, models = [], {}
    r2A = None
    model_list = list(MODELS) + ([("F  A + Voronoi volume", BASE + ["V_vor"]), ("G  all descriptors + Voronoi volume", BASE + EXTRA_E + ["V_vor"])] if vol is not None else [])
    for name, cols in model_list:
        gb, r2 = fit_eval(df, cols, tr, te); models[name] = (gb, cols); r2A = r2 if r2A is None else r2A
        rows.append(dict(analysis="model", name=name, n_features=len(cols), r2_test=r2, delta_vs_A=r2 - r2A))
        print(f"  {name:48s} R2(test) = {r2:.3f}   change vs A = {r2 - r2A:+.3f}   [{time.time()-t0:.0f} s]")
    if expect_r2A is not None:
        flag = "MATCHES" if abs(rows[0]["r2_test"] - expect_r2A) < 0.003 else "DOES NOT MATCH"
        print(f"  reproduction check: model A R2 {rows[0]['r2_test']:.3f} vs earlier analysis {expect_r2A:.3f} -> {flag}")

    # ---- how much of the Ca-Zn premium does local structure alone predict?
    y = df.y.values
    gbF = GradientBoostingRegressor(random_state=0); r2F = r2_score(y[te], gbF.fit(df.loc[tr, STRUCT_ONLY], y[tr]).predict(df.loc[te, STRUCT_ONLY]))
    oof = cross_val_predict(GradientBoostingRegressor(random_state=0), df[STRUCT_ONLY], y, cv=3)
    d_raw = y[ca].mean() - y[zn].mean(); d_pred = oof[ca].mean() - oof[zn].mean(); d_res = (y - oof)[ca].mean() - (y - oof)[zn].mean()
    print(f"\n  STRUCTURE-ONLY model (no identity, no alpha): R2(test) = {r2F:.3f}   [{time.time()-t0:.0f} s]")
    print(f"     Ca-Zn log-stress difference: observed {d_raw:+.3f}; predicted from local structure alone {d_pred:+.3f} ({100*d_pred/d_raw:.0f}% of it); residual {d_res:+.3f}")
    if vol is not None:
        so2 = STRUCT_ONLY + ["V_vor"]
        oof2 = cross_val_predict(GradientBoostingRegressor(random_state=0), df[so2], y, cv=3)
        d_pred2 = oof2[ca].mean() - oof2[zn].mean(); d_res2 = (y - oof2)[ca].mean() - (y - oof2)[zn].mean()
        print(f"  STRUCTURE-ONLY + Voronoi volume: predicted {d_pred2:+.3f} ({100*d_pred2/d_raw:.0f}% of the observed {d_raw:+.3f}); residual {d_res2:+.3f}")
        rows.append(dict(analysis="structure_only_with_volume", name="structure-only + V_vor (out-of-fold)", n_features=len(so2), diff_observed=d_raw, diff_predicted=d_pred2, diff_residual=d_res2, explained_fraction=d_pred2 / d_raw))
    rows.append(dict(analysis="structure_only", name="structure-only (out-of-fold)", n_features=len(STRUCT_ONLY), r2_test=r2F, diff_observed=d_raw, diff_predicted=d_pred, diff_residual=d_res, explained_fraction=d_pred / d_raw))

    # ---- inside each species
    wcols = [c for c in BASE + EXTRA_E if c != "own_type_Ca"]; wbase = [c for c in BASE if c != "own_type_Ca"]
    print("\n  WITHIN-SPECIES (identity held fixed):")
    for lab, mask in (("Ca atoms", ca), ("Zn atoms", zn)):
        sub = df[mask].reset_index(drop=True); a, b = train_test_split(np.arange(len(sub)), test_size=0.25, random_state=0)
        r_b = fit_eval(sub, wbase, a, b)[1]; r_all = fit_eval(sub, wcols, a, b)[1]
        rows.append(dict(analysis="within_species", name=lab, n_features=len(wcols), r2_test=r_all, r2_base_counts_only=r_b))
        print(f"     {lab}: R2 with 3.25 A counts only = {r_b:.3f}; with all structural descriptors = {r_all:.3f}   [{time.time()-t0:.0f} s]")

    # ---- SHAP importance of the full model
    try:
        import shap
        gbE, colsE = models[model_list[-1][0]]; samp = df.loc[te, colsE].sample(min(shap_rows, len(te)), random_state=0)
        sv = np.abs(shap.TreeExplainer(gbE)(samp).values).mean(axis=0)
        print("\n  mean |SHAP| in the last model above (5000 held-out atoms):")
        for c, v in sorted(zip(colsE, sv), key=lambda t: -t[1]): print(f"     {c:14s} {v:.4f}"); rows.append(dict(analysis="shap_E", name=c, mean_abs_shap=v))
    except Exception as e: print("  (SHAP skipped:", e, ")")
    if vol is not None:
        print("\n  NOTE (volume run): the Voronoi volume differs systematically between Ca and Zn (see TARGET line), so as a descriptor it"
              "\n  can itself act as a proxy for chemical identity; read 'explained fraction' with that in mind.")
    print("\n  INTERPRETATION CAVEAT: local structure and chemical identity are strongly collinear here (Ca atoms sit in measurably different"
          "\n  environments), so a large 'explained fraction' shows that structure is a PROXY for Ca, not that it is the mechanism; a small one"
          "\n  would mean identity matters beyond the measured structure.  No Voronoi volume was available, so size effects are only proxied.")
    if out: pd.DataFrame(rows).to_csv(out, index=False); print(f"\nwrote {out}")
    return rows

def _selftest():
    print("SELF-TEST: synthetic glass-like system with a KNOWN stress law (identity effect + q6 effect + noise)")
    rng = np.random.default_rng(3); n = 9000; L = (n / 0.0439) ** (1 / 3); origin = np.array([11.0, -4.0, 60.0])
    pos = rng.uniform(0, L, (n, 3)) + origin; box = np.column_stack([origin, origin + L])
    types = rng.choice([MG_TYPE, ZN_TYPE, CA_TYPE], n, p=[0.60, 0.35, 0.05])
    F, _ = structural_features(types, pos, box); q6 = np.nan_to_num(F.q6.values, nan=np.nanmean(F.q6.values))
    log_vm = 0.7 * (types == CA_TYPE) + 4.0 * (q6 - q6.mean()) + rng.normal(0, 0.3, n); vm = np.exp(log_vm)
    tmp = tempfile.mkdtemp(); path = os.path.join(tmp, "syn.dump")
    with open(path, "w") as f:
        for ts in (20000, 25000):
            f.write(f"ITEM: TIMESTEP\n{ts}\nITEM: NUMBER OF ATOMS\n{n}\nITEM: BOX BOUNDS pp pp pp\n")
            f.write(f"{origin[0]} {origin[0]+L}\n{origin[1]} {origin[1]+L}\n{origin[2]} {origin[2]+L}\nITEM: ATOMS id type x y z " + " ".join(STRESS) + "\n")
            for k in rng.permutation(n):
                f.write(f"{k+1} {types[k]} {pos[k,0]:.5f} {pos[k,1]:.5f} {pos[k,2]:.5f} 0 0 0 {vm[k]/np.sqrt(3):.6f} 0 0\n")   # von Mises = sqrt(3)*|sxy|
    idx = fl.index_dump(path); info, s, _ = fl.frame_at_strain(idx, 0.005); fr = fl.read_frame(path, info)
    rows = analyze(fr, shap_rows=2000)
    R = {r["name"]: r for r in rows if r["analysis"] == "model"}
    rA, rB, rE = [R[m[0]]["r2_test"] for m in (MODELS[0], MODELS[1], MODELS[-1])]
    so = [r for r in rows if r["analysis"] == "structure_only"][0]; shap_rows_ = sorted([r for r in rows if r["analysis"] == "shap_E"], key=lambda r: -r["mean_abs_shap"])
    print(f"\nchecks: R2 A={rA:.3f} B={rB:.3f} E={rE:.3f}; explained fraction {so['explained_fraction']:.2f}; top SHAP feature {shap_rows_[0]['name']}")
    assert rB - rA > 0.25, "the planted q6 effect must raise R2 clearly"
    assert rE > rB - 0.03
    assert so["explained_fraction"] < 0.35, "identity was planted independently of structure, so structure should explain little of the Ca-Zn gap"
    assert shap_rows_[0]["name"] == "q6", "the planted q6 effect should be the top SHAP feature"
    print("SELF-TEST PASSED")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("dump", nargs="?"); ap.add_argument("strain", nargs="?", type=float)
    ap.add_argument("--erate_per_ps", type=float, default=0.001); ap.add_argument("--out"); ap.add_argument("--expect_r2A", type=float)
    ap.add_argument("--volumes", help="volumes.dump from voronoi_py.py: divide stress by Voronoi volume and add volume as a descriptor")
    ap.add_argument("--selftest", action="store_true"); a = ap.parse_args()
    if a.selftest: _selftest(); sys.exit(0)
    if not a.dump or a.strain is None: ap.error("give a dump path and a target strain, or --selftest")
    idx = fl.index_dump(a.dump); info, s, _ = fl.frame_at_strain(idx, a.strain, a.erate_per_ps)
    print(f"{a.dump}: {len(idx)} frames; target strain {a.strain:.3f} -> frame at {s:.3f} (timestep {info['timestep']})")
    fr = fl.read_frame(a.dump, info)
    vol = load_volumes(a.volumes, info["timestep"], fr["id"]) if a.volumes else None
    analyze(fr, expect_r2A=None if a.volumes else a.expect_r2A, out=a.out, vol=vol)
