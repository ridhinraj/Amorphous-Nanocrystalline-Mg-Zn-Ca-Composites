#!/usr/bin/env python3
"""
d2min_sensitivity.py -- are the D2min cluster-population conclusions robust to the analysis parameters?
(manuscript Sections 3.5 / 4.3 / 4.6: "sensitivity to cutoff, threshold and clustering parameters was not tested")

For each configuration (one dump per case) the INCREMENTAL D2min cluster count is recomputed at several plateau frames
while ONE parameter at a time is varied from the baseline used in the paper:
    cutoff      3.5 / 4.0 / [4.5] / 5.0 / 5.5 A      (neighbor radius of the affine fit)
    mult        3 / [4] / 5                           (threshold = probe mean + mult x probe std)
    eps         4 / [6] / 8 A                         (DBSCAN radius)
    min_samples 10 / [15] / 25                        (DBSCAN core size)
The BASELINE row must reproduce the numbers already in the paper (printed alongside, from the earlier CSVs).  When the three
volume-fraction cases (vf_1, vf_2, vf_3) are run together, the script also reports, for every setting, whether the monotone
ordering  vf_1 > vf_2 > vf_3  survives and the ratio vf_3/vf_1; for the composition cases (run_5, run_1, run_6) it reports the spread.

Usage:  python3 d2min_sensitivity.py vf_1=vf_1/tension_trajectory.dump vf_2=vf_2/tension_trajectory.dump vf_3=vf_3/tension_trajectory.dump --out sens_vf.csv
        python3 d2min_sensitivity.py --selftest
"""
import argparse, sys, os, tempfile, time
import numpy as np, pandas as pd
from sklearn.cluster import DBSCAN
import fastlib as fl

BASE = dict(cutoff=4.5, mult=4.0, eps=6.0, min_samples=15)
OFAT = dict(cutoff=[3.5, 4.0, 4.5, 5.0, 5.5], mult=[3.0, 4.0, 5.0], eps=[4.0, 6.0, 8.0], min_samples=[10, 15, 25])
PLATEAU = (0.09, 0.10, 0.11, 0.12, 0.13, 0.14)
# mean incremental cluster count over exactly these six strains, from the earlier loop-based analysis (baseline settings)
EXPECTED = {"run_1": 115.5, "run_5": 131.0, "run_6": 126.17, "vf_1": 87.17, "vf_2": 75.17, "vf_3": 55.0}

def configs():
    seen, out = set(), []
    for name, vals in OFAT.items():
        for v in vals:
            c = dict(BASE); c[name] = v; key = tuple(sorted(c.items()))
            if key not in seen: seen.add(key); out.append((("baseline" if c == BASE else f"{name}={v:g}"), c))
    return sorted(out, key=lambda t: (t[0] != "baseline", t[0]))

def count_clusters(coords, eps, min_samples):
    if len(coords) < min_samples: return 0, 0
    lab = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(coords); real = lab[lab >= 0]
    return (len(set(real)), int(np.bincount(real).max())) if len(real) else (0, 0)

def run_case(label, path, erate, plateau, quiet=False):
    t0 = time.time(); idx = fl.index_dump(path)
    pairs = [(plateau_s, *fl.frame_at_strain(idx, plateau_s, erate)[::2]) for plateau_s in plateau]       # (strain, k) per plateau frame
    need = sorted({0, 1} | {k for _, k in [(s, k) for s, _, k in pairs]} | {k - 1 for _, _, k in pairs})
    fr = {k: fl.read_frame(path, idx[k]) for k in need}
    print(f"[{label}] {len(idx)} frames; reading {len(need)} frames took {time.time()-t0:.0f} s", flush=True)
    rows = []
    cache = {}
    for cutoff in sorted({c["cutoff"] for _, c in configs()}):
        probe = fl.d2min(fr[0]["pos"], fr[0]["box"], fr[1]["pos"], fr[1]["box"], cutoff); pv = probe[~np.isnan(probe)]
        d2 = {k: fl.d2min(fr[k - 1]["pos"], fr[k - 1]["box"], fr[k]["pos"], fr[k]["box"], cutoff) for _, _, k in pairs}
        cache[cutoff] = (pv.mean(), pv.std(), d2)
    print(f"[{label}] D2min computed for {len(cache)} cutoffs x {len(pairs)+1} frame pairs ({time.time()-t0:.0f} s)", flush=True)
    for name, c in configs():
        m, s, d2 = cache[c["cutoff"]]; thr = m + c["mult"] * s; per = []
        for _, _, k in pairs:
            flagged = (d2[k] > thr) & ~np.isnan(d2[k])
            per.append(count_clusters(fr[k]["pos"][flagged], c["eps"], c["min_samples"]) + (int(flagged.sum()),))
        per = np.array(per, float)
        rows.append(dict(case=label, setting=name, **c, mean_clusters=per[:, 0].mean(), sd_clusters=per[:, 0].std(ddof=1) if len(per) > 1 else 0.0,
                         mean_largest=per[:, 1].mean(), mean_flagged=per[:, 2].mean()))
    print(f"[{label}] done in {time.time()-t0:.0f} s", flush=True)
    return rows

def summarize(df, check_expected=True):
    print("\n" + "=" * 100 + "\nMEAN SIMULTANEOUS CLUSTERS over the plateau frames, by setting  (sd = frame-to-frame)\n" + "=" * 100)
    piv = df.pivot_table(index="setting", columns="case", values="mean_clusters", sort=False).round(1)
    print(piv.to_string())
    base = df[df.setting == "baseline"].set_index("case").mean_clusters
    print("\nBASELINE reproduction of the numbers already in the paper:")
    for c in base.index:
        if check_expected and c in EXPECTED:
            ok = abs(base[c] - EXPECTED[c]) / EXPECTED[c] < 0.05
            print(f"   {c:6s} now {base[c]:7.2f}   earlier {EXPECTED[c]:7.2f}   {'OK (within 5%)' if ok else 'DIFFERS by %.0f%%' % (100*(base[c]/EXPECTED[c]-1))}")
    cases = set(df.case)
    if {"vf_1", "vf_2", "vf_3"} <= cases:
        print("\nVOLUME-FRACTION CONCLUSION (vf_1 > vf_2 > vf_3?) under each setting:")
        n_ok = 0
        for s in piv.index:
            a, b, c = piv.loc[s, "vf_1"], piv.loc[s, "vf_2"], piv.loc[s, "vf_3"]; mono = a > b > c; n_ok += mono
            print(f"   {s:16s} {a:6.1f} {b:6.1f} {c:6.1f}   monotone: {'YES' if mono else 'NO '}   vf_3/vf_1 = {c/a:.2f}")
        print(f"   -> ordering survives in {n_ok} of {len(piv)} settings")
    comp = [c for c in ("run_5", "run_1", "run_6") if c in cases]
    if len(comp) == 3:
        print("\nCOMPOSITION CONCLUSION (similar counts across Ca content?): spread (max/min) under each setting")
        for s in piv.index: v = piv.loc[s, comp]; print(f"   {s:16s} {v.round(1).tolist()}   max/min = {v.max()/v.min():.2f}")

def _make_synthetic(path, n_clusters, seed):
    rng = np.random.default_rng(seed); n = 6000; L = (n / 0.0439) ** (1 / 3); origin = np.array([21.0, -7.0, 33.0])
    base = rng.uniform(0, L, (n, 3)) + origin
    grid = np.array([[x, y, z] for x in (0.2, 0.7) for y in (0.2, 0.7) for z in (0.2, 0.7)]) * L + origin     # 8 well-separated sites
    def kick(p, k):
        p = p.copy()
        for c in grid[rng.permutation(8)[:k]]:
            v = rng.normal(size=3); v *= 1.3 / np.linalg.norm(v); p[np.linalg.norm(p - c, axis=1) < 4.0] += v
        return p
    frames = {20000: base}; frames[25000] = base + rng.normal(0, 0.03, (n, 3))
    frames[115000] = base + rng.normal(0, 0.03, (n, 3)); frames[120000] = kick(frames[115000] + rng.normal(0, 0.03, (n, 3)), n_clusters)
    frames[135000] = frames[120000] + rng.normal(0, 0.03, (n, 3)); frames[140000] = kick(frames[135000] + rng.normal(0, 0.03, (n, 3)), n_clusters)
    with open(path, "w") as f:
        for ts in sorted(frames):
            f.write(f"ITEM: TIMESTEP\n{ts}\nITEM: NUMBER OF ATOMS\n{n}\nITEM: BOX BOUNDS pp pp pp\n{origin[0]} {origin[0]+L}\n{origin[1]} {origin[1]+L}\n{origin[2]} {origin[2]+L}\nITEM: ATOMS id type x y z\n")
            for k in range(n): f.write(f"{k+1} 1 {frames[ts][k,0]:.5f} {frames[ts][k,1]:.5f} {frames[ts][k,2]:.5f}\n")

def _selftest():
    print("SELF-TEST: three synthetic cases with 6 / 4 / 2 planted rearrangement clusters per frame pair (mimicking vf_1 > vf_2 > vf_3)")
    tmp = tempfile.mkdtemp(); rows = []
    for lab, k in (("vf_1", 6), ("vf_2", 4), ("vf_3", 2)):
        p = os.path.join(tmp, lab + ".dump"); _make_synthetic(p, k, seed=k); rows += run_case(lab, p, 0.001, (0.10, 0.12))
    df = pd.DataFrame(rows); summarize(df, check_expected=False)
    b = df[df.setting == "baseline"].set_index("case").mean_clusters
    print(f"\nchecks: baseline counts {b.round(1).to_dict()} (planted 6 / 4 / 2)")
    assert abs(b["vf_1"] - 6) <= 1 and abs(b["vf_2"] - 4) <= 1 and abs(b["vf_3"] - 2) <= 1, "planted cluster counts not recovered"
    mono = [(r.setting, r.mean_clusters) for r in df.itertuples()]; piv = df.pivot_table(index="setting", columns="case", values="mean_clusters")
    assert (piv.vf_1 >= piv.vf_2).all() and (piv.vf_2 >= piv.vf_3).all(), "planted ordering should survive every setting"
    print("SELF-TEST PASSED")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cases", nargs="*", help="LABEL=dump_path"); ap.add_argument("--erate_per_ps", type=float, default=0.001)
    ap.add_argument("--plateau_strains", type=float, nargs="+", default=list(PLATEAU)); ap.add_argument("--out", default="d2min_sensitivity.csv"); ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest: _selftest(); sys.exit(0)
    if not a.cases: ap.error("give LABEL=dump cases or --selftest")
    rows = []
    for item in a.cases:
        lab, path = item.split("=", 1) if "=" in item else (os.path.basename(os.path.dirname(item)) or item, item)
        rows += run_case(lab, path, a.erate_per_ps, a.plateau_strains)
        pd.DataFrame(rows).to_csv(a.out, index=False)
    df = pd.DataFrame(rows); summarize(df); print(f"\nwrote {a.out}")
