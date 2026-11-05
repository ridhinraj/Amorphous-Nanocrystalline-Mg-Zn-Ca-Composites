#!/usr/bin/env python3
"""
interphase_all_seeds.py -- the weak-interphase test for EVERY seed of a volume-fraction composite, at several strains, with paired
pre-yield / yield comparisons.   (manuscript Section 3.4 addendum: 10 of 59 seeds at yield, 2 seeds before yield)

Pipeline (identical conventions to find_seed_centers.py + interphase_shell_analysis.py, vectorized):
  1. seeds are located in the FIRST frame (Mg-type atoms with q6 > 0.40 and Mg-neighbor fraction > 0.95, DBSCAN eps 6 A, min_samples 15);
  2. for each requested strain the incremental D2min (frame k-1 -> k, cutoff 4.5 A) is computed ONCE for all atoms;
  3. for each seed, mean D2min in the shell |r - R| <= 3 A is divided by the bulk mean.  Two bulk definitions are reported:
       ratio_orig  bulk = atoms farther than 7.5 A from THIS seed's surface (what the earlier single-seed runs used);
       ratio_clean bulk = atoms farther than 7.5 A from EVERY seed's surface (excludes other seeds' interphases from the reference).
  4. across seeds: mean, 95% CI, t-test and Wilcoxon test against 1; and PAIRED tests between strains (does the excess develop at yield?).
Radii: --radii 12 found   (12 = nominal build radius; 'found' = each seed's own recovered radius, which underestimates the true one).

Usage: python3 interphase_all_seeds.py vf_1/tension_trajectory.dump --strains 0.02 0.04 0.07 0.10 --radii 12 found --out vf1_interphase.csv
       python3 interphase_all_seeds.py --selftest
"""
import argparse, sys, os, tempfile, time
import numpy as np, pandas as pd
from scipy import stats
from sklearn.cluster import DBSCAN
import fastlib as fl

MG = 1
# per-seed ratios already obtained one seed at a time (radius 12 A, frame at 0.070, ORIGINAL bulk definition): used as a reproduction check
EARLIER = {(88.54, 41.90, 57.48): 1.29, (21.82, 58.02, 68.88): 1.18, (90.81, 52.32, 98.80): 1.16, (74.29, 56.75, 126.70): 1.05, (60.72, 99.51, 18.85): 1.32,
           (106.09, 73.57, 138.13): 1.54, (46.11, 106.23, 182.76): 1.70, (76.56, 146.95, 145.95): 1.40, (188.51, 42.26, 119.49): 1.23, (182.47, 91.50, 75.23): 0.99}

def find_seeds(fr0, q6_thr=0.40, purity_thr=0.95, cutoff=4.0, eps=6.0, min_samples=15):
    types, pos, box = fr0["type"], fr0["pos"], fr0["box"]; n = len(types)
    i, j, d, r = fl.neighbor_pairs(pos, box, cutoff)
    q6 = fl.steinhardt(n, i, d, r, ls=(6,))[6]
    ntot = np.bincount(i, minlength=n); nmg = np.bincount(i, weights=(types[j] == MG).astype(float), minlength=n)
    cand = (types == MG) & (q6 > q6_thr) & (nmg / np.maximum(ntot, 1) > purity_thr) & (ntot > 0)
    if cand.sum() < min_samples: return []
    coords = pos[cand]; lab = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(coords); out = []
    for k in sorted(set(lab) - {-1}):
        c = coords[lab == k]; ctr = c.mean(axis=0)
        out.append(dict(center=ctr, radius_found=float(np.percentile(np.linalg.norm(c - ctr, axis=1), 90)), n_candidates=len(c)))
    return out

def seed_ratios(d2, cur_pos, cur_box, seeds, radii, shell=3.0, cutoff=4.5):
    """per-seed interphase/bulk ratios for one D2min array"""
    L = cur_box[:, 1] - cur_box[:, 0]; n = len(cur_pos); valid = ~np.isnan(d2)
    signed = np.empty((len(seeds), n), np.float32)
    for s, (sd, R) in enumerate(zip(seeds, radii)):
        dv = cur_pos - sd["center"]; dv -= L * np.round(dv / L); signed[s] = np.linalg.norm(dv, axis=1) - R
    rmin = signed.min(axis=0); bulk_clean = valid & (rmin > shell + cutoff); mb = d2[bulk_clean].mean()
    rows = []; union = np.zeros(n, bool)
    for s in range(len(seeds)):
        sh = valid & (signed[s] >= -shell) & (signed[s] <= shell); union |= sh
        bo = valid & (signed[s] > shell + cutoff)
        rows.append(dict(seed=s, n_shell=int(sh.sum()), shell_mean=float(d2[sh].mean()) if sh.any() else np.nan,
                         ratio_clean=float(d2[sh].mean() / mb) if sh.any() else np.nan, ratio_orig=float(d2[sh].mean() / d2[bo].mean()) if sh.any() else np.nan))
    return rows, dict(bulk_clean_mean=float(mb), pooled_ratio=float(d2[union].mean() / mb), n_bulk_clean=int(bulk_clean.sum()))

def describe(x, label):
    x = np.asarray(x, float); x = x[~np.isnan(x)]; m, sem = x.mean(), stats.sem(x); ci = stats.t.interval(0.95, len(x) - 1, loc=m, scale=sem)
    t, p = stats.ttest_1samp(x, 1.0); w = stats.wilcoxon(x - 1.0).pvalue
    print(f"   {label:34s} n={len(x):3d} mean {m:.3f} (95% CI {ci[0]:.3f}-{ci[1]:.3f})  median {np.median(x):.3f}  {int((x > 1).sum())}/{len(x)} above 1   t-test p={p:.2g}  Wilcoxon p={w:.2g}")

def run(path, strains, radii_tokens, erate, shell=3.0, cutoff=4.5, out=None):
    t0 = time.time(); idx = fl.index_dump(path); fr0 = fl.read_frame(path, idx[0])
    seeds = find_seeds(fr0); print(f"{path}: {len(idx)} frames; {len(seeds)} seeds found ({time.time()-t0:.0f} s)", flush=True)
    if not seeds: sys.exit("no seeds found")
    print(f"   recovered radii: mean {np.mean([s['radius_found'] for s in seeds]):.2f} A (min {min(s['radius_found'] for s in seeds):.2f}, max {max(s['radius_found'] for s in seeds):.2f})")
    recs = []
    for st in strains:
        info, actual, k = fl.frame_at_strain(idx, st, erate)
        prev, cur = fl.read_frame(path, idx[k - 1]), fl.read_frame(path, idx[k]); d2 = fl.d2min(prev["pos"], prev["box"], cur["pos"], cur["box"], cutoff)
        for tok in radii_tokens:
            radii = [s["radius_found"] for s in seeds] if tok == "found" else [float(tok)] * len(seeds)
            rows, glob = seed_ratios(d2, cur["pos"], cur["box"], seeds, radii, shell, cutoff)
            for r_, sd in zip(rows, seeds):
                recs.append(dict(strain_target=st, strain_frame=actual, radius=tok, **r_, cx=sd["center"][0], cy=sd["center"][1], cz=sd["center"][2], radius_found=sd["radius_found"],
                                 pooled_ratio=glob["pooled_ratio"], bulk_clean_mean=glob["bulk_clean_mean"]))
        print(f"   strain {st:.3f} (frame {actual:.3f}) done ({time.time()-t0:.0f} s)", flush=True)
    df = pd.DataFrame(recs)
    if out: df.to_csv(out, index=False); print(f"wrote {out}")
    for tok in radii_tokens:
        print("\n" + "=" * 104 + f"\nRADIUS = {tok}\n" + "=" * 104)
        sub = df[df.radius == tok]
        for st in strains:
            g = sub[sub.strain_target == st]
            print(f" strain {st:.3f}  (frame {g.strain_frame.iloc[0]:.3f}),  pooled shell/bulk = {g.pooled_ratio.iloc[0]:.3f}")
            describe(g.ratio_clean, "per-seed ratio, clean bulk"); describe(g.ratio_orig, "per-seed ratio, original bulk")
        if len(strains) > 1:
            print(" PAIRED comparisons between strains (same seeds, clean-bulk ratio):")
            for a in range(len(strains)):
                for b in range(a + 1, len(strains)):
                    x = sub[sub.strain_target == strains[a]].ratio_clean.values; y = sub[sub.strain_target == strains[b]].ratio_clean.values
                    diff = y - x; print(f"   {strains[a]:.3f} -> {strains[b]:.3f}: mean change {diff.mean():+.3f}  paired t p={stats.ttest_rel(y, x).pvalue:.2g}  Wilcoxon p={stats.wilcoxon(diff).pvalue:.2g}  ({int((diff > 0).sum())}/{len(diff)} seeds increase)")
    # reproduction of the earlier single-seed runs
    chk = df[(df.radius == "12") & (np.isclose(df.strain_frame, 0.07, atol=0.0026))]
    if len(chk):
        print("\nREPRODUCTION of earlier single-seed results (radius 12, ~0.07 strain, original bulk definition):")
        for (cx, cy, cz), want in EARLIER.items():
            m = chk[(abs(chk.cx - cx) < 0.05) & (abs(chk.cy - cy) < 0.05) & (abs(chk.cz - cz) < 0.05)]
            if len(m): got = m.ratio_orig.iloc[0]; print(f"   seed ({cx}, {cy}, {cz}): now {got:.2f}   earlier {want:.2f}   {'OK' if abs(got - want) <= 0.015 else 'DIFFERS'}")
    return df

def _selftest():
    print("SELF-TEST: synthetic composite, 3 HCP seeds; shells of seeds 1-2 rearrange at the later ('yield') frame only")
    rng = np.random.default_rng(23); L = 72.0; origin = np.array([15.0, -20.0, 5.0]); a0, ca_ratio = 3.209, 1.619
    a1, a2, a3 = np.array([1, 0, 0]) * a0, np.array([0, np.sqrt(3), 0]) * a0, np.array([0, 0, ca_ratio]) * a0
    basis = [(0, 0, 0), (.5, .5, 0), (.5, 5/6, .5), (0, 1/3, .5)]; centers = origin + L * np.array([[.25, .25, .25], [.75, .3, .7], [.4, .75, .3]]); R = 12.0
    seedpos = []
    for c in centers:
        for i in range(-5, 6):
            for j in range(-5, 6):
                for k in range(-5, 6):
                    for f in basis:
                        p = i*a1 + j*a2 + k*a3 + f[0]*a1 + f[1]*a2 + f[2]*a3
                        if np.linalg.norm(p) < R: seedpos.append(p + c)
    seedpos = np.array(seedpos); n_m = int(0.0439 * L ** 3); mat = rng.uniform(0, L, (n_m, 3)) + origin
    from scipy.spatial import cKDTree
    mat = mat[cKDTree(seedpos).query(mat)[0] > 2.0]; mt = rng.choice([1, 2, 3], len(mat), p=[.6, .35, .05])
    pos0 = np.vstack([mat, seedpos]); types = np.concatenate([mt, np.full(len(seedpos), 1)]); n = len(pos0)
    f1 = pos0 + rng.normal(0, .03, (n, 3)); f2 = f1 + rng.normal(0, .03, (n, 3)); f3 = f2 + rng.normal(0, .03, (n, 3))
    for c in centers[:2]:
        sel = np.abs(np.linalg.norm(f3 - c, axis=1) - R) <= 3.0; v = rng.normal(size=3); f3[sel] += v / np.linalg.norm(v) * 1.0
    tmp = tempfile.mkdtemp(); path = os.path.join(tmp, "syn.dump"); frames = {20000: pos0, 25000: f1, 30000: f2, 35000: f3}
    with open(path, "w") as f:
        for ts in sorted(frames):
            f.write(f"ITEM: TIMESTEP\n{ts}\nITEM: NUMBER OF ATOMS\n{n}\nITEM: BOX BOUNDS pp pp pp\n{origin[0]} {origin[0]+L}\n{origin[1]} {origin[1]+L}\n{origin[2]} {origin[2]+L}\nITEM: ATOMS id type x y z\n")
            for k in range(n): f.write(f"{k+1} {types[k]} {frames[ts][k,0]:.5f} {frames[ts][k,1]:.5f} {frames[ts][k,2]:.5f}\n")
    df = run(path, [0.010, 0.015], ["12"], 0.001)
    hi = df[np.isclose(df.strain_frame, 0.015)].sort_values("seed"); lo = df[np.isclose(df.strain_frame, 0.010)].sort_values("seed")
    print("\nchecks: seeds found", len(hi), "| yield ratios", hi.ratio_clean.round(1).tolist(), "| pre-yield ratios", lo.ratio_clean.round(2).tolist())
    assert len(hi) == 3
    order = np.argsort([np.linalg.norm(hi[["cx", "cy", "cz"]].values[s] - centers[0]) for s in range(3)])        # identify seeds 1,2,3 by position
    kicked = [int(np.argmin([np.linalg.norm(hi[["cx", "cy", "cz"]].values[s] - centers[t]) for s in range(3)])) for t in range(3)]
    assert hi.ratio_clean.values[kicked[0]] > 5 and hi.ratio_clean.values[kicked[1]] > 5, "kicked shells must stand out"
    assert 0.5 < hi.ratio_clean.values[kicked[2]] < 1.6 and (lo.ratio_clean.values > 0.5).all() and (lo.ratio_clean.values < 1.6).all()
    # equivalence with the earlier single-seed tool (loop D2min) on the same data
    import interphase_shell_analysis as old
    cs = hi[["cx", "cy", "cz"]].values[kicked[0]]; res = old.analyze_shell(path, 0.015, tuple(cs), 12.0, shell_width=3.0, cutoff=4.5, erate_per_ps=0.001, mode="incremental")
    r_old = res["regions"]["interphase_shell"]["mean_d2min"] / res["regions"]["bulk_matrix"]["mean_d2min"]; r_new = hi.ratio_orig.values[kicked[0]]
    print(f"earlier single-seed tool: {r_old:.6f}   this script (original bulk definition): {r_new:.6f}   identical: {abs(r_old - r_new) < 1e-6 * max(1, abs(r_old))}")
    assert abs(r_old - r_new) < 1e-6 * max(1, abs(r_old))
    print("SELF-TEST PASSED")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("dump", nargs="?"); ap.add_argument("--strains", type=float, nargs="+", default=[0.02, 0.04, 0.07, 0.10])
    ap.add_argument("--radii", nargs="+", default=["12"]); ap.add_argument("--erate_per_ps", type=float, default=0.001); ap.add_argument("--out"); ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest: _selftest(); sys.exit(0)
    if not a.dump: ap.error("give a dump path or --selftest")
    run(a.dump, a.strains, [("12" if r in ("12", "12.0") else r) for r in a.radii], a.erate_per_ps, out=a.out)
