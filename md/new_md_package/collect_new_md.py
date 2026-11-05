#!/usr/bin/env python3
"""
collect_new_md.py -- yield strength, modulus and softening for every finished new case, with the SAME definitions as
extract_props.py (peak of the Savitzky-Golay curve, 0.010 strain width, cubic), plus group means / SDs and the
matched-geometry calcium comparison.

    python3 collect_new_md.py [stress_strain_properties_all25.csv]
Needs cases_manifest.csv (from setup_new_md.py) and <case>/tension_stress_strain.txt. Writes new_md_results.csv and
new_md_groups.csv.  The optional csv (the 25 existing runs) lets existing cases join their groups.
"""
import os, sys, numpy as np, pandas as pd
from scipy.signal import savgol_filter

def smooth(e, s, width=0.010, order=3):
    dt = float(np.median(np.diff(e))); win = int(round(width / dt)); win += (win % 2 == 0)
    return savgol_filter(s, max(win, order + 2 + (order + 2) % 2 + 1), order)

def props(e, s, width=0.010):
    sm = smooth(e, s, width); i = int(np.argmax(sm))
    m = (e >= 0.001) & (e <= 0.010); E = np.polyfit(e[m], s[m], 1)[0]
    w = e >= e[-1] - 0.010; sl, ic = np.polyfit(e[w], s[w], 1); flow = sl * e[-1] + ic
    return dict(E_GPa=E, yield_strength_GPa=sm[i], yield_strain=e[i], flow_stress_15pct_GPa=flow,
                softening_pct=100 * (sm[i] - flow) / sm[i], final_strain=e[-1])

COL = {"x": "Stress_X_GPa", "y": "Stress_Y_GPa", "z": "Stress_Z_GPa"}
# existing runs that belong to a new group (matched geometry: 12 A seeds)
EXISTING = {"r1_vf5": "vf_1", "r1_vf10": "vf_2", "r1_vf15": "vf_3", "r6_vf5": "vf_r6_1", "r6_vf10": "vf_r6_2",
            "r6_vf15": "vf_r6_3"}

man = pd.read_csv("cases_manifest.csv"); rows = []
for _, r in man.iterrows():
    p = f"{r.case}/tension_stress_strain.txt"
    if not os.path.exists(p): continue
    df = pd.read_csv(p, sep=r"\s+")
    if len(df) < 50: continue
    d = props(df["Strain"].values, df[COL[r.axis]].values)
    rows.append(dict(r.to_dict(), **d, complete=bool(d["final_strain"] >= 0.1495)))
res = pd.DataFrame(rows)
if res.empty: sys.exit("no finished cases found")
res.round(4).to_csv("new_md_results.csv", index=False)
pd.set_option("display.width", 220)
print(res[["case", "group", "xCa", "vf", "yield_strength_GPa", "E_GPa", "softening_pct", "complete"]].round(4).to_string(index=False))
if not res["complete"].all(): print("\nWARNING: incomplete curves:", ", ".join(res[~res.complete].case))

ex = pd.read_csv(sys.argv[1]).set_index("case") if len(sys.argv) > 1 else None
out = []
for g, d in res.groupby("group"):
    vals = list(d.yield_strength_GPa)
    if ex is not None and g in EXISTING and EXISTING[g] in ex.index: vals.append(float(ex.loc[EXISTING[g], "yield_strength_GPa"]))
    out.append(dict(group=g, n=len(vals), mean_GPa=np.mean(vals), sd_GPa=np.std(vals, ddof=1) if len(vals) > 1 else np.nan,
                    cv_pct=100*np.std(vals, ddof=1)/np.mean(vals) if len(vals) > 1 else np.nan, cases=" ".join(d.case)))
G = pd.DataFrame(out).round(4); G.to_csv("new_md_groups.csv", index=False)
print("\nGroup means (existing matched-geometry run included; r*_vf10 mixes placement and independent cases -- see the split below):"); print(G.drop(columns="cases").to_string(index=False))

# matched-geometry calcium comparison at 12 A / 10 %.
# Placement set  : seed-placement replicates on ONE glass per composition (rep_*, plus the matched-series case and the existing run).
# Independent set: independent quenches (q*), with the original-glass 10 % case counted once as realization 1.
from scipy import stats
ORIG10 = {5: "vf_r5_2", 1: "vf_2", 6: "vf_r6_2"}          # original-glass 12 A / 10 % case of each composition
def yv(case):
    if case in set(res.case): return [float(res.loc[res.case == case, "yield_strength_GPa"].iloc[0])]
    if ex is not None and case in ex.index: return [float(ex.loc[case, "yield_strength_GPa"])]
    return []
def sets(run, indep):
    g = res[res.group == f"r{run}_vf10"]
    if indep: return yv(ORIG10[run]) + [float(x) for x in g[g.case.str.startswith("q")].yield_strength_GPa]
    return yv(ORIG10[run]) + [float(x) for x in g[g.case.str.startswith("rep_")].yield_strength_GPa]
for label, indep in (("PLACEMENT replicates (one glass per composition)", False), ("INDEPENDENT glasses (original + independent quenches)", True)):
    print(f"\n--- Matched-geometry Ca comparison, 12 A / 10 %: {label}")
    V = {run: sets(run, indep) for run in (5, 1, 6)}
    for run, v in V.items():
        print(f"  Composition {run} (xCa={ {5:0.03,1:0.05,6:0.07}[run] }): n={len(v)}" + (f"  mean={np.mean(v):.4f} GPa  sd={np.std(v, ddof=1):.4f}" if len(v) > 1 else ""))
    for a, b in ((5, 1), (1, 6), (5, 6)):
        if len(V[a]) > 1 and len(V[b]) > 1:
            t = stats.ttest_ind(V[b], V[a], equal_var=False)
            print(f"  Composition {a} -> {b} (Welch): {100*(np.mean(V[b])/np.mean(V[a])-1):+.1f} %  p={t.pvalue:.3g}  (n={len(V[a])},{len(V[b])})")
print("\nNote: with n of 3-5 per group report effect size and SD; p-values are indicative only.")
