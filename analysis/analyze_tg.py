"""
analyze_tg.py -- glass-transition temperature (Tg), fictive temperature (Tf) and excess volume from the
thermo output of a LAMMPS melt-quench log (columns: Step Temp PotEng Press Volume Density).

Method (standard dilatometric / energetic construction):
  * the cooling ramp's specific volume v = V/N and potential energy per atom e = PE/N are plotted against T;
  * a CONTINUOUS TWO-LINE model (glass line below the break, liquid line above) is fitted by least squares,
    with the break temperature searched on a 1 K grid; the break is Tg;
  * because the liquid branch is curved over a wide range, the fit is repeated for several upper bounds of
    the fitting window and the spread is reported as the uncertainty (a single hand-picked window would hide it);
  * fictive temperature of a glass with property P at T0:  Tf = (P(T0) - b_g*T0 - a_l) / (b_l - b_g), where
    a_l,b_l are the liquid-line intercept/slope and b_g the glass-line slope (= where the glass line,
    extended through the glass's own value at T0, meets the equilibrium-liquid line).
Blocks are classified from the data themselves (ramp up / hold / ramp down), not assumed from their order.
"""
import re, sys, numpy as np, pandas as pd

def parse_thermo(path):
    blocks, cur, cols, n_atoms = [], None, None, None
    for line in open(path, errors="replace"):
        s = line.strip()
        if s.startswith("Step") and "Temp" in s and "Volume" in s and "PotEng" in s:
            cols, cur = s.split(), []
            continue
        if cur is not None:
            if s.startswith("Loop time") or s.startswith("ERROR"):
                m = re.search(r"with (\d+) atoms", s); n_atoms = int(m.group(1)) if m else n_atoms
                if cur: blocks.append(dict(df=pd.DataFrame(cur, columns=cols), n=n_atoms))
                cur = None; continue
            p = s.split()
            if len(p) == len(cols):
                try: cur.append([float(x) for x in p])
                except ValueError: pass
    if cur: blocks.append(dict(df=pd.DataFrame(cur, columns=cols), n=n_atoms))   # truncated final block
    N = next((b["n"] for b in blocks if b["n"]), None)
    for b in blocks:
        b["n"] = b["n"] or N
        d = b["df"]; k = max(5, len(d) // 20)
        t0, t1 = d.Temp.iloc[k:2*k].median(), d.Temp.iloc[-k:].median()   # skip the first rows (velocity set-up)
        b["T0"], b["T1"] = t0, t1
        b["kind"] = "ramp_up" if t1 - t0 > 100 else "ramp_down" if t0 - t1 > 100 else "hold"
        b["v"] = d.Volume.values / b["n"]; b["e"] = d.PotEng.values / b["n"]; b["T"] = d.Temp.values; b["step"] = d.Step.values
    return blocks

def bilinear(T, y, tb_grid):
    """continuous piecewise-linear fit with a free break; returns dict(Tg, c, s_glass, s_liquid, sse)"""
    best = None
    for tb in tb_grid:
        X = np.column_stack([np.ones_like(T), np.minimum(T - tb, 0.0), np.maximum(T - tb, 0.0)])
        coef, res, *_ = np.linalg.lstsq(X, y, rcond=None)
        sse = float(np.sum((X @ coef - y) ** 2))
        if best is None or sse < best["sse"]:
            best = dict(Tg=tb, c=coef[0], s_glass=coef[1], s_liquid=coef[2], sse=sse)
    return best

def tg_with_sensitivity(T, y, t_lo=300.0, t_hi_list=(700, 800, 900, 1000, 1100), tb_range=(330, 650)):
    out = []
    for t_hi in t_hi_list:
        m = (T >= t_lo) & (T <= t_hi)
        if m.sum() < 50: continue
        fit = bilinear(T[m], y[m], np.arange(*tb_range, 1.0)); fit["t_hi"] = t_hi
        # a break pinned to the edge of the search range means "no interior transition found"
        fit["edge"] = fit["Tg"] <= tb_range[0] + 2 or fit["Tg"] >= tb_range[1] - 2
        out.append(fit)
    return out

def summarize(fits):
    tg = np.array([f["Tg"] for f in fits]); return tg.mean(), tg.min(), tg.max()

def fictive(P_glass_T0, T0, fit):
    """Tf from the glass's own value at T0 and one two-line fit (liquid intercept/slope, glass slope)."""
    s_l, s_g, c, tb = fit["s_liquid"], fit["s_glass"], fit["c"], fit["Tg"]
    a_l = c - s_l * tb
    return (P_glass_T0 - s_g * T0 - a_l) / (s_l - s_g)


COMPOSITIONS = {1:(.35,.05), 2:(.30,.04), 3:(.40,.06), 4:(.25,.05), 5:(.35,.03), 6:(.35,.07), 7:(.30,.06), 8:(.40,.04), 9:(.32,.05), 10:(.38,.05)}
_AMU = dict(Mg=24.305, Zn=65.38, Ca=40.078)
def infer_composition(blocks):
    """mean atomic mass = density*volume/(N*1.66054) from the last 300 K settle, matched to the study's ten compositions.
    Lets every log be checked against the composition it is supposed to belong to."""
    holds = [b for b in blocks if b["kind"] == "hold" and b["T0"] < 350]
    if not holds: return None
    d = holds[-1]["df"]; N = holds[-1]["n"]
    M = d.Density.values[-100:].mean() * d.Volume.values[-100:].mean() / (N * 1.66054)
    mm = lambda z, c: (1 - z - c) * _AMU["Mg"] + z * _AMU["Zn"] + c * _AMU["Ca"]
    k = min(COMPOSITIONS, key=lambda j: abs(mm(*COMPOSITIONS[j]) - M))
    return f"#{k} (xZn={COMPOSITIONS[k][0]}, xCa={COMPOSITIONS[k][1]}); mean mass {M:.2f} vs {mm(*COMPOSITIONS[k]):.2f} amu"

def _role_blocks(B):
    """pick the stages out of the parsed blocks from their temperatures, not their order"""
    r = dict(quench=None, settle1=None, heat=None, hold600=None, recool=None, settle2=None)
    for i, b in enumerate(B):
        if b["kind"] == "ramp_down" and b["T0"] > 1000 and r["quench"] is None: r["quench"] = i
    q = r["quench"]
    for i, b in enumerate(B):
        if q is not None and i <= q: continue
        if b["kind"] == "hold" and b["T0"] < 350 and r["settle1"] is None: r["settle1"] = i
        elif b["kind"] == "ramp_up" and 500 < b["T1"] < 700 and r["heat"] is None: r["heat"] = i
        elif b["kind"] == "hold" and 550 < b["T0"] < 650 and r["hold600"] is None: r["hold600"] = i
        elif b["kind"] == "ramp_down" and b["T0"] < 700 and r["recool"] is None: r["recool"] = i
    holds300 = [i for i, b in enumerate(B) if b["kind"] == "hold" and b["T0"] < 350 and (q is not None and i > q)]
    if len(holds300) > 1: r["settle2"] = holds300[-1]
    return r

def _rate(b, lo, hi, dt_ps):
    m = (b["T"] > lo) & (b["T"] < hi)
    return -np.polyfit(b["step"][m], b["T"][m], 1)[0] / dt_ps * 1e12          # K/s (positive = cooling)

def analyze_log(path, label=None):
    label = label or path
    B = parse_thermo(path); R = _role_blocks(B)
    dt_ps = 0.0005
    for line in open(path, errors="replace"):
        m = re.match(r"\s*timestep\s+([0-9.eE+-]+)", line)
        if m: dt_ps = float(m.group(1)); break
    out = dict(label=label, n_atoms=B[0]["n"], blocks_found=len(B), inferred_composition=infer_composition(B))
    nan = float("nan")
    q = B[R["quench"]] if R["quench"] is not None else None
    if q is None: out["note"] = "no 1500->300 K quench block found"; return out, B, R
    out["quench_rate_Ks"] = _rate(q, 400, 1200, dt_ps)
    for key, nm in (("v", "vol"), ("e", "PE")):
        f = tg_with_sensitivity(q["T"], q[key], 300.0, (700, 800, 1000), (330, 650)); m, lo, hi = summarize(f)
        out[f"Tg_quench_{nm}"], out[f"Tg_quench_{nm}_min"], out[f"Tg_quench_{nm}_max"] = m, lo, hi
        for hi_w in (560, 580, 600):
            out[f"Tg_quench_{nm}_w{hi_w}"] = tg_with_sensitivity(q["T"], q[key], 300.0, (hi_w,), (330, 540))[0]["Tg"]
    if R["heat"] is not None:
        h = B[R["heat"]]
        for key, nm in (("v", "vol"), ("e", "PE")):
            f = tg_with_sensitivity(h["T"], h[key], 300.0, (580, 595), (340, 560)); out[f"Tg_heat_{nm}"] = summarize(f)[0] if f else nan
    if R["recool"] is not None:
        c2 = B[R["recool"]]; out["recool_rate_Ks"] = _rate(c2, 350, 550, dt_ps)
        for key, nm in (("v", "vol"), ("e", "PE")):
            for hi_w in (560, 580, 600):
                out[f"Tg_recool_{nm}_w{hi_w}"] = tg_with_sensitivity(c2["T"], c2[key], 300.0, (hi_w,), (330, 540))[0]["Tg"]
    if R["hold600"] is not None:
        hd = B[R["hold600"]]; e = hd["e"]; t = np.arange(len(e)) * 1000 * dt_ps
        out["hold600_PE_plateau"] = e[-200:].mean()
        out["hold600_PE_drift_meV_per_ns"] = np.polyfit(t[len(e)//2:], e[len(e)//2:], 1)[0] * 1000 * 1000
        out["hold600_max_row_jump_meV"] = np.abs(np.diff(e)).max() * 1000
        f = tg_with_sensitivity(q["T"], q["e"], 300.0, (700,), (330, 650))[0]
        out["hold600_below_ramp_liquid_line_meV"] = (e[-200:].mean() - (f["c"] + f["s_liquid"] * (600 - f["Tg"]))) * 1000
    if R["settle1"] is not None and R["settle2"] is not None:
        s1, s2 = B[R["settle1"]], B[R["settle2"]]
        for key, nm in (("e", "PE"), ("v", "vol")):
            a, b_ = s1[key][-100:], s2[key][-100:]
            out[f"d{nm}_final_minus_first"] = b_.mean() - a.mean()
            out[f"d{nm}_se"] = float(np.hypot(a.std(ddof=1), b_.std(ddof=1)) / np.sqrt(100))
            tf1, tf2 = [], []
            for hi_w in (600, 700, 800, 1000):
                f = tg_with_sensitivity(q["T"], q[key], 300.0, (hi_w,), (330, 650))[0]
                tf1.append(fictive(a.mean(), 300.0, f)); tf2.append(fictive(b_.mean(), 300.0, f))
            out[f"Tf_first_{nm}"], out[f"Tf_final_{nm}"] = np.mean(tf1), np.mean(tf2)
            out[f"Tf_first_{nm}_range"] = f"{min(tf1):.0f}-{max(tf1):.0f}"; out[f"Tf_final_{nm}_range"] = f"{min(tf2):.0f}-{max(tf2):.0f}"
        out["density_first"] = s1["df"].Density.values[-100:].mean(); out["density_final"] = s2["df"].Density.values[-100:].mean()
    return out, B, R

def print_summary(o):
    g = lambda k, f="{:.0f}": (f.format(o[k]) if k in o and o[k] == o[k] else "n/a")
    print(f"\n=== {o['label']}  ({o.get('n_atoms')} atoms) ===\n  composition inferred from density x volume: {o.get('inferred_composition')}")
    if "note" in o: print("  ", o["note"]); return
    print(f"  effective quench rate through the ramp : {o['quench_rate_Ks']:.3e} K/s" + (f"   (re-cool {o['recool_rate_Ks']:.3e} K/s)" if 'recool_rate_Ks' in o else ""))
    print(f"  Tg, first quench   volume {g('Tg_quench_vol')} K [{g('Tg_quench_vol_min')}-{g('Tg_quench_vol_max')}]   PE {g('Tg_quench_PE')} K [{g('Tg_quench_PE_min')}-{g('Tg_quench_PE_max')}]   (fit windows 300-700/800/1000 K)")
    print(f"  Tg, same 300-600 K window: quench {g('Tg_quench_vol_w600')}/{g('Tg_quench_PE_w600')}  re-cool {g('Tg_recool_vol_w600')}/{g('Tg_recool_PE_w600')}  reheat(3e12 K/s) {g('Tg_heat_vol')}/{g('Tg_heat_PE')}   [vol/PE]")
    print(f"  600 K hold: PE plateau {g('hold600_PE_plateau','{:.5f}')} eV/atom, drift {g('hold600_PE_drift_meV_per_ns','{:+.2f}')} meV/atom/ns, largest one-row PE jump {g('hold600_max_row_jump_meV','{:.2f}')} meV/atom, {g('hold600_below_ramp_liquid_line_meV','{:+.1f}')} meV/atom vs the quench's liquid line at 600 K")
    print(f"  300 K glass, final minus first: dPE {g('dPE_final_minus_first','{:+.2f}') if False else ('{:+.2f}'.format(1000*o['dPE_final_minus_first']) if 'dPE_final_minus_first' in o else 'n/a')} +/- {('%.2f' % (1000*o['dPE_se'])) if 'dPE_se' in o else 'n/a'} meV/atom;  dV {('{:+.4f}'.format(o['dvol_final_minus_first'])) if 'dvol_final_minus_first' in o else 'n/a'} A^3/atom")
    print(f"  fictive temperature (PE): first {g('Tf_first_PE')} K [{o.get('Tf_first_PE_range','n/a')}] -> final {g('Tf_final_PE')} K [{o.get('Tf_final_PE_range','n/a')}];   (volume): first {g('Tf_first_vol')} K -> final {g('Tf_final_vol')} K")

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("logs", nargs="+", help="LAMMPS log files; optionally LABEL=PATH")
    ap.add_argument("--out", default="tg_summary.csv")
    a = ap.parse_args()
    rows = []
    for item in a.logs:
        label, path = item.split("=", 1) if "=" in item else (item, item)
        o, B, R = analyze_log(path, label); print_summary(o); rows.append(o)
    pd.DataFrame(rows).to_csv(a.out, index=False); print(f"\nwrote {a.out}")
