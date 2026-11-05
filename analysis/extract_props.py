import sys, os
import numpy as np, pandas as pd
from scipy.signal import savgol_filter

def load(path):
    df = pd.read_csv(path, sep=r"\s+")
    return df["Strain"].values, df["Stress_X_GPa"].values

def smooth(e, s, width=0.010, order=3):
    """Savitzky-Golay smoothing over a fixed STRAIN width (not a fixed number of points), so that runs
    at different strain rates (different row spacing) are smoothed identically."""
    dt = float(np.median(np.diff(e)))
    win = int(round(width / dt))
    win += (win % 2 == 0)
    return savgol_filter(s, max(win, order + 2 + (order + 2) % 2 + 1), order)

def props(e, s, width=0.010, E_window=(0.001, 0.010), flow_window=0.010):
    """Mechanical properties from one tensile stress-strain curve.
    yield strength = maximum of the Savitzky-Golay-smoothed curve (fixed width of 0.01 strain, cubic)
    E              = slope of a least-squares line over 0.1-1.0 % strain
    flow stress    = value at the final strain of a linear fit to the last 0.01 strain
    softening      = (yield strength - flow stress) / yield strength
    """
    sm = smooth(e, s, width)
    i_pk = int(np.argmax(sm))
    m = (e >= E_window[0]) & (e <= E_window[1])
    E = np.polyfit(e[m], s[m], 1)[0]
    w = e >= e[-1] - flow_window
    slope, icpt = np.polyfit(e[w], s[w], 1)
    flow = slope * e[-1] + icpt
    return dict(
        E_GPa=E, peak_GPa=sm[i_pk], peak_strain=e[i_pk],
        s15_fit=flow, softening_pct=100 * (sm[i_pk] - flow) / sm[i_pk],
        s15_last3=s[-3:].mean(),
        slope_end=np.polyfit(e[e >= e[-1] - 0.03], sm[e >= e[-1] - 0.03], 1)[0],   # over the last 0.03 strain
    )

CASES = (["run_%d" % i for i in range(1, 11)] +
         ["vf_1", "vf_2", "vf_3", "vf_20_1", "vf_25_1", "vf_r6_1", "vf_r6_2", "vf_r6_3",
          "vf_rad15_1", "vf_rad20_1", "t310_1", "rate1e8_1", "rate5e8_1", "rate2e9_1",
          "rate2e9_run6"])   # 25th simulation: Composition 6 at 2e9 /s (Section 3.7 addendum)

def find_files(root):
    """Accepts either layout:
       <root>/<case>/tension_stress_strain.txt   (project directory)
       <root>/<case>__tension_stress_strain.txt  (flat upload folder, e.g. ss_upload)"""
    found = {}
    for c in CASES:
        for cand in (f"{root}/{c}/tension_stress_strain.txt", f"{root}/{c}__tension_stress_strain.txt"):
            if os.path.exists(cand):
                found[c] = cand
                break
    return found

if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) > 1 else "."
    found = find_files(root)
    missing = [c for c in CASES if c not in found]
    rows = []
    for c, p in found.items():
        e, s = load(p)
        rows.append(dict(case=c, **props(e, s)))
    if not rows:
        sys.exit(f"No tension_stress_strain files found under '{root}'.")
    out = pd.DataFrame(rows).rename(columns={
        "E_GPa": "E_GPa", "peak_GPa": "yield_strength_GPa", "peak_strain": "yield_strain",
        "s15_fit": "flow_stress_15pct_GPa", "softening_pct": "softening_pct",
        "s15_last3": "mean_last3_points_GPa_legacy", "slope_end": "end_slope_GPa_per_strain"})
    out.round(4).to_csv("stress_strain_properties.csv", index=False)
    pd.set_option("display.width", 200)
    print(out.round(4).to_string(index=False))
    print(f"\n{len(rows)} of {len(CASES)} cases processed; wrote stress_strain_properties.csv")
    if missing:
        print("Missing:", ", ".join(missing))
