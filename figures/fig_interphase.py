import sys; sys.path.insert(0, "/tmp/w/fig15")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats
import style as S
S.apply()

def load(f):
    d = pd.read_csv("/mnt/user-data/uploads/" + f)
    d["radius"] = d["radius"].astype(str)
    return d

B5, B15 = load("vf1_interphase_all.csv"), load("vf3_interphase_all.csv")

def summ(d, rad):
    s = d[d.radius == rad]; out = []
    for st in sorted(s.strain_target.unique()):
        v = s[s.strain_target == st].ratio_clean.values
        m = v.mean(); h = stats.t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
        out.append((st, m, h))
    return np.array(out)

fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 3.2), gridspec_kw={"width_ratios": [1.15, 1]})
c5, c15 = S.VF5[1], S.VF5[3]
a.axvspan(0.067, 0.069, color="#E6E6E6", lw=0, zorder=0)
a.text(0.0695, 0.93, "peak", color="#555", fontsize=7.5, va="bottom")
a.axhline(1, color="#222", lw=0.8)
for d, rad, col, mk, ls, lab, ms in ((B5, "12", c5, "o", "-", "B5 (59 seeds), nominal shell", 5),
                                      (B15, "12", c15, "s", "-", "B15 (177 seeds), nominal shell", 5),
                                      (B5, "found", c5, "o", "--", "B5, recovered radius (9.7 Å)", 3.5)):
    r = summ(d, rad)
    a.errorbar(r[:, 0], r[:, 1], yerr=r[:, 2], color=col, marker=mk, ls=ls, ms=ms, capsize=2.5, lw=1.3,
               alpha=0.55 if rad == "found" else 1, label=lab)
a.set_xlim(0.01, 0.11); a.set_ylim(0.9, 1.95); a.set_xticks([0.02, 0.04, 0.06, 0.08, 0.10])
a.set_xlabel("Engineering strain"); a.set_ylabel("Interphase / bulk D$^2_\\mathregular{min}$ ratio")
a.legend(loc="upper right", handlelength=1.8)
S.panel(a, "a", "All seeds, every strain")

s = B5[B5.radius == "12"]
p4 = s[s.strain_target == 0.04].sort_values("seed").ratio_clean.values
p7 = s[s.strain_target == 0.07].sort_values("seed").ratio_clean.values
for u, v in zip(p4, p7):
    b.plot([0, 1], [u, v], color="#9E9E9E", lw=0.6, alpha=0.7)
b.plot([0, 1], [p4.mean(), p7.mean()], color=c5, lw=3, marker="o", ms=7, zorder=5)
b.axhline(1, color="#222", lw=0.8)
pt = stats.ttest_rel(p4, p7).pvalue
b.text(0.05, 0.97, "mean change %s\n(paired t, p = %.3f)\n%d of %d seeds decrease"
       % (("%+.2f" % (p7 - p4).mean()).replace("-", "\u2212"), pt, int((p7 < p4).sum()), len(p4)),
       transform=b.transAxes, va="top", fontsize=7.5)
b.set_xticks([0, 1]); b.set_xticklabels(["4%\n(before peak)", "7%\n(at peak)"]); b.set_xlim(-0.15, 1.15)
b.set_ylim(0.6, 2.4); b.set_ylabel("Per-seed ratio")
S.panel(b, "b", "The 59 seeds of B5")
fig.tight_layout(w_pad=2.0)
S.save(fig, "fig10_interphase")
print("ok", (p7 - p4).mean(), pt)
