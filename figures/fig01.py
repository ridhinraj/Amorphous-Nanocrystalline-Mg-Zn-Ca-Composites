import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
"""Figure 1 - study workflow. 6.5 in wide, Arial (Liberation Sans fallback), all text >= 8 pt, 600 dpi."""
import logging; logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams["font.family"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
FS = 8
W, H = 6.5, 6.55
fig = plt.figure(figsize=(W, H)); ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
BLUE, GREY, GREEN, ORANGE, PURPLE = "#444444", "#8C8C8C", "#222222", "#555555", "#000000"
FILL = {BLUE: "#E9E9E9", GREY: "#F5F5F5", GREEN: "#DADADA", ORANGE: "#F0F0F0", PURPLE: "#FFFFFF"}

def box(x, y, w, h, col, title, sec, body, lw=1.4):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.06",
                                fc=FILL[col], ec=col, lw=lw))
    ax.text(x + w/2, y + h - 0.12, title, ha="center", va="top", fontsize=FS, fontweight="bold")
    ax.text(x + w/2, y + h - 0.33, sec, ha="center", va="top", fontsize=FS, style="italic", color="#444")
    ax.text(x + 0.07, y + h - 0.54, body, ha="left", va="top", fontsize=FS, linespacing=1.25)

def arrow(p, q, col="#333", ls="-"):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=9, lw=1.2, color=col, linestyle=ls,
                                 shrinkA=0, shrinkB=0))

# ---- (a) model construction
ax.text(0.05, 6.35, "(a) Model construction", fontsize=FS + 1, fontweight="bold", va="center")
bw, gap, x0, ya, ha = 1.50, 0.133, 0.05, 4.75, 1.40
a = [("Interatomic\npotential", "Sec. 2.1", "2NN MEAM for\nMg-Zn-Ca with ZBL\ncore repulsion\n(LAMMPS)"),
     ("Amorphous\nmatrix", "Sec. 2.2", "about 390,000 atoms\n10 compositions\nMelt, quench at\n10$^{12}$ K/s, 600 K\nreheat and re-cool"),
     ("Nanocrystal\ncomposites", "Sec. 2.3", "HCP $\\alpha$-Mg seeds:\nnone (pure glass),\none seed (15-35 Å)\nor many 12 Å seeds\n(5-18.9 vol%)"),
     ("Relaxation\nand tension", "Sec. 2.4", "Minimization, NPT\nat 300 K\nTension at 10$^9$ s$^{-1}$\nto 15% strain")]
for i, (t, s, b) in enumerate(a):
    x = x0 + i*(bw + gap)
    ax.add_patch(FancyBboxPatch((x, ya), bw, ha, boxstyle="round,pad=0,rounding_size=0.06", fc=FILL[BLUE], ec=BLUE, lw=1.4))
    ax.text(x + bw/2, ya + ha - 0.10, t, ha="center", va="top", fontsize=FS, fontweight="bold", linespacing=1.15)
    ax.text(x + bw/2, ya + ha - 0.46, s, ha="center", va="top", fontsize=FS, style="italic", color="#444")
    ax.text(x + 0.07, ya + ha - 0.65, b, ha="left", va="top", fontsize=FS, linespacing=1.2)
    if i < 3: arrow((x + bw + 0.01, ya + ha/2), (x + bw + gap - 0.01, ya + ha/2))

# ---- validation / thermal history / (b) simulation set
yb, hb = 3.0, 1.60
for i, (t, s, b) in enumerate([("Validation", "Sec. 3.1", "Lattice, cohesive and\nstacking-fault energy,\nelastic constants,\nliquid mixing enthalpy,\ncrystal tension"),
                               ("Thermal history", "Sec. 3.2", "Tg, fictive temperature,\neffective quench rate;\nRDF and density\nof the matrices")]):
    x = x0 + i*(bw + gap)
    ax.add_patch(FancyBboxPatch((x, yb), bw, hb, boxstyle="round,pad=0,rounding_size=0.06", fc=FILL[GREY], ec=GREY, lw=1.4))
    ax.text(x + bw/2, yb + hb - 0.10, t, ha="center", va="top", fontsize=FS, fontweight="bold")
    ax.text(x + bw/2, yb + hb - 0.30, s, ha="center", va="top", fontsize=FS, style="italic", color="#444")
    ax.text(x + 0.07, yb + hb - 0.50, b, ha="left", va="top", fontsize=FS, linespacing=1.2)
    arrow((x + bw/2, ya - 0.02), (x + bw/2, yb + hb + 0.02), GREY)
sx = x0 + 2*(bw + gap); sw = 2*bw + gap
ax.add_patch(FancyBboxPatch((sx, yb), sw, hb, boxstyle="round,pad=0,rounding_size=0.06", fc=FILL[GREEN], ec=GREEN, lw=1.6))
ax.text(sx + 0.10, yb + hb - 0.10, "(b) Simulation set", fontsize=FS + 1, fontweight="bold", va="top")
ax.text(sx + sw - 0.10, yb + hb - 0.12, "48 tensile simulations", fontsize=FS, fontweight="bold", va="top", ha="right")
cells = [("10 compositions", "Sec. 3.3"), ("3 pure-glass controls,\n2 crystal references", "Sec. 3.1, 3.3"),
         ("22 volume-fraction and\nseed-geometry cases", "Sec. 3.4"), ("11 rate and temperature\ncases, 4 compositions", "Sec. 3.7")]
cw, ch = (sw - 0.30)/2, 0.56
for k, (t, s) in enumerate(cells):
    cx = sx + 0.10 + (k % 2)*(cw + 0.10); cy = yb + 0.72 - (k//2)*(ch + 0.06)
    ax.add_patch(FancyBboxPatch((cx, cy - 0.0), cw, ch, boxstyle="round,pad=0,rounding_size=0.04", fc="white", ec=GREEN, lw=1.1))
    ax.text(cx + 0.06, cy + ch - 0.05, t, fontsize=FS, va="top", ha="left", linespacing=1.15)
    ax.text(cx + cw - 0.05, cy + 0.04, s, fontsize=FS, va="bottom", ha="right", style="italic", color="#444")
arrow((sx + bw/2, ya - 0.02), (sx + bw/2, yb + hb + 0.02), GREEN)
arrow((sx + bw + gap + bw/2, ya - 0.02), (sx + bw + gap + bw/2, yb + hb + 0.02), GREEN)

# ---- (c) analyses
ax.text(0.05, 2.78, "(c) Analyses of the 48 simulations", fontsize=FS + 1, fontweight="bold", va="center")
cwid, cg = 1.17, 0.10; yc, hc = 0.08, 2.00
ax.plot([sx + sw/2, sx + sw/2], [yb - 0.01, 2.45], color="#333", lw=1.4)
cx_list = [0.05 + i*(cwid + cg) for i in range(5)]
ax.plot([cx_list[0] + cwid/2, cx_list[4] + cwid/2], [2.45, 2.45], color="#333", lw=1.4)
C = [(ORANGE, "Mechanical\nresponse", "Sec. 2.4, 3.3,\n3.4, 3.7", "Yield strength\n(peak stress),\nmodulus, softening\nvs composition,\nvolume fraction,\nrate, temperature"),
     (ORANGE, "Deformation\ncharacter", "Sec. 3.5", "Shear-band\ndetection (von\nMises, DBSCAN);\nD$^2_{min}$ clusters;\nsensitivity to\nsettings"),
     (ORANGE, "Atomic-scale\norigin", "Sec. 3.6", "Volume-normalized\nstress regression,\nSHAP, q4/q6,\nWarren-Cowley order;\npure-glass control"),
     (ORANGE, "Interphase", "Sec. 3.5", "Shell D$^2_{min}$ vs\nbulk for every\nseed at 2-10%\nstrain (5% and\n15% volume\nfractions)"),
     (PURPLE, "Surrogate\nmodel", "Sec. 3.8", "Gaussian process\non 13-25 configs;\nretrospective\nheld-out tests;\nfollow-up ranked\nby uncertainty")]
for (col, t, s, b), x in zip(C, cx_list):
    ax.add_patch(FancyBboxPatch((x, yc), cwid, hc, boxstyle="round,pad=0,rounding_size=0.06", fc=FILL[col], ec=col, lw=1.4))
    ax.text(x + cwid/2, yc + hc - 0.10, t, ha="center", va="top", fontsize=FS, fontweight="bold", linespacing=1.15)
    ax.text(x + cwid/2, yc + hc - 0.50, s, ha="center", va="top", fontsize=FS, style="italic", color="#444", linespacing=1.15)
    ax.text(x + 0.06, yc + hc - 0.88, b, ha="left", va="top", fontsize=FS, linespacing=1.2)
    arrow((x + cwid/2, 2.45), (x + cwid/2, yc + hc + 0.02))
# follow-up loop: GP box back up to the simulation set
gx = cx_list[4] + cwid + 0.08
ax.add_patch(FancyArrowPatch((gx, yc + hc - 0.1), (gx, yb - 0.01), arrowstyle="-|>", mutation_scale=9, lw=1.4,
                             color=PURPLE, linestyle=(0, (4, 2)), connectionstyle="arc3,rad=0"))
for t in ax.texts:
    assert t.get_fontsize() >= 8, t.get_text()
import sys
for ext, kw in (("png", {}), ("tif", {"pil_kwargs": {"compression": "tiff_lzw"}}), ("pdf", {})):
    fig.savefig(os.path.join(HERE, "rendered", "fig01_workflow." + ext), dpi=600, **kw)
