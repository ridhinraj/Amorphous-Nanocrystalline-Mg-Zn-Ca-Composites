#!/usr/bin/env python3
"""voronoi_py.py -- per-atom Voronoi volumes without LAMMPS' VORONOI package (numpy + scipy only).

  python3 voronoi_py.py frames_selected.dump -o volumes.dump

Reads the positions and box of every frame (orthogonal periodic box, columns id type x y z ...), and writes a dump with
columns  id type c_vor[1] c_vor[2]  (volume in A^3, number of Voronoi neighbours) in the same layout that
voronoi_merge.py expects, so the rest of the pipeline is unchanged.

Method: the box is cut into slabs along x; each slab is padded by --pad A on both sides, periodic images are added, scipy's
Qhull Voronoi is run on the padded slab, and only the volumes of the slab's own atoms are kept.  A cell is exact when the
padding exceeds the cell's circumradius (about 4-5 A in these alloys); the default 8 A is generous.  Cells that touch the
padding edge are reported (should be zero).
"""
import argparse, sys, time
import numpy as np
from scipy.spatial import Voronoi, ConvexHull

def read_frames(path):
    with open(path) as f:
        while True:
            line = f.readline()
            if not line: return
            if not line.startswith("ITEM: TIMESTEP"): continue
            step = int(f.readline()); f.readline(); n = int(f.readline()); f.readline()
            box = np.array([[float(x) for x in f.readline().split()[:2]] for _ in range(3)])
            hdr = f.readline().split()[2:]
            data = np.loadtxt([f.readline() for _ in range(n)], ndmin=2)
            yield step, box, hdr, data

def volumes(pos, box, pad=8.0, slab=40.0):
    lo, L = box[:, 0], box[:, 1] - box[:, 0]
    if pad >= L.min() / 2: sys.exit("pad too large for this box")
    p = np.mod(pos - lo, L)                                   # wrap into [0, L)
    n = len(p)
    shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], float) * L
    img_p = (p[None, :, :] + shifts[:, None, :]).reshape(-1, 3)
    img_id = np.tile(np.arange(n), len(shifts))
    keep = np.all((img_p > -pad) & (img_p < L + pad), axis=1)  # drop images farther than pad from the box
    img_p, img_id = img_p[keep], img_id[keep]
    nslab = max(1, int(L[0] // slab)); edges = np.linspace(0, L[0], nslab + 1)
    vol = np.full(n, np.nan); nfaces = np.zeros(n, int); bad = 0
    for s in range(nslab):
        x0, x1 = edges[s], edges[s + 1]
        sel = (img_p[:, 0] > x0 - pad) & (img_p[:, 0] < x1 + pad)
        pts, ids = img_p[sel], img_id[sel]
        core = np.where((pts[:, 0] >= x0) & (pts[:, 0] < x1) &
                        np.all((pts >= 0) & (pts < L), axis=1))[0]
        vor = Voronoi(pts)
        counts = np.bincount(vor.ridge_points.ravel(), minlength=len(pts))
        for i in core:
            reg = vor.regions[vor.point_region[i]]
            if -1 in reg or len(reg) < 4: bad += 1; continue
            v = vor.vertices[reg]
            if (v[:, 0].min() < x0 - pad * 0.9) or (v[:, 0].max() > x1 + pad * 0.9): bad += 1; continue
            vol[ids[i]] = ConvexHull(v).volume; nfaces[ids[i]] = counts[i]
    return vol, nfaces, bad

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump"); ap.add_argument("-o", default="volumes.dump")
    ap.add_argument("--pad", type=float, default=8.0); ap.add_argument("--slab", type=float, default=40.0)
    a = ap.parse_args()
    with open(a.o, "w") as out:
        for step, box, hdr, d in read_frames(a.dump):
            t0 = time.time()
            ids = d[:, hdr.index("id")].astype(int); order = np.argsort(ids); d = d[order]
            pos = d[:, [hdr.index("x"), hdr.index("y"), hdr.index("z")]]
            vol, nf, bad = volumes(pos, box, a.pad, a.slab)
            vbox = float(np.prod(box[:, 1] - box[:, 0]))
            print(f"step {step}: N={len(d)}  sum(V)/V_box={np.nansum(vol) / vbox:.5f}  unresolved cells={bad + int(np.isnan(vol).sum() - bad)}  ({time.time() - t0:.0f} s)")
            if np.isnan(vol).any(): sys.exit("some cells unresolved -- increase --pad")
            out.write(f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n{len(d)}\nITEM: BOX BOUNDS pp pp pp\n")
            for k in range(3): out.write(f"{box[k, 0]:.16e} {box[k, 1]:.16e}\n")
            out.write("ITEM: ATOMS id type c_vor[1] c_vor[2]\n")
            for i in range(len(d)):
                out.write(f"{int(d[i, hdr.index('id')])} {int(d[i, hdr.index('type')])} {vol[i]:.6f} {nf[i]}\n")
    print("wrote", a.o)
main()
