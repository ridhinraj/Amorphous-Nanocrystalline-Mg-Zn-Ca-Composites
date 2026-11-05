#!/usr/bin/env python3
"""voronoi_merge.py -- convert the project's raw per-atom stress (stress x volume) into stress, using
Voronoi volumes from voronoi_rerun.in, and re-do the Ca-vs-Zn comparison.

  python3 voronoi_merge.py frames_selected.dump volumes.dump --label run_5 \
        --start-step 20000 --erate-ps 1e-3 --dt-ps 0.001 \
        --stress-file tension_stress_strain.txt --save-atoms run_5_atoms.npz

Input columns of frames_selected.dump:  id type x y z c_peratom_stress[1..6]  (sxx syy szz sxy sxz syz)
Types: 1 = Mg, 2 = Zn, 3 = Ca.   Output: a table per frame and species, plus self-checks.

Units:  LAMMPS compute stress/atom (metal units) is in bar*Angstrom^3 and equals -(pressure tensor) x volume,
so s_i / V_i is already the tension-positive stress used in the project (stress_x = -pxx).  GPa = bar * 1e-4.
The von Mises value does not depend on the sign.  Self-check: sum(s_i)/V_box must equal the thermo stress.
"""
import argparse, sys
import numpy as np

NAMES = {1: "Mg", 2: "Zn", 3: "Ca"}
BAR_A3_TO_GPA_A3 = 1.0e-4

def read_frames(path, ncols_expected=None):
    """yield dict(step, box(3,2), cols(list), data(N,ncol))"""
    with open(path) as f:
        while True:
            line = f.readline()
            if not line: return
            if not line.startswith("ITEM: TIMESTEP"): continue
            step = int(f.readline())
            f.readline(); n = int(f.readline())
            f.readline()
            box = np.array([[float(x) for x in f.readline().split()[:2]] for _ in range(3)])
            hdr = f.readline().split()[2:]
            data = np.loadtxt([f.readline() for _ in range(n)], ndmin=2)
            order = np.argsort(data[:, hdr.index("id")])
            yield dict(step=step, box=box, cols=hdr, data=data[order])

def vm(s):
    sxx, syy, szz, sxy, sxz, syz = s.T
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
                   + 3.0 * (sxy ** 2 + sxz ** 2 + syz ** 2))

def global_stress_at(path, strain):
    """nearest row of tension_stress_strain.txt (Strain Stress_X Stress_Y Stress_Z), or None"""
    try:
        arr = np.loadtxt(path, skiprows=1)
    except Exception:
        return None
    i = int(np.argmin(np.abs(arr[:, 0] - strain)))
    return arr[i, 0], arr[i, 1:4]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stress_dump"); ap.add_argument("volume_dump")
    ap.add_argument("--label", default="case")
    ap.add_argument("--start-step", type=int, default=20000)
    ap.add_argument("--erate-ps", type=float, default=1.0e-3)
    ap.add_argument("--dt-ps", type=float, default=0.001)
    ap.add_argument("--stress-file", default=None)
    ap.add_argument("--save-atoms", default=None)
    a = ap.parse_args()
    per_step = a.erate_ps * a.dt_ps

    vols = {f["step"]: f for f in read_frames(a.volume_dump)}
    rows, atoms_out, warn = [], {}, 0
    first = None            # (step, ids, sig) of the first frame, for load-induced changes
    part = []              # load-partition rows
    print(f"=== {a.label} ===")
    for fr in read_frames(a.stress_dump):
        step, d, cols = fr["step"], fr["data"], fr["cols"]
        if step not in vols: sys.exit(f"step {step} missing from volume dump")
        vf = vols[step]
        if not np.array_equal(d[:, cols.index("id")], vf["data"][:, vf["cols"].index("id")]):
            sys.exit("atom ids differ between the two dumps")
        typ = d[:, cols.index("type")].astype(int)
        sc = [c for c in cols if c.startswith("c_") and "[" in c]
        if len(sc) != 6: sys.exit(f"expected six c_*[1..6] stress columns, found {sc}")
        s_raw = d[:, [cols.index(c) for c in sc]]                    # bar*A^3
        vol_i = vf["data"][:, vf["cols"].index("c_vor[1]")]          # A^3
        box = fr["box"]; vbox = float(np.prod(box[:, 1] - box[:, 0]))
        strain = (step - a.start_step) * per_step
        sig = s_raw / vol_i[:, None] * BAR_A3_TO_GPA_A3              # GPa, tension positive (s = -P*V)
        vm_raw, vm_gpa = vm(s_raw), vm(sig)
        # --- checks ---
        vsum = vol_i.sum()
        sum_stress = s_raw[:, :3].sum(0) / vbox * BAR_A3_TO_GPA_A3   # GPa; equals the tension-positive global stress
        sum_vor = (sig[:, :3] * vol_i[:, None]).sum(0) / vsum        # volume-weighted mean
        msg = f"step {step}  strain {strain:.4f}  N={len(d)}  sum(V_i)/V_box = {vsum / vbox:.5f}"
        if abs(vsum / vbox - 1) > 1e-3: msg += "   <-- WARNING: Voronoi volumes do not tile the box"; warn += 1
        print(msg)
        print(f"   mean stress xx,yy,zz from atoms = {sum_stress[0]:+.3f} {sum_stress[1]:+.3f} {sum_stress[2]:+.3f} GPa")
        g = global_stress_at(a.stress_file, strain) if a.stress_file else None
        if g is not None:
            print(f"   thermo stress (strain {g[0]:.4f})  = {g[1][0]:+.3f} {g[1][1]:+.3f} {g[1][2]:+.3f} GPa"
                  f"   (should agree to a few %, residual = noise + thermostat KE term)")
        for t in (1, 2, 3):
            m = typ == t
            if not m.any(): continue
            rows.append(dict(label=a.label, step=step, strain=strain, species=NAMES[t], n=int(m.sum()),
                             vol_A3=vol_i[m].mean(), vm_raw_mean=vm_raw[m].mean(), vm_GPa_mean=vm_gpa[m].mean(),
                             vm_raw_logmean=np.log(vm_raw[m][vm_raw[m] > 0]).mean(),
                             vm_GPa_logmean=np.log(vm_gpa[m][vm_gpa[m] > 0]).mean()))
        if first is None:
            first = (step, d[:, cols.index("id")].copy(), sig.copy())
        else:
            if not np.array_equal(first[1], d[:, cols.index("id")]): sys.exit("atom ids differ from first frame")
            dsig = sig[:, :3] - first[2][:, :3]                      # GPa, per atom, since the first frame
            gl = (dsig * vol_i[:, None]).sum(0) / vsum               # volume-weighted mean change = global change
            for t_ in (1, 2, 3):
                m_ = typ == t_
                if not m_.any(): continue
                sp = (dsig[m_] * vol_i[m_, None]).sum(0) / vol_i[m_].sum()
                part.append(dict(label=a.label, step=step, strain=strain, species=NAMES[t_],
                                 vol_fraction=vol_i[m_].sum() / vsum, dsxx=sp[0], dsyy=sp[1], dszz=sp[2],
                                 global_dsxx=gl[0], factor_xx=sp[0] / gl[0] if abs(gl[0]) > 1e-9 else float("nan"),
                                 se_factor=(dsig[m_, 0].std() / np.sqrt(m_.sum())) / abs(gl[0]) if abs(gl[0]) > 1e-9 else float("nan")))
        if a.save_atoms:
            atoms_out[f"id_{step}"] = d[:, cols.index("id")].astype(np.int64)
            atoms_out[f"type_{step}"] = typ; atoms_out[f"vol_{step}"] = vol_i
            atoms_out[f"vm_gpa_{step}"] = vm_gpa; atoms_out[f"vm_raw_{step}"] = vm_raw
    # --- summary ---
    print("\nspecies-resolved means:")
    print(f"{'step':>8} {'strain':>7} {'sp':>3} {'n':>7} {'V (A^3)':>8} {'vm raw':>10} {'vm GPa':>8} {'ln raw':>8} {'ln GPa':>8}")
    for r in rows:
        print(f"{r['step']:8d} {r['strain']:7.4f} {r['species']:>3} {r['n']:7d} {r['vol_A3']:8.2f} "
              f"{r['vm_raw_mean']:10.1f} {r['vm_GPa_mean']:8.3f} {r['vm_raw_logmean']:8.3f} {r['vm_GPa_logmean']:8.3f}")
    print("\nCa vs Zn  (the quantity behind 'Ca carries ~2x the stress of Zn'):")
    for step in sorted({r["step"] for r in rows}):
        zn = next((r for r in rows if r["step"] == step and r["species"] == "Zn"), None)
        ca = next((r for r in rows if r["step"] == step and r["species"] == "Ca"), None)
        if zn and ca:
            print(f"  step {step}: volume Ca/Zn = {ca['vol_A3'] / zn['vol_A3']:.2f} | "
                  f"raw (stress*V) ratio = {ca['vm_raw_mean'] / zn['vm_raw_mean']:.2f}, d ln = {ca['vm_raw_logmean'] - zn['vm_raw_logmean']:+.3f} | "
                  f"stress ratio = {ca['vm_GPa_mean'] / zn['vm_GPa_mean']:.2f}, d ln = {ca['vm_GPa_logmean'] - zn['vm_GPa_logmean']:+.3f}")
    if part:
        print("\nLoad-induced stress change since the first frame (volume-weighted mean per species, GPa, tension positive):")
        print("  factor = species mean change / global change along the loading axis; 1.00 = carries its fair share")
        print(f"{'strain':>7} {'sp':>3} {'vol frac':>9} {'d sxx':>8} {'d syy':>8} {'d szz':>8} {'global dsxx':>12} {'factor':>7} {'+/-':>6}")
        for r in part:
            print(f"{r['strain']:7.4f} {r['species']:>3} {r['vol_fraction']:9.4f} {r['dsxx']:8.3f} {r['dsyy']:8.3f} "
                  f"{r['dszz']:8.3f} {r['global_dsxx']:12.3f} {r['factor_xx']:7.2f} {r['se_factor']:6.2f}")
        with open(f"{a.label}_load_partition.csv", "w", newline="") as f_:
            w_ = __import__("csv").DictWriter(f_, fieldnames=list(part[0])); w_.writeheader(); w_.writerows(part)
        print(f"wrote {a.label}_load_partition.csv")
    import csv
    out = f"{a.label}_species_summary.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print("wrote", out)
    if a.save_atoms:
        np.savez_compressed(a.save_atoms, **atoms_out); print("wrote", a.save_atoms)
    if warn: print(f"\n{warn} frame(s) failed the volume-tiling check -- do not use these volumes")
main()
