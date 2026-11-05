#!/usr/bin/env python3
"""select_frames.py -- copy only chosen frames out of a (large) multi-frame LAMMPS dump.

  python3 select_frames.py tension_trajectory.dump --list
  python3 select_frames.py tension_trajectory.dump --strains 0 0.0585 0.15 \
          --start-step 20000 --erate-ps 1e-3 --dt-ps 0.001 -o frames_selected.dump
  python3 select_frames.py tension_trajectory.dump --steps 20000 78500 170000 -o frames_selected.dump

Strain -> step:  step = start_step + strain / (erate_ps * dt_ps).   erate_ps is the `fix deform ... erate`
value exactly as written in the tension input (metal units: 1/ps; 1e9 1/s = 1e-3 1/ps).  The requested
strain is snapped to the nearest frame actually present in the dump, and the strain really used is printed.
"""
import argparse, sys

def frames(path):
    """yield (timestep, [lines of the frame])"""
    buf, step = [], None
    with open(path) as f:
        for line in f:
            if line.startswith("ITEM: TIMESTEP"):
                if buf:
                    yield step, buf
                buf = [line]
                step = int(next(f)); buf.append(f"{step}\n")
            else:
                buf.append(line)
        if buf:
            yield step, buf

def scan_steps(path):
    steps = []
    with open(path) as f:
        for line in f:
            if line.startswith("ITEM: TIMESTEP"):
                steps.append(int(next(f)))
    return steps

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump"); ap.add_argument("-o", default="frames_selected.dump")
    ap.add_argument("--list", action="store_true", help="only list the timesteps present")
    ap.add_argument("--steps", type=int, nargs="*"); ap.add_argument("--strains", type=float, nargs="*")
    ap.add_argument("--start-step", type=int, default=20000)
    ap.add_argument("--erate-ps", type=float, default=1.0e-3)
    ap.add_argument("--dt-ps", type=float, default=0.001)
    a = ap.parse_args()
    steps = scan_steps(a.dump)
    per_step = a.erate_ps * a.dt_ps
    if a.list:
        print(f"{len(steps)} frames")
        for s in steps: print(s, f"strain={(s - a.start_step) * per_step:.5f}")
        return
    want = list(a.steps or [])
    for e in (a.strains or []):
        target = a.start_step + e / per_step
        want.append(min(steps, key=lambda s: abs(s - target)))
    want = sorted(set(want))
    missing = [s for s in want if s not in steps]
    if missing: sys.exit(f"steps not in dump: {missing}")
    keep = set(want)
    with open(a.o, "w") as out:
        for s, lines in frames(a.dump):
            if s in keep:
                out.writelines(lines)
    for s in want: print(f"kept step {s}  strain={(s - a.start_step) * per_step:.5f}")
    print("wrote", a.o)
main()
