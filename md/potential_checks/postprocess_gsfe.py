#!/usr/bin/env python3
"""Usage: python3 postprocess_gsfe.py gsfe_partial_z.log [gsfe_partial_zx.log gsfe_straight_z.log]
Prints gamma(n) in mJ/m2 and the summary numbers (USF = max over the first third of the <10-10> path, ISF at 1/3)."""
import sys, re
EV_A2_TO_MJ_M2 = 16021.766
for fn in sys.argv[1:]:
    rows=[]
    for l in open(fn):
        if l.startswith("GSFE_POINT"):
            p=l.split(); rows.append((int(p[2]),float(p[3]),float(p[4]),float(p[5])))
    rows.sort()
    e0=rows[0][2]; 
    g=[(n,u,(e-e0)/A*EV_A2_TO_MJ_M2) for n,u,e,A in rows]
    print(f"\n== {fn}  ({len(g)} points)")
    for n,u,gg in g: print(f"  n={n:2d} u={u:7.4f} A  gamma={gg:9.3f} mJ/m2")
    last=g[-1][2]
    print(f"  periodicity: gamma(last) = {last:.4f} mJ/m2 (should be ~0)")
    if len(g)==37:
        third=g[:13]; usf=max(third,key=lambda r:r[2])
        print(f"  USF (max on 0..1/3 of <10-10>) = {usf[2]:.2f} mJ/m2 at n={usf[0]};  ISF (n=12) = {g[12][2]:.2f} mJ/m2; n=24: {g[24][2]:.2f}")
    elif len(g)==21:
        m=max(g,key=lambda r:r[2]); print(f"  max along <11-20> = {m[2]:.2f} mJ/m2 at n={m[0]}")
