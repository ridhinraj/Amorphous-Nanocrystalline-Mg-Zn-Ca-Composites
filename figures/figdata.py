import glob, numpy as np, pandas as pd
import os, sys
ROOT = os.environ.get("MGZNCA_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
sys.path.insert(0, os.path.join(ROOT, "analysis"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_props import smooth, props
def _find(case):
    cands=[os.path.join(DATA, "stress_strain", f"{case}__tension_stress_strain.txt")]
    cands.append(os.path.join(ROOT, case, "tension_stress_strain.txt"))   # or a case directory in the project root
    for c in cands:
        if os.path.exists(c) and os.path.getsize(c)>10000: return c
    return None
def curve(case):
    p=_find(case)
    if p is None: raise FileNotFoundError(case)
    d=pd.read_csv(p,sep=r"\s+"); col="Stress_Z_GPa" if case=="hcp_z" else "Stress_X_GPa"
    e=d["Strain"].values; s=d[col].values; return e,s,smooth(e,s)
