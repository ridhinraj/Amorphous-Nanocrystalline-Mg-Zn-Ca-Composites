import os, sys
ROOT = os.environ.get("MGZNCA_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
sys.path.insert(0, os.path.join(ROOT, "analysis"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd, numpy as np, matplotlib.pyplot as plt, style as S
S.apply()
o=pd.read_csv(os.path.join(DATA,"stress_strain_properties_all25.csv")).set_index("case")
n=pd.read_csv(os.path.join(DATA,"new_md_properties_final.csv")).set_index("case")
pk=lambda c:o.loc[c,"yield_strength_GPa"] if c in o.index else n.loc[c,"peak_stress_GPa"]
fl=lambda c:o.loc[c,"flow_stress_15pct_GPa"] if c in o.index else n.loc[c,"flow_stress_15pct_GPa"]
fig,(a,b)=plt.subplots(1,2,figsize=(6.5,2.9))
# (a) Zn
S.panel(a,"a","Zinc content (x$_\\mathregular{Ca}$ ≈ 0.05)")
zs=[("run_4",.25),("run_9",.32),("run_1",.35),("run_10",.38)]
a.plot([x for _,x in zs],[pk(c) for c,_ in zs],"-o",color=S.ZN["mid"],label="Peak strength")
a.plot([x for _,x in zs],[fl(c) for c,_ in zs],"--o",color=S.ZN["mid"],mfc="white",label="Flow stress, 15% strain")
a.set_xlabel("Zn content, x$_\\mathregular{Zn}$"); a.set_ylabel("Stress (GPa)"); a.set_ylim(0.5,1.35); a.legend(loc="lower right")
# (b) Ca
S.panel(b,"b","Calcium content")
ser={0.30:([("run_2",.04),("run_7",.06)],S.ZN["light"],"s"),0.35:([("run_5",.03),("run_1",.05),("run_6",.07)],S.ZN["mid"],"o"),0.40:([("run_8",.04),("run_3",.06)],S.ZN["dark"],"^")}
for z,(pts,col,m) in ser.items():
    b.plot([x for _,x in pts],[pk(c) for c,_ in pts],"-",marker=m,color=col,label=f"x$_\\mathregular{{Zn}}$ = {z:.2f}")
    b.plot([x for _,x in pts],[fl(c) for c,_ in pts],"--",marker=m,color=col,mfc="white")
g=[("glass0_r5",.03),("glass0_r1",.05),("glass0_r6",.07)]
b.plot([x for _,x in g],[pk(c) for c,_ in g],"-D",color=S.GLASS,ms=4.2,label="Pure glass, x$_\\mathregular{Zn}$ = 0.35")
b.plot([x for _,x in g],[fl(c) for c,_ in g],"--D",color=S.GLASS,ms=4.2,mfc="white")
b.set_xlabel("Ca content, x$_\\mathregular{Ca}$"); b.set_ylim(0.5,1.35); b.set_xticks([.03,.04,.05,.06,.07]); b.legend(loc="lower right",ncol=1)
b.text(0.03,0.97,"solid: peak strength\ndashed, open: flow stress, 15% strain",transform=b.transAxes,va="top",fontsize=6.8,color="#444")
fig.tight_layout(w_pad=1.5); S.save(fig,"fig06_composition")
