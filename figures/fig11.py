import os, sys
ROOT = os.environ.get("MGZNCA_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
sys.path.insert(0, os.path.join(ROOT, "analysis"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, matplotlib.pyplot as plt, style as S
S.apply()
fig,(a,b)=plt.subplots(1,2,figsize=(6.5,2.8))
x=np.arange(3); w=0.36; lab=["Ca 0.03","Ca 0.05","Ca 0.07"]
zq=[.4916,.4917,.4921]; cq=[.3120,.3143,.3165]
a.bar(x-w/2,zq,w,color=S.ZN["mid"],label="Zn atoms",ec="white"); a.bar(x+w/2,cq,w,color=S.CA["mid"],label="Ca atoms",ec="white")
a.axhline(.484,color="#555",ls="--",lw=0.9); a.text(2.55,.484+0.012,"bulk HCP 0.484",ha="right",fontsize=6.5,color="#444")
a.axhline(.575,color="#555",ls=":",lw=1.0); a.text(2.55,.575+0.012,"bulk FCC 0.575",ha="right",fontsize=6.5,color="#444")
a.set_ylim(0,0.8); a.set_xticks(x); a.set_xticklabels(lab); a.set_ylabel("Mean q$_6$"); S.panel(a,"a","Local packing order"); a.legend(loc="upper left",ncol=2,columnspacing=1.0)
za=[.222,.220,.211]; ca=[.511,.498,.463]
b.bar(x-w/2,za,w,color=S.ZN["mid"],label="Zn–Zn",ec="white"); b.bar(x+w/2,ca,w,color=S.CA["mid"],label="Ca–Ca",ec="white")
for xi,(zv,cv) in enumerate(zip(za,ca)): b.text(xi+w/2,cv+0.012,f"{cv/zv:.2f}×",ha="center",fontsize=7,color="#444")
b.set_ylim(0,0.7); b.set_xticks(x); b.set_xticklabels(lab); b.set_ylabel("Warren–Cowley α (own type)"); S.panel(b,"b","Chemical short-range order"); b.legend(loc="upper right",ncol=2,columnspacing=1.0)
fig.tight_layout(w_pad=1.5); S.save(fig,"fig11_structure")
