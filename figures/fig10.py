import os, sys
ROOT = os.environ.get("MGZNCA_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
sys.path.insert(0, os.path.join(ROOT, "analysis"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np, matplotlib.pyplot as plt, style as S
S.apply()
fig,(a,b)=plt.subplots(1,2,figsize=(6.5,2.9),gridspec_kw={"width_ratios":[1,1.15]})
# (a) mean |SHAP| glass0_r1, 15% strain, volume-normalized, all descriptors + Voronoi volume
feat=[("Voronoi volume",.1713,"vol"),("Own type = Ca",.0992,"id"),("Mg neighbors (3.25 Å)",.0693,"o"),("Ca neighbors (4.0 Å)",.0628,"o"),("Ca neighbors (3.25 Å)",.0429,"o"),("Zn neighbors (3.25 Å)",.0418,"o"),("Warren–Cowley α",.0259,"o"),("q$_6$",.0125,"o")]
col={"vol":"#111111","id":S.CA["mid"],"o":"#9A9A9A"}
y=np.arange(len(feat))[::-1]
a.barh(y,[f[1] for f in feat],color=[col[f[2]] for f in feat],height=0.62)
for yi,f in zip(y,feat): a.text(f[1]+0.003,yi,f"{f[1]:.3f}",va="center",fontsize=6.8,color="#444")
a.set_yticks(y); a.set_yticklabels([f[0] for f in feat]); a.set_xlim(0,0.21); a.set_xlabel("Mean |SHAP value|"); S.panel(a,"a","Feature importance")
a.text(0.97,0.06,"pure glass, Ca 0.05\n15% strain",transform=a.transAxes,ha="right",fontsize=6.5,color="#555")
# (b) raw vs volume-normalized Ca-Zn log difference at 15 % strain
cases=["Ca 0.03","Ca 0.05","Ca 0.07"]
raw_g=[0.689,0.670,0.644]; nrm_g=[0.473,0.458,0.435]; raw_c=[0.692,0.672,0.644]; nrm_c=[0.476,0.460,0.435]
x=np.arange(3); w=0.19
b.bar(x-1.5*w,raw_g,w,color="#BDBDBD",label="stress × volume (raw)",ec="white")
b.bar(x-0.5*w,nrm_g,w,color=S.CA["mid"],label="volume-normalized stress",ec="white")
b.bar(x+0.5*w,raw_c,w,color="#BDBDBD",ec="white",hatch="////")
b.bar(x+1.5*w,nrm_c,w,color=S.CA["mid"],ec="white",hatch="////")
b.set_xticks(x); b.set_xticklabels(cases); b.set_ylim(0,1.05); b.axhline(np.log(1.24),color="#555",ls=":",lw=0.9,zorder=0); b.set_ylabel("Ca−Zn log-stress difference"); S.panel(b,"b","Raw versus volume-normalized")
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
b.legend(handles=[Patch(fc="#BDBDBD",label="stress × volume (raw)"),Patch(fc=S.CA["mid"],label="volume-normalized stress"),Patch(fc="white",ec="#444",label="solid: pure glass"),Patch(fc="white",ec="#444",hatch="////",label="hatched: single-seed composite"),Line2D([0],[0],color="#555",ls=":",lw=0.9,label="volume part, ln(V$_\\mathregular{Ca}$/V$_\\mathregular{Zn}$) = 0.21")],loc="upper right",ncol=1,fontsize=6.5,handlelength=1.3)

fig.tight_layout(w_pad=1.5); fig.savefig(os.path.join(HERE,"rendered","fig10_stress_asymmetry.png"),bbox_inches="tight",pad_inches=0.03,dpi=600); fig.savefig(os.path.join(HERE,"rendered","fig10_stress_asymmetry.pdf"),bbox_inches="tight",pad_inches=0.03)
