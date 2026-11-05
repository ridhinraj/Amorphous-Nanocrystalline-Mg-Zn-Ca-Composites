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
fig,(a,b,c)=plt.subplots(1,3,figsize=(6.5,2.6),gridspec_kw={"width_ratios":[1.1,1,0.8]})
S.panel(a,"a","Strain rate, 300 K")
# small multiplicative offsets along the (log) rate axis so that compositions whose
# yield strengths coincide (Comp. 3 and Comp. 6 at 2x10^9) remain separately visible
comp={"Comp. 1":(["rate1e8_1","rate5e8_1","run_1","rate2e9_1"],[1e8,5e8,1e9,2e9],"o","#111111",1.00),
      "Comp. 6":(["run_6","rate2e9_run6"],[1e9,2e9],"s","#444444",0.93),
      "Comp. 4":(["run_4","rate2e9_run4"],[1e9,2e9],"^","#8A8A8A",1.00),
      "Comp. 3":(["run_3","rate2e9_run3"],[1e9,2e9],"D","#CFCFCF",1.08)}
for k,(cs,r,m,col,jit) in list(comp.items())[::-1]:
    a.plot([x*jit for x in r],[pk(x) for x in cs],"-",marker=m,color=col,ms=4.6,mec="#111",mew=0.6,label=k)
a.set_xscale("log"); a.set_xlabel("Strain rate (s$^{-1}$)"); a.set_ylabel("Peak strength, σ$_\\mathregular{max}$ (GPa)"); a.legend(loc="lower right",handlelength=1.4)
a.set_xticks([1e8,1e9]); a.set_xticklabels(["10$^8$","10$^9$"]); a.set_xlim(7e7,3e9)
S.panel(b,"b","Temperature, 10$^9$ s$^{-1}$")
T=np.array([250,300,350,400]); y=np.array([pk(f"temp_{t}") for t in T]); f=np.array([fl(f"temp_{t}") for t in T])
for t,yy,ff in zip(T,y,f): b.plot(t,yy,"o",color=S.TEMP[t],ms=5.2,mec="#222",mew=0.5); b.plot(t,ff,"o",color=S.TEMP[t],ms=5.2,mfc="white",mec=S.TEMP[t],mew=1.1)
p=np.polyfit(T,y,1); xx=np.array([240,410]); b.plot(xx,np.polyval(p,xx),"-",color="#555",lw=0.9,zorder=0)
b.plot(310,pk("t310_1"),"X",color="#111",ms=5)
from matplotlib.lines import Line2D
b.legend(handles=[Line2D([],[],marker="o",ls="",color="#B87700",mec="#222",mew=0.5,ms=5,label="peak strength"),Line2D([],[],marker="o",ls="",color="#B87700",mfc="white",mew=1.1,ms=5,label="flow stress, 15% strain"),Line2D([],[],marker="X",ls="",color="#111",ms=5,label="310 K (earlier run)")],loc="lower left",handlelength=1.0,borderaxespad=0.1)
b.text(0.97,0.97,f"slope {p[0]*1000:.2f} MPa/K\n(−0.25% per K at 300 K)".replace("-","−"),transform=b.transAxes,ha="right",va="top",fontsize=6.6,color="#444")

b.set_xlabel("Temperature (K)"); b.set_ylabel("Stress (GPa)"); b.set_xlim(235,415)
S.panel(c,"c","Doubling the rate")
names=["Comp. 1","Comp. 6","Comp. 4","Comp. 3"]; pairs=[("run_1","rate2e9_1"),("run_6","rate2e9_run6"),("run_4","rate2e9_run4"),("run_3","rate2e9_run3")]
r=[100*(pk(b_)/pk(a_)-1) for a_,b_ in pairs]; yy=np.arange(4)[::-1]
c.barh(yy,r,color=["#111111","#555555","#888888","#BBBBBB"],height=0.55)
for yi,v in zip(yy,r): c.text(v+0.2,yi,f"+{v:.1f}%",va="center",fontsize=6.8)
c.set_yticks(yy); c.set_yticklabels(names); c.set_xlim(0,14); c.set_xlabel("Increase in peak strength (%)")
fig.tight_layout(w_pad=1.1); S.save(fig,"fig12_rate_temperature")
