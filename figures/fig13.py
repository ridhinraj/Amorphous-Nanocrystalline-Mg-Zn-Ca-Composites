import os, sys
ROOT = os.environ.get("MGZNCA_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
sys.path.insert(0, os.path.join(ROOT, "analysis"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, matplotlib.pyplot as plt, style as S, warnings
warnings.filterwarnings("ignore")
import gp_heldout as H
S.apply(); np.random.seed(0)
d16=H.old.loc[H.S16]; rows=[]
for i in range(len(d16)):
    tr=d16.drop(d16.index[i]); gp,sc=H.fit(tr); m,s=H.pred(gp,sc,d16.iloc[[i]]); rows.append((d16.index[i],d16.y.iloc[i],m[0],s[0]))
loo=pd.DataFrame(rows,columns=["case","sim","pred","sd"]); loo["z"]=(loo.sim-loo.pred)/loo.sd
loo.to_csv(os.path.join(DATA,"loo_16point.csv"),index=False); print(loo.round(3).to_string())
inf=loo[loo.case!="t310_1"].copy(); inf["out"]=abs(inf.z)>1.96
print("within 95%:",(~inf.out).sum(),"of",len(inf),"; z std all16 %.2f"%loo.z.std(ddof=0))
ho=pd.read_csv(os.path.join(DATA,"heldout_predictions.csv"))
fig,(a,b)=plt.subplots(1,2,figsize=(6.5,2.95))
lim=(0.88,1.32)
a.plot(lim,lim,"-",color="#999",lw=0.8,zorder=0)
for _,r in inf.iterrows():
    a.errorbar(r.sim,r.pred,yerr=1.96*r.sd,fmt="o",color="#222",mfc="white" if r.out else "#222",mew=1.1,ms=4.2,capsize=2,lw=0.9)
a.set_xlim(lim); a.set_ylim(lim); a.set_xlabel("Simulated peak strength (GPa)"); a.set_ylabel("Leave-one-out prediction (GPa)")
S.panel(a,"a",f"Leave-one-out ({(~inf.out).sum()}/{len(inf)} within 95%)")
from matplotlib.lines import Line2D
a.legend(handles=[Line2D([],[],marker="o",ls="",color="#222",ms=4.2,label="within 95% interval"),Line2D([],[],marker="o",ls="",color="#222",mfc="white",mew=1.1,ms=4.2,label="outside 95% interval")],loc="upper left")
# (b) frozen-surrogate tests
grp={"VF = 0 (pure glass)":("glass0",S.GLASS,"D"),"5–15% VF, Ca 0.03":("vf_r5",S.VF5[2],"o"),"10% packings, Ca 0.03":("rep_r5",S.VF5[1],"s"),"10% packings, Ca 0.07":("rep_r6",S.VF5[3],"^"),"2×10$^9$ s$^{-1}$, two new compositions":("rate2e9",S.RATE["5e8"],"P")}
b.plot(lim,lim,"-",color="#999",lw=0.8,zorder=0)
for name,(key,col,m) in grp.items():
    sub=ho[ho.case.str.startswith(key)]
    for j,(_,r) in enumerate(sub.iterrows()):
        out=abs(r.z)>1.96
        wide=1.96*r.sigma>0.06
        b.errorbar(r.simulated,r.pred,yerr=1.96*r.sigma,fmt=m,color=col,mfc="white" if out else col,mec=col if out else "#222",mew=1.1 if out else 0.4,ms=4.6,capsize=2,elinewidth=1.0 if wide else 1.4,alpha=0.55 if wide else 1.0,lw=0.9,label=name if j==0 else None)
for sim,pred,sd in ((0.928,0.926,0.011),(0.917,0.931,0.010)):
    b.errorbar(sim,pred,yerr=1.96*sd,fmt="*",color="#444",ms=6.5,capsize=2,lw=0.9,label="18.5%, 18.9% VF (13-point model)" if sim==0.928 else None)
b.set_xlim(lim[0],1.42); b.set_ylim(lim[0],1.42); b.plot((lim[0],1.42),(lim[0],1.42),"-",color="#999",lw=0.8,zorder=0)
b.set_xlabel("Simulated peak strength (GPa)"); b.set_ylabel("Predicted (GPa)"); S.panel(b,"b","Tests of the frozen surrogate")
hs=[Line2D([],[],marker=m,ls="",color=col,mec="#222",mew=0.4,ms=4.6,label=name) for name,(key,col,m) in grp.items()]+[Line2D([],[],marker="*",ls="",color="#444",ms=6.5,label="18.5%, 18.9% VF (13-point model)"),Line2D([],[],marker="o",ls="",color="#222",mfc="white",mew=1.1,ms=4.2,label="open: outside 95% interval")]
b.legend(handles=hs,loc="upper left",handlelength=1.0,labelspacing=0.3,fontsize=6.4)
fig.tight_layout(w_pad=1.5); S.save(fig,"fig13_surrogate")
