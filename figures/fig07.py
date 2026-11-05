import os, sys
ROOT = os.environ.get("MGZNCA_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
sys.path.insert(0, os.path.join(ROOT, "analysis"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd, numpy as np, matplotlib.pyplot as plt, style as S
S.apply()
o=pd.read_csv(os.path.join(DATA,"stress_strain_properties_all25.csv")).set_index("case")
n=pd.read_csv(os.path.join(DATA,"new_md_properties_final.csv")).set_index("case")
g=lambda c,k_o,k_n:(o.loc[c,k_o] if c in o.index else n.loc[c,k_n])
pk=lambda c:g(c,"yield_strength_GPa","peak_stress_GPa"); so=lambda c:g(c,"softening_pct","softening_pct")
VF={"vf_r5_1":4.97,"vf_r5_2":10.02,"vf_r5_3":14.98,"vf_1":5.00,"vf_2":10.01,"vf_3":15.01,"vf_25_1":18.49,"vf_20_1":18.91,"vf_r6_1":5.04,"vf_r6_2":10.00,"vf_r6_3":15.04,
    "rep_r5_a":10.0,"rep_r5_b":10.0,"rep_r5_c":10.0,"rep_r1_a":10.0,"rep_r1_b":10.0,"rep_r1_c":10.0,"rep_r6_a":10.0,"rep_r6_b":10.0,"rep_r6_c":10.0}
M={"r5":dict(ca=.03,glass="glass0_r5",s=["vf_r5_1","vf_r5_2","vf_r5_3"],rep=["vf_r5_2","rep_r5_a","rep_r5_b","rep_r5_c"],single="run_5",svf=.16,m="o"),
   "r1":dict(ca=.05,glass="glass0_r1",s=["vf_1","vf_2","vf_3","vf_25_1","vf_20_1"],rep=["vf_2","rep_r1_a","rep_r1_b","rep_r1_c"],single="run_1",svf=.74,m="s"),
   "r6":dict(ca=.07,glass="glass0_r6",s=["vf_r6_1","vf_r6_2","vf_r6_3"],rep=["vf_r6_2","rep_r6_a","rep_r6_b","rep_r6_c"],single="run_6",svf=2.02,m="^")}
fig,axs=plt.subplots(2,2,figsize=(6.5,5.3)); a,b,c,d=axs.ravel()
for k,m in M.items():
    col=S.CA_LEVEL[m["ca"]]; lab=f"Ca {m['ca']:.2f}"
    x=[0]+[VF[s] for s in m["s"]]; y=[pk(m["glass"])]+[pk(s) for s in m["s"]]
    a.plot(x[:4],y[:4],"-",color=col,lw=1.3); a.plot(x[4:],y[4:],"o",color=col,ms=4) if len(x)>4 else None
    a.plot(x[1:],y[1:],m["m"],color=col,ms=4.6,label=lab,mec="white",mew=0.5); a.plot(0,pk(m["glass"]),m["m"],color=col,mfc="white",ms=5.2,mew=1.2)
    r=np.array([pk(s) for s in m["rep"]]); a.errorbar(10.0,r.mean(),yerr=r.std(ddof=1),color="black",capsize=2.5,lw=0.9,zorder=5)
    a.plot(m["svf"],pk(m["single"]),m["m"],color=col,ms=3.6,mfc="none",mew=0.9)
    # relative change
    rel=[100*(pk(s)/pk(m["glass"])-1) for s in m["s"]]; b.plot([0]+[VF[s] for s in m["s"]],[0]+rel,"-",color=col,lw=1.3)
    b.plot([VF[s] for s in m["s"]],rel,m["m"],color=col,ms=4.6,mec="white",mew=0.5,label=lab)
    b.plot(m["svf"],100*(pk(m["single"])/pk(m["glass"])-1),m["m"],color=col,ms=3.6,mfc="none",mew=0.9)
    # softening
    sx=[0]+[VF[s] for s in m["s"]]; sy=[so(m["glass"])]+[so(s) for s in m["s"]]
    c.plot(sx,sy,"-",color=col,lw=1.3); c.plot(sx[1:],sy[1:],m["m"],color=col,ms=4.6,mec="white",mew=0.5,label=lab); c.plot(0,sy[0],m["m"],color=col,mfc="white",ms=5.2,mew=1.2)
# fit line
allx=[];ally=[]
for k,m in M.items():
    for s in m["s"][:3]: allx.append(VF[s]); ally.append(100*(pk(s)/pk(m["glass"])-1))
sl=np.sum(np.array(allx)*np.array(ally))/np.sum(np.array(allx)**2); xx=np.array([0,19]); b.plot(xx,sl*xx,":",color="#555",lw=1.0); b.text(12.8,sl*13.5+3.2,f"\u2212{abs(sl):.2f}% per 1% VF",fontsize=6.8,color="#444",rotation=0)
print("pooled slope through origin (5-15 %):",sl)
b.axhline(0,color="#999",lw=0.6)
for ax,l,t in ((a,"a","Peak strength"),(b,"b","Change relative to the matched glass"),(c,"c","Post-peak softening")): S.panel(ax,l,t)
a.set_ylabel("Peak strength, σ$_\\mathregular{max}$ (GPa)"); b.set_ylabel("Relative change (%)"); c.set_ylabel("Softening (%)")
for ax in (a,b,c): ax.set_xlabel("Nanocrystal volume fraction (%)"); ax.set_xlim(-1,20)
a.legend(loc="lower left",title="matrix",title_fontsize=7,handlelength=1.2); a.text(0.97,0.97,"hollow, VF = 0: pure glass\nbars: s.d. of 4 packings\nsmall open: single-seed",transform=a.transAxes,ha="right",va="top",fontsize=6.3,color="#444")
# seed radius
S.panel(d,"d","Seed radius, 10% VF, Ca 0.05")
r=[12,15,20]; yr=[pk("vf_2"),pk("vf_rad15_1"),pk("vf_rad20_1")]
d.plot(r,yr,"-",color=S.SEED["15"],lw=1.3)
for rr_,yy,k in zip(r,yr,("12","15","20")): d.plot(rr_,yy,"D",color=S.SEED[k],ms=5,mec="white",mew=0.5)
rep=np.array([pk(s) for s in M["r1"]["rep"]]); d.axhspan(rep.mean()-rep.std(ddof=1),rep.mean()+rep.std(ddof=1),color="#DDDDDD",zorder=0)
d.text(20.8,rep.mean(),"packing scatter (±1 s.d., 12 Å)",fontsize=6.3,color="#444",ha="right",va="center")
d.set_xlabel("Seed radius (Å)"); d.set_ylabel("Peak strength, σ$_\\mathregular{max}$ (GPa)"); d.set_xticks(r); d.set_xlim(11,21)
fig.tight_layout(h_pad=1.1,w_pad=1.4); S.save(fig,"fig07_volume_fraction")
