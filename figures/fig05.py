import os, sys
ROOT = os.environ.get("MGZNCA_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
sys.path.insert(0, os.path.join(ROOT, "analysis"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, matplotlib.pyplot as plt, style as S, figdata as D
from extract_props import props
S.apply()
fig,axs=plt.subplots(4,2,figsize=(6.5,8.6)); A=axs.ravel()
def draw(ax,case,color,label,ls="-",lw=1.4,dots=True,zorder=2):
    e,s,sm=D.curve(case); ax.plot(e,s,color=color,lw=0.4,alpha=0.25,zorder=1); ax.plot(e,sm,color=color,lw=lw,ls=ls,label=label,zorder=zorder)
    if dots:
        i=np.argmax(sm); ax.plot(e[i],sm[i],"o",ms=3.6,mfc=color,mec="white",mew=0.6,zorder=3)
# (a) pure glass, Ca series
S.panel(A[0],"a","Pure glass, three Ca contents")
for c,k in (("glass0_r5",0.03),("glass0_r1",0.05),("glass0_r6",0.07)): draw(A[0],c,S.CA_LEVEL[k],f"x$_\\mathregular{{Ca}}$ = {k:.2f}")
# (b-d) VF series per matrix with glass
vfs={"b":("Ca 0.03 matrix",["vf_r5_1","vf_r5_2","vf_r5_3"],"glass0_r5"),
     "c":("Ca 0.05 matrix",["vf_1","vf_2","vf_3","vf_25_1","vf_20_1"],"glass0_r1"),
     "d":("Ca 0.07 matrix",["vf_r6_1","vf_r6_2","vf_r6_3"],"glass0_r6")}
labs={"vf_r5_1":"5%","vf_r5_2":"10%","vf_r5_3":"15%","vf_1":"5%","vf_2":"10%","vf_3":"15%","vf_25_1":"18.5%","vf_20_1":"18.9%","vf_r6_1":"5%","vf_r6_2":"10%","vf_r6_3":"15%"}
for ax,(k,(t,cs,gl)) in zip(A[1:4],vfs.items()):
    S.panel(ax,k,f"Volume fraction, {t}")
    draw(ax,gl,S.GLASS,"0% (glass)",ls="--",lw=1.1)
    cols=S.VF5 if len(cs)==5 else S.VF5[:3]
    for c,col in zip(cs,cols): draw(ax,c,col,labs[c])
# (e) seed radius at 10 %
S.panel(A[4],"e","Seed radius, 10% volume fraction")
for c,r in (("vf_2","12"),("vf_rad15_1","15"),("vf_rad20_1","20")): draw(A[4],c,S.SEED[r],f"{r} Å")
# (f) rate
S.panel(A[5],"f","Strain rate, 300 K")
for c,k,l in (("rate1e8_1","1e8","1×10$^8$ s$^{-1}$"),("rate5e8_1","5e8","5×10$^8$ s$^{-1}$"),("run_1","1e9","1×10$^9$ s$^{-1}$"),("rate2e9_1","2e9","2×10$^9$ s$^{-1}$")): draw(A[5],c,S.RATE[k],l)
# (g) temperature
S.panel(A[6],"g","Temperature, 10$^9$ s$^{-1}$")
for c,T in (("temp_250",250),("temp_300",300),("temp_350",350),("temp_400",400)): draw(A[6],c,S.TEMP[T],f"{T} K")
# (h) perfect Mg crystal
S.panel(A[7],"h","Perfect Mg crystal, 300 K")
draw(A[7],"hcp_z","#333333","[0001]"); draw(A[7],"hcp_x","#888888","[2$\\bar{1}\\bar{1}$0]")
for i,ax in enumerate(A):
    ax.set_xlim(0,0.15); ax.set_xticks([0,0.05,0.10,0.15]); ax.set_ylim(-0.02 if i!=7 else 0,None)
    ax.set_ylabel("Tensile stress (GPa)" if i%2==0 else "")
    if i>=6: ax.set_xlabel("Engineering strain")
    ax.legend(loc="upper left" if i==7 else "lower right",handlelength=1.6,labelspacing=0.25,borderaxespad=0.3,ncol=1)
for ax in A[:6]: ax.set_xlabel("Engineering strain")
fig.tight_layout(h_pad=1.0,w_pad=1.2)
S.save(fig,"fig05_stress_strain")
