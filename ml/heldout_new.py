import re, warnings, numpy as np, pandas as pd
import gp_heldout as H
np.random.seed(0)
fin=pd.read_csv("/mnt/user-data/outputs/new_md_properties_final.csv").set_index("case")
def vf(case):
    for r in ("/tmp/f1","/tmp/f2"):
        try:
            t=open(f"{r}/{case}/log.2_build").read(); return float(re.search(r"achieved VF=([0-9.]+)",t).group(1))
        except Exception: pass
    return None
COMP={"r5":(.35,.03),"r1":(.35,.05),"r6":(.35,.07)}
rows=[]
for c in fin.index:
    if c.startswith("hcp"): continue
    y=fin.loc[c,"peak_stress_GPa"]
    if c.startswith("glass0"):  m=c[-2:]; rows.append((c,"M16","VF=0 (no seed)",*COMP[m],0.0,0.0,1e9,300,y))
    elif c.startswith("vf_r5"): rows.append((c,"M16","Ca 0.03, 12 A",*COMP["r5"],vf(c),12,1e9,300,y))
    elif c.startswith("rep_"):  m=c[4:6]; rows.append((c,"M16","10 pct VF, 12 A, Ca {:.2f}".format(COMP[m][1]),*COMP[m],vf(c),12,1e9,300,y))
    elif c.startswith("temp_"): T=float(c[5:]); rows.append((c,"M25","T series",.35,.05,.0074,25,1e9,T,y))
    elif c=="rate2e9_run4": rows.append((c,"M25","2e9 at comp 4",.25,.05,.0074,25,2e9,300,y))
    elif c=="rate2e9_run3": rows.append((c,"M25","2e9 at comp 3",.40,.06,.0128,30,2e9,300,y))
df=pd.DataFrame(rows,columns=["case","model","group","xZn","xCa","volume_fraction","seed_radius","strain_rate","temperature_K","y"])
df["volume_fraction"]=df.volume_fraction.astype(float)
df.loc[df.model!="M16",["volume_fraction"]]=df.loc[df.model!="M16","volume_fraction"]  # already fractions
df["loading_mode_shear"]=0.0
# M16 / M25 volume fraction for old-style rows are fractions already (e.g. .0074); new rows from logs are fractions too
out=[]
gp16,sc16=H.fit(H.old.loc[H.S16]); gp25,sc25=H.fit(H.old.loc[H.S25])
for r in df.itertuples():
    gp,sc=(gp16,sc16) if r.model=="M16" else (gp25,sc25)
    m,s=H.pred(gp,sc,df.loc[[r.Index]]); z=(r.y-m[0])/s[0]
    out.append(dict(case=r.case,group=r.group,model="16-run" if r.model=="M16" else "25-run",pred=round(m[0],3),sigma=round(s[0],3),simulated=round(r.y,3),err_pct=round(100*(r.y-m[0])/m[0],1),z=round(z,2),inside95=abs(z)<=1.96))
o=pd.DataFrame(out); pd.set_option("display.width",200); print(o.to_string(index=False))
o.to_csv("/tmp/a1/heldout_predictions.csv",index=False)
print("\nsummary by group (mean abs % error, coverage 95%):")
for g,x in o.groupby("group"): print(" %-24s n=%d  MAE %.1f%%  inside95 %d/%d  mean sigma %.3f"%(g,len(x),x.err_pct.abs().mean(),x.inside95.sum(),len(x),x.sigma.mean()))
