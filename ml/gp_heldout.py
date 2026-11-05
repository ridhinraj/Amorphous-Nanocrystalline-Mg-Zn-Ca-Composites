"""gp_heldout.py -- frozen-surrogate held-out test for the new MD cases.
Reconstructs the training sets of the draft (13-point, 16-point, 25-run) from Table S1, checks the draft's own
calibration numbers, then predicts the new cases. Retrospective: predictions are made after the new runs finished,
from models fitted only on earlier data; nothing is tuned on the new cases."""
import warnings, numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.exceptions import ConvergenceWarning
import gp_active_learning as G
warnings.simplefilter("ignore", category=ConvergenceWarning)
OLD=[ # case, xZn, xCa, VF%, r, rate, T, sigma_y
("run_1",.35,.05,.74,25,1e9,300,1.165),("run_2",.30,.04,.38,20,1e9,300,1.084),("run_3",.40,.06,1.28,30,1e9,300,1.253),
("run_4",.25,.05,.74,25,1e9,300,1.105),("run_5",.35,.03,.16,15,1e9,300,1.074),("run_6",.35,.07,2.02,35,1e9,300,1.249),
("run_7",.30,.06,.74,25,1e9,300,1.178),("run_8",.40,.04,.38,20,1e9,300,1.158),("run_9",.32,.05,1.04,28,1e9,300,1.150),
("run_10",.38,.05,.50,22,1e9,300,1.174),("vf_1",.35,.05,5.00,12,1e9,300,1.094),("vf_2",.35,.05,10.01,12,1e9,300,1.035),
("vf_3",.35,.05,15.01,12,1e9,300,0.970),("vf_20_1",.35,.05,18.91,12,1e9,300,0.928),("vf_25_1",.35,.05,18.49,12,1e9,300,0.917),
("vf_r6_1",.35,.07,5.04,12,1e9,300,1.179),("vf_r6_2",.35,.07,10.00,12,1e9,300,1.102),("vf_r6_3",.35,.07,15.04,12,1e9,300,1.035),
("vf_rad15_1",.35,.05,9.94,15,1e9,300,1.056),("vf_rad20_1",.35,.05,9.81,20,1e9,300,1.097),("t310_1",.35,.05,.74,25,1e9,310,1.137),
("rate1e8_1",.35,.05,.74,25,1e8,300,0.864),("rate5e8_1",.35,.05,.74,25,5e8,300,1.059),("rate2e9_1",.35,.05,.74,25,2e9,300,1.286),
("rate2e9_run6",.35,.07,2.02,35,2e9,300,1.377)]
def mk(rows):
    d=pd.DataFrame(rows,columns=["case","xZn","xCa","vf","seed_radius","strain_rate","temperature_K","y"])
    d["volume_fraction"]=d.vf/100; d["loading_mode_shear"]=0.0; return d.set_index("case",drop=False)
old=mk(OLD)
S13=["run_%d"%i for i in range(1,11)]+["vf_1","vf_2","vf_3"]
S16=S13+["vf_20_1","vf_25_1","t310_1"]
S25=list(old.index)
def fit(d):
    X=d[G.FEATURE_COLUMNS].values; sc=StandardScaler().fit(X)
    gp=G.build_gp(); gp.fit(sc.transform(X),d.y.values); return gp,sc
def pred(gp,sc,d):
    m,s=gp.predict(sc.transform(d[G.FEATURE_COLUMNS].values),return_std=True); return m,s
if __name__=="__main__":
    np.random.seed(0)
    # 1. reproduce the draft's calibration
    d16=old.loc[S16]; n=len(d16); z=[]
    for i in range(n):
        tr=d16.drop(d16.index[i]); gp,sc=fit(tr); m,s=pred(gp,sc,d16.iloc[[i]]); z.append((d16.y.iloc[i]-m[0])/max(s[0],1e-8))
    z=np.array(z); print("LOO 16-point: std z=%.2f, within 1 sigma %d/16, within 1.96 sigma %d/16"%(z.std(),(abs(z)<=1).sum(),(abs(z)<=1.96).sum()))
    print("  largest |z|:",[(d16.index[i],round(z[i],2)) for i in np.argsort(-abs(z))[:3]])
    gp,sc=fit(old.loc[S13]); m,s=pred(gp,sc,old.loc[["vf_20_1","vf_25_1"]]); print("13-pt -> vf_20_1 %.3f±%.3f, vf_25_1 %.3f±%.3f (draft 0.926±0.011, 0.931±0.010)"%(m[0],s[0],m[1],s[1]))
    gp,sc=fit(old.loc[S16]); m,s=pred(gp,sc,old.loc[["rate2e9_1"]]); print("16-pt -> rate2e9_1 %.3f±%.3f (draft 1.104±1.569)"%(m[0],s[0]))
