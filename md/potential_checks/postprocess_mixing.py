#!/usr/bin/env python3
"""Run in the folder holding mix_<TAG>_s<SEED>.log and mixblocks_<TAG>_s<SEED>_<T>K.dat for TAG in Mg, Zn, Ca, MgZnCa
(and optionally Ca_oldprotocol).  Prints, for each temperature,
  * per-reference liquid check (volume per atom and mean-squared displacement after 20 ps; liquid if D = MSD/(6t) > 0.1 A^2/ps, undercooled liquid 0.03-0.1, otherwise not liquid; volume > 90 A^3/atom = gas)
  * dH_mix = E_ternary - sum_i x_i E_i  (x_i = ACTUAL atom fractions of the ternary run), kJ/mol,
    with the error from the block means (standard error of 4 blocks) and, if several seeds exist, the spread of seeds."""
import glob, re, math, statistics as st
EV_TO_KJMOL = 96.4853
T_PROD_PS = 20.0
def load(tag, seed, T):
    log=open(f"mix_{tag}_s{seed}.log").read()
    cnt=re.search(r"^COUNTS .*nMg=(\d+) nZn=(\d+) nCa=(\d+) N=(\d+)",log,re.M)
    n=[int(cnt.group(i)) for i in (1,2,3,4)]
    pt=[l for l in log.splitlines() if l.startswith("MIXPT") and f"T={T:.1f}" in l][-1]
    vol=float(re.search(r"vol_per_atom=([\d.]+)",pt).group(1)); msd=float(re.search(r"msd_A2=([\d.]+)",pt).group(1))
    rows=[l.split() for l in open(f"mixblocks_{tag}_s{seed}_{T:.1f}K.dat") if not l.startswith("#")]
    pe=[float(r[1]) for r in rows]
    if not pe: raise ValueError(f'no block means in mixblocks_{tag}_s{seed}_{T:.1f}K.dat')
    return n,vol,msd,pe
def mean_se(x): 
    m=sum(x)/len(x); se=(st.stdev(x)/math.sqrt(len(x))) if len(x)>1 else float("nan"); return m,se
seeds=sorted({int(re.search(r"_s(\d+)\.log",f).group(1)) for f in glob.glob("mix_MgZnCa_s*.log")})
for T in (1200,1000):
    print(f"\n===== T = {T} K =====")
    res=[]
    for s in seeds:
        try:
            ref={}
            for tag in ("Mg","Zn","Ca"):
                n,vol,msd,pe=load(tag,s,T); ref[tag]=(mean_se(pe),vol,msd)
                D=msd/(6*T_PROD_PS)
                if vol>90: status="GAS (evaporated): INVALID"
                elif D>0.1: status="LIQUID"
                elif D>0.03: status="undercooled liquid (slow; a crystal gives MSD ~ 0.5 A^2)"
                else: status="NOT LIQUID (crystalline/glassy)"
                bad = vol>90 or D<=0.03
                print(f" seed {s} pure {tag:2s}: E = {ref[tag][0][0]:.5f} eV/atom, V = {vol:.2f} A^3/atom, MSD(20 ps) = {msd:.1f} A^2, D = {D:.3f} A^2/ps -> {status}")
                if bad: raise ValueError(f"reference {tag} is not a liquid at {T} K")
            n,vol,msd,pe=load("MgZnCa",s,T); (Et,Ese)=mean_se(pe)
            x=[n[0]/n[3],n[1]/n[3],n[2]/n[3]]
            print(f" seed {s} ternary: x = Mg {x[0]:.4f} Zn {x[1]:.4f} Ca {x[2]:.4f}, E = {Et:.5f} eV/atom, MSD = {msd:.1f} A^2")
            E=[ref[t][0][0] for t in ("Mg","Zn","Ca")]; SE=[ref[t][0][1] for t in ("Mg","Zn","Ca")]
            dH=(Et-sum(xi*e for xi,e in zip(x,E)))*EV_TO_KJMOL
            se=math.sqrt(Ese**2+sum((xi*e)**2 for xi,e in zip(x,SE)))*EV_TO_KJMOL
            print(f" seed {s}: dH_mix = {dH:.2f} +- {se:.2f} kJ/mol (block standard error)")
            res.append(dH)
        except (FileNotFoundError,IndexError,AttributeError,ValueError) as e:
            print(f" seed {s}: incomplete ({e})")
    if len(res)>1: print(f" ALL SEEDS: mean {st.mean(res):.2f} kJ/mol, sd {st.stdev(res):.2f}, n={len(res)}")
if glob.glob("mix_Ca_oldprotocol_s*.log"):
    print("\n===== old protocol for pure Ca (fcc start, 1000 K, no melt) =====")
    for f in sorted(glob.glob("mix_Ca_oldprotocol_s*.log")):
        s=int(re.search(r"_s(\d+)\.log",f).group(1)); n,vol,msd,pe=load("Ca_oldprotocol",s,1000)
        print(f" seed {s}: E = {mean_se(pe)[0]:.5f} eV/atom, V = {vol:.2f}, MSD(20 ps) = {msd:.1f} A^2 -> {'liquid' if msd/(6*T_PROD_PS)>0.1 else 'CRYSTALLINE (superheated): old Ca reference was not a liquid'}")
