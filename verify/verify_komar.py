#!/usr/bin/env python3
"""③ Komar: M_V = 2 Omega_H J_V.  Re-arithmetic from residuals/komar_residuals.csv
(the from-B-field Komar integral is BLOCKED: no Python integrand module exists yet --
equations/ contains only .md, no eom.py / komar_integrand.py). Output verify/out/komar_check.csv"""
import csv,json,os,numpy as np
R=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT=f'{R}/verify/out'
X=list(csv.DictReader(open(f'{R}/residuals/komar_residuals.csv')))
out=[]
for x in X:
    MV,JV,Om=[float(x[k]) for k in ('MV','JV','OmegaH')]
    two=2*Om*JV; rel=abs(MV-two)/abs(MV) if MV else float('nan')
    out.append(dict(family=x['family'],idx=x['idx'],MV=MV,JV=JV,OmegaH=Om,twoOmJV=two,
                    rel_recomputed=rel,rel_stored=float(x['rel_residual']),
                    agrees=abs(rel-float(x['rel_residual']))<1e-5))
with open(f'{OUT}/komar_check.csv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(out[0].keys()));w.writeheader()
    for x in out: w.writerow({k:(f"{v:.10g}" if isinstance(v,float) else v) for k,v in x.items()})
rel=np.array([x['rel_recomputed'] for x in out])
fam={}
for x in out: fam.setdefault(x['family'],[]).append(x['rel_recomputed'])
print("③ M_V = 2 Omega_H J_V : %d rows; median rel=%.3e  p90=%.3e  max=%.3e"%(len(rel),np.median(rel),np.percentile(rel,90),rel.max()))
for k,v in fam.items(): print("   family %-28s n=%3d median=%.3e max=%.3e"%(k,len(v),np.median(v),max(v)))
print("   stored rel_residual agreement: %.1f%% rows"%(100*np.mean([x['agrees'] for x in out])))
print("   NOTE: from-B-field Komar integral NOT recomputed (no equations/eom.py or integrand module in this workspace)")
json.dump(dict(n=len(rel),median=float(np.median(rel)),p90=float(np.percentile(rel,90)),max=float(rel.max()),
               family_medians={k:float(np.median(v)) for k,v in fam.items()},
               block="no equations/eom.py; from-field Komar integral unavailable"),
          open(f'{OUT}/komar_check.json','w'),indent=1)
