#!/usr/bin/env python3
"""② Smarr: M = 2 T_H S + 2 Omega_H J + Q Phi, from INDEPENDENTLY extracted quantities (script 1).
Cross-checks residuals/smarr_residuals.csv. Output verify/out/smarr_check.csv / .json"""
import csv,json,os,numpy as np, sys
R=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT=f'{R}/verify/out'
rows=list(csv.DictReader(open(f'{OUT}/extracted.csv')))
res=[];have=0
for r in rows:
    f=lambda k: float(r[k]) if r.get(k) not in (None,'') else None
    M,TH,AH,S,Om,J=f('M'),f('T_H'),f('A_H'),f('S'),f('Omega_H'),f('J')
    Q,Phi=f('Q'),f('Phi')
    if None in (M,TH,S,Om,J,Q,Phi): continue
    have+=1
    lhs=M; rhs=2*TH*S+2*Om*J+Q*Phi; rel=abs(lhs-rhs)/abs(M)
    res.append(dict(source=r['source'],M=M,TH_S_term=2*TH*S,OmJ_term=2*Om*J,QPhi_term=Q*Phi,rhs=rhs,rel=rel))
with open(f'{OUT}/smarr_check.csv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(res[0].keys()));w.writeheader()
    for x in res: w.writerow({k:(f"{v:.10g}" if isinstance(v,float) else v) for k,v in x.items()})
# cross-check with residuals/smarr_residuals.csv
X=list(csv.DictReader(open(f'{R}/residuals/smarr_residuals.csv')))
ref=np.array([float(x['Smarr_reldiff_recomputed']) for x in X])
med=float(np.median([x['rel'] for x in res]))
print("② Smarr check: %d solutions with gauge fields"%have)
for x in res: print("   %-34s rel=%.3e"%(x['source'],x['rel']))
print("   median over these = %.3e   (criterion ~1e-11; BRIEF target)"%med)
print("   cross-check residuals/smarr_residuals.csv: n=%d median=%.3e p90=%.3e"%(len(ref),np.median(ref),np.percentile(ref,90)))
json.dump(dict(n_with_gauge=have,median_rel=med,n_ref=len(ref),ref_median=float(np.median(ref)),ref_p90=float(np.percentile(ref,90))),open(f'{OUT}/smarr_check.json','w'),indent=1)
