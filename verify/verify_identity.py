#!/usr/bin/env python3
"""④ T_H A_H = 2 r_H sqrt(g_H) (=4 r_H for r_H=1, g_H=4), exact for the paper's ansatz.
Uses only independently extracted quantities. Output verify/out/identity_check.csv"""
import csv,json,os,numpy as np
R=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT=f'{R}/verify/out'
rows=list(csv.DictReader(open(f'{OUT}/extracted.csv')))
res=[]
for r in rows:
    rel=float(r['TH_AH_identity_rel']); gH=float(r['g_H'])
    res.append(dict(source=r['source'],g_H=gH,rel_dev=rel,
                    verdict=('qualified' if rel<1e-6 else ('grey' if rel<1e-2 else 'FAIL'))))
with open(f'{OUT}/identity_check.csv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(res[0].keys()));w.writeheader()
    for x in res: w.writerow({k:(f"{v:.10g}" if isinstance(v,float) else v) for k,v in x.items()})
q=[x['rel_dev'] for x in res if x['verdict']=='qualified']
print("④ T_H A_H = 4 r_H identity:")
for x in res: print("   %-34s g_H=%.10f rel_dev=%.3e  %s"%(x['source'],x['g_H'],x['rel_dev'],x['verdict']))
print("   qualified median = %.3e (criterion ~7e-9); FAIL count = %d"%(np.median(q) if q else float('nan'),sum(1 for x in res if x['verdict']=='FAIL')))
ref=list(csv.DictReader(open(f'{R}/solutions/physical_quantities.csv')))
col=[float(r['TH_AH_rel_dev']) for r in ref if r.get('TH_AH_rel_dev')]
print("   cross-check solutions/physical_quantities.csv TH_AH_rel_dev: n=%d median=%.3e max=%.3e"%(len(col),np.median(col),max(col)))
json.dump(dict(n=len(res),qualified_median=float(np.median(q)) if q else None,
               fails=[x['source'] for x in res if x['verdict']=='FAIL'],
               ref_n=len(col),ref_median=float(np.median(col)),ref_max=float(max(col))),
          open(f'{OUT}/identity_check.json','w'),indent=1)
