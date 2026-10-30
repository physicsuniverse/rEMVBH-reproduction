#!/usr/bin/env python3
"""① extract M,J,chi,Omega_H,T_H,A_H,S,Phi,Q from the raw fields -- NOT from stored header values.
Outputs verify/out/extracted.csv"""
import sys,os,csv,glob,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import _common as C
OUT=f'{C.ROOT}/verify/out'; os.makedirs(OUT,exist_ok=True)
x,th=C.read_grids(); rows=[]
for p in sorted(glob.glob(f'{C.ROOT}/solutions/field_profiles/sol_*.txt')):
    par,F=C.read_profile(p); rh=par['rh']
    q=C.quantities_full(x,th,rh,F['f'],F['g'],F['h'],F['W'],F['At'],F['Bt'])
    rows.append(dict(source=os.path.basename(p),kind='profile',**{k:q[k] for k in
        ('M','J','chi','Omega_H','T_H','A_H','S','Phi','Q','Q_over_M','g_H','TH_AH_identity_rel')}))
for p in sorted(glob.glob(f'{C.ROOT}/solutions/metric_files/*.dat')):
    par,F,xm,thm=C.read_metric(p)
    q=C.quantities_from_fghW(xm,thm,1.0,F['f'],F['g'],F['h'],F['W'])
    q['chi']=q['J']/q['M']**2
    rows.append(dict(source=os.path.basename(p),kind='metric',**{k:q[k] for k in
        ('M','J','chi','Omega_H','T_H','A_H','S','g_H','TH_AH_identity_rel')}))
keys=sorted({k for r in rows for k in r},key=lambda k:(k!='source',k!='kind',k))
with open(f'{OUT}/extracted.csv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=keys);w.writeheader()
    for r in rows: w.writerow({k:(f"{r[k]:.12g}" if isinstance(r.get(k),float) else r.get(k,'')) for k in keys})
print(f"wrote {OUT}/extracted.csv  ({len(rows)} solutions)")
for r in rows:
    th_=f" TH_AH_id={r['TH_AH_identity_rel']:.2e}" if 'TH_AH_identity_rel' in r else ''
    print(f"  {r['source']:34s} M={r['M']:.9g} J={r['J']:.9g} chi={r['chi']:.6f} OmH={r['Omega_H']:.7g} TH={r['T_H']:.7g} AH={r['A_H']:.8g}{th_}")
