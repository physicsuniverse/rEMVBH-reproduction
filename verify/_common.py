"""Shared readers/formulas. Reads DATA (grids, fields) but NOT the stored physical-quantity values."""
import numpy as np, os, glob
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def read_grids():
    """xGrid (Nx) and thetaGrid (Ny) from solutions/field_profiles/grids.txt"""
    L=[l for l in open(f'{ROOT}/solutions/field_profiles/grids.txt') if l.strip() and not l.startswith('#')]
    x=np.array([float(v) for v in L[0].split()]); th=np.array([float(v) for v in L[1].split()])
    return x,th

def read_profile(path):
    """returns (params Nx,Ny,rh ; fields dict Bt,Bp,At,Ap,f,g,h,W) -- skips the stored physical values (line 4)"""
    L=open(path).read().split('\n')
    p=[float(v) for v in L[1].split()]      # q chi alpha rh Nx Ny (structural only)
    nx,ny,rh=int(p[4]),int(p[5]),p[3]
    names=['Bt','Bp','At','Ap','f','g','h','W']
    F={}
    for k,nm in enumerate(names):
        v=np.array([float(t) for t in L[5+k].split()])
        assert v.size==nx*ny,(nm,v.size,nx*ny)
        F[nm]=v.reshape(nx,ny)              # x slowest per format line
    return dict(nx=nx,ny=ny,rh=rh),F

def read_metric(path):
    """metric .dat: line1 = {q,chi,M} (ignored); line2 = one brace group with
       Nx + Ny + 4*Nx*Ny values: xGrid, thetaGrid, f, g, h, W  (x slowest)."""
    txt=open(path).read(); g=[t for t in txt.split('}') if '{' in t]
    v=np.array([float(x) for x in g[-1].split('{')[-1].split(',') if x.strip()])
    import re
    m=re.search(r'_(\d+)x(\d+)',path)
    if m: nx,ny=int(m.group(1)),int(m.group(2))
    else:
        nx,ny=100,50                       # known grid for the two metric files in this workspace
    assert v.size==nx+ny+4*nx*ny,(v.size,nx,ny)
    x=v[:nx]; th=v[nx:nx+ny]; rest=v[nx+ny:].reshape(4,nx,ny)
    F={nm:rest[k] for k,nm in enumerate(['f','g','h','W'])}
    return dict(nx=nx,ny=ny),F,x,th

def horizon(iv,rh,y):
    """6-point polynomial extrapolation of a horizon quantity in delta=1+x"""
    V=np.vander(iv[:6],6,increasing=True)
    return np.array([np.linalg.solve(V,col[:6])[0] for col in y.T])

def thetavec(ny):
    return np.pi*(2*np.arange(ny)+1)/(4*ny)

def quantities_from_fghW(x,th,rh,f,g,h,W):
    r=2.0*rh/(1.0-x); t=1.0/r; s2=np.sin(th)[None,:]**2
    N=1.0-rh/r[:,None]; om=(W/r[:,None])*(rh/r[:,None])
    gtt=-f*N**2+(g/f)*r[:,None]**2*s2*om**2
    gtp=-(g/f)*rh*W*s2
    NP=8; out={}
    def fit1(v):
        A=np.vander(t[-NP:],NP,increasing=True); c=np.linalg.solve(A,v[-NP:]); return c[0],c[1]
    M=[];J=[]
    for j in range(len(th)):
        a0,a1=fit1((gtt[:,j]+1.0)*r/2.0); M.append(a0)
        b0,b1=fit1(gtp[:,j]*r/(2.0*s2[0,j])); J.append(-b0)
    out['M']=float(np.mean(M)); out['J']=float(np.mean(J))
    dh=1.0+x
    fH=horizon(dh,rh,f); gH=horizon(dh,rh,g); hH=horizon(dh,rh,h); WH=horizon(dh,rh,W)
    gH=float(np.mean(gH)); out['Omega_H']=float(np.mean(WH)/rh); out['g_H']=gH
    out['T_H']=float(np.mean(fH/(2*np.pi*rh*np.sqrt(gH*hH))))
    # spectral theta quadrature of sin(th) g sqrt(h)/f at the horizon
    z=2*th/(np.pi/2)-1; integ=np.sin(th)*gH*np.sqrt(hH)/fH; qs=[]
    for deg in (8,10,12,15,18):
        if deg>len(th)-2: continue
        c=np.polynomial.chebyshev.chebfit(z,integ,deg); s=0.0
        for k,ck in enumerate(c):
            Ik=2.0 if k==0 else (0.0 if k%2 else 2.0/(1.0-k*k)); s+=ck*Ik*np.pi/4.0
        qs.append(s)
    out['A_H']=float(2*np.pi*rh**2*2.0*np.mean(qs)); out['S']=out['A_H']/4.0
    out['TH_AH_identity_rel']=abs(out['T_H']*out['A_H']/(2*rh*np.sqrt(gH))-1.0)
    return out

def quantities_full(x,th,rh,f,g,h,W,At,Bt):
    o=quantities_from_fghW(x,th,rh,f,g,h,W)
    r=2.0*rh/(1.0-x); t=1.0/r; NP=8
    def fit(v):
        A=np.vander(t[-NP:],NP,increasing=True); return np.linalg.solve(A,v[-NP:])
    ph=[];Q=[]
    for j in range(len(th)):
        c=fit(At[:,j]); ph.append(c[0]); Q.append(-c[1])
    o['Phi']=float(np.mean(ph)); o['Q']=float(np.mean(Q)); o['chi']=o['J']/o['M']**2
    o['Q_over_M']=o['Q']/o['M']
    return o
