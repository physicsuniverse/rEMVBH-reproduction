#!/usr/bin/env python3
r"""
equations/eom.py  --  field equations of the EMV model on the rotating background,
                      in pure Python (sympy for the algebra, numpy for evaluation).

Action (paper Eq. action; 16 pi G = 1):
    S = \int sqrt(-g) [ R - U(|B|^2) F_ab F^ab - V_ab V^ab ],
    U = exp(alpha |B|^2),   V_ab = d_a B_b - d_b B_a,   F_ab = d_a A_b - d_b A_a.

Field equations
    Einstein :  G_ab = 2 ( T^F_ab + T^V_ab )
    Maxwell  :  grad_a ( U F^{ab} ) = 0
    vector   :  grad_a V^{ab} = 1/2 U' F_cd F^cd B^b

    T^F_ab = U ( F_a{}^c F_bc - 1/4 F_cd F^cd g_ab )
    T^V_ab = V_a{}^c V_bc - 1/4 V_cd V^cd g_ab + 1/2 U' F_cd F^cd B_a B_b

Background ansatz (paper Eqs. metricansatz / field ansatz):
    ds^2 = -f N^2 dt^2 + (g/f)[ h(dr^2 + r^2 dth^2)
                               + r^2 sin^2(th) (dphi - (W/r)(1-N) dt)^2 ]
    A_mu = ( A_t - (W/r)(1-N) A_phi sin^2 th, 0, 0, A_phi sin^2 th )
    B_mu = ( B_t - (W/r)(1-N) B_phi sin^2 th, 0, 0, B_phi sin^2 th )
    x = 1 - 2 r_H/r  =>  r = 2 r_H/(1-x),  N = (1+x)/2,  1-N = (1-x)/2

Design
------
* `symbolic_field_strengths()` returns sympy expressions for F_ab, V_ab, F^2, B^2,
  U, U', T^F_ab, T^V_ab -- all low-order, cheap to lambdify.
* curvature is assembled **numerically** from the metric array with spectral
  derivatives (exact for the spectral data), via the standard formulas.
* `residuals(profile)` returns the three residual vectors, normalised.

Usage
    python3 equations/eom.py <profile.txt>          # one solution
    python3 equations/eom.py --all                  # every solutions/field_profiles/sol_*.txt
"""
import argparse
import glob
import json
import os

import numpy as np
import sympy as sp

# --------------------------------------------------------------------------
# symbolic part
# --------------------------------------------------------------------------
t, x, th, ph = sp.symbols("t x theta phi", real=True)
alpha = sp.symbols("alpha", positive=True)
rh = sp.symbols("r_H", positive=True)

At = sp.Function("At")(x, th)
Ap = sp.Function("Ap")(x, th)
Bt = sp.Function("Bt")(x, th)
Bp = sp.Function("Bp")(x, th)
f = sp.Function("f")(x, th)
g = sp.Function("g")(x, th)
h = sp.Function("h")(x, th)
W = sp.Function("W")(x, th)

COORDS = (t, x, th, ph)


def _setup():
    r = 2 * rh / (1 - x)
    N = sp.Rational(1, 2) * (1 + x)
    omN = sp.Rational(1, 2) * (1 - x)          # 1 - N
    s2 = sp.sin(th) ** 2
    shift = (W / r) * omN

    g = sp.zeros(4, 4)
    g[0, 0] = -f * N ** 2 + (g_ / f) * r * r * s2 * shift ** 2
    return None


def background_metric():
    """covariant metric g_{ab} for the ansatz (see module docstring)"""
    r = 2 * rh / (1 - x)
    N = sp.Rational(1, 2) * (1 + x)
    omN = sp.Rational(1, 2) * (1 - x)
    s2 = sp.sin(th) ** 2
    drdx = sp.diff(r, x)
    M = sp.zeros(4, 4)
    M[0, 0] = -f * N ** 2 + (g / f) * W ** 2 * s2 * omN ** 2
    M[0, 3] = M[3, 0] = -(g / f) * W * r * s2 * omN
    M[1, 1] = (g / f) * h * drdx ** 2
    M[2, 2] = (g / f) * h * r ** 2
    M[3, 3] = (g / f) * r ** 2 * s2
    return M


def vector_potentials():
    """A_mu and B_mu one-form components (0=t, 1=x, 2=th, 3=phi)"""
    r = 2 * rh / (1 - x)
    omN = sp.Rational(1, 2) * (1 - x)
    s2 = sp.sin(th) ** 2
    A = [At - (W / r) * omN * Ap * s2, sp.Integer(0), sp.Integer(0), Ap * s2]
    B = [Bt - (W / r) * omN * Bp * s2, sp.Integer(0), sp.Integer(0), Bp * s2]
    return A, B


def _d(expr, mu):
    return sp.diff(expr, COORDS[mu])


def field_strengths():
    """Return F_ab, V_ab (antisymmetric 4x4 matrices) symbolically."""
    A, B = vector_potentials()
    F = sp.zeros(4, 4)
    V = sp.zeros(4, 4)
    for a in range(4):
        for b in range(a + 1, 4):
            F[a, b] = _d(A[b], a) - _d(A[a], b)
            F[b, a] = -F[a, b]
            V[a, b] = _d(B[b], a) - _d(B[a], b)
            V[b, a] = -V[a, b]
    return F, V


def invariants():
    """F^2 = F_ab F^ab, B^2 = B_a B^a, U, U', 1/2 U' F^2"""
    g = background_metric()
    ginv = g.inv()
    F, V = field_strengths()
    A, B = vector_potentials()
    # raise indices
    Fup = sp.Matrix(4, 4, lambda a, b: sum(
        ginv[a, c] * ginv[b, d] * F[c, d] for c in range(4) for d in range(4)))
    F2 = sp.simplify(sum(F[a, b] * Fup[a, b] for a in range(4) for b in range(4)))
    B2 = sp.simplify(sum(ginv[a, b] * B[a] * B[b] for a in range(4) for b in range(4)))
    U = sp.exp(alpha * B2)
    Up = sp.diff(U, sp.Symbol("__dummy")) if False else alpha * U
    return F, V, Fup, F2, B2, U, Up


def stress_energy():
    """T^F_ab and T^V_ab symbolically (low order -- cheap to lambdify)."""
    g = background_metric()
    ginv = g.inv()
    F, V = field_strengths()
    A, B = vector_potentials()
    n = 4
    Fup = sp.Matrix(n, n, lambda a, b: sum(
        ginv[a, c] * ginv[b, d] * F[c, d] for c in range(n) for d in range(n)))
    Vup = sp.Matrix(n, n, lambda a, b: sum(
        ginv[a, c] * ginv[b, d] * V[c, d] for c in range(n) for d in range(n)))
    F2 = sum(F[a, b] * Fup[a, b] for a in range(n) for b in range(n))
    V2 = sum(V[a, b] * Vup[a, b] for a in range(n) for b in range(n))
    B2 = sum(ginv[a, b] * B[a] * B[b] for a in range(n) for b in range(n))
    U = sp.exp(alpha * B2)
    Up = alpha * U
    TF = sp.Matrix(n, n, lambda a, b: sp.simplify(
        U * (sum(F[a, c] * Fup[b, c] if False else F[a, c] * ginv[b, d] * F[c, d]
                 for c in range(n) for d in range(n))
             - sp.Rational(1, 4) * F2 * g[a, b])))
    TV = sp.Matrix(n, n, lambda a, b: sp.simplify(
        sum(V[a, c] * ginv[b, d] * V[c, d] for c in range(n) for d in range(n))
        - sp.Rational(1, 4) * V2 * g[a, b]
        + sp.Rational(1, 2) * Up * F2 * B[a] * B[b]))
    return TF, TV, dict(F2=F2, V2=V2, B2=B2, U=U)


# --------------------------------------------------------------------------
# numerical machinery
# --------------------------------------------------------------------------
def _cheb_diff(n):
    """Chebyshev differentiation matrix on the Chebyshev-Gauss points"""
    xg = -np.cos((2 * np.arange(n) + 1) * np.pi / (2 * n))
    B = np.cos(np.arange(n)[None, :] * np.arccos(np.clip(xg, -1, 1))[:, None])
    dB = np.zeros_like(B)
    s = np.sqrt(np.maximum(1 - xg ** 2, 1e-300))
    for k in range(1, n):
        dB[:, k] = k * np.sin(k * np.arccos(np.clip(xg, -1, 1))) / s
    return xg, dB @ np.linalg.inv(B)


def _theta_diff(n, thg):
    B = np.cos(2.0 * np.arange(n)[None, :] * thg[:, None])
    dB = -2.0 * np.arange(n)[None, :] * np.sin(2.0 * np.arange(n)[None, :] * thg[:, None])
    return dB @ np.linalg.inv(B)


class Grid:
    """spectral differentiation on the (x, theta) tensor grid"""

    def __init__(self, xg, thg):
        self.xg = np.asarray(xg, float)
        self.thg = np.asarray(thg, float)
        _, self.Dx = _cheb_diff(len(xg))
        self.Dt = _theta_diff(len(thg), self.thg)

    def d(self, A, axis):
        if axis == 0:
            return self.Dx @ A
        return A @ self.Dt.T


def load_profile(path):
    """read sol_*.txt -> dict with meta, header quantities, 8 fields on (Nx,Ny)"""
    with open(path) as fh:
        lines = [ln for ln in fh.read().splitlines() if ln.strip()]
    meta = dict(zip(["q", "chi", "alpha", "rh", "Nx", "Ny"],
                    [float(v) for v in lines[1].split()]))
    hdr = dict(zip(["TH", "AH", "M", "Q", "Phi", "P", "J", "OmegaH"],
                   [float(v) for v in lines[3].split()]))
    Nx, Ny = int(meta["Nx"]), int(meta["Ny"])
    names = ["Bt", "Bp", "At", "Ap", "f", "g", "h", "W"]
    fields = {}
    for i, nm in enumerate(names):
        arr = np.array([float(v) for v in lines[5 + i].split()])
        if arr.size != Nx * Ny:
            raise ValueError(f"{path}: field {nm} has {arr.size} values, expected {Nx*Ny}")
        fields[nm] = arr.reshape(Nx, Ny, order="C")   # x slowest
    return meta, hdr, fields


def build_metric_numeric(grid, F, rh_v=1.0):
    """metric and its first/second derivatives as arrays (Nx,Ny,4,4)"""
    xg = grid.xg[:, None]
    th = grid.thg[None, :]
    ff, gg, hh, WW = F["f"], F["g"], F["h"], F["W"]
    r = 2 * rh_v / (1 - xg) * np.ones_like(ff)
    drdx = 2 * rh_v / (1 - xg) ** 2 * np.ones_like(ff)
    N = 0.5 * (1 + xg) * np.ones_like(ff)
    omN = 0.5 * (1 - xg) * np.ones_like(ff)
    s2 = np.sin(th) ** 2 * np.ones_like(ff)

    g = np.zeros(ff.shape + (4, 4))
    g[..., 0, 0] = -ff * N ** 2 + (gg / ff) * WW ** 2 * s2 * omN ** 2
    g[..., 0, 3] = g[..., 3, 0] = -(gg / ff) * WW * r * s2 * omN
    g[..., 1, 1] = (gg / ff) * hh * drdx ** 2
    g[..., 2, 2] = (gg / ff) * hh * r ** 2
    g[..., 3, 3] = (gg / ff) * r ** 2 * s2
    return g


def christoffel_numeric(grid, g):
    """Gamma[a,b,c] = 1/2 g^{ad}(d_b g_dc + d_c g_db - d_d g_bc)"""
    ginv = np.linalg.inv(g)
    # only x (1) and theta (2) derivatives are nonzero
    dg = np.zeros(g.shape + (4,))
    for i in range(4):
        for j in range(4):
            dg[..., i, j, 1] = grid.d(g[..., i, j], 0)
            dg[..., i, j, 2] = grid.d(g[..., i, j], 1)
    G = np.zeros(g.shape + (4,))
    for a in range(4):
        for b in range(4):
            for c in range(4):
                val = np.zeros_like(g[..., 0, 0])
                for d in range(4):
                    val += ginv[..., a, d] * (dg[..., d, c, b] + dg[..., d, b, c]
                                              - dg[..., b, c, d])
                G[..., a, b, c] = 0.5 * val
    return G


def riemann_numeric(grid, G):
    """R[a,b,c,d] = d_c G[a,b,d] - d_d G[a,b,c] + G[a,c,e]G[e,b,d] - G[a,d,e]G[e,b,c]"""
    n = 4
    dG = np.zeros(G.shape + (4,))
    for a in range(n):
        for b in range(n):
            for c in range(n):
                dG[..., a, b, c, 1] = grid.d(G[..., a, b, c], 0)
                dG[..., a, b, c, 2] = grid.d(G[..., a, b, c], 1)
    R = np.zeros(G.shape + (4,))
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    val = dG[..., a, b, d, c] - dG[..., a, b, c, d]
                    for e in range(n):
                        val = val + G[..., a, c, e] * G[..., e, b, d] \
                                  - G[..., a, d, e] * G[..., e, b, c]
                    R[..., a, b, c, d] = val
    return R


def curvature_numeric(grid, g):
    G = christoffel_numeric(grid, g)
    R = riemann_numeric(grid, G)
    Ric = np.zeros(g.shape)
    for b in range(4):
        for d in range(4):
            Ric[..., b, d] = sum(R[..., a, b, a, d] for a in range(4))
    ginv = np.linalg.inv(g)
    Rs = np.einsum("...ab,...ab->...", ginv, Ric)
    E = Ric - 0.5 * Rs[..., None, None] * g
    return G, Ric, Rs, E


# --------------------------------------------------------------------------
# residual assembly
# --------------------------------------------------------------------------
def residuals(profile_path, verbose=False):
    meta, hdr, F = load_profile(profile_path)
    Nx, Ny = int(meta["Nx"]), int(meta["Ny"])
    grid = Grid(_cheb_diff(Nx)[0], _theta_nodes(profile_path, Ny))
    rh_v = meta["rh"]

    g = build_metric_numeric(grid, F, rh_v)
    G, Ric, Rs, E = curvature_numeric(grid, g)
    ginv = np.linalg.inv(g)

    # --- field strengths (numeric, from the same ansatz) ---
    A_mu, B_mu = _potentials_numeric(grid, F, rh_v)
    dA = _deriv_1form(grid, A_mu)
    dB = _deriv_1form(grid, B_mu)
    Fmn = np.zeros(A_mu.shape + (4,))
    Vmn = np.zeros(A_mu.shape + (4,))
    for a in range(4):
        for b in range(4):
            Fmn[..., a, b] = dA[..., b, a] - dA[..., a, b]
            Vmn[..., a, b] = dB[..., b, a] - dB[..., a, b]
    # index raising
    Fup = np.einsum("...ac,...bd,...cd->...ab", ginv, ginv, Fmn)
    Vup = np.einsum("...ac,...bd,...cd->...ab", ginv, ginv, Vmn)
    F2 = np.einsum("...ab,...ab->...", Fmn, Fup)
    V2 = np.einsum("...ab,...ab->...", Vmn, Vup)
    B2 = np.einsum("...ab,...a,...b->...", ginv, B_mu, B_mu)
    U = np.exp(np.clip(meta["alpha"] * B2, -700, 700))
    Up = meta["alpha"] * U
    # stress tensors
    TF = np.zeros(Fmn.shape)
    for a in range(4):
        for b in range(4):
            TF[..., a, b] = U * (sum(Fmn[..., a, c] * Fup[..., b, c] for c in range(4))
                                 - 0.25 * F2 * g[..., a, b])
    TV = np.zeros(Fmn.shape)
    for a in range(4):
        for b in range(4):
            TV[..., a, b] = (sum(Vmn[..., a, c] * Vup[..., b, c] for c in range(4))
                             - 0.25 * V2 * g[..., a, b]
                             + 0.5 * Up * F2 * B_mu[..., a] * B_mu[..., b])
    T = TF + TV

    # --- residual 1: Einstein  G_ab - 2 T_ab ---
    Eres = E - 2 * T
    # --- residual 2: Maxwell  grad_a (U F^{ab}) ---
    UFup = U[..., None, None] * Fup
    divUF = _div2_numeric(grid, g, G, UFup)
    # --- residual 3: vector  grad_a V^{ab} - 1/2 U' F^2 B^b ---
    divV = _div2_numeric(grid, g, G, Vup)
    B_up = np.einsum("...ab,...b->...a", ginv, B_mu)
    Vres = divV - 0.5 * (Up * F2)[..., None] * B_up

    return dict(meta=meta, hdr=hdr, grid=grid, E=Eres, M=divUF, Vv=Vres,
                T=T, F2=F2, U=U, scale=_scale(E, T, divUF)), _indices(Nx, Ny)


def _indices(Nx, Ny):
    return np.arange(Nx), np.arange(Ny)


def _theta_nodes(profile_path, Ny):
    """the theta grid used by the profiles: cosine-type midpoints on [0, pi/2]"""
    return (2 * np.arange(Ny) + 1) * np.pi / (4 * Ny)


def _potentials_numeric(grid, F, rh_v):
    xg = grid.xg[:, None]
    th = grid.thg[None, :]
    ff, gg, hh, WW = F["f"], F["g"], F["h"], F["W"]
    r = 2 * rh_v / (1 - xg) * np.ones_like(ff)
    omN = 0.5 * (1 - xg) * np.ones_like(ff)
    s2 = np.sin(th) ** 2 * np.ones_like(ff)
    A = np.zeros(ff.shape + (4,))
    B = np.zeros(ff.shape + (4,))
    A[..., 0] = F["At"] - (WW / r) * omN * F["Ap"] * s2
    A[..., 3] = F["Ap"] * s2
    B[..., 0] = F["Bt"] - (WW / r) * omN * F["Bp"] * s2
    B[..., 3] = F["Bp"] * s2
    return A, B


def _deriv_1form(grid, A):
    dA = np.zeros(A.shape + (4,))
    for a in range(4):
        dA[..., a, 1] = grid.d(A[..., a], 0)
        dA[..., a, 2] = grid.d(A[..., a], 1)
    return dA


def _div2_numeric(grid, g, G, Sup):
    """div_a S^{ab} for an antisymmetric S^{ab} (contravariant)"""
    n = 4
    dS = np.zeros(Sup.shape + (4,))
    for a in range(n):
        for b in range(n):
            dS[..., a, b, 1] = grid.d(Sup[..., a, b], 0)
            dS[..., a, b, 2] = grid.d(Sup[..., a, b], 1)
    out = np.zeros(Sup.shape[:2] + (4,))
    for b in range(n):
        val = np.zeros(Sup.shape[:2])
        for a in range(n):
            val += dS[..., a, b, a]
            for c in range(n):
                val += G[..., a, a, c] * Sup[..., c, b] + G[..., b, a, c] * Sup[..., a, c]
        out[..., b] = val
    return out


def _scale(E, T, M):
    return dict(E_absmax=float(np.max(np.abs(E))),
                T_absmax=float(np.max(np.abs(T))),
                M_absmax=float(np.max(np.abs(M))),
                rel_E=float(np.max(np.abs(E)) / max(np.max(np.abs(T)), 1e-300)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("profile", nargs="?")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()
    root = os.environ.get(
        "EMVBH_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    files = ([args.profile] if args.profile else
             sorted(glob.glob(os.path.join(root, "solutions/field_profiles/sol_*.txt"))))
    if not files:
        print(f"[FAIL] no profiles found under {root}/solutions/field_profiles/")
        return 1
    rows = []
    for p in files:
        try:
            r, _ = residuals(p)
        except Exception as exc:            # noqa: BLE001
            print(f"[FAIL] {os.path.basename(p)}: {type(exc).__name__}: {exc}")
            continue
        s = r["scale"]
        print(f"[ok ] {os.path.basename(p):28s} |G-2T|max={s['E_absmax']:.3e} "
              f"|2T|max={s['T_absmax']:.3e} rel={s['rel_E']:.3e} "
              f"|div(UF)|max={s['M_absmax']:.3e}")
        rows.append(dict(file=os.path.basename(p), **{k: float(v) for k, v in r["meta"].items()},
                         **s))
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data",
                       "eom_residuals.csv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if rows:
        import csv
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"### wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
