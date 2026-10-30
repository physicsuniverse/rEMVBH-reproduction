"""Shared helpers for the known-limit unit tests.

Pure Python: numpy / scipy / sympy.  No Mathematica, no dependency on the
other teammates' scripts (they may not have landed yet) -- this module is
self-contained and only reads files under ``solutions/`` and ``data/``.

Metric ansatz (verbatim from the paper, Eq. ``metricansatz``):

    ds^2 = -f N^2 dt^2 + (g/f) [ h (dr^2 + r^2 dtheta^2)
                                 + r^2 sin^2(theta) (dphi - W/r (1-N) dt)^2 ]

with ``N = 1 - r_H/r``.  Coordinates are quasi-isotropic.

The **exact vacuum solution** of this ansatz (verified symbolically in
``test_schwarzschild_limit.py``) is

    g = (1 + r_H/r)^2,   f = 1/g,   h = 1,   W = 0

for which the areal radius is ``R = r*g`` and ``g_tt * g_RR = -1`` exactly.
"""

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "tests", "_cache")
os.makedirs(CACHE, exist_ok=True)

# --------------------------------------------------------------------------
# 1. the EMV metric ansatz
# --------------------------------------------------------------------------


def emv_metric(f, g, h, W, r, theta, rh):
    """4x4 metric of the paper's ansatz at a single (r, theta).  (-,+,+,+) sign."""
    N = 1.0 - rh / r
    Phi = (g / f) * h
    Psi = f * N**2
    # dphi - (W/r)(1-N) dt   with (1-N) = rh/r
    om = (W / r) * (rh / r)
    s2 = np.sin(theta) ** 2
    gtt = -Psi + Phi * r**2 * s2 * om**2
    gtp = -Phi * r**2 * s2 * om
    gpp = Phi * r**2 * s2
    grr = Phi
    gthth = Phi * r**2
    return np.array(
        [
            [gtt, 0.0, 0.0, gtp],
            [0.0, grr, 0.0, 0.0],
            [0.0, 0.0, gthth, 0.0],
            [gtp, 0.0, 0.0, gpp],
        ]
    )


def schwarzschild_in_ansatz(r, theta, rh):
    """Exact vacuum solution written in the paper's ansatz: g=(1+r_H/r)^2, f=1/g, h=1."""
    g = (1.0 + rh / r) ** 2
    return emv_metric(f=1.0 / g, g=g, h=1.0, W=0.0, r=r, theta=theta, rh=rh)


def areal_radius(r, rh):
    """R = r sqrt(Phi) with sqrt(Phi) = g = (1+r_H/r)^2 for the vacuum solution."""
    return r * (1.0 + rh / r) ** 2


# --------------------------------------------------------------------------
# 2. symbolic Ricci (sympy), cached to disk
# --------------------------------------------------------------------------

_SYMBOLS = {}


def _symbols():
    if not _SYMBOLS:
        import sympy as sp

        t, r, th, ph = sp.symbols("t r theta phi", real=True)
        _SYMBOLS.update(
            sp=sp,
            t=t,
            r=r,
            th=th,
            ph=ph,
            X=[t, r, th, ph],
        )
    return _SYMBOLS


def ricci_symbolic(metric):
    """Ricci tensor of a 4x4 sympy Matrix.  Returns a 4x4 sympy Matrix."""
    S = _symbols()
    sp, X = S["sp"], S["X"]
    g = metric
    gi = g.inv()
    n = 4
    Gam = [[[0] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                Gam[a][b][c] = sp.simplify(
                    sum(
                        gi[a, d]
                        * (
                            sp.diff(g[d, b], X[c])
                            + sp.diff(g[d, c], X[b])
                            - sp.diff(g[b, c], X[d])
                        )
                        for d in range(n)
                    )
                    / 2
                )
    Ric = sp.zeros(n, n)
    for a in range(n):
        for b in range(n):
            s = 0
            for k in range(n):
                s += sp.diff(Gam[k][a][b], X[k]) - sp.diff(Gam[k][a][k], X[b])
                for l in range(n):
                    s += Gam[k][l][k] * Gam[l][a][b] - Gam[k][l][b] * Gam[l][a][k]
            Ric[a, b] = sp.simplify(s)
    return Ric


def cached_ricci(kind, builder):
    """Run ``builder()`` (returning a 4x4 sympy Matrix) and cache the diagonalised
    result as JSON: list of [i, j, str(component)] for the non-zero components."""
    path = os.path.join(CACHE, "ricci_%s.json" % kind)
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    t0 = time.time()
    Ric = builder()
    nz = []
    for i in range(4):
        for j in range(4):
            v = _symbols()["sp"].simplify(Ric[i, j])
            if v != 0:
                nz.append([i, j, str(v)])
    out = dict(
        kind=kind,
        n_nonzero=len(nz),
        nonzero=nz,
        seconds=round(time.time() - t0, 2),
    )
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    return out


# --------------------------------------------------------------------------
# 3. Kerr (Boyer-Lindquist) + shadow
# --------------------------------------------------------------------------


def kerr_bl_metric(r, theta, a, M=1.0):
    """Kerr in Boyer-Lindquist coordinates, (-,+,+,+)."""
    Sigma = r**2 + a**2 * np.cos(theta) ** 2
    Delta = r**2 - 2 * M * r + a**2
    s2 = np.sin(theta) ** 2
    g = np.zeros((4, 4))
    g[0, 0] = -(1 - 2 * M * r / Sigma)
    g[1, 1] = Sigma / Delta
    g[2, 2] = Sigma
    g[3, 3] = (r**2 + a**2 + 2 * M * r * a**2 * s2 / Sigma) * s2
    g[0, 3] = g[3, 0] = -2 * M * r * a * s2 / Sigma
    return g


def kerr_shadow_edge(a, M=1.0, n=20001):
    """Bardeen's exact shadow rim of Kerr for a distant equatorial observer (i=pi/2).

    Returns (alpha, beta) in units of M, the celestial coordinates of the rim.
    """
    from scipy.optimize import brentq

    def eta(r):
        return -r**3 * (r**3 - 6 * M * r**2 + 9 * M**2 * r - 4 * a**2 * M) / (
            a**2 * (r - M) ** 2
        )

    def xi(r):
        return (r**3 - 3 * M * r**2 + a**2 * r + a**2 * M) / (a * (r - M))

    r1 = 2 * M * (1 + np.cos(2.0 / 3.0 * np.arccos(-a / M)))
    r2 = 2 * M * (1 + np.cos(2.0 / 3.0 * np.arccos(a / M)))
    rlo, rhi = min(r1, r2), max(r1, r2)
    rs = np.linspace(rlo + 1e-9, rhi - 1e-9, 4001)
    ok = eta(rs) >= 0
    idx = np.where(ok)[0]
    lo, hi = rs[idx[0]], rs[idx[-1]]
    r = np.linspace(lo, hi, n)
    alpha = -xi(r)
    beta = np.sqrt(np.maximum(eta(r), 0.0))
    return alpha, beta


def ellipticity(alpha, beta):
    """Paper's definition: eps = H/W - 1, W = equatorial diameter, H = polar diameter."""
    W = float(np.max(alpha) - np.min(alpha))
    H = float(2.0 * np.max(np.abs(beta)))
    return H / W - 1.0


def circle_edge(R=5.196152422706632, n=4001):
    """A perfect circle: the reference case that must give eps == 0 exactly."""
    th = np.linspace(0.0, 2.0 * np.pi, n)
    return R * np.cos(th), R * np.sin(th)


def pixel_contour(mask):
    """Boundary pixels of a boolean mask, returned as (col, row) float arrays."""
    from scipy import ndimage
    m = np.asarray(mask, dtype=bool)
    er = ndimage.binary_erosion(m, border_value=0)
    bnd = m & ~er
    rows, cols = np.nonzero(bnd)
    return cols.astype(float), rows.astype(float)


# --------------------------------------------------------------------------
# 4. analytic Smarr checks
# --------------------------------------------------------------------------


def smarr_schwarzschild(M):
    rH = 2 * M
    TH = 1.0 / (8 * np.pi * M)
    AH = 4 * np.pi * rH**2
    S = AH / 4.0
    return dict(M=M, J=0.0, Q=0.0, Omega=0.0, Phi=0.0, TH=TH, AH=AH, S=S,
                lhs=M, rhs=2 * TH * S)


def smarr_rn(M, Q):
    """Reissner-Nordstrom (|Q| <= M).  r_+ = M + sqrt(M^2-Q^2)."""
    rp = M + np.sqrt(max(M * M - Q * Q, 0.0))
    TH = (rp - M) / (2 * np.pi * rp**2)
    AH = 4 * np.pi * rp**2
    S = AH / 4.0
    Phi = Q / rp
    return dict(M=M, Q=Q, J=0.0, Omega=0.0, Phi=Phi, TH=TH, AH=AH, S=S,
                lhs=M, rhs=2 * TH * S + Q * Phi)


def smarr_kerr(M, a):
    """Kerr with spin parameter a = J/M (|a| <= M); Omega_H = a/(r_+^2+a^2)."""
    rp = M + np.sqrt(max(M * M - a * a, 0.0))
    Delta = rp**2 - 2 * M * rp + a**2
    TH = (rp - M) / (2 * np.pi * (rp**2 + a**2))
    AH = 4 * np.pi * (rp**2 + a**2)
    S = AH / 4.0
    Omega = a / (rp**2 + a**2)
    J = a * M
    return dict(M=M, J=J, Q=0.0, Omega=Omega, Phi=0.0, TH=TH, AH=AH, S=S,
                lhs=M, rhs=2 * TH * S + 2 * Omega * J)


# --------------------------------------------------------------------------
# 5. data readers (only files that already exist in the repo)
# --------------------------------------------------------------------------


def read_physical_quantities():
    import csv

    p = os.path.join(ROOT, "solutions", "physical_quantities.csv")
    with open(p) as fh:
        rows = list(csv.DictReader(fh))
    return rows


def read_fig1_massofV():
    """alpha, chi, M, M_V, M_V/M in percent (data/fig1_massofV.csv)."""
    import csv

    p = os.path.join(ROOT, "data", "fig1_massofV.csv")
    with open(p) as fh:
        rows = list(csv.DictReader(fh))
    out = {}
    for r in rows:
        a = int(round(float(r["alpha"])))
        out.setdefault(a, []).append(
            (float(r["chi"]), float(r["MV"]), float(r["M"]))
        )
    for a in out:
        out[a].sort()
    return out


# --------------------------------------------------------------------------
# 6. tiny test harness -- works with or without pytest
# --------------------------------------------------------------------------


def run_tests(namespace, prefix="test_"):
    """Run every ``test_*`` function in ``namespace``; return list of result dicts."""
    results = []
    for name in sorted(k for k in namespace if k.startswith(prefix)):
        fn = namespace[name]
        t0 = time.time()
        try:
            fn()
            results.append(dict(name=name, status="PASS", seconds=round(time.time() - t0, 3), error=""))
        except AssertionError as exc:
            results.append(dict(name=name, status="FAIL", seconds=round(time.time() - t0, 3), error=str(exc)))
        except Exception as exc:  # noqa: BLE001
            results.append(dict(name=name, status="ERROR", seconds=round(time.time() - t0, 3),
                                error="%s: %s" % (type(exc).__name__, exc)))
    return results


def report(module_name, results):
    """Print a human table plus one machine-readable RESULT line per test."""
    print("\n=== %s ===" % module_name)
    for r in results:
        print("  [%-5s] %-34s %6.2fs  %s" % (r["status"], r["name"], r["seconds"], r["error"][:110]))
    for r in results:
        print("RESULT " + json.dumps(dict(name=r["name"], status=r["status"],
                                          seconds=r["seconds"], error=r["error"][:200])))
    n_fail = sum(1 for r in results if r["status"] != "PASS")
    print("  %d/%d passed" % (len(results) - n_fail, len(results)))
    return 0 if n_fail == 0 else 1


# --------------------------------------------------------------------------
# 7. independent numerical Ricci (nested finite differences)
# --------------------------------------------------------------------------


def _christoffel_numeric(metric_fn, x, h):
    n = 4
    g0 = metric_fn(x)
    gi = np.linalg.inv(g0)
    dg = np.zeros((2, n, n))
    for j in range(2):
        e = np.zeros(2); e[j] = h
        dg[j] = (metric_fn(x + e) - metric_fn(x - e)) / (2.0 * h)
    DIR = {1: 0, 2: 1}  # coordinate index -> finite-difference slot
    Gam = np.zeros((n, n, n))
    for a in range(n):
        for b in range(n):
            for c in range(n):
                acc = 0.0
                for d in range(n):
                    c1 = dg[DIR[c]][d, b] if c in DIR else 0.0   # d_c g_{db}
                    c2 = dg[DIR[b]][d, c] if b in DIR else 0.0   # d_b g_{dc}
                    c3 = dg[DIR[d]][b, c] if d in DIR else 0.0   # d_d g_{bc}
                    acc += gi[a, d] * (c1 + c2 - c3)
                Gam[a, b, c] = 0.5 * acc
    return Gam


def ricci_numeric(metric_fn, x, h=1e-4):
    """Ricci tensor from nested finite differences of ``metric_fn(x)``.

    Independent of the sympy path: shares no algebra with it.  Coordinates are
    ordered (t, r, theta, phi); the metric is assumed to depend only on (r, theta),
    so the only non-zero derivative directions are r (index 1) and theta (index 2).
    """
    n = 4
    DIR = {1: 0, 2: 1}  # coordinate index -> finite-difference slot
    x = np.asarray(x, float)
    Gam0 = _christoffel_numeric(metric_fn, x, h)
    dG = np.zeros((2, n, n, n))
    for ci, slot in DIR.items():
        e = np.zeros(2); e[slot] = h
        Gp = _christoffel_numeric(metric_fn, x + e, h)
        Gm = _christoffel_numeric(metric_fn, x - e, h)
        dG[slot] = (Gp - Gm) / (2.0 * h)
    Ric = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            s = 0.0
            # + d_c Gamma^c_{ab}
            for c in range(n):
                if c in DIR:
                    s += dG[DIR[c]][c, a, b]
            # - d_b Gamma^c_{ac}
            if b in DIR:
                for c in range(n):
                    s -= dG[DIR[b]][c, a, c]
            # + Gamma^c_{cd} Gamma^d_{ab} - Gamma^c_{ad} Gamma^d_{bc}
            for c in range(n):
                for d in range(n):
                    s += Gam0[c, c, d] * Gam0[d, a, b] - Gam0[c, a, d] * Gam0[d, b, c]
            Ric[a, b] = s
    return Ric
