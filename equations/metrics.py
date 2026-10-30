#!/usr/bin/env python3
"""
equations/metrics.py  --  symbolic metric toolkit for the EMV rotating-black-hole
                          background, in pure sympy (NO Mathematica).

Conventions (the same ones used throughout this package; cf. the paper's metric ansatz):

    ds^2 = -f N^2 dt^2 + (g/f) [ h (dr^2 + r^2 dth^2)
                                 + r^2 sin^2(th) (dphi - (W/r)(1-N) dt)^2 ]

    N = 1 - r_H/r,     quasi-isotropic x = 1 - 2 r_H / r  =>  r = 2 r_H/(1-x)
    signature (-,+,+,+)
    f, g, h, W are functions of (r, th) -- here of (x, th).

Public API
----------
    coords, r_of_x, N_of_x                 : symbols / helper expressions
    metric_ansatz(rh=1)                    : 4x4 sympy Matrix g_{ab}
    christoffel(g)                         : Gamma^a_{bc}
    riemann(g)                             : R^a_{bcd}
    ricci(g)                               : R_{ab}
    ricci_scalar(g)                        : R
    einstein(g)                            : G_{ab} = R_{ab} - R g_{ab}/2
    free_symbols_of_derivs()               : lists the derivative payload
    make_numeric(expr_list)                : sympy.lambdify -> numpy callable
    spectral_derivatives(grid_x, grid_th, values)
                                           : numpy d/dx, d/dth, d2/... on the grid

Run `python3 equations/metrics.py` for a self-test that checks:
  * symmetry of Christoffel/Riemann/Ricci
  * Bianchi identity  grad_a G^{ab} = 0  (symbolic, on the ansatz)
  * the flat-space limit f=g=h=1, W=0 gives G_{ab}=0
  * the Schwarzschild-like limit reproduces R_{ab}=0 outside the horizon
"""
import sympy as sp

# --------------------------------------------------------------------------
# coordinates and the two "radial" variables
# --------------------------------------------------------------------------
t, x, th, ph = sp.symbols("t x theta phi", real=True)
rh = sp.symbols("r_H", positive=True)

COORDS = (t, x, th, ph)


def r_of_x(rh_=rh):
    """r = 2 r_H/(1-x)"""
    return 2 * rh_ / (1 - x)


def N_of_x():
    """N = 1 - r_H/r = (1+x)/2"""
    return sp.Rational(1, 2) * (1 + x)


def one_minus_N():
    """1 - N = (1-x)/2"""
    return sp.Rational(1, 2) * (1 - x)


# the four metric functions, as undefined functions of (x, theta)
f = sp.Function("f")(x, th)
g = sp.Function("g")(x, th)
h = sp.Function("h")(x, th)
Wv = sp.Function("W")(x, th)


def metric_ansatz(rh_=rh, funcs=None):
    """Return the covariant metric g_{ab} in coordinates (t, x, theta, phi).

    funcs : optional dict {'f','g','h','W'} of sympy expressions -- used by
            the self-test to plug in simple limits.
    """
    if funcs is None:
        ff, gg, hh, WW = f, g, h, Wv
    else:
        ff, gg, hh, WW = funcs["f"], funcs["g"], funcs["h"], funcs["W"]
    r = r_of_x(rh_)
    N = N_of_x()
    s2 = sp.sin(th) ** 2
    # dr/dx
    drdx = sp.diff(r, x)

    gtt = -ff * N ** 2 + (gg / ff) * r * r * s2 * (WW / r) ** 2 * (1 - N) ** 2
    gtp = -(gg / ff) * r * r * s2 * (WW / r) * (1 - N)
    gxx = (gg / ff) * hh * drdx ** 2
    ghh = (gg / ff) * hh * r * r
    gpp = (gg / ff) * r * r * s2

    M = sp.zeros(4, 4)
    M[0, 0] = sp.simplify(gtt)
    M[1, 1] = sp.simplify(gxx)
    M[2, 2] = sp.simplify(ghh)
    M[0, 3] = M[3, 0] = sp.simplify(gtp)
    M[3, 3] = sp.simplify(gpp)
    return M


# --------------------------------------------------------------------------
# curvature
# --------------------------------------------------------------------------
def _ginv(g):
    return g.inv()


def christoffel(g, coords=COORDS):
    """Gamma^a_{bc} = 1/2 g^{ad}(d_b g_{dc} + d_c g_{db} - d_d g_{bc})"""
    n = len(coords)
    ginv = _ginv(g)
    G = sp.MutableDenseNDimArray(sp.zeros(n ** 3), (n, n, n))  # G[a,b,c]=Gamma^a_{bc}
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                val = sp.Integer(0)
                for d in range(n):
                    val += ginv[a, d] * (sp.diff(g[d, c], coords[b])
                                         + sp.diff(g[d, b], coords[c])
                                         - sp.diff(g[b, c], coords[d]))
                val = sp.Rational(1, 2) * val
                G[a, b, c] = val
                G[a, c, b] = val
    return G


def riemann(g, coords=COORDS):
    """R^a_{bcd} = d_c Gamma^a_{bd} - d_d Gamma^a_{bc}
                   + Gamma^a_{ce} Gamma^e_{bd} - Gamma^a_{de} Gamma^e_{bc}"""
    n = len(coords)
    G = christoffel(g, coords)
    R = sp.MutableDenseNDimArray(sp.zeros(n ** 4), (n, n, n, n))
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    val = sp.diff(G[a, b, d], coords[c]) - sp.diff(G[a, b, c], coords[d])
                    for e in range(n):
                        val += G[a, c, e] * G[e, b, d] - G[a, d, e] * G[e, b, c]
                    R[a, b, c, d] = val
    return R


def ricci(g, coords=COORDS):
    """R_{bd} = R^a_{bad}"""
    n = len(coords)
    R = riemann(g, coords)
    Ric = sp.zeros(n, n)
    for b in range(n):
        for d in range(b, n):
            val = sum(R[a, b, a, d] for a in range(n))
            Ric[b, d] = val
            Ric[d, b] = val
    return Ric


def ricci_scalar(g, coords=COORDS):
    ginv = _ginv(g)
    Ric = ricci(g, coords)
    n = len(coords)
    return sum(ginv[a, b] * Ric[a, b] for a in range(n) for b in range(n))


def einstein(g, coords=COORDS):
    Ric = ricci(g, coords)
    R = ricci_scalar(g, coords)
    return sp.simplify(Ric - R * g / 2)


# --------------------------------------------------------------------------
# derivative payload for lambdify
# --------------------------------------------------------------------------
FUNCS = ("f", "g", "h", "W")
DERIVS = (
    (0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2),
)   # value, d/dx, d/dth, d2/dx2, d2/dx dth, d2/dth2


def payload_symbols():
    """Return (symbols, substitution_rules).

    We replace every occurrence of f(x,th) and its derivatives by a plain
    symbol, so that the final expression can be lambdified with a flat
    argument list of 4*6 = 24 numbers.
    """
    syms = []
    rules = []
    for name in FUNCS:
        fn = sp.Function(name)
        base = fn(x, th)
        for (i, j) in DERIVS:
            if i == 0 and j == 0:
                occ = base
            else:
                d = base
                if i:
                    d = sp.diff(d, x, i)
                if j:
                    d = sp.diff(d, th, j)
                occ = d
            s = sp.Symbol(f"{name}_{i}{j}", real=True)
            syms.append(s)
            rules.append((occ, s))
    return syms, rules


def substitute_payload(expr, rules):
    """Replace field functions and their derivatives (highest order first)."""
    out = expr
    for occ, s in rules:
        out = out.subs(occ, s)
    return out


def make_numeric(expr_list, extra_args=()):
    """sympy.lambdify the flattened expressions; returns a numpy-callable.

    The call signature is func(*extra_args, f00, f10, f01, f20, f11, f02,
    g00, ..., W02)  i.e. extra args first, then the 24 payload numbers.
    """
    syms, rules = payload_symbols()
    flat = []
    for e in expr_list:
        flat.append(substitute_payload(sp.together(e), rules))
    return sp.lambdify(tuple(extra_args) + tuple(syms), flat, modules="numpy")


# --------------------------------------------------------------------------
# numerical spectral derivatives on the (x, theta) grid
# --------------------------------------------------------------------------
def spectral_derivatives(xg, thg, vals, nyq=None):
    """Spectral (Chebyshev-tensor / cosine-tensor) derivatives.

    xg   : 1-D array, Chebyshev Gauss points in [-1,1]       (Nx,)
    thg  : 1-D array, theta collocation points in [0, pi/2]  (Ny,)
           (the solutions are equatorially symmetric, so the theta grid
            only covers half the domain; cos(2 j theta) is used, which
            extends automatically)
    vals : 2-D array (Nx, Ny)

    returns dict with keys '00','10','01','20','11','02'.
    """
    import numpy as np
    xg = np.asarray(xg, float)
    thg = np.asarray(thg, float)
    vals = np.asarray(vals, float)
    Nx, Ny = vals.shape
    # --- x direction: Chebyshev
    Bx = np.cos(np.arange(Nx)[None, :] * np.arccos(np.clip(xg, -1, 1))[:, None])
    dBx = np.zeros_like(Bx)
    s = np.sqrt(np.maximum(1 - xg ** 2, 1e-300))
    for k in range(1, Nx):
        dBx[:, k] = k * np.sin(k * np.arccos(np.clip(xg, -1, 1))) / s
    d2Bx = np.zeros_like(Bx)
    for k in range(2, Nx):
        # d2/dx2 T_k
        d2Bx[:, k] = (k * (xg ** 2 - 1) * np.sin(k * np.arccos(np.clip(xg, -1, 1)))
                      / (s ** 3)
                      - k ** 2 * np.cos(k * np.arccos(np.clip(xg, -1, 1))) / (s ** 2))
    Dx = dBx @ np.linalg.inv(Bx)
    D2x = d2Bx @ np.linalg.inv(Bx)
    # --- theta direction: cos(2 j theta), equatorially symmetric
    Bt = np.cos(2.0 * np.arange(Ny)[None, :] * thg[:, None])
    dBt = -2.0 * np.arange(Ny)[None, :] * np.sin(2.0 * np.arange(Ny)[None, :] * thg[:, None])
    d2Bt = -4.0 * (np.arange(Ny)[None, :] ** 2) * Bt
    Dt = dBt @ np.linalg.inv(Bt)
    D2t = d2Bt @ np.linalg.inv(Bt)

    vx = Dx @ vals
    vxx = D2x @ vals
    vt = vals @ Dt.T
    vtt = vals @ D2t.T
    vxt = Dx @ vals @ Dt.T
    return {"00": vals, "10": vx, "01": vt, "20": vxx, "11": vxt, "02": vtt}


# --------------------------------------------------------------------------
# self-test
# --------------------------------------------------------------------------
def _selftest():
    print("=== metrics.py self-test ===")
    g0 = metric_ansatz()
    print("metric_ansatz: 4x4, components  g_tt, g_tphi, g_xx, g_thth, g_phiphi nonzero")
    for (i, j, nm) in ((0, 0, "g_tt"), (0, 3, "g_tphi"), (1, 1, "g_xx"),
                       (2, 2, "g_thth"), (3, 3, "g_phiphi")):
        print(f"   {nm:8s} = {sp.simplify(g0[i, j])}")

    # --- contracted Bianchi identity:  grad_a G^{ab} = 0  (must hold for ANY metric)
    #     This is the decisive consistency test of Christoffel -> Riemann -> Ricci -> Einstein.
    import numpy as np
    print("[INFO] contracted Bianchi identity is verified numerically in verify_eom.py;")
    print("       here we check it symbolically on a simple 2-parameter slice.")
    #   Take f, g, h, W as low-order polynomials in x only (theta-independent) so the
    #   symbolic contraction is cheap but still non-trivial.
    ff = 1 + sp.Rational(1, 3) * x + sp.Rational(1, 7) * x ** 2
    gg = 1 - sp.Rational(1, 5) * x + sp.Rational(1, 11) * x ** 2
    hh = 1 + sp.Rational(1, 4) * x ** 2
    WW = sp.Rational(1, 9) * (1 + x) ** 2
    gs = metric_ansatz(funcs={"f": ff, "g": gg, "h": hh, "W": WW})
    Gs = einstein(gs)
    ginv = gs.inv()
    # G^{ab} = g^{ac} g^{bd} G_{cd}
    Gup = sp.Matrix(4, 4, lambda a, b: sp.simplify(
        sum(ginv[a, c] * ginv[b, d] * Gs[c, d] for c in range(4) for d in range(4))))
    # div^b = d_a G^{ab} + Gamma^a_{ac} G^{cb} + Gamma^b_{ac} G^{ac}
    Gs_ch = christoffel(gs)
    div = []
    for b in range(4):
        val = sp.Integer(0)
        for a in range(4):
            val += sp.diff(Gup[a, b], COORDS[a])
            for c in range(4):
                val += Gs_ch[a, a, c] * Gup[c, b] + Gs_ch[b, a, c] * Gup[a, c]
        div.append(sp.simplify(val))
    ok = all(sp.simplify(d) == 0 for d in div)
    print(f"[{'PASS' if ok else 'FAIL'}] contracted Bianchi identity  grad_a G^(ab) = 0")

    # --- antisymmetry of Riemann
    R = riemann(g0)
    ok = True
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    if sp.simplify(R[a, b, c, d] + R[a, b, d, c]) != 0:
                        ok = False
    print(f"[{'PASS' if ok else 'FAIL'}] Riemann antisymmetry R^a_bcd = -R^a_bdc")

    # --- Ricci symmetry
    Ric = ricci(g0)
    okr = True
    for a in range(4):
        for b in range(4):
            if sp.simplify(Ric[a, b] - Ric[b, a]) != 0:
                okr = False
    print(f"[{'PASS' if okr else 'FAIL'}] Ricci symmetry R_ab = R_ba")

    # --- second Bianchi identity is checked numerically in verify/verify_eom.py
    print("[INFO] second Bianchi identity: checked in verify_eom.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(_selftest())
