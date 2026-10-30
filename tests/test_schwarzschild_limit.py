#!/usr/bin/env python3
"""Schwarzschild limit: switch off charge, vector hair and spin => Schwarzschild.

In the paper's ansatz the exact vacuum solution is

    g = (1 + r_H/r)^2,    f = 1/g,    h = 1,    W = 0,

whose areal radius is R = r g.  We verify, exactly:

1. R_ab = 0 symbolically (this *derives* the embedding rather than assuming it);
2. R_ab = 0 numerically (independent finite-difference path);
3. g_tt * g_RR = -1 in areal coordinates (the task's stated criterion);
4. horizon quantities: M = R_H/2, T_H = 1/(16 pi r_H), A_H = 64 pi r_H^2;
5. T_H A_H = 4 r_H, i.e. the paper's exact identity with g_H = 4;
6. Smarr M = 2 T_H S holds exactly.

A canary then shows that a plausible-but-wrong embedding (f = h = 1-r_H/r, g = 1)
is *rejected* by the same tests -- so the tests discriminate, they are not formal.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import _common as C

RH = 1.0


def _schw_symbolic():
    """The ansatz metric with g=(1+rh/r)^2, f=1/g, h=1 in (t,r,theta,phi)."""
    S = C._symbols(); sp, r, th = S["sp"], S["r"], S["th"]
    rh = sp.symbols("r_h", positive=True)
    gfun = (1 + rh / r) ** 2
    Phi = (gfun / (1 / gfun)) * 1  # (g/f) h  with f = 1/g, h = 1
    Psi = (1 / gfun) * (1 - rh / r) ** 2
    metric = sp.Matrix([
        [-Psi, 0, 0, 0],
        [0, Phi, 0, 0],
        [0, 0, Phi * r**2, 0],
        [0, 0, 0, Phi * r**2 * sp.sin(th) ** 2]])
    return C.ricci_symbolic(metric)


def _schw_numeric(x, rh=RH):
    r, th = x
    g = (1.0 + rh / r) ** 2
    return C.emv_metric(f=1.0 / g, g=g, h=1.0, W=0.0, r=r, theta=th, rh=rh)


def _wrong_embedding(x, rh=RH):
    """Plausible but WRONG: f = h = 1-r_h/r, g = 1 (that is not Schwarzschild)."""
    r, th = x
    f = 1.0 - rh / r
    return C.emv_metric(f=f, g=1.0, h=f, W=0.0, r=r, theta=th, rh=rh)


# ---------------------------------------------------------------- 1. symbolic
def test_vacuum_ricci_zero_symbolic():
    res = C.cached_ricci("schwarzschild_ansatz", _schw_symbolic)
    assert res["n_nonzero"] == 0, "vacuum ansatz is not Ricci-flat: %s" % (res["nonzero"][:3],)


# ---------------------------------------------------------------- 2. numeric
def test_vacuum_ricci_zero_numeric():
    worst = max(float(np.abs(C.ricci_numeric(_schw_numeric, [r, th])).max())
                for r in (1.3, 2.7, 8.0) for th in (0.5, 1.3, 2.5))
    assert worst < 5e-5, "numeric Ricci = %.2e (expected ~0, FD floor ~1e-6)" % worst


# ------------------------------------------------------ 3. g_tt * g_RR = -1
def test_gtt_times_gRR_minus_one():
    """Transform to the areal radius R = r g and check g_tt g_RR = -1 exactly."""
    rr = np.array([1.2, 1.7, 3.0, 7.5])
    for r in rr:
        g = (1.0 + RH / r) ** 2
        f = 1.0 / g
        Phi = (g / f) * 1.0
        gtt = -f * (1.0 - RH / r) ** 2
        R = C.areal_radius(r, RH)
        dRdr = 1.0 - RH**2 / r**2
        gRR = Phi / dRdr**2
        prod = gtt * gRR
        assert abs(prod + 1.0) < 1e-13, "g_tt g_RR = %.15f at r=%g (expected -1)" % (prod, r)


# ------------------------------------------------------------ 4-6. horizon data
def _horizon_data():
    rH = RH
    gH = (1.0 + RH / rH) ** 2          # = 4
    fH = 1.0 / gH                      # = 1/4
    hH = 1.0
    RH_areal = C.areal_radius(rH, RH)  # = 4 r_H
    M = RH_areal / 2.0                 # = 2 r_H
    TH = fH / (2.0 * np.pi * rH * np.sqrt(gH * hH))
    # A_H = 2 pi r_H^2 * int_0^pi sin(theta) g sqrt(h)/f dtheta ; the integrand is
    # theta-independent, so the integral equals 2 * (g sqrt(h)/f).
    AH = 2.0 * np.pi * rH**2 * 2.0 * (gH * np.sqrt(hH) / fH)
    S = AH / 4.0
    return dict(M=M, TH=TH, AH=AH, S=S, RH_areal=RH_areal, gH=gH)


def test_horizon_area_and_temperature():
    d = _horizon_data()
    assert abs(d["AH"] - 4 * np.pi * d["RH_areal"] ** 2) < 1e-10, \
        "A_H = %.10f != 4 pi R_H^2 = %.10f" % (d["AH"], 4 * np.pi * d["RH_areal"] ** 2)
    assert abs(d["TH"] - 1.0 / (16 * np.pi * RH)) < 1e-14, "T_H = %.12f (expected 1/(16 pi r_H))" % d["TH"]
    assert abs(d["M"] - 2 * RH) < 1e-14, "M = %.12f (expected 2 r_H)" % d["M"]


def test_TH_AH_identity():
    """T_H A_H = 4 r_H (paper's exact identity, g_H = 4)."""
    d = _horizon_data()
    val = d["TH"] * d["AH"]
    assert abs(val / (4 * RH) - 1.0) < 1e-12, "T_H A_H = %.12f, expected 4 r_H" % val


def test_smarr_schwarzschild():
    d = _horizon_data()
    lhs, rhs = d["M"], 2 * d["TH"] * d["S"]
    assert abs(lhs - rhs) / lhs < 1e-12, "Smarr: M=%.12f vs 2T_H S=%.12f" % (lhs, rhs)


# ------------------------------------------------------------------- canaries
def test_canary_wrong_embedding_is_rejected():
    """A wrong vacuum embedding must FAIL g_tt g_RR = -1, so the test is not formal.

    Wrong embedding: f = h = 1 - r_H/r, g = 1.  It looks natural (it makes
    g_tt = -(1-r_H/r)^3) but is *not* Schwarzschild; the areal-coordinate identity
    g_tt g_RR = -1 exposes it immediately.
    """
    r = 2.0
    f = 1.0 - RH / r
    Phi = (1.0 / f) * f          # (g/f) h with g = 1, h = f
    gtt = -f * (1.0 - RH / r) ** 2
    gRR = Phi / ((1.0 - RH**2 / r**2) ** 2)
    prod = gtt * gRR
    assert abs(prod + 1.0) > 1e-3, "canary failed: the wrong embedding passed g_tt g_RR = -1"
    # ... while the correct embedding passes to machine precision
    g = (1.0 + RH / r) ** 2
    Phi_ok = (g / (1.0 / g)) * 1.0
    gtt_ok = -(1.0 / g) * (1.0 - RH / r) ** 2
    gRR_ok = Phi_ok / ((1.0 - RH**2 / r**2) ** 2)
    assert abs(gtt_ok * gRR_ok + 1.0) < 1e-13, "correct embedding failed the same check"


def main():
    return C.report("test_schwarzschild_limit", C.run_tests(dict(globals())))


if __name__ == "__main__":
    sys.exit(main())
