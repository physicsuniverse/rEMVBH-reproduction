#!/usr/bin/env python3
"""Kerr limit: charge = 0 and vector hair = 0 => Kerr; and a -> 0 => Schwarzschild.

Tests
-----
1. Kerr (Boyer-Lindquist) is Ricci flat -- verified numerically (independent path);
2. the horizon is at Delta(r_+) = 0; a = M is extremal and T_H -> 0 there;
3. frame dragging: g_t phi != 0 for a != 0 and g_t phi -> 0 as a -> 0;
4. the a -> 0 limit reproduces Schwarzschild (g_tt(R) = -(1-2M/R) in areal coordinates);
5. canary: a Kerr metric with a sign error in g_t phi is REJECTED by the Ricci test.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import _common as C

M = 1.0


def test_kerr_ricci_flat():
    worst = 0.0
    for a in (0.0, 0.3, 0.7, 0.9):
        for r in (3.0, 6.0, 20.0):
            for th in (0.5, 1.2, 2.4):
                if a == 0.0:
                    met = lambda x, a=a: C.kerr_bl_metric(x[0], x[1], 1e-8)
                else:
                    met = lambda x, a=a: C.kerr_bl_metric(x[0], x[1], a)
                worst = max(worst, float(np.abs(C.ricci_numeric(met, [r, th], h=1e-3)).max()))
    assert worst < 1e-4, "Kerr Ricci = %.2e (expected 0, FD floor ~1e-6)" % worst


def test_horizon_and_extremality():
    for a in (0.0, 0.5, 0.9, M):
        rp = M + np.sqrt(max(M * M - a * a, 0.0))
        Delta = rp**2 - 2 * M * rp + a**2
        assert abs(Delta) < 1e-12, "Delta(r_+) = %.3e at a=%g" % (Delta, a)
    # extremal: T_H -> 0
    a = M
    rp = M + np.sqrt(max(M * M - a * a, 0.0))
    TH = (rp - M) / (2 * np.pi * (rp**2 + a**2))
    assert abs(TH) < 1e-12, "T_H at extremality = %.3e (expected 0)" % TH


def test_frame_dragging():
    g = C.kerr_bl_metric(4.0, np.pi / 2, 0.7)
    assert abs(g[0, 3]) > 1e-2, "g_tphi vanishes for a=0.7 -- frame dragging missing"
    g0 = C.kerr_bl_metric(4.0, np.pi / 2, 1e-10)
    assert abs(g0[0, 3]) < 1e-9, "g_tphi does not vanish as a->0: %.2e" % abs(g0[0, 3])


def test_a_to_zero_is_schwarzschild():
    """The a->0 Kerr limit, written in areal coordinates, must equal Schwarzschild."""
    a = 1e-6
    for r in (2.5, 5.0, 12.0):
        g = C.kerr_bl_metric(r, np.pi / 2, a)
        gtt = g[0, 0]
        sch = -(1.0 - 2 * M / r)
        assert abs(gtt - sch) < 1e-9, "g_tt(a->0) = %.12f vs Schwarzschild %.12f" % (gtt, sch)


def test_ansatz_vacuum_matches_schwarzschild():
    """The paper-ansatz vacuum solution has g_tt(R) = -(1-2M/R) at R = r(1+r_H/r)^2."""
    rh = 1.0
    for r in (1.4, 2.2, 6.0):
        g = C.schwarzschild_in_ansatz(r, np.pi / 2, rh)
        R = C.areal_radius(r, rh)
        M_areal = 2 * rh
        gtt = g[0, 0]
        ref = -(1.0 - 2 * M_areal / R)
        assert abs(gtt - ref) < 1e-12, "ansatz g_tt=%.12f vs Schwarzschild %.12f at R=%.4f" % (gtt, ref, R)


def test_canary_bad_kerr_sign_rejected():
    """A wrong normalisation of g_tphi (a common slip) must be caught by the Ricci test."""
    def bad(x):
        r, th = x
        g = C.kerr_bl_metric(r, th, 0.7)
        # NOTE: flipping the sign of g_tphi is just a -> -a, i.e. still Kerr.
        # A genuine error is a wrong *magnitude* (a wrong normalisation of J).
        g[0, 3] *= 1.25
        g[3, 0] = g[0, 3]
        return g
    val = float(np.abs(C.ricci_numeric(bad, [6.0, 1.2], h=1e-3)).max())
    assert val > 1e-4, "canary failed: a wrong-normalisation Kerr passed the Ricci test (%.2e)" % val


def main():
    return C.report("test_kerr_limit", C.run_tests(dict(globals())))


if __name__ == "__main__":
    sys.exit(main())
