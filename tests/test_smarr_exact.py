#!/usr/bin/env python3
"""The Smarr relation itself, tested on exact analytic metrics (machine precision).

For a stationary, axisymmetric, asymptotically flat black hole,
    M = 2 T_H S + 2 Omega_H J + Q Phi
must hold exactly (it is an identity once the horizon quantities are defined
consistently).  We check it on Schwarzschild, Reissner-Nordstrom and Kerr, where
every ingredient is known in closed form, so the test isolates the *formula* from
any numerical solution.

A canary shows that dropping the 2 Omega_H J term is detected for Kerr -- the test
has teeth rather than being a tautology.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import _common as C

TOL = 1e-13


def test_schwarzschild():
    for M in (0.5, 1.0, 3.7):
        s = C.smarr_schwarzschild(M)
        assert abs(s["lhs"] - s["rhs"]) / abs(s["lhs"]) < TOL, "Schwarzschild M=%g: %.3e" % (M, s["lhs"] - s["rhs"])


def test_reissner_nordstrom():
    for M, Q in ((1.0, 0.0), (1.0, 0.3), (1.0, 0.8), (2.0, 1.9)):
        s = C.smarr_rn(M, Q)
        assert abs(s["lhs"] - s["rhs"]) / abs(s["lhs"]) < TOL, "RN M=%g Q=%g: %.3e" % (M, Q, s["lhs"] - s["rhs"])


def test_kerr():
    for M, a in ((1.0, 0.0), (1.0, 0.5), (1.0, 0.9), (2.0, 1.7)):
        s = C.smarr_kerr(M, a)
        assert abs(s["lhs"] - s["rhs"]) / abs(s["lhs"]) < TOL, "Kerr M=%g a=%g: %.3e" % (M, a, s["lhs"] - s["rhs"])


def test_extremal_kerr_and_rn():
    """Extremal limits: T_H -> 0, and the Smarr relation still holds."""
    s = C.smarr_kerr(1.0, 1.0)
    assert abs(s["TH"]) < 1e-14, "extremal Kerr T_H = %.2e" % s["TH"]
    assert abs(s["lhs"] - s["rhs"]) / s["lhs"] < TOL
    s = C.smarr_rn(1.0, 1.0)
    assert abs(s["TH"]) < 1e-14, "extremal RN T_H = %.2e" % s["TH"]
    assert abs(s["lhs"] - s["rhs"]) / s["lhs"] < TOL


def test_identity_is_non_trivial():
    """The three contributions are individually nonzero and none dominates trivially."""
    s = C.smarr_kerr(1.0, 0.9)
    th_term = 2 * s["TH"] * s["S"]
    rot_term = 2 * s["Omega"] * s["J"]
    assert th_term > 0.1 * s["M"] and rot_term > 0.1 * s["M"], \
        "Kerr a=0.9: thermal=%.3f rotational=%.3f (both must be sizeable)" % (th_term, rot_term)


def test_canary_missing_rotation_term_rejected():
    """Dropping 2 Omega_H J must FAIL for Kerr -- otherwise the test is vacuous."""
    s = C.smarr_kerr(1.0, 0.9)
    wrong = 2 * s["TH"] * s["S"] + s["Q"] * s["Phi"]
    assert abs(s["lhs"] - wrong) / s["lhs"] > 0.1, \
        "canary failed: the no-rotation formula still matches Kerr (err=%.3e)" % (s["lhs"] - wrong)


def main():
    return C.report("test_smarr_exact", C.run_tests(dict(globals())))


if __name__ == "__main__":
    sys.exit(main())
