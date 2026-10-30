#!/usr/bin/env python3
"""Flat limit: every field off and no mass => Minkowski, R_ab = 0.

Tests
-----
1. symbolic: R_ab(Minkowski) = 0 identically (cached sympy result);
2. the paper's ansatz with f=g=h=1, W=0, r_H=0 reproduces the Minkowski line
   element at random points;
3. independent numeric Ricci (nested finite differences of the metric) vanishes.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import _common as C


def _flat_symbolic():
    sp = C._symbols()["sp"]
    r, th = sp.symbols("r theta", positive=True)
    return C.ricci_symbolic(sp.diag(-1, 1, r**2, r**2 * sp.sin(th) ** 2))


def test_flat_symbolic_ricci_zero():
    """R_ab = 0 identically for the flat metric (symbolic, cached)."""
    res = C.cached_ricci("minkowski", _flat_symbolic)
    assert res["n_nonzero"] == 0, "flat metric has non-zero Ricci: %s" % (res["nonzero"][:3],)


def test_ansatz_reduces_to_minkowski():
    """f=g=h=1, W=0, r_H=0 in the paper's ansatz => diag(-1,1,r^2,r^2 sin^2)."""
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(25):
        r = rng.uniform(1.5, 12.0)
        th = rng.uniform(0.05, np.pi - 0.05)
        g = C.emv_metric(f=1.0, g=1.0, h=1.0, W=0.0, r=r, theta=th, rh=0.0)
        ref = np.diag([-1.0, 1.0, r**2, r**2 * np.sin(th) ** 2])
        worst = max(worst, float(np.abs(g - ref).max()))
    assert worst < 1e-14, "ansatz deviates from Minkowski by %.2e" % worst


def test_flat_numeric_ricci_zero():
    """Independent numeric Ricci: |R_ab| < 1e-6 for flat space (finite differences)."""
    def flat(x):
        r, th = x
        return np.diag([-1.0, 1.0, r**2, r**2 * np.sin(th) ** 2])

    worst = max(float(np.abs(C.ricci_numeric(flat, [r, th])).max())
                for r in (1.7, 3.3, 9.1) for th in (0.4, 1.2, 2.6))
    assert worst < 1e-6, "numeric Ricci of flat space = %.2e (expected ~0)" % worst


def test_flat_curvature_scalar_zero():
    """R = g^{ab}R_ab = 0 exactly (symbolic)."""
    sp = C._symbols()["sp"]
    r, th = sp.symbols("r theta", positive=True)
    g = sp.diag(-1, 1, r**2, r**2 * sp.sin(th) ** 2)
    Ric = C.ricci_symbolic(g)
    assert sp.simplify(sp.trace(Ric * g.inv())) == 0


def main():
    return C.report("test_flat_limit", C.run_tests(dict(globals())))


if __name__ == "__main__":
    sys.exit(main())
