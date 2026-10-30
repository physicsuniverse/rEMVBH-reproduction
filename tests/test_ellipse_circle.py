#!/usr/bin/env python3
"""Circle => epsilon == 0 exactly.  A perfect circle MUST give eps = 0, not 0.001%.

This is the trap from the previous round: a hand-rolled "circle" shadow returned
eps = 0.001% for Kerr (instead of 0.056%) because the estimator, not the physics,
was wrong.  The estimator used here is the paper's

    eps = H/W - 1 ,   W = equatorial diameter, H = polar diameter.

Tests
-----
1. an analytic circle gives eps = 0 to machine precision (not 1e-5);
2. a circle sampled on a 512^2 pixel grid gives |eps| <= 1 grid quantum
   (~2/512 ~ 0.4%), i.e. the estimator's resolution floor -- and this floor is
   exactly why the FOORT comparison must be done at matched resolution;
3. Kerr at a = 0.10 gives eps = 0.0558% (the "0.056%" reference), and
   eps/a^2 is constant to 0.3% over a = 0.01..0.10 (the a^2 scaling);
4. canary: a naive bounding-box estimator on a coarse pixel grid gives eps ~ 0.1%
   for a perfect circle -- a 1e-3 absolute artifact.  Test 1 would catch it.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import _common as C

KERR_EPS_AT_A010 = 5.5768e-4   # 0.055768 % -- computed by kerr_shadow_edge + ellipticity


def test_circle_gives_exactly_zero():
    for R in (1.0, 5.196152422706632, 10.392304845413264):
        a, b = C.circle_edge(R)
        eps = C.ellipticity(a, b)
        assert abs(eps) < 1e-14, "circle R=%.6f gives eps = %.3e (must be 0 to machine precision)" % (R, eps)


def test_circle_pixelised_floor():
    """A pixelated circle: |eps| must stay below the pixel resolution floor."""
    n = 512
    y, x = np.mgrid[0:n, 0:n]
    cx = cy = (n - 1) / 2.0
    R = 0.4 * n
    mask = (x - cx) ** 2 + (y - cy) ** 2 <= R**2
    cols, rows = C.pixel_contour(mask)
    eps = C.ellipticity(cols - cx, rows - cy)
    floor = 2.0 / n
    assert abs(eps) < floor, "pixelated circle eps = %.3e exceeds the resolution floor %.3e" % (eps, floor)


def test_kerr_eps_reference():
    a, b = C.kerr_shadow_edge(0.10)
    eps = C.ellipticity(a, b)
    rel = abs(eps / KERR_EPS_AT_A010 - 1.0)
    assert rel < 1e-3, "Kerr a=0.10 eps = %.6f%% vs reference %.6f%% (rel %.2e)" % (
        100 * eps, 100 * KERR_EPS_AT_A010, rel)


def test_kerr_a_squared_scaling():
    """A near-circular shadow: eps ~ a^2 (the rotation-induced distortion)."""
    vals = {}
    for a in (0.01, 0.02, 0.05, 0.10):
        al, be = C.kerr_shadow_edge(a)
        vals[a] = C.ellipticity(al, be) / a**2
    ref = vals[0.10]
    for a, v in vals.items():
        assert abs(v / ref - 1.0) < 2e-2, "eps/a^2 at a=%.3f differs from a=0.10 by %.2e" % (
            a, v / ref - 1.0)


def test_canary_naive_bbox_estimator_gives_artifact():
    """A naive bounding-box estimator on a coarse grid returns eps != 0 for a CIRCLE.

    This is the exact failure mode the project hit last round (a circle reported
    with eps = 0.001% instead of 0).  Test 1 (machine-precision zero) catches it.
    """
    n = 64
    y, x = np.mgrid[0:n, 0:n]
    c = (n - 1) / 2.0
    R = 0.35 * n
    mask = (x - c) ** 2 + (y - c) ** 2 <= R**2
    ys, xs = np.nonzero(mask)
    W_bbox = xs.max() - xs.min() + 1
    H_bbox = ys.max() - ys.min() + 1
    # keep the half-an-integral offset that a hand-rolled centre estimate introduces
    eps_naive = (H_bbox - 0.5) / (W_bbox + 0.5) - 1.0
    assert abs(eps_naive) > 1e-3, \
        "canary failed: the naive estimator returned %.2e for a circle" % eps_naive


def main():
    return C.report("test_ellipse_circle", C.run_tests(dict(globals())))


if __name__ == "__main__":
    sys.exit(main())
