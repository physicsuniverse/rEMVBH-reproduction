#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig4_shadows.py -- Python redraw of Fig. 4 "Shadows" (4 panels, 2x2)
================================================================================

solid   = rEMVBH shadow edge, from **our own FOORT ray tracing**
          (this package ships the extracted boundaries in `data/shadow_boundaries/`).
dashed  = analytic Kerr-Newman shadow edge (Bardeen asymptotic impact parameters).

WARNING: the two pipelines are not exactly equivalent (methodology note that **must** accompany the figure)
------------------------------------------------------------------------------
FOORT constructs the photon directions at the observer position r_obs = 1000 M in a **local orthonormal
frame**, and its output (alpha/M, beta/M) is the tangential projection on the screen; the Bardeen (alpha, beta)
are **asymptotic** impact parameters. The two differ by O(M/r_obs). Measured in this repository (native KN
through FOORT vs analytic KN; at 1024^2 |dW/W| = 0.04%, |dH/H| = 0.10%, max contour deviation 0.0118 M; at
512^2 the effect is buried by pixel quantisation (0.0293 M/px) as 0.2-0.35%). Data: `data/fig4_shadows_pipeline_check.csv`,
discussion in `REPORT_fig4.md` section 3 -- which explains that the 0.2-0.35% quoted for C1 last round was a resolution artefact.
The KN dashed curve is therefore labelled "asymptotic Bardeen impact parameters"; the O(0.02 M) offset between
the solid and dashed curves **must not** be read as a physical effect.

Panel parameters (verified against the legend of the original figure and `01_paper_analysis/per-figure analysis.md` section 4;
------------------------------------------------------------------------------
  (a) top left     chi = 0.10       scan q = 0.2347 / 0.4887 / 0.7407 / 0.9947 / 1.2527 (5 EMV curves)
  (b) top right    q = 0.6887       scan chi = 0.010 / 0.100 / 0.300 / 0.500 / 0.725 (5 EMV curves)
  (c) bottom left  chi = 0.10       EMV q = 0.9947, EMV q = 1.2527  vs  KN q = 0.9947
  (d) bottom right q = 0.6887, chi = 0.725   EMV  vs  KN

Usage
----
    export PATH=/usr/local/texlive/2026/bin/x86_64-linux:$PATH
    export PYTHONPATH=<workspace>/tools/pylibs
    python3 07_figures/fig4_shadows.py
"""

import os
import sys
import csv

# matplotlib writes its LaTeX cache to ~/.cache by default, which is not writable in this sandbox; the cache
# must be pointed at a writable directory before importing figstyle (which imports matplotlib). Override with MPLCONFIGDIR.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(_HERE)          # reproduction package root
sys.path.insert(0, _HERE)
from figstyle import (apply_style, panel_grid, panel_label, style_axes, save,
                      CURVE_COLORS, LW_CURVE, LW_CONTOUR, FS_LEG, FS_TITLE)

BND = os.path.join(ROOT, "data", "shadow_boundaries")
OUT = os.path.join(_HERE, "out")


# ════════════════════════════════════════════════════════════════════════════
# 1. analytic Kerr-Newman shadow (Bardeen asymptotic impact parameters), self-implemented
# ════════════════════════════════════════════════════════════════════════════
def kn_shadow(q, chi, theta_obs=np.pi / 2.0, n=6001):
    """Closed KN shadow contour (alpha, beta), in units of M = 1.

    In M = 1 units a = chi and Q = q; Delta = r^2 - 2r + a^2 + Q^2.
    The spherical photon orbits are given by R(r) = 0 and dR/dr = 0:
        R/E^2 = [r^2 + a^2 - a xi]^2 - Delta[(xi - a)^2 + eta]
        xi(r)  = (r^2 + a^2 - 4 r Delta / Delta') / a
        eta(r) = 16 r^2 Delta / Delta'^2 - (xi - a)^2
    Celestial coordinates: alpha = -xi / sin(theta_obs),
               beta  = +- sqrt(eta + a^2 cos^2(theta_obs) - xi^2 cot^2(theta_obs))
    (consistent with Tsukamoto:2017fxq / Shaikh:2019fpu; this implementation was validated independently
     to 0.02 M against FOORT's native KerrNewmanMetric in
     `plotting/kn_reference.py` in the full analysis tree, see REPORT_fig4.md.)
    """
    a, Q = float(chi), float(q)
    if a * a + Q * Q >= 1.0:
        raise ValueError(f"no Kerr-Newman black hole for q={q}, chi={chi}")

    def Delta(r):
        return r * r - 2.0 * r + a * a + Q * Q

    def Dp(r):
        return 2.0 * r - 2.0

    if a == 0.0:                                    # Schwarzschild / RN limit
        r = (3.0 + np.sqrt(9.0 - 8.0 * Q * Q)) / 2.0
        b2 = 4.0 * r ** 3 / (2.0 * r - 2.0)
        xi = np.linspace(-np.sqrt(b2), np.sqrt(b2), n)
        eta = b2 - xi ** 2
    else:
        def xi_of_r(r):
            return (r * r + a * a - 4.0 * r * Delta(r) / Dp(r)) / a

        def eta_of_r(r):
            return (16.0 * r * r * Delta(r) / Dp(r) ** 2) - (xi_of_r(r) - a) ** 2

        rs = np.linspace(1.0 + 1e-9, 14.0, 400001)
        with np.errstate(divide="ignore", invalid="ignore"):
            et = eta_of_r(rs)
        ok = np.isfinite(et) & (et >= 0.0)
        if not ok.any():
            raise ValueError(f"no spherical photon orbits for q={q}, chi={chi}")
        i0 = int(np.argmax(ok))
        i1 = len(ok) - 1 - int(np.argmax(ok[::-1]))

        def bisect(f, lo, hi, it=200):
            flo = f(lo)
            for _ in range(it):
                mid = 0.5 * (lo + hi)
                fm = f(mid)
                if (fm > 0) == (flo > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return 0.5 * (lo + hi)

        r1 = bisect(eta_of_r, rs[max(i0 - 1, 0)], rs[i0])
        r2 = bisect(eta_of_r, rs[i1], rs[min(i1 + 1, len(rs) - 1)])
        s = 0.5 * (1.0 - np.cos(np.pi * np.linspace(0.0, 1.0, n)))
        rr = r1 + (r2 - r1) * s
        xi = xi_of_r(rr)
        eta = eta_of_r(rr)

    st, ct = np.sin(theta_obs), np.cos(theta_obs)
    alpha = -xi / st
    disc = np.clip(eta + a * a * ct * ct - xi * xi * (ct / st) ** 2, 0.0, None)
    beta = np.sqrt(disc)
    return (np.concatenate([alpha, alpha[::-1], alpha[:1]]),
            np.concatenate([beta, -beta[::-1], beta[:1]]))


# ════════════════════════════════════════════════════════════════════════════
# 2. data: FOORT-generated edges
# ════════════════════════════════════════════════════════════════════════════
# run name -> (q, chi, M, metric file (relative to FOORT/), grid)
EMV_RUNS = {
    "q2347_1024":         (0.2347, 0.100, 2.457682, "<archive>/chi100difq/metricdatas/q2347chi100a20.dat", 1024),
    "q4887_1024":         (0.4887, 0.100, 2.620235, "<archive>/chi100difq/metricdatas/q4887chi100a20.dat", 1024),
    "q7407_1024":         (0.7407, 0.100, 2.952962, "<archive>/chi100difq/metricdatas/q7407chi100a20.dat", 1024),
    "q9947_1024":         (0.9947, 0.100, 3.852004, "<archive>/chi100difq/metricdatas/q9947chi100a20.dat", 1024),
    "q12527_2048":        (1.2527, 0.100, 5.780872, "<archive>/chi100difq/metricdatas/q12527chi100a20.dat", 2048),
    "q6887chi10_1024":    (0.6887, 0.010, 2.775376, "<archive>/q6887difchi/metricdatas/q6887chi10a20.dat", 1024),
    "q6887chi100_1024":   (0.6887, 0.100, 2.809558, "<archive>/q6887difchi/metricdatas/q6887chi100a20.dat", 1024),
    "q6887chi300_1024":   (0.6887, 0.300, 3.038972, "<archive>/q6887difchi/metricdatas/q6887chi300a20.dat", 1024),
    "q6887chi500_1024":   (0.6887, 0.500, 3.642166, "<archive>/q6887difchi/metricdatas/q6887chi500a20.dat", 1024),
    "q6887chi725_1024":   (0.6887, 0.725, 7.481893, "<archive>/G/metricdata/q6887chi725.dat", 1024),
}


def load_boundary(run, close=True):
    p = os.path.join(BND, run + "_boundary.csv")
    d = np.loadtxt(p, delimiter=",", skiprows=1)
    a, b = d[:, 0].copy(), d[:, 1].copy()
    if close:
        a = np.append(a, a[0])
        b = np.append(b, b[0])
    return a, b


def write_csv(name, rows):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name)
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["curve", "kind", "q", "chi", "alpha_over_M", "beta_over_M"])
        w.writerows(rows)
    print(f"  wrote {p}  ({len(rows)} rows)")
    return p


# ════════════════════════════════════════════════════════════════════════════
# 3. plotting
# ════════════════════════════════════════════════════════════════════════════
FS_NOTE = 6.0

# Legends must never sit on top of the curves: in (a)(b) the curves fill the panel -> move below the axes;
# in (c)(d) the shadow interior is largely empty -> place the legend in the centre of the panel.
LEG_KW = dict(frameon=False, borderpad=0.0, labelspacing=0.30,
              handlelength=1.9, handletextpad=0.5)
LEG_BELOW = dict(ncol=3, fontsize=6.6, columnspacing=1.0,
                 handlelength=1.6, frameon=False)


def _limits(ax, xlim, ylim):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$\alpha/M$")
    ax.set_ylabel(r"$\beta/M$")
    ax.axhline(0.0, color="0.82", lw=0.6, ls=(0, (5, 4)), zorder=0)
    ax.axvline(0.0, color="0.82", lw=0.6, ls=(0, (5, 4)), zorder=0)
    style_axes(ax)


def main():
    apply_style()
    # reserve a band at the bottom of the canvas for the footnote; width exactly = COL_DOUBLE (7.00 in)
    fig, axes = panel_grid(2, 2, width="double", height=7.30,
                           hspace=0.34, wspace=0.26,
                           left=0.082, right=0.985, top=0.965, bottom=0.190)
    (axa, axb), (axc, axd) = axes
    all_rows = {}

    # ---------------- (a) χ = 0.10, scan q ----------------
    qs = [0.2347, 0.4887, 0.7407, 0.9947, 1.2527]
    runs = ["q2347_1024", "q4887_1024", "q7407_1024", "q9947_1024", "q12527_2048"]
    rows = []
    for q, run, c in zip(qs, runs, CURVE_COLORS):
        a, b = load_boundary(run)
        axa.plot(a, b, color=c, lw=LW_CURVE, solid_capstyle="round",
                 label=rf"$q={q:.4f}$")
        rows += [("EMV", "EMV", q, 0.10, x, y) for x, y in zip(a, b)]
    axa.set_title(r"$\chi = 0.10$,  rEMVBH: $q = 0.2347 \ldots 1.2527$",
                  fontsize=FS_TITLE - 1.5, pad=4.0)
    axa.legend(loc="upper center", bbox_to_anchor=(0.5, -0.155), **LEG_BELOW)
    _limits(axa, (-5.35, 5.75), (-5.50, 5.50))
    panel_label(axa, "a")
    all_rows["a"] = rows

    # ---------------- (b) q = 0.6887, scan χ ----------------
    chis = [0.010, 0.100, 0.300, 0.500, 0.725]
    runs = ["q6887chi10_1024", "q6887chi100_1024", "q6887chi300_1024",
            "q6887chi500_1024", "q6887chi725_1024"]
    rows = []
    for chi, run, c in zip(chis, runs, CURVE_COLORS):
        a, b = load_boundary(run)
        axb.plot(a, b, color=c, lw=LW_CURVE, solid_capstyle="round",
                 label=rf"$\chi={chi:.3f}$")
        rows += [("EMV", "EMV", 0.6887, chi, x, y) for x, y in zip(a, b)]
    axb.set_title(r"$q = 0.6887$,  rEMVBH: $\chi = 0.010 \ldots 0.725$",
                  fontsize=FS_TITLE - 1.5, pad=4.0)
    axb.legend(loc="upper center", bbox_to_anchor=(0.5, -0.155), **LEG_BELOW)
    _limits(axb, (-5.35, 6.95), (-5.50, 5.50))
    panel_label(axb, "b")
    all_rows["b"] = rows

    # ---------------- (c) χ = 0.10, EMV vs KN ----------------
    rows = []
    a, b = load_boundary("q9947_1024")
    axc.plot(a, b, color=CURVE_COLORS[0], lw=LW_CURVE, solid_capstyle="round",
             label=r"rEMVBH: $q=0.9947$")
    rows += [("EMV_q0.9947", "EMV", 0.9947, 0.10, x, y) for x, y in zip(a, b)]
    a, b = load_boundary("q12527_2048")
    axc.plot(a, b, color=CURVE_COLORS[2], lw=LW_CURVE, solid_capstyle="round",
             label=r"rEMVBH: $q=1.2527$")
    rows += [("EMV_q1.2527", "EMV", 1.2527, 0.10, x, y) for x, y in zip(a, b)]
    ka, kb = kn_shadow(0.9947, 0.10)
    axc.plot(ka, kb, color="k", lw=LW_CONTOUR + 0.5, ls=(0, (5, 3.2)),
             label=r"KN: $q=0.9947$")
    rows += [("KN_q0.9947", "KN", 0.9947, 0.10, x, y) for x, y in zip(ka, kb)]
    axc.set_title(r"$\chi = 0.10$:  rEMVBH $q = 0.9947,\,1.2527$;  KN $q = 0.9947$",
                  fontsize=FS_TITLE - 1.5, pad=4.0)
    axc.legend(loc="center", fontsize=FS_LEG - 0.5, **LEG_KW)
    _limits(axc, (-4.85, 4.85), (-4.85, 4.85))
    panel_label(axc, "c")
    all_rows["c"] = rows

    # ---------------- (d) q = 0.6887, χ = 0.725, EMV vs KN ----------------
    rows = []
    a, b = load_boundary("q6887chi725_1024")
    axd.plot(a, b, color=CURVE_COLORS[0], lw=LW_CURVE, solid_capstyle="round",
             label="rEMVBH")
    rows += [("EMV", "EMV", 0.6887, 0.725, x, y) for x, y in zip(a, b)]
    ka, kb = kn_shadow(0.6887, 0.725)
    axd.plot(ka, kb, color="k", lw=LW_CONTOUR + 0.5, ls=(0, (5, 3.2)), label="KN")
    rows += [("KN", "KN", 0.6887, 0.725, x, y) for x, y in zip(ka, kb)]
    axd.set_title(r"$q = 0.6887$,  $\chi = 0.725$:  rEMVBH vs KN",
                  fontsize=FS_TITLE - 1.5, pad=4.0)
    axd.legend(loc="center", fontsize=FS_LEG - 0.5, **LEG_KW)
    _limits(axd, (-3.95, 6.85), (-5.35, 5.35))
    panel_label(axd, "d")
    all_rows["d"] = rows

    # ---------------- in-figure methodology footnote (two pipelines) ----------------
    # placed **inside** the canvas (y in [0,1]). Reason: for long usetex strings matplotlib's
    # get_window_extent returns a grossly overestimated width, and with savefig.bbox="tight" that
    # inflates the page from 7.00 in to 13.9 in. Putting the content on the canvas with bbox=None is the
    # only way to keep the page exactly equal to figsize (= COL_DOUBLE wide).
    fig.text(0.5, 0.104,
             r"Solid curves: rEMVBH shadow, FOORT ray tracing in the observer's local "
             r"orthonormal frame at $r_{\rm obs}=1000\,M$.",
             ha="center", va="bottom", fontsize=FS_NOTE + 0.5)
    fig.text(0.5, 0.062,
             r"Dashed curves: analytic Kerr--Newman shadow, asymptotic Bardeen "
             r"impact parameters $(\alpha,\beta)$.",
             ha="center", va="bottom", fontsize=FS_NOTE + 0.5)
    fig.text(0.5, 0.020,
             r"Pipeline offset (FOORT native KN vs analytic KN): "
             r"$|\Delta W/W|,\,|\Delta H/H|\lesssim0.1\%$, contour $\lesssim0.02\,M$.",
             ha="center", va="bottom", fontsize=FS_NOTE + 0.5)

    for k in ("a", "b", "c", "d"):
        write_csv(f"fig4_shadows_{k}.csv", all_rows[k])

    import matplotlib.pyplot as _plt
    _plt.rcParams["savefig.bbox"] = None      # see the note above
    save(fig, "fig4_shadows", fmt=("pdf", "png"))


if __name__ == "__main__":
    main()
