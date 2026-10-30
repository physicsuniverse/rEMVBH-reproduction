#!/usr/bin/env python3
"""
fig10_spectral.py -- Fig. 10: spectral evidence for the smoothness of the numerical solutions (2x2)
================================================================================
Redrawn with the unified style module 07_figures/figstyle.py, replacing the older version in 05_solution_plots.

Four panels
--------
(a) angular spectrum: max_i |c_ij| vs j (8 fields, geometric decay)
(b) radial spectrum: max_j |c_ij| vs i (Chebyshev coefficients, geometric decay)
(c) deviation from the analytic g = (1 + r_H/r)^2 (horizontal axis x = 1 - 2 r_H/r)
(d) derivative spectrum i^2 max_j |c_ij| (shows that second derivatives are fully resolved)

Data: solutions/field_profiles/sol_*.txt (8 fields on a 60x20 grid, x slowest)
Coordinates: x = 1 - 2 r_H/r in [-1, 1], theta in [0, pi/2]
"""
import os
import sys
import glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "tools", "pylibs"))
from figstyle import *                       # noqa: E402
import matplotlib.pyplot as plt              # noqa: E402

DATA = os.path.join(ROOT, "solutions", "field_profiles")
NAMES = ["Bt", "Bp", "At", "Ap", "f", "g", "h", "W"]
LABELS = [r"$B_t$", r"$B_\varphi$", r"$A_t$", r"$A_\varphi$",
          r"$f$", r"$g$", r"$h$", r"$W$"]


# ----------------------------------------------------------- data loading
def load_grids():
    with open(os.path.join(DATA, "grids.txt")) as fh:
        lines = [l for l in fh if not l.startswith("#") and l.strip()]
    x = np.array([float(v) for v in lines[0].split()])
    th = np.array([float(v) for v in lines[1].split()])
    return x, th


def load_solution(path):
    with open(path) as fh:
        lines = fh.readlines()
    # L0: "# q chi alpha rh Nx Ny"
    par = np.array([float(v) for v in lines[1].split()])
    # L2: "# TH AH M Q Phi P J OmegaH"
    phy = np.array([float(v) for v in lines[3].split()])
    # L4: the 8 lines after "# fields: ..." are the 8 fields
    F = np.array([[float(v) for v in lines[5 + k].split()] for k in range(8)])
    Nx, Ny = int(par[4]), int(par[5])
    d = {n: F[k].reshape(Nx, Ny) for k, n in enumerate(NAMES)}
    d.update(q=par[0], chi=par[1], alpha=par[2], rh=par[3], Nx=Nx, Ny=Ny,
             TH=phy[0], AH=phy[1], M=phy[2], P=phy[5], OmH=phy[7],
             path=os.path.basename(path))
    return d


def cheb_coeffs(f):
    """Chebyshev coefficients of f sampled at x_i = cos(pi(2i+1)/(2N)).
    (consistent with the validated implementation in 05_solution_plots/plot_solutions.py)"""
    N = len(f)
    k = np.arange(N)
    M = np.cos(np.pi * np.outer(k, 2 * np.arange(N) + 1) / (2 * N))
    c = (2.0 / N) * (M @ f)
    c[0] /= 2.0
    return c


def cos_coeffs(g):
    """cos(2 j theta) coefficients of g sampled at theta_j = pi(2j+1)/(4Ny)."""
    N = len(g)
    j = np.arange(N)
    M = np.cos(np.pi * np.outer(j, 2 * np.arange(N) + 1) / (2 * N))
    c = (2.0 / N) * (M @ g)
    c[0] /= 2.0
    return c


def spectrum2d(F):
    """2-D spectral coefficients c[i,j] of F(x,theta) on the solver's grid。
    Transform in theta for each radial point first, then transform radially for each angular mode."""
    Nx, Ny = F.shape
    A = np.array([cos_coeffs(F[i, :]) for i in range(Nx)])       # (Nx, Ny)
    C = np.array([cheb_coeffs(A[:, j]) for j in range(Ny)]).T    # (Nx, Ny)
    return C


# ----------------------------------------------------------- main figure
def main():
    apply_style()
    xg, thg = load_grids()
    files = sorted(glob.glob(os.path.join(DATA, "sol_*.txt")))
    sols = [load_solution(p) for p in files]
    ref = next(s for s in sols if abs(s["chi"] - 0.5) < 1e-6 and abs(s["q"] - 0.6937) < 1e-6)

    fig, axes = panel_grid(2, 2, width=COL_DOUBLE, height=6.05,
                           hspace=0.30, wspace=0.24)

    # -- (a) angular spectrum
    ax = axes[0, 0]
    for k, n in enumerate(NAMES):
        c = spectrum2d(ref[n])
        ang = np.abs(c).max(axis=0)
        ax.semilogy(np.arange(len(ang)), np.maximum(ang, 1e-18),
                    color=CURVE_COLORS[k % 6], lw=1.2 if k < 6 else 1.2,
                    ls="-" if k < 6 else "--", label=LABELS[k], alpha=0.9)
    ax.set_xlabel(r"angular mode $j$   [basis $\cos 2j\theta$]")
    ax.set_ylabel(r"$\max_i\,|c_{ij}|$")
    ax.set_ylim(1e-17, 3e0)
    ax.set_xlim(-0.6, 25)
    ax.legend(ncol=3, loc="upper right", fontsize=6.2, columnspacing=0.7,
              handlelength=1.2, labelspacing=0.22)
    panel_label(ax, "a")

    # -- (b) radial spectrum
    ax = axes[0, 1]
    for k, n in enumerate(NAMES):
        c = spectrum2d(ref[n])
        rad = np.abs(c).max(axis=1)
        ax.semilogy(np.arange(len(rad)), np.maximum(rad, 1e-18),
                    color=CURVE_COLORS[k % 6], lw=1.2,
                    ls="-" if k < 6 else "--", label=LABELS[k], alpha=0.9)
    ax.set_xlabel(r"radial mode $i$   [basis $T_i(x)$]")
    ax.set_ylabel(r"$\max_j\,|c_{ij}|$")
    ax.set_ylim(1e-17, 3e0)
    ax.set_xlim(-1.5, 61)
    panel_label(ax, "b")

    # -- (c) deviation of g from the analytic expression
    ax = axes[1, 0]
    keep = [(0.70, 0.11), (1.00, 0.11), (0.6937, 0.50), (0.6937, 0.73)]
    for k, (qq, cc) in enumerate(keep):
        sol = next((x for x in sols
                    if abs(x["chi"] - cc) < 5e-3 and abs(x["q"] - qq) < 5e-3), None)
        if sol is None:
            print(f"  [warn] panel(c): no solution for q={qq}, chi={cc}")
            continue
        gnum = sol["g"]
        gana = (1.0 + sol["rh"] / (2.0 * sol["rh"] / (1.0 - xg))[:, None]) ** 2
        dev = np.abs(gnum - gana).max(axis=1)
        lbl = rf"$q={sol['q']:.4f},\ \chi={sol['chi']:.2f}$"
        ax.semilogy(xg, np.maximum(dev, 1e-17), lw=1.4,
                    color=CURVE_COLORS[k], label=lbl)
    ax.set_yscale("log")
    ax.set_xlabel(r"$x = 1 - 2r_H/r$")
    ax.set_ylabel(r"$\max_\theta\,|g - (1+r_H/r)^2|$")
    ax.set_xlim(-1, 1)
    ax.set_ylim(3e-16, 3e-2)
    ax.axhline(1e-4, color="k", ls=":", lw=0.9)
    ax.text(0.99, 1.5e-4, r"threshold $10^{-4}$",
            fontsize=6.2, ha="right")
    ax.legend(loc="lower left", fontsize=6.0, ncol=2,
              handlelength=1.3, labelspacing=0.22, columnspacing=0.9)
    panel_label(ax, "c")

    # -- (d) derivative spectrum
    ax = axes[1, 1]
    for k, n in enumerate(NAMES):
        c = spectrum2d(ref[n])
        d = np.arange(c.shape[0]) ** 2 * np.abs(c).max(axis=1)
        ax.semilogy(np.arange(len(d)), np.maximum(d, 1e-18),
                    color=CURVE_COLORS[k % 6], lw=1.2,
                    ls="-" if k < 6 else "--", label=LABELS[k], alpha=0.9)
    ax.set_xlabel(r"radial mode $i$")
    ax.set_ylabel(r"$i^2\,\max_j|c_{ij}|$   $[\partial_x^2 \mathcal{F}]$")
    ax.set_ylim(1e-17, 1e4)
    panel_label(ax, "d")

    for a in axes.ravel():
        style_axes(a, grid=True)

    fig.suptitle(rf"Smoothness of the numerical solution "
                 rf"($q={ref['q']:.4f}$, $\chi={ref['chi']:.3f}$, $\alpha={ref['alpha']:.0f}$)",
                 fontsize=FS_TITLE + 0.5, y=0.995)
    save(fig, "fig10_spectral")

    # -- report numbers
    cg = spectrum2d(ref["g"])
    ang_g = np.abs(cg).max(axis=0)
    print(f"  [check] angular coefficients of g, j=0: {ang_g[0]:.4f}   max for j>=1: {ang_g[1:].max():.3e}")
    for n in ("Bt", "f"):
        c = spectrum2d(ref[n]); r = np.abs(c).max(axis=1)
        print(f"  [check] {n}: radial spectrum {r[0]:.3e} -> {r[-1]:.3e}  "
              f"({np.log10(r[0]/max(r[-1],1e-300)):.1f} decades)")


if __name__ == "__main__":
    main()
