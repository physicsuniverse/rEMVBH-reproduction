#!/usr/bin/env python3
"""
fig2_lowq.py -- Fig. 2 `lowq`: colour map of the domain of existence at low charge, in (q, chi)
=========================================================================
Data: `06_team_run2/domain_scan/validation/student_grid_*.csv` (full export of the student's archive,
47,512 points, alpha = 20, q in [0.0507, 0.6937] with step 0.001, chi in [0.01, 0.92]).
T~_H = 8 Pi M T_H,   S~ = A_H / (16 Pi M^2)。

Each of the three colours on the plot **has a definite meaning** (see REPORT_fig23.md section 6):
  * colour -- a solution exists in the archive **and** its Smarr relative residual is < 1e-4
              (the paper's second criterion, the g deviation, is absent: measured ~1e-11, never triggered)
  * grey   -- a solution exists **but** the Smarr relative residual is > 1e-4 ("computed but not reliable")
  * white  -- **no** archived solution at this (q, chi) (chi > chi_max(q), beyond the computed domain)

The presence mask uses the exact chi_max(q) boundary, not convex-hull interpolation -- otherwise cells
inside the hull but outside the data would be painted grey by mistake.
"""
import os, sys, glob, csv
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import *                                    # noqa: F401,F403
from figstyle import (apply_style, panel_grid, panel_label, style_axes, save,
                      cbar, FS_LEG, FS_TITLE, LW_REF, LW_CONTOUR, GRAY_NA)
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.interpolate import griddata, NearestNDInterpolator

apply_style()

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(_HERE)          # reproduction package root
# NOTE: this script reads the raw parameter-scan archive, which is NOT shipped with
# the package (too large). Set EMVBH_ARCHIVE to the directory that holds it, e.g.
#   export EMVBH_ARCHIVE=/path/to/archive
ROOT = os.environ.get("EMVBH_ARCHIVE", ROOT)
DS = ROOT + "/06_team_run2/domain_scan"
HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = os.path.join(HERE, "data")
os.makedirs(OUTD, exist_ok=True)

# Colours: same structure as ColorData["TemperatureMap"] of the original figure (blue-white-red),
# but perceptually more uniform and clearer in print. Candidate comparison in preview/cmap_candidates.png.
CMAP = os.environ.get("FIG_CMAP", "RdYlBu_r")

SMRR = 1e-4
CONT_A = [0.66, 0.77, 0.85, 0.90]
CONT_B = [1.129, 0.936, 0.805, 0.704]
CB_A = (0.439, 0.938)          # paper Fig.2(a) colorbar
CB_B = (0.532, 1.368)          # paper Fig.2(b) colorbar
XLIM, YLIM = (-0.02, 0.66), (0.0, 1.0)
STU_HDR = "q,chi,alpha,maxBt,maxBp,M,TH,AH,Q,Phi,P,J,OmH".split(",")


def load_student():
    out = []
    for f in sorted(glob.glob(DS + "/validation/student_grid_*.csv")):
        with open(f) as fh:
            for row in csv.reader(fh):
                if len(row) != len(STU_HDR):
                    continue
                try:
                    out.append(dict(zip(STU_HDR, [float(x) for x in row])))
                except ValueError:
                    pass
    return out


def add_derived(rows):
    for r in rows:
        r["tth"] = 8 * np.pi * r["M"] * r["TH"]
        r["st"] = r["AH"] / (16 * np.pi * r["M"] ** 2)
        r["smarr"] = abs(r["TH"] * r["AH"] / 2 + 2 * r["J"] * r["OmH"]
                         + r["Q"] * r["Phi"] - r["M"]) / r["M"]
        r["rel"] = int(r["smarr"] < SMRR)
    return rows


def regular_grid(rows, dq=0.0010, dc=0.0020):
    """interpolate onto a regular lattice; the DATA mask comes from the exact
    chi_max(q) edge, not from the convex hull of griddata."""
    q = np.array([r["q"] for r in rows]); c = np.array([r["chi"] for r in rows])
    rel = np.array([r["rel"] for r in rows], float)
    qg = np.arange(q.min(), q.max() + dq, dq)
    cg = np.arange(max(c.min(), 1e-6), c.max() + dc, dc)
    uq = np.unique(np.round(q, 6))
    cmax = np.array([c[np.abs(q - v) < 1e-9].max() for v in uq])
    Q, C = np.meshgrid(qg, cg, indexing="ij")
    cmax_g = np.interp(qg, uq, cmax)
    inside = C <= (cmax_g[:, None] + dc / 2.0)
    pts = np.column_stack([q, c])
    M = NearestNDInterpolator(pts, rel)(Q, C)
    vals = {}
    for key in ("tth", "st"):
        v = np.array([r[key] for r in rows])
        vals[key] = griddata(pts, v, (Q, C), method="linear")
    return qg, cg, vals["tth"], vals["st"], M, inside


def main():
    rows = add_derived(load_student())
    nrel = sum(r["rel"] for r in rows)
    print(f"### lowq archive: {len(rows)} points; reliable {nrel} "
          f"({100*nrel/len(rows):.0f}%); unreliable {len(rows)-nrel}")
    print(f"### colormap: {CMAP}")

    with open(OUTD + "/fig2_lowq.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["q", "chi", "TH_tilde", "S_tilde", "reliable"])
        for r in sorted(rows, key=lambda r: (r["q"], r["chi"])):
            w.writerow([f"{r['q']:.4f}", f"{r['chi']:.4f}",
                        f"{r['tth']:.6g}", f"{r['st']:.6g}", r["rel"]])
    print("### wrote data/fig2_lowq.csv")

    qg, cg, A, B, M, inside = regular_grid(rows)
    n_no = int((~inside).sum()); n_unrel = int((inside & (M < 0.5)).sum())
    print(f"### cells: reliable {int((inside & (M>=0.5)).sum())}, "
          f"unreliable(grey) {n_unrel}, no-archive-solution(white) {n_no}")

    fig, axes = panel_grid(1, 2, width="double", height=2.85, wspace=0.34)
    grey_cmap = ListedColormap([GRAY_NA])
    tt = np.linspace(0, 1, 500)
    kn_q, kn_c = np.sqrt(np.clip(1 - tt ** 2, 0, None)), tt

    mappables = []
    for ax, V, clev, cbr, title, letter in zip(
            axes[0], (A, B), (CONT_A, CONT_B), (CB_A, CB_B),
            (r"$\widetilde T_H=8\pi M\,T_H$", r"$\widetilde S=A_H/(16\pi M^2)$"),
            ("a", "b")):
        good = np.where(inside & (M >= 0.5), V, np.nan)
        bad = np.where(inside & (M < 0.5), 1.0, np.nan)
        ax.pcolormesh(qg, cg, bad.T, cmap=grey_cmap, shading="nearest",
                      vmin=0, vmax=1, rasterized=True)
        pc = ax.pcolormesh(qg, cg, good.T, cmap=CMAP, shading="nearest",
                           vmin=cbr[0], vmax=cbr[1], rasterized=True)
        mappables.append(pc)
        cs = ax.contour(qg, cg, good.T, levels=sorted(clev), colors="k",
                        linewidths=LW_CONTOUR, linestyles="--")
        ax.clabel(cs, inline=True, fontsize=FS_LEG - 1.0, fmt="%.3g",
                  inline_spacing=2)
        ax.plot(kn_q, kn_c, "k-", lw=LW_REF + 0.4, zorder=6)
        ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
        ax.set_xlabel(r"$q=Q/M$"); ax.set_ylabel(r"$\chi=J/M^2$")
        ax.set_title(title, fontsize=FS_TITLE, pad=5)
        style_axes(ax)
        panel_label(ax, letter, dx=0.025, dy=0.975)

    handles = [Line2D([], [], color="k", lw=LW_REF + 0.4,
                      label=r"$\chi^2+q^2=1$"),
               Line2D([], [], color="k", lw=LW_CONTOUR, ls="--",
                      label=r"$\widetilde T_H,\widetilde S=\mathrm{const}$"),
               Patch(facecolor=GRAY_NA, edgecolor="none",
                     label=r"solution exists, Smarr $>10^{-4}$"),
               Patch(facecolor="white", edgecolor="0.6",
                     label="no archived solution")]
    axes[0][0].legend(handles=handles, loc="upper center",
                      bbox_to_anchor=(1.10, -0.20), ncol=2, frameon=False,
                      fontsize=FS_LEG - 0.5)
    for ax, pc in zip(axes[0], mappables):
        cbar(fig, pc, ax, label="", fraction=0.048, pad=0.025)
    save(fig, "fig2_lowq")


if __name__ == "__main__":
    main()
