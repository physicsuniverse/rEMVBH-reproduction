#!/usr/bin/env python3
"""
fig3_highq.py -- Fig. 3 `highq`: colour map of the domain of existence at high charge (q > 0.65)
=========================================================================
**Data source (established this round)**:
  <archive>/getallaxisols/axisolsofQ/
      uplimitofq/uplimitofqwithLRandphyp/upofqforchi<NN>.mx      NN = chi*100
  **77 files** in total, chi = 0.02 ... 0.745 (step 0.01); each file **pushes q upwards from q = 0.6947**
  up to the upper edge q_max(chi) of the domain of existence at that chi. This is the two-dimensional
  grid behind Fig. 3 of the paper: exported to `06_team_run2/domain_scan/validation/highq_grid_*.csv`,
  **24,968 points**, each record carrying {T_H, A_H, M, Q, Phi, P, J, OmH}, so T~_H and S~ follow directly.

  The narrow strip at q in [0.65, 0.6947], chi in [0.01, 0.9] in the lower-left corner is filled in from
  the **lowq archive** (`validation/student_grid_*.csv`) so that the panel matches q in [0.65, 1.37].

Reliability criterion as in Fig. 2: Smarr relative residual < 1e-4 (the archive carries no g-deviation
column; measured to be ~1e-11 throughout). Colorbar ranges are those measured from the paper:
"""
import os, sys, glob, csv
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import *                                    # noqa: F401,F403
from figstyle import (apply_style, panel_grid, panel_label, style_axes, save,
                      cbar, GRAY_NA, FS_LEG, FS_TITLE,
                      LW_REF, LW_CONTOUR)
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

# Colours: same structure as ColorData["TemperatureMap"] of the original figure (blue-white-red).
# Candidate comparison in preview/cmap_candidates.png.
CMAP = os.environ.get("FIG_CMAP", "RdYlBu_r")

SMRR = 1e-4
CB_A = (0.458, 1.327)          # paper Fig.3(a) colorbar
CB_B = (0.261, 0.777)          # paper Fig.3(b) colorbar
CONT_A = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
CONT_B = [0.3, 0.4, 0.5, 0.6, 0.7]
XLIM, YLIM = (0.645, 1.40), (0.0, 0.80)

HQ_HDR = "q,chi,alpha,maxBt,maxBp,TH,AH,M,Q,Phi,P,J,OmH,e9,e10".split(",")
LQ_HDR = "q,chi,alpha,maxBt,maxBp,M,TH,AH,Q,Phi,P,J,OmH".split(",")


def read_rows(paths, hdr):
    out = []
    for f in paths:
        with open(f) as fh:
            for row in csv.reader(fh):
                if len(row) != len(hdr):
                    continue
                try:
                    out.append(dict(zip(hdr, [float(x) for x in row])))
                except ValueError:
                    pass
    return out


def add_derived(rows):
    good = []
    for r in rows:
        if not all(np.isfinite(r[k]) for k in ("M", "TH", "AH", "Q", "Phi", "J", "OmH")):
            continue
        if r["M"] <= 0:
            continue
        r = dict(r)
        r["tth"] = 8 * np.pi * r["M"] * r["TH"]
        r["st"] = r["AH"] / (16 * np.pi * r["M"] ** 2)
        r["smarr"] = abs(r["TH"] * r["AH"] / 2 + 2 * r["J"] * r["OmH"]
                         + r["Q"] * r["Phi"] - r["M"]) / r["M"]
        r["rel"] = int(r["smarr"] < SMRR)
        good.append(r)
    return good


def regular_grid(rows, dq=0.0010, dc=0.0025):
    """interpolate onto a regular lattice; the DATA mask is the exact
    q_max(chi) edge of the archive, NOT the convex hull of griddata."""
    q = np.array([r["q"] for r in rows]); c = np.array([r["chi"] for r in rows])
    rel = np.array([r["rel"] for r in rows], float)
    qg = np.arange(q.min(), q.max() + dq, dq)
    cg = np.arange(max(c.min() - dc / 2, 1e-6), c.max() + dc, dc)
    uc = np.unique(np.round(c, 6))
    qmax = np.array([q[np.abs(c - v) < 1e-9].max() for v in uc])
    Q, C = np.meshgrid(qg, cg, indexing="ij")
    qmax_g = np.interp(cg, uc, qmax)
    inside = Q <= (qmax_g[None, :] + dq / 2.0)
    pts = np.column_stack([q, c])
    M = NearestNDInterpolator(pts, rel)(Q, C)
    vals = {}
    for key in ("tth", "st"):
        v = np.array([r[key] for r in rows])
        vals[key] = griddata(pts, v, (Q, C), method="linear")
    return qg, cg, vals["tth"], vals["st"], M, inside


def main():
    hq = add_derived(read_rows(sorted(glob.glob(DS + "/validation/highq_grid_*.csv")),
                               HQ_HDR))
    lq = add_derived(read_rows(sorted(glob.glob(DS + "/validation/student_grid_*.csv")),
                               LQ_HDR))
    lq = [r for r in lq if r["q"] >= XLIM[0]]          # only the left strip
    rows = hq + lq
    print(f"### highq archive: {len(hq)} points;  +lowq left strip: {len(lq)}")
    print(f"### q range {min(r['q'] for r in rows):.4f} .. {max(r['q'] for r in rows):.4f}"
          f"   chi range {min(r['chi'] for r in rows):.3f} .. {max(r['chi'] for r in rows):.3f}")
    print(f"### reliable {sum(r['rel'] for r in rows)} / {len(rows)}")

    rows_sorted = sorted(rows, key=lambda r: (r["q"], r["chi"]))
    with open(OUTD + "/fig3_highq.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["q", "chi", "TH_tilde", "S_tilde", "reliable"])
        for r in rows_sorted:
            w.writerow([f"{r['q']:.4f}", f"{r['chi']:.4f}",
                        f"{r['tth']:.6g}", f"{r['st']:.6g}", r["rel"]])
    print("### wrote data/fig3_highq.csv")

    qg, cg, A, B, M, inside = regular_grid(rows)
    print(f"### cells: reliable {int((inside & (M>=0.5)).sum())}, "
          f"unreliable(grey) {int((inside & (M<0.5)).sum())}, "
          f"no-archive-solution(white) {int((~inside).sum())}")

    def draw(swap, fname):
        """swap=False -> x=q, y=chi (task spec);  swap=True -> the paper's own
        Fig.3 frame (chi horizontal, q vertical)."""
        fig, axes = panel_grid(1, 2, width="double", height=2.85, wspace=0.34)
        grey_cmap = ListedColormap([GRAY_NA])
        tt = np.linspace(0, 1, 800)
        kn_q = np.sqrt(np.clip(1 - tt ** 2, 0, None))
        mappables = []
        for ax, V, clev, cbr, title, letter in zip(
                axes[0], (A, B), (CONT_A, CONT_B), (CB_A, CB_B),
                (r"$\widetilde T_H=8\pi M\,T_H$", r"$\widetilde S=A_H/(16\pi M^2)$"),
                ("a", "b")):
            good = np.where(inside & (M >= 0.5), V, np.nan)
            bad = np.where(inside & (M < 0.5), 1.0, np.nan)
            if swap:
                Xg, Yg, G, Bd = cg, qg, good, bad
                xlim, ylim = YLIM, XLIM
                xlab, ylab = r"$\chi=J/M^2$", r"$q=Q/M$"
            else:
                Xg, Yg, G, Bd = qg, cg, good.T, bad.T
                xlim, ylim = XLIM, YLIM
                xlab, ylab = r"$q=Q/M$", r"$\chi=J/M^2$"
            ax.pcolormesh(Xg, Yg, Bd, cmap=grey_cmap, shading="nearest",
                          vmin=0, vmax=1, rasterized=True)
            pc = ax.pcolormesh(Xg, Yg, G, cmap=CMAP, shading="nearest",
                               vmin=cbr[0], vmax=cbr[1], rasterized=True)
            mappables.append(pc)
            cs = ax.contour(Xg, Yg, G, levels=sorted(clev), colors="k",
                            linewidths=LW_CONTOUR, linestyles="--")
            ax.clabel(cs, inline=True, fontsize=FS_LEG - 1.0, fmt="%.2g",
                      inline_spacing=2)
            if swap:
                ax.plot(tt, kn_q, "k-", lw=LW_REF + 0.4, zorder=6)
                ax.axvline(0.1, color="k", lw=1.0, ls="-.", zorder=7)
            else:
                sel = kn_q >= XLIM[0]
                ax.plot(kn_q[sel], tt[sel], "k-", lw=LW_REF + 0.4, zorder=6)
            ax.set_xlim(*xlim); ax.set_ylim(*ylim)
            ax.set_xlabel(xlab); ax.set_ylabel(ylab)
            ax.set_title(title, fontsize=FS_TITLE, pad=5)
            style_axes(ax)
            panel_label(ax, letter, dx=0.025, dy=0.975)

        handles = [Line2D([], [], color="k", lw=LW_REF + 0.4,
                          label=r"$\chi^2+q^2=1$ (extremal KN)"),
                   Line2D([], [], color="k", lw=LW_CONTOUR, ls="--",
                          label=r"$\widetilde T_H,\widetilde S=\mathrm{const}$"),
                   Patch(facecolor=GRAY_NA, edgecolor="none",
                         label=r"solution exists, Smarr $>10^{-4}$"),
                   Patch(facecolor="white", edgecolor="0.6",
                         label="no archived solution")]
        if swap:
            handles.insert(1, Line2D([], [], color="k", lw=1.0, ls="-.",
                                     label=r"$\chi=0.1$"))
        axes[0][0].legend(handles=handles, loc="upper center",
                          bbox_to_anchor=(1.10, -0.20), ncol=2,
                          frameon=False, fontsize=FS_LEG - 0.5)
        for ax, pc in zip(axes[0], mappables):
            cbar(fig, pc, ax, label="", fraction=0.048, pad=0.025)
        save(fig, fname)

    draw(False, "fig3_highq")
    draw(True, "fig3_highq_paperframe")

    # ---- q_max(chi): data upper edge vs the paper-reader's measurements ----
    paper = {0.10: 1.355, 0.15: 1.320, 0.20: 1.284, 0.25: 1.252, 0.30: 1.213,
             0.35: 1.178, 0.40: 1.130, 0.45: 1.070, 0.50: 1.016, 0.55: 0.960,
             0.60: 0.901, 0.65: 0.845, 0.70: 0.780, 0.74: 0.712}
    print("\n### q_max(chi): archive vs paper-reader (01_paper_analysis)")
    print("  chi   q_max(data)  q_max(paper)   diff   q_max_reliable  KN bound")
    for c in sorted(paper):
        a = [r for r in rows if abs(r["chi"] - c) <= 0.0051]
        ar = [r for r in a if r["rel"] == 1]
        if not a:
            continue
        qm = max(r["q"] for r in a)
        qmr = max(r["q"] for r in ar) if ar else float("nan")
        print(f"  {c:4.2f}   {qm:8.4f}     {paper[c]:8.3f}   {qm-paper[c]:+7.4f}"
              f"      {qmr:8.4f}      {np.sqrt(max(0,1-c*c)):.4f}")


if __name__ == "__main__":
    main()
