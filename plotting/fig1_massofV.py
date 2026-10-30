#!/usr/bin/env python3
"""
fig1_massofV.py —— Fig.1 ``MassofV``:  M_V/M (%) vs chi, six alpha curves.

Data
----
The curves are the Komar volume integral M_V recomputed by *us* for every chi of
the student's own fixed-alpha chi-continuation archives

    06_team_run2/crosscheck/student_alpha_archives/axisolsofchialpha-{5,8,10,12,15,25,30,35,40}.mx
    (alpha = 20 comes from numerics/axisols-example/axisolsofchiNofq-1.mx)

with the author's recipe (script 06_team_run2/alpha_scan/scripts/fig1_archive_mv.wls):

    GetNumExpr[{intofMVx,intofJVx} * 2 rh/(-1+x)^2 /. rules, {x,th}, fun, {xn,yn}, fields]
    MV = -1/(4 Pi) 2 Pi SPIntegrate2D[GetSPFun[..], {x,-1,1},{th,0,Pi}, 120, 60]

Overlaid: the points we obtained *independently* by our own alpha-continuation
Newton solves (06_team_run2/alpha_scan/alpha_waypoints.csv) and the vector-digitized
paper curve (02_numerics_repro/physics/data/MassofV_digitized.json).

Run:
    export PATH=/usr/local/texlive/2026/bin/x86_64-linux:$PATH
    export PYTHONPATH=/path/to/tools/pylibs   # only needed for the LaTeX rendering helpers
    python3 07_figures/fig1_massofV.py
"""
import csv
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(_HERE)          # reproduction package root
# NOTE: this script reads the raw parameter-scan archive, which is NOT shipped with
# the package (too large). Set EMVBH_ARCHIVE to the directory that holds it, e.g.
#   export EMVBH_ARCHIVE=/path/to/archive
ROOT = os.environ.get("EMVBH_ARCHIVE", ROOT)
sys.path.insert(0, os.path.join(ROOT, "07_figures"))
from figstyle import *          # noqa: F401,F403
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.join(ROOT, "07_figures")
RAW = os.path.join(ROOT, "06_team_run2/alpha_scan/fig1")
WAYPOINTS = os.path.join(ROOT, "06_team_run2/alpha_scan/alpha_waypoints.csv")
DIGITIZED = os.path.join(ROOT, "02_numerics_repro/physics/data/MassofV_digitized.json")

def wf(x):
    """Wolfram InputForm float -> python float ('1.23*^-6' -> 1.23e-6)."""
    if isinstance(x, (int, float)):
        return float(x)
    x = str(x).strip().replace("*^", "e")
    try:
        return float(x)
    except ValueError:
        return None


MAIN_ALPHAS = [5, 10, 15, 20, 30, 40]           # the six curves of the paper's Fig. 1
EXTRA_ALPHAS = [8, 12, 25, 35]                  # further student archives (shown faintly)


# ─────────────────────────────────────────────────────────────── data
def load_archive(a):
    """-> (chi[], MV/M[] , M[], MV[]) from our recomputation for archive alpha=a"""
    f = os.path.join(RAW, f"alpha{a}.csv")
    if not os.path.exists(f):
        raise FileNotFoundError(f)
    chi, mvm, M, MV = [], [], [], []
    with open(f) as fh:
        for r in csv.DictReader(fh):
            c, mv, m, v = (wf(r["chi"]), wf(r["MV_over_M_percent"]),
                           wf(r["M"]), wf(r["MV"]))
            if None in (c, mv, m, v):
                continue
            chi.append(c); mvm.append(mv); M.append(m); MV.append(v)
    o = np.argsort(chi)
    return (np.array(chi)[o], np.array(mvm)[o], np.array(M)[o], np.array(MV)[o])


def load_waypoints():
    """our independent alpha-continuation points"""
    out = []
    with open(WAYPOINTS) as fh:
        for r in csv.DictReader(fh):
            try:
                out.append(dict(alpha=wf(r["alpha"]), chi=wf(r["chi"]),
                                M=wf(r["M"]), mvm=wf(r["MV_over_M_percent"]),
                                maxBt=wf(r["maxBt"]),
                                res=wf(r["maxResidual_fromJson"])))
            except KeyError:
                continue
            if None in (out[-1]["alpha"], out[-1]["chi"], out[-1]["mvm"]):
                out.pop()
    return out


def load_digitized():
    with open(DIGITIZED) as fh:
        d = json.load(fh)
    return {int(k): sorted(v) for k, v in d.items()}


# ─────────────────────────────────────────────────────────────── assemble CSV
def write_data_csv(archives):
    d = os.path.join(HERE, "data")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "fig1_massofV.csv")
    n = 0
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["alpha", "chi", "M", "MV", "MV_over_M_percent"])
        for a in MAIN_ALPHAS + EXTRA_ALPHAS:
            chi, mvm, M, MV = archives[a]
            for k in range(len(chi)):
                w.writerow([f"{a:g}", f"{chi[k]:.10g}", f"{M[k]:.12g}",
                            f"{MV[k]:.12g}", f"{mvm[k]:.10g}"])
                n += 1
    print(f"  wrote {path}  ({n} rows, {len(archives)} alpha chains)")
    return path


# ─────────────────────────────────────────────────────────────── figure
def main():
    archives = {}
    for a in MAIN_ALPHAS + EXTRA_ALPHAS:
        archives[a] = load_archive(a)
    write_data_csv(archives)

    wp = load_waypoints()
    dig = load_digitized()
    col = {a: CURVE_COLORS[i % len(CURVE_COLORS)] for i, a in enumerate(MAIN_ALPHAS)}

    apply_style()
    fig, axes = panel_grid(1, 1, width="single", height=2.9)
    ax = axes[0, 0]

    # --- faint extra archives (8, 12, 25, 35) --------------------------------
    for i, a in enumerate(EXTRA_ALPHAS):
        chi, mvm, _, _ = archives[a]
        ax.plot(chi, mvm, color=GRAY_NA, lw=0.9, alpha=0.85, zorder=1,
                label=r"other $\alpha$" if i == 0 else None)

    # --- digitized paper curves (dashed targets) ----------------------------
    for i, a in enumerate(MAIN_ALPHAS):
        if a not in dig:
            continue
        p = np.array(dig[a])
        # The digitized target is drawn UNDER the solid curve with the same width.
        # Because the two agree to ~1e-3 percentage points (far below the line
        # width) the dashes are completely hidden wherever the ranges overlap --
        # the only visible part is where the digitized curve runs slightly longer
        # in chi.  That is the point: the agreement is at the line-width level.
        ax.plot(p[:, 0], p[:, 1], color=col[a], lw=0.75, ls=(0, (3.5, 2.0)),
                alpha=0.8, zorder=2,
                label="digitized paper curve" if i == 0 else None)

    # --- our recomputed Komar curves ----------------------------------------
    for a in MAIN_ALPHAS:
        chi, mvm, _, _ = archives[a]
        ax.plot(chi, mvm, color=col[a], lw=LW_CURVE, zorder=3,
                label=rf"$\alpha={a}$")

    # --- our independent alpha-continuation points --------------------------
    on_curve = [p for p in wp if any(abs(p["alpha"] - a) < 1e-6 for a in MAIN_ALPHAS)]
    off_curve = [p for p in wp if p not in on_curve]
    for i, p in enumerate(on_curve):
        ax.plot([p["chi"]], [p["mvm"]], marker="o", ms=4.2, mew=0.7, ls="none",
                mfc=col.get(p["alpha"], "k"), mec="white", zorder=6,
                label="independent Newton solve" if i == 0 else None)
    if off_curve:
        ax.plot([p["chi"] for p in off_curve], [p["mvm"] for p in off_curve],
                marker="s", ms=3.2, mew=0.7, ls="none", mfc="none", mec="#6b6b6b",
                zorder=5, label=r"intermediate $\alpha$ (this work)")

    ax.set_xlabel(r"$\chi$")
    ax.set_ylabel(r"$M_V/M\ (\%)$")
    ax.set_xlim(0.0, 0.78)
    ax.set_ylim(0.0, 27.0)
    style_axes(ax)

    # two-column legend, upper-left (empty region: curves start near zero)
    h, lab = ax.get_legend_handles_labels()
    order = [lab.index(rf"$\alpha={a}$") for a in MAIN_ALPHAS]
    order += [lab.index("digitized paper curve")]
    order += [lab.index("independent Newton solve")]
    if r"intermediate $\alpha$" in lab:
        order += [lab.index(r"intermediate $\alpha$")]
    if r"other $\alpha$" in lab:
        order += [lab.index(r"other $\alpha$")]
    leg = ax.legend([h[i] for i in order], [lab[i] for i in order],
                    ncol=2, loc="upper left", fontsize=6.0,
                    handlelength=1.5, columnspacing=0.9, labelspacing=0.26,
                    borderpad=0.22, borderaxespad=0.3,
                    title=r"curves: student archives $\;\;|\;\;$ markers: this work")
    leg.get_title().set_fontsize(6.0)
    leg._legend_box.align = "left"

    save(fig, "fig1_massofV")

    # ───────────────────────────────────────────── numbers for the report
    print("\n  chi = 0.5 readings (%).  archive-recomputed | digitized | this work")
    for a in MAIN_ALPHAS:
        chi, mvm, _, _ = archives[a]
        j = int(np.argmin(np.abs(chi - 0.5)))
        arc = mvm[j]
        pdig = np.array(dig[a])
        dg = np.interp(0.5, pdig[:, 0], pdig[:, 1])
        mine = [p["mvm"] for p in wp
                if abs(p["chi"] - 0.5) < 1e-9 and abs(p["alpha"] - a) < 1e-6]
        s = f"    alpha={a:>3}  {arc:8.4f}  {dg:8.4f}"
        if mine:
            s += f"  {mine[0]:8.4f}   rel.diff vs paper = {abs(mine[0]-dg)/dg*100:.4f} %"
        print(s)
    return 0


if __name__ == "__main__":
    sys.exit(main())
