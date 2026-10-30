"""
figstyle.py -- unified style module for all paper figures
================================================================================
Goal: clean, consistent and readable; every formula goes through LaTeX (matplotlib usetex).

Design principles
--------
1. **One set of sizes**: every figure is designed at its final typeset width
     single column 3.4 in (86 mm) / double column 7.0 in (178 mm)
2. **One set of font sizes**: axis labels 9 pt, ticks 8 pt, panel labels 10 pt bold
   (draw at the final width -> nothing is rescaled on embedding, so the sizes are the true sizes)
3. **One set of colours**:
     - signed fields     -> `RdBu_r` (diverging, perceptually uniform, colour-blind friendly)
     - non-negative scalars -> `inferno` (perceptually uniform, clearly readable in print)
     - curves             -> 6-colour qualitative cycle (a `TAB10` subset with maximal separation)
4. **One set of line styles**: curves 1.6 pt, contours 0.9 pt, reference lines 1.2 pt dashed
5. **Rasterisation trap**: use `rasterized=True` for 3D surfaces, otherwise the PDF becomes huge

Usage
----
    import sys; sys.path.insert(0, '<workspace>/07_figures')
    from figstyle import *
    apply_style()
    fig, axes = panel_grid(1, 2, width='double', height=2.6)
    ...
    save(fig, 'lowq_repro')
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ------------------------------------ sizes
IN = 1.0                                   # inches
COL_SINGLE = 3.40                          # PRD single-column text width
COL_DOUBLE = 7.00                          # PRD double-column text width

# ------------------------ font sizes (pt; drawn at the final width, never rescaled)
FS_LABEL = 9.0      # axis labels
FS_TICK  = 8.0      # ticks
FS_LEG   = 8.0      # legend
FS_PANEL = 10.0     # (a)(b) panel labels
FS_TITLE = 9.5      # subplot titles
FS_CBAR  = 8.0      # colorbar ticks

# ------------------------------------ lines
LW_CURVE = 1.6
LW_CONTOUR = 0.9
LW_REF = 1.2
MS_MARK = 3.6

# ------------------------------------ colours
CMAP_SIGNED = "RdBu_r"      # signed quantities (fields, $M_V$, ...)
CMAP_SCALAR = "RdYlBu_r"    # non-negative scalars ($\tilde T_H$, $\tilde S$): same lineage as
                            # the TemperatureMap of the original figures (blue-white-yellow-red),
                            # but more saturated, perceptually uniform and colour-blind friendly
CMAP_DENS   = "magma"       # spectral densities etc. (dark base, light top; clear on white)
CURVE_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c",
                "#d62728", "#9467bd", "#8c564b"]
GRAY_NA = "#b8b8b8"         # not computed / not reliable
GRID_C = "#d9d9d9"

__all__ = [
    "IN", "COL_SINGLE", "COL_DOUBLE",
    "FS_LABEL", "FS_TICK", "FS_LEG", "FS_PANEL", "FS_TITLE", "FS_CBAR",
    "LW_CURVE", "LW_CONTOUR", "LW_REF", "MS_MARK",
    "CMAP_SIGNED", "CMAP_SCALAR", "CMAP_DENS", "CURVE_COLORS", "GRAY_NA",
    "apply_style", "panel_grid", "panel_label", "save", "cbar",
    "style_axes", "figdir",
]

FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def figdir():
    os.makedirs(FIGDIR, exist_ok=True)
    return FIGDIR


# ════════════════════════════════════════════════════════════════════════════
def apply_style(usetex=True, tex_preamble_extra=""):
    """Global style. Must be called once before creating any figure."""
    pre = (r"\usepackage{amsmath}\usepackage{amssymb}\usepackage{bm}"
           r"\usepackage{siunitx}" + tex_preamble_extra)
    rcParams.update({
        # ---- LaTeX ----
        "text.usetex": bool(usetex),
        "text.latex.preamble": pre,
        "font.family": "serif" if not usetex else "serif",
        "mathtext.fontset": "cm",
        # ---- font sizes ----
        "font.size": FS_LABEL,
        "axes.labelsize": FS_LABEL,
        "axes.titlesize": FS_TITLE,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEG,
        # ---- line widths ----
        "axes.linewidth": 0.8,
        "lines.linewidth": LW_CURVE,
        "lines.markersize": MS_MARK,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6,
        "ytick.minor.width": 0.6,
        "xtick.major.size": 3.2,
        "ytick.major.size": 3.2,
        "xtick.minor.size": 1.8,
        "ytick.minor.size": 1.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        # ---- legend ----
        "legend.frameon": False,
        "legend.handlelength": 1.6,
        "legend.handletextpad": 0.5,
        "legend.labelspacing": 0.30,
        "legend.borderaxespad": 0.4,
        "legend.columnspacing": 1.1,
        # ---- miscellaneous ----
        "figure.dpi": 140,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
        "axes.unicode_minus": False,
    })


def style_axes(ax, grid=False, minor=True):
    """Unified axis appearance: thin frame, inward ticks, optional minor ticks and light grid."""
    for s in ("top", "right", "bottom", "left"):
        ax.spines[s].set_linewidth(0.8)
    if minor:
        ax.minorticks_on()
        ax.tick_params(which="minor", direction="in", top=True, right=True)
    if grid:
        ax.grid(True, which="major", color=GRID_C, lw=0.5, alpha=0.7, zorder=0)
        ax.set_axisbelow(True)
    return ax


def panel_grid(nrow, ncol, width="double", height=None, aspect=None,
               hspace=0.30, wspace=0.28, left=None, right=None,
               top=None, bottom=None):
    """Create a panel grid. width: 'single' | 'double' | a number in inches. height = total height (inch)."""
    W = {"single": COL_SINGLE, "double": COL_DOUBLE}.get(width, width)
    if height is None:
        height = W * nrow / max(ncol, 1) * 0.78
    fig, axes = plt.subplots(nrow, ncol, figsize=(W, height), squeeze=False,
                             constrained_layout=False)
    fig.subplots_adjust(hspace=hspace, wspace=wspace,
                        left=left, right=right, top=top, bottom=bottom)
    return fig, axes


def panel_label(ax, letter, loc="upper left", dx=0.02, dy=0.97, color="k",
                boxed=True):
    """Bold (a)/(b) panel label.

    boxed=True (the default) puts a **white background box** behind the label so that it stays
    readable even when curves or a legend overlap it -- the standard fix for a label hidden by data.
    """
    x = dx if "left" in loc else (0.98 if "right" in loc else 0.5)
    y = dy if "upper" in loc else 0.03
    ha = "left" if "left" in loc else ("right" if "right" in loc else "center")
    va = "top" if "upper" in loc else "bottom"
    bb = (dict(facecolor="white", edgecolor="none", pad=1.4, alpha=0.85)
          if boxed else None)
    ax.text(x, y, rf"\textbf{{({letter})}}", transform=ax.transAxes,
            ha=ha, va=va, fontsize=FS_PANEL, color=color, zorder=20, bbox=bb)


def cbar(fig, mappable, ax, label="", pad=0.015, fraction=0.046, aspect=18,
         extend="neither", ticks=None):
    """Colorbar in the unified style."""
    cb = fig.colorbar(mappable, ax=ax, pad=pad, fraction=fraction,
                      aspect=aspect, extend=extend)
    cb.ax.tick_params(labelsize=FS_CBAR, width=0.8, length=2.8, direction="in")
    cb.outline.set_linewidth(0.8)
    if label:
        cb.set_label(label, fontsize=FS_LABEL)
    if ticks is not None:
        cb.set_ticks(ticks)
    return cb


def save(fig, name, fmt=("pdf",), outdir=None, verbose=True):
    """Save to 07_figures/out/<name>.pdf (vector)."""
    d = outdir or figdir()
    os.makedirs(d, exist_ok=True)
    paths = []
    for f in fmt:
        p = os.path.join(d, f"{name}.{f}")
        fig.savefig(p)
        paths.append(p)
    plt.close(fig)
    if verbose:
        for p in paths:
            print(f"  wrote {p}  ({os.path.getsize(p)/1024:.0f} KB)")
    return paths
