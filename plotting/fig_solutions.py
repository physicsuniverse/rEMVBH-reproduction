#!/usr/bin/env python3
"""
fig_solutions.py --- redraw the paper's Figs. 5-8 ("sol1".."sol4"):
the eight solved functions {B_t, B_phi, A_t, A_phi, f, g, h, W} of the four
solutions I-IV as 3-D surfaces on the numerical coordinates (x, theta).

Differences from the student's original Mathematica Plot3D figure
-----------------------------------------------------------------
* perceptually-uniform colormaps instead of the rainbow sweep:
    signed  fields -> CMAP_SIGNED (RdBu_r), symmetric about zero
    non-neg fields -> viridis (both ends far from white)  [--cmap ...]
* every surface rasterised (keeps the PDF small while axes/labels stay vector)
* one fixed view for all four figures, so they can be compared side by side
* large panel labels, minimal in-panel text, true final-size fonts

Data (all four are the STUDENT'S OWN 100x50 solutions = the paper's source):
    data/sol1_fields.dat ... data/sol4_fields.dat   (written by
    data/export_fig_data.wls from
      <archive>/gpt/tensorNoG100x50-qIndex*.mx
      <archive>/gpt/highres_upq_3points/chi10-q1p2527-100x50.mx)

Usage
-----
    export PATH=/usr/local/texlive/2026/bin/x86_64-linux:$PATH
    export PYTHONPATH=<root>/tools/pylibs
    python3 07_figures/fig_solutions.py [--cmap family|inferno] [--preview]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from figstyle import (COL_SINGLE, COL_DOUBLE, FS_LABEL, FS_TICK, FS_TITLE,
                      CMAP_SIGNED, apply_style, save, figdir, GRID_C)

# Colormap for the NON-NEGATIVE fields, chosen so that BOTH ends are far from
# white: 'viridis' runs dark purple -> yellow.  (figstyle's CMAP_SCALAR is
# 'inferno', whose black low end hides the large near-zero plateaus; the
# previous white->red ramp hid them just as badly on a white page.)
CMAP_POS = "viridis"
CMAP_POS_CHOICES = ["viridis", "cividis", "RdYlBu_r", "inferno", "family"]

FIELDS = ["Bt", "Bp", "At", "Ap", "f", "g", "h", "W"]
LABELS = {
    "Bt": r"$B_t$", "Bp": r"$B_\varphi$", "At": r"$A_t$", "Ap": r"$A_\varphi$",
    "f":  r"$f$",   "g":  r"$g$",         "h":  r"$h$",     "W":  r"$W$",
}
def is_signed(raw, rel_tol=1e-6):
    """True iff the solution's own grid values really span both signs.

    The decision is made PER PANEL from the RAW grid data (not the spectrally
    interpolated surface, which can ring by ~1e-6 and fake a sign change), and
    not per field name.  A diverging map is only meaningful for a quantity that
    crosses zero: using one for a field that merely touches zero would paint
    its lower part in the pale middle of the map -- the "the surface looks half
    missing" effect reported by the user.
    """
    lo, hi = float(np.min(raw)), float(np.max(raw))
    scale = max(1.0, abs(lo), abs(hi))
    return lo < -rel_tol * scale and hi > rel_tol * scale


def colour_range(raw, Z, signed):
    """(vmin, vmax) for the colour scale of one panel.

    `signed` panels are symmetric about zero.  For a non-negative panel whose
    lower bound is less than 1% of its upper bound the lower end is snapped to
    exactly 0 -- it is physically zero and the colour bar would otherwise print
    meaningless labels such as 3.8e-4.
    """
    if signed:
        m = max(abs(float(np.min(raw))), abs(float(np.max(raw))))
        return -m, m
    lo, hi = float(np.min(raw)), float(np.max(raw))
    if lo < 0.01 * hi:
        lo = 0.0
    return lo, hi

# design width: the paper is revtex4-2 "preprint" (single column, ~6.5 in text
# width) and includes these at 0.72\linewidth  ->  4.68 in final size.
# Landscape 2 rows x 4 columns, drawn at exactly the text width of the journal.
# History: the original 4x2 portrait layout was 4.55 x 8.45 in (aspect 1.86) and
# overflowed the revtex page -> LaTeX reported "Float too large for page"
# (7.4 pt for sol1-4, 74.8 pt for Fig. 9) and SILENTLY TRUNCATED the caption.
# The fix was to go landscape instead of shrinking the artwork.
# JHEP is single-column with generous margins: its text width is 5.95 in (A4,
# 595x842 pt, with the text block spanning x = 83.3 .. 511.9 pt).  Making the
# figure exactly that wide means includegraphics[width=textwidth] applies NO
# scaling, so the 9 pt type is really 9 pt.  (At the previous 7.00 in the JHEP
# text width would have shrunk the artwork by 14%, i.e. 9 pt -> 7.8 pt.)
FIG_W = 5.95                # = JHEP text width
FIG_H = 3.62 * 5.95 / 7.00  # ~3.077 in, keeps the panel proportions

VIEW = dict(elev=24, azim=-55)
VIEW_C = dict(elev=22, azim=-62)

# ---- presentation style for the surface ("long strip" mitigation) -------
# The fields are nearly x-only (median theta-share of the variation is 0.004),
# so a bare surface reads as a single ridge.  A faint lattice makes the flat
# parts read as a surface.  A shorter theta box also helps.
# ---- "aesthetics" looks (task-48) ---------------------------------------
# The user asked for the 3-D panels to look more "elegant / generous".  Four
# things were hurting that:
#   (1) the grey lattice on the three box walls is pure chartjunk -- it covers a
#       lot of area and carries no information;
#   (2) eight colour bars per figure (one per panel) fragment the layout;
#   (3) the large flat near-zero plateaus visually "press down" the surface;
#   (4) the panels' z ranges differ a lot, which reads as ragged.
# Each LOOK below answers those differently.
LOOK = {
    # current delivered look: grey walls + grid + per-panel bar
    "G":  dict(panes="walls", grid=True,  cbar="panel", zticks=False, barw=0.038),
    # P1 minimal: no walls, no grid, only a whisper-light floor; hairline
    #            frameless colour bar
    "P1": dict(panes="floor", grid=False, cbar="panel", zticks=False, barw=0.024),
    # P2 informative: as P1 but the absolute z scale is printed on the panel again
    #            (two ticks) and the floor keeps a barely-there grid for depth
    "P2": dict(panes="floor", grid="floor", cbar="panel", zticks=True, barw=0.024),
    # P3 low-camera: as P1 but a much lower elevation so the flat plateaus read as
    #            receding planes instead of pressing the surface down
    # P3 tall-z: as P1 but a noticeably taller z box, so the relief of the
    #          surface gets more of the panel instead of the flat plateaus
    "P3": dict(panes="floor", grid=False, cbar="panel", zticks=False, barw=0.024),
    # ---- task-49: the surface must look like a SMOOTH analytic solution ----
    # The sparse quad lattice of "style G" made the surface look like assorted
    # tiles ("rough, not trustworthy"), even though the solution is a spectral
    # series whose coefficients decay by 15 orders.  Both Q looks draw the
    # surface at full resolution with NO edges at all.
    #   Q1: no panes whatsoever -- the surface floats on pure white
    #   Q2: same, but a whisper-light floor is kept for depth/contrast
    "Q1": dict(panes="none",  grid=False, cbar="panel", zticks=False, barw=0.024),
    "Q2": dict(panes="floor", grid=False, cbar="panel", zticks=False, barw=0.024),
    # Q3: Q2 + the diverging map is CUT at the centre, so zero maps to a light
    #     blue / light red instead of white.  Keeps sign-changing fields fully
    #     visible even where they sit at zero (the white band of RdBu_r is
    #     invisible on a white page).
    "Q3": dict(panes="floor", grid=False, cbar="panel", zticks=False, barw=0.024,
               divcut=True),
    # Q4: Q3 with antialiased=True (polygon edges blended instead of hard)
    "Q4": dict(panes="floor", grid=False, cbar="panel", zticks=False, barw=0.024,
               divcut=True, aa=True),
}

STYLE = "A"
STYLE_CFG = {
    #        box_aspect            wire   view
    "A": dict(box=(1.20, 1.00, 0.92), wire=None, view=None),
    "B": dict(box=(1.20, 1.00, 0.92), wire=(6, 2), view=None),
    "C": dict(box=(1.32, 0.78, 0.95), wire=None, view=VIEW_C),
    "D": dict(box=(1.32, 0.78, 0.95), wire=(6, 2), view=VIEW_C),
    # E: quad-edge lattice drawn by plot_surface itself (lead's suggestion) --
    # a coarser rcount/ccount keeps the surface smooth while showing the mesh
    "E": dict(box=(1.32, 0.78, 0.95), wire=None, view=VIEW_C, edge=(46, 18)),
    "F": dict(box=(1.20, 1.00, 0.92), wire=None, view=None,    edge=(46, 18)),
    # G: the lattice spaced for print (lines every ~1/17 inch on the page)
    "G": dict(box=(1.30, 0.82, 0.95), wire=None, view=VIEW_C, edge=(30, 12)),
    # --- aesthetics candidates: same data treatment, different presentation
    "P1": dict(box=(1.32, 0.80, 1.00), wire=None, view=VIEW_C, edge=(30, 12)),
    "P2": dict(box=(1.32, 0.80, 1.00), wire=None, view=VIEW_C, edge=(30, 12)),
    "P3": dict(box=(1.30, 0.74, 1.34), wire=None, view=VIEW_C, edge=(30, 12)),
    "Q1": dict(box=(1.30, 0.84, 1.10), wire=None, view=VIEW_C, edge=None),
    "Q2": dict(box=(1.30, 0.84, 1.10), wire=None, view=VIEW_C, edge=None),
    "Q3": dict(box=(1.30, 0.84, 1.10), wire=None, view=VIEW_C, edge=None, shade=False),
    "Q4": dict(box=(1.30, 0.84, 1.10), wire=None, view=VIEW_C, edge=None, shade=False,
               aa=True),
}
BOX  = (1.20, 1.0, 0.92)
# Evaluation grid for the SURFACE.  matplotlib's plot_surface colours each quad
# FLAT (one colour per quad, no Gouraud interpolation), so on a steep gradient a
# coarse grid shows visible colour steps -- at 400 dpi that reads as banding and
# makes an analytic solution look like assembled tiles.  The archived grid is
# 100x50, but the solution is a spectral series F = sum c_ij T_i(x) cos(2 j th)
# whose coefficients decay by ~15 orders, so evaluating it on a finer grid is
# exact, not an interpolation artefact.
XD = np.linspace(-1.0, 1.0, 360)
# theta in [0, pi/2]: this IS the solved domain.  The equator theta = pi/2 is a
# Neumann boundary (reflection symmetry), so [0, pi/2] is the fundamental domain.
# (The cos(2 j theta) spectral basis is exactly even about pi/2, verified to
# 4e-16, so the full meridian could be displayed -- but the user prefers the
# honest solved domain, and the fields are overwhelmingly x-dependent anyway:
# the median theta-share of the total variation is 0.004.)
THD = np.linspace(0.0, np.pi / 2, 151)


# ─────────────────────────────────────────────────────────────── data
def cheb_gauss_x(n):
    return -np.cos((2 * np.arange(n) + 1) * np.pi / (2 * n))


def theta_grid(n):
    return (2 * np.arange(n) + 1) * np.pi / (4 * n)


def load_solution(tag_or_path):
    """Read a solution file in the shared layout used by
    `05_solution_plots/data/sol_*.txt` (which fig_solutions.py writes too):

        L1  # q chi alpha rh Nx Ny
        L2   <values>
        L3  # TH AH M Q Phi P J OmegaH
        L4   <values>
        L5  # fields: Bt Bp At Ap f g h W (1 line per field, x slowest)
        L6..L13  8 field lines of Nx*Ny values

    Optional trailing comments (written by 07_figures/data/export_fig_data.wls)
    supply the grids and newSmarr; if absent the Chebyshev-Gauss / midpoint
    grids are reconstructed from Nx and Ny (they are fixed by the solver).
    `tag_or_path` may be a tag ("sol1") or a path to any such file.
    """
    path = tag_or_path if os.sep in tag_or_path else \
        os.path.join(HERE, "data", f"{tag_or_path}_fields.dat")
    lines = [l.rstrip("\n") for l in open(path)]
    meta = np.array([float(v) for v in lines[1].split()])
    phys = np.array([float(v) for v in lines[3].split()])
    Nx, Ny = int(meta[4]), int(meta[5])
    raw = np.array([[float(v) for v in lines[5 + k].split()] for k in range(8)])

    x = th = None
    new_smarr = float("nan")
    for i, l in enumerate(lines):
        if l.startswith("# grid x") and i + 1 < len(lines):
            x = np.array([float(v) for v in lines[i + 1].split()])
        elif l.startswith("# grid theta") and i + 1 < len(lines):
            th = np.array([float(v) for v in lines[i + 1].split()])
        elif l.startswith("# newSmarr"):
            try:
                new_smarr = float(l.split()[2])
            except (IndexError, ValueError):
                pass
    if x is None or x.size != Nx:
        x = cheb_gauss_x(Nx)
    if th is None or th.size != Ny:
        th = theta_grid(Ny)

    F = {nm: raw[k].reshape(Nx, Ny) for k, nm in enumerate(FIELDS)}
    return dict(tag=os.path.basename(path).split("_fields")[0],
                path=path, x=x, th=th, F=F,
                q=float(meta[0]), chi=float(meta[1]), alpha=float(meta[2]),
                rh=float(meta[3]), Nx=Nx, Ny=Ny,
                phys=phys, newSmarr=new_smarr)


def interp_field(values, x, th, xd, thd):
    """clean two-step spectral interpolation (rewritten for clarity)."""
    Nx, Ny = x.size, th.size
    # --- x direction: values = Bx @ C  with Bx[i,n] = T_n(x_i)
    Bx = np.cos(np.arange(Nx)[None, :] * np.arccos(np.clip(x, -1, 1))[:, None])
    C = np.linalg.solve(Bx, values)                    # [Nx, Ny]
    Bxd = np.cos(np.arange(Nx)[None, :] * np.arccos(np.clip(xd, -1, 1))[:, None])
    tmp = Bxd @ C                                      # [len(xd), Ny]
    # --- theta direction: tmp = Ct @ Bth with Bth[j,m] = cos(2 m th_j)
    Bth = np.cos(2.0 * np.arange(Ny)[None, :] * th[:, None])       # [Ny, Ny]
    Ct = np.linalg.solve(Bth, tmp.T)                   # [Ny, len(xd)]
    Bthd = np.cos(2.0 * np.arange(Ny)[None, :] * thd[:, None])     # [len(thd), Ny]
    return (Bthd @ Ct).T                               # [len(xd), len(thd)]


# ─────────────────────────────────────────────────────────────── figure
def make_fig(sol, cmap_mode=CMAP_POS):
    x, th, F = sol["x"], sol["th"], sol["F"]
    X, T = np.meshgrid(XD, THD / np.pi, indexing="ij")
    look = LOOK.get(STYLE, LOOK["G"])

    apply_style()
    fig, axes = plt.subplots(2, 4, figsize=(FIG_W, FIG_H), squeeze=False,
                             subplot_kw={"projection": "3d"})
    # At the JHEP text width each of the 4 columns is only ~1.5 in, so the
    # columns need a big relative gap: each panel must leave room for its own
    # colour bar plus tick labels without touching the next panel.  `right` is
    # pulled in so the 4th column's bar stays INSIDE the nominal 5.95 in box --
    # otherwise bbox=tight widens the PDF and includegraphics[width=textwidth]
    # would shrink the artwork after all.
    if STYLE == "G":
        fig.subplots_adjust(left=0.004, right=0.930, top=0.905, bottom=0.004,
                            wspace=0.62, hspace=0.10)
    else:
        # hairline bars need much less room -> the panels themselves get bigger,
        # which is most of what makes the figure look "generous".  P2 also prints
        # the z ticks on the LEFT of each panel, which needs a wider gutter.
        ws = 0.52 if look["zticks"] else 0.30
        # right pulled in so the 4th column's bar + labels stay inside the
        # nominal 5.95 in box (otherwise bbox=tight widens the PDF to ~6.06 in
        # and includegraphics[width=textwidth] shrinks the artwork by 1.8%)
        fig.subplots_adjust(left=0.004, right=0.937, top=0.930, bottom=0.004,
                            wspace=ws, hspace=0.06)

    rd = plt.get_cmap(CMAP_SIGNED)
    # Colormap rule (unchanged): non-negative fields -> viridis, whose LOW end is
    # dark purple and HIGH end yellow so neither end is near white (a white
    # surface on a white page is invisible -- the original "looks half empty"
    # complaint).  Genuinely sign-changing fields -> RdBu_r about zero.
    cmap_pos = (matplotlib.colors.LinearSegmentedColormap
                .from_list("rdbu_pos", rd(np.linspace(0.5, 0.92, 256)))) \
        if cmap_mode == "family" else plt.get_cmap(cmap_mode)
    # RdBu_r with the central (near-white) band removed: zero becomes a light
    # blue / light red, so a sign-changing field stays visible at zero too.
    cmap_div_cut = matplotlib.colors.LinearSegmentedColormap.from_list(
        "rdbu_divcut",
        rd(np.concatenate([np.linspace(0.00, 0.42, 160),
                           np.linspace(0.58, 1.00, 160)])))

    for k, nm in enumerate(FIELDS):
        ax = axes.ravel()[k]
        Z = interp_field(F[nm], x, th, XD, THD)
        raw = F[nm]
        signed = is_signed(raw)
        vmin, vmax = colour_range(raw, Z, signed)
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = rd if signed else cmap_pos
        if signed and look.get("divcut"):
            cmap = cmap_div_cut
        cfg = STYLE_CFG[STYLE]

        ecfg = cfg.get("edge")
        if ecfg:
            surf = ax.plot_surface(X, T, Z, cmap=cmap, norm=norm,
                                   linewidth=0.11, edgecolor=(0.18, 0.18, 0.18, 0.40),
                                   antialiased=True, shade=True,
                                   rcount=ecfg[0], ccount=ecfg[1],
                                   rasterized=True)
        else:
            surf = ax.plot_surface(X, T, Z, cmap=cmap, norm=norm,
                                   linewidth=0,
                                   antialiased=bool(cfg.get("aa", False)),
                                   # shade=False: matplotlib's lighting multiplies the
                                   # colormap colour, which (a) makes the colour->value
                                   # mapping inexact (the colour bar then lies) and
                                   # (b) at full resolution produces a moire hatch on
                                   # steep faces.  Without it the colour is the value.
                                   shade=bool(cfg.get("shade", True)),
                                   rcount=Z.shape[0], ccount=Z.shape[1],
                                   rasterized=True)

        # ---- hairline colour bar, still carrying 3 readable magnitudes
        cax = ax.inset_axes([1.020, 0.34, look["barw"], 0.42])
        cb = fig.colorbar(surf, cax=cax)
        cb.outline.set_visible(False)
        cb.ax.tick_params(labelsize=5.6, width=0.0, length=0.0, pad=0.9)
        ticks = np.linspace(norm.vmin, norm.vmax, 3)
        scale = max(1.0, abs(norm.vmax))
        cb.set_ticks(ticks)
        cb.ax.set_yticklabels(
            ["0" if abs(t) < 1e-9 * scale else f"{t:.3g}" for t in ticks])
        cb.ax.yaxis.set_ticks_position("right")

        # ---- z limits: 12% of the range under the surface so it never lies flat
        # ON the floor pane (A_t, f, h sit at their lower bound over big areas).
        zlo, zhi = float(np.min(Z)), float(np.max(Z))
        zr = zhi - zlo if zhi > zlo else max(abs(zhi), 1.0)
        ax.set_zlim(zlo - 0.12 * zr, zhi + 0.02 * zr)

        view = dict(cfg["view"] or VIEW)
        if look.get("elev") is not None:
            view["elev"] = look["elev"]
        ax.view_init(**view)
        ax.set_box_aspect(cfg["box"])

        # ---- axes: only what carries information
        ax.set_xlabel(r"$x$", labelpad=-6, fontsize=FS_LABEL)
        ax.set_ylabel("")
        ax.set_zlabel("")
        ax.set_xticks([-1, 0, 1])
        # NOTE: the theta plotting coordinate is T = theta/pi (see the meshgrid
        # above), so the ticks must be placed at theta/pi too.  Placing them at
        # the raw radians [0, pi/4, pi/2] stretched the axis to 1.571 while the
        # data only reaches 0.5, leaving ~2/3 of the axis empty -- and mislabelled
        # the tick positions.
        ax.set_yticks(np.array([0, np.pi / 4, np.pi / 2]) / np.pi)
        ax.set_yticklabels([r"$0$", r"$\pi/4$", r"$\pi/2$"])
        ax.xaxis.set_tick_params(pad=-1.5, labelsize=FS_TICK - 0.5)
        ax.yaxis.set_tick_params(pad=0.5, labelsize=FS_TICK - 1.5)
        if look["zticks"]:
            # Two ticks only: the panel's own min and max, restoring the absolute
            # vertical scale the hairline colour bar cannot convey.  matplotlib
            # draws the z ticks at the RIGHT edge, i.e. straight into the colour
            # bar, so juggle the z axis over to the left.
            ax.zaxis._axinfo["juggled"] = (1, 2, 0)
            ax.set_zticks([zlo, zhi])
            zscale = max(1.0, abs(zlo), abs(zhi))
            fmt = lambda v: "0" if abs(v) < 1e-6 * zscale else f"{v:.3g}"
            ax.set_zticklabels([fmt(zlo), fmt(zhi)])
            ax.zaxis.set_tick_params(pad=1.0, labelsize=5.6)
        else:
            ax.set_zticklabels([])
            ax.zaxis.set_tick_params(pad=-0.5, labelsize=0)

        # ---- panes: walls only if the look asks for them
        if look["panes"] == "walls":
            for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
                axis.pane.set_facecolor((0.921, 0.921, 0.941, 1.0))
                axis.pane.set_alpha(1.0)
                axis.pane.set_edgecolor("0.72")
        else:
            for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
                axis.pane.set_alpha(0.0)
                axis.pane.set_edgecolor("none")
            if look["panes"] == "floor":
                # a whisper-light floor: anchors the surface in depth and keeps
                # contrast for pale surfaces, without the "cage"
                ax.zaxis.pane.set_facecolor((0.945, 0.945, 0.960, 1.0))
                ax.zaxis.pane.set_alpha(0.85)
        if look["grid"] is True:
            ax.grid(True, linewidth=0.3, color="0.82", zorder=0)
        elif look["grid"] == "floor":
            ax.grid(True, linewidth=0.22, color="0.90", zorder=0)
            for axis in (ax.xaxis, ax.yaxis):
                axis._axinfo["grid"]["linewidth"] = 0.0
        else:
            ax.grid(False)

        ax.set_title(LABELS[nm], fontsize=FS_TITLE + 2.5, fontweight="bold",
                     pad=1.0)
    return fig

    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cmap", default=CMAP_POS, choices=CMAP_POS_CHOICES,
                    help="colormap for the non-negative fields (default viridis); "
                         "'family' = the old white->red half of RdBu_r")
    ap.add_argument("--preview", action="store_true",
                    help="also write a PNG preview for visual inspection")
    ap.add_argument("--only", default="", help="comma list, e.g. sol1,sol4")
    ap.add_argument("--style", default="G", choices=sorted(STYLE_CFG),
                    help="presentation style; see STYLE_CFG")
    ap.add_argument("--fig9", action="store_true",
                    help="additionally render Fig. 9 for the (q,chi)=(0.6937,0.500) solution")
    args = ap.parse_args()
    global STYLE
    STYLE = args.style

    tags = [t.strip() for t in args.only.split(",") if t.strip()] or \
           ["sol1", "sol2", "sol3", "sol4"]
    for tag in tags:
        sol = load_solution(tag)
        print(f"### {tag}: q={sol['q']:.4f} chi={sol['chi']:.3f} "
              f"grid {sol['Nx']}x{sol['Ny']}  newSmarr={sol['newSmarr']:.3e}")
        fig = make_fig(sol, cmap_mode=args.cmap)
        save(fig, f"fig_{tag}")
        if args.preview:
            fig2 = make_fig(sol, cmap_mode=args.cmap)
            p = os.path.join(HERE, "preview", f"{tag}_{args.cmap}.png")
            os.makedirs(os.path.dirname(p), exist_ok=True)
            fig2.savefig(p, dpi=110)
            plt.close(fig2)
            print("   preview", p)

    if args.fig9:
        # Fig. 9 of the paper: the (q, chi, alpha) = (0.6937, 0.500, 20) solution.
        # That solution only exists on the 60x20 archive grid.
        p9 = os.path.join(HERE, "data", "sol9_q6937_chi0p50_fields.dat")
        sol = load_solution(p9)
        print(f"### fig9: q={sol['q']:.4f} chi={sol['chi']:.3f} "
              f"grid {sol['Nx']}x{sol['Ny']}")
        fig = make_fig(sol, cmap_mode=args.cmap)
        save(fig, "fig9_profiles")
        if args.preview:
            fig2 = make_fig(sol, cmap_mode=args.cmap)
            q = os.path.join(HERE, "preview", f"fig9_{args.cmap}.png")
            fig2.savefig(q, dpi=110)
            plt.close(fig2)
            print("   preview", q)


if __name__ == "__main__":
    main()
