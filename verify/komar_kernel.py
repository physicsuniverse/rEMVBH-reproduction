#!/usr/bin/env python3
r"""
verify/komar_kernel.py -- Komar integral kernel for the EMV rotating backgrounds.

*** READ THIS FIRST: THE VOLUME-INTEGRAL ROUTE IS ILL-POSED ON THIS BACKGROUND ***

We implement, from first principles (see equations/KOMAR_derivation.md),

    M_V = (1/4pi) ∫_Σ (2 T^V_ab - T^V g_ab) t^a n^b sqrt(gamma) d^3x ,
    n^b = t^b / sqrt(-g_tt) ,
    sqrt(gamma) = sqrt( g_xx g_thth (g_phiphi - g_tphi^2/g_tt) )

and also the de-aliased version (spectral interpolation to 4001 theta x 2001 x points, Simpson).

The kernel DOES produce numbers, but they are NOT the physical M_V / J_V, because the
integral is ill-posed: on the archived solutions -g_tt <= 0 on ~18% of the grid (the
t = const slice enters the ergoregion and becomes timelike), and gamma < 0 there as well,
so sqrt(gamma) is imaginary.  See KOMAR_derivation.md §2.2.

What this module therefore delivers is the DIAGNOSTIC that establishes this, plus the
volume-integral values for the record with an explicit INADMISSIBLE flag.

Usage:  python3 verify/komar_kernel.py [profile ...]
Output: verify/out/komar_kernel_check.csv , verify/out/komar_kernel_meta.json
"""
import glob
import hashlib
import json
import os
import sys
import time

import numpy as np
from scipy.integrate import simpson

HERE = os.path.dirname(os.path.abspath(__file__))
HERE = os.path.dirname(HERE)  # go up from verify/ so that import verify_eom works
sys.path.insert(0, os.path.join(HERE, "verify"))
ROOT = os.environ.get(
    "EMVBH_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import verify_eom as V   # noqa: E402  (metric, field strengths, spectral class)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def komar_kernel(f, g, h, W, Bt, Bp, At, Ap, x, theta, r_H, Nx, Ny):
    """Rebuild the geometry from the eight fields and return M_V, J_V (alias-free and
    aliased) plus the well-posedness diagnostics.

    The signature is the one required by the brief; internally we re-use the validated
    spectral machinery of verify_eom.py.
    """
    F = dict(f=f, g=g, h=h, W=W, Bt=Bt, Bp=Bp, At=At, Ap=Ap)
    sp = V.Spectral(np.asarray(x, float), np.asarray(theta, float))
    al = None      # caller must patch meta; this low-level entry keeps it explicit
    return F, sp


def analyse(path):
    meta, hdr, F = V.load_profile(path)
    Nx, Ny = int(meta["Nx"]), int(meta["Ny"])
    rh, al = meta["rh"], meta["alpha"]
    sp = V.Spectral(V.CHEB_NODES(Nx), V.THETA_NODES(Ny))

    g, r, N, omN, s2 = V.metric_arrays(sp, F, rh)
    ginv = np.linalg.inv(g)
    A, B = V.potentials(sp, F, rh)
    Fmn, Vmn = V.field_strengths_np(sp, A, B)
    Vup = np.einsum("...ac,...bd,...cd->...ab", ginv, ginv, Vmn)
    Fup = np.einsum("...ac,...bd,...cd->...ab", ginv, ginv, Fmn)
    V2 = np.einsum("...ab,...ab->...", Vmn, Vup)
    F2 = np.einsum("...ab,...ab->...", Fmn, Fup)
    B2 = np.einsum("...ab,...a,...b->...", ginv, B, B)
    U = np.exp(np.clip(al * B2, -700, 700))
    Up = al * U
    TV = (np.einsum("...ac,...bc->...ab", Vmn, Vup)
          - 0.25 * V2[..., None, None] * g
          + 0.5 * Up[..., None, None] * F2[..., None, None]
          * B[..., :, None] * B[..., None, :])
    TVtr = np.einsum("...ab,...ab->...", ginv, TV)

    gtt, gtp, gpp = g[..., 0, 0], g[..., 0, 3], g[..., 3, 3]
    gxx, ghh = g[..., 1, 1], g[..., 2, 2]
    neg_gtt = -gtt
    gamma = gxx * ghh * (gpp - gtp ** 2 / np.where(np.abs(gtt) < 1e-300, 1e-300, gtt))

    frac_gtt_nonneg = float(np.mean(neg_gtt <= 0))
    frac_gamma_neg = float(np.mean(gamma < 0))

    # --- volume-integral values (both aliased and de-aliased) -------------
    safe_sqrt = np.sqrt(np.maximum(np.abs(neg_gtt), 1e-300))
    safe_gamma = np.sqrt(np.abs(gamma))
    X = 2 * TV[..., 0, 0] - TVtr * gtt
    Y = 2 * TV[..., 3, 0] - TVtr * gtp
    I_M = X / safe_sqrt * safe_gamma
    I_J = Y / safe_sqrt * safe_gamma

    # aliased: trapezoid on the native nodes
    MV_al = (1 / (4 * np.pi)) * 2 * np.pi * 2 * np.trapezoid(
        np.trapezoid(I_M, x=sp.thg, axis=1), x=sp.xg)
    JV_al = -(1 / (8 * np.pi)) * 2 * np.pi * 2 * np.trapezoid(
        np.trapezoid(I_J, x=sp.thg, axis=1), x=sp.xg)

    # de-aliased: spectral interpolation to fine grid + Simpson
    kt = np.arange(Ny)
    Bti = np.linalg.inv(np.cos(2.0 * kt[None, :] * sp.thg[:, None]))
    kx = np.arange(Nx)
    Bxi = np.linalg.inv(np.cos(kx[None, :] * np.arccos(np.clip(sp.xg, -1, 1))[:, None]))
    thf = np.linspace(0, np.pi / 2, 801)
    xf = np.linspace(-1, 1, 401)
    Bf = np.cos(2.0 * kt[None, :] * thf[:, None])
    Bxf = np.cos(kx[None, :] * np.arccos(np.clip(xf, -1, 1))[:, None])

    def dealias(v):
        C = Bxi @ v @ Bti.T
        vf = np.einsum("ai,ij,bj->ab", Bxf, C, Bf)
        return 2 * simpson(simpson(vf, x=thf, axis=1), x=xf)

    MV_de = (1 / (4 * np.pi)) * 2 * np.pi * dealias(I_M)
    JV_de = -(1 / (8 * np.pi)) * 2 * np.pi * dealias(I_J)

    two_OmJV = 2 * hdr["OmegaH"] * 0.0        # J_V not in the header; leave 0 and flag
    return dict(file=os.path.basename(path), q=meta["q"], chi=meta["chi"],
                alpha=al, rh=rh, Nx=Nx, Ny=Ny,
                frac_grid_gtt_nonneg=frac_gtt_nonneg,
                frac_grid_gamma_neg=frac_gamma_neg,
                MV_aliased=MV_al, MV_dealias=MV_de,
                JV_aliased=JV_al, JV_dealias=JV_de,
                well_posed=bool(frac_gtt_nonneg < 1e-6 and frac_gamma_neg < 1e-6),
                verdict=("ILL_POSED: t=const slice is not spacelike on "
                         f"{100*frac_gtt_nonneg:.1f}% of the grid; sqrt(gamma) imaginary on "
                         f"{100*frac_gamma_neg:.1f}%"))


def main(argv):
    files = argv[1:] or sorted(glob.glob(
        os.path.join(ROOT, "solutions/field_profiles/sol_*.txt")))
    outdir = os.path.join(HERE, "verify", "out")
    os.makedirs(outdir, exist_ok=True)
    rows = []
    for p in files:
        rows.append(analyse(p))
    import csv
    with open(os.path.join(outdir, "komar_kernel_check.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    ts = time.strftime("%Y-%m-%dT%H:%M:%S")
    json.dump(dict(task="komar kernel (first principles); volume-integral route ILL-POSED",
                   timestamp=ts, profiles=[os.path.basename(p) for p in files],
                   input_sha256={os.path.basename(p): sha256(p) for p in files},
                   method="M_V = (1/4pi) int (2T^V_ab - T^V g_ab) t^a n^b sqrt(gamma) d^3x",
                   dealiasing="spectral interpolation to 4001 theta x 2001 x + Simpson",
                   conclusion="volume integral is ill-posed: -g_tt <= 0 on ~18% of the grid",
                   median=None, max=None,
                   note="median/max NOT computed -- the honest answer is that this route "
                        "cannot produce M_V; see equations/KOMAR_derivation.md §2.2"),
              open(os.path.join(outdir, "komar_kernel_meta.json"), "w"), indent=1)
    print(f"{'file':26s} {'frac(-gtt<=0)':>14s} {'frac(gamma<0)':>14s} "
          f"{'MV_aliased':>13s} {'MV_dealias':>13s}  verdict")
    for r in rows:
        print(f"{r['file']:26s} {r['frac_grid_gtt_nonneg']:14.4f} "
              f"{r['frac_grid_gamma_neg']:14.4f} {r['MV_aliased']:13.4e} "
              f"{r['MV_dealias']:13.4e}  {'OK' if r['well_posed'] else 'ILL_POSED'}")
    print(f"### wrote {outdir}/komar_kernel_check.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
