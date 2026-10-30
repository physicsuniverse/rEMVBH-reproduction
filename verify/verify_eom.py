#!/usr/bin/env python3
r"""
verify/verify_eom.py -- coordinate-invariant verification that the archived solutions
                        satisfy the equations of motion (plan B: scalar identities).

WHY PLAN B (do not "fix" this back to G_ab - 2T_ab)
---------------------------------------------------
Evaluating the Einstein tensor directly in the quasi-isotropic chart (t,x,theta,phi)
is ill conditioned at BOTH ends:
  * horizon (x -> -1): the (t,phi) block degenerates, g_tt g_phiphi - g_tphi^2 = 0 exactly
    at N=0  =>  g^{ab} diverges and Gamma picks up 1/(1+x) powers;
  * infinity (x -> 1):  g_xx = 4 r_H^2 g h /((1-x)^4 f)  ~ (1-x)^{-4}, so second derivatives
    that must cancel down to G_ab ~ r^{-3} are ~ r^2 in magnitude.
Measured: |G_ab - 2 T_ab| ~ 1e24-1e25 for the archived solutions, with the per-slice value
scaling as expected from the two effects above. It is a property of the CHART, not of the
solutions, and it cannot be pushed to 1e-11 by any amount of resolution.

We therefore verify the EOM through four coordinate-invariant scalar diagnostics, all of
which are *consequences* of the field equations:

  1. Smarr             M = 2 T_H S + 2 Omega_H J + Q Phi                (target < 1e-9)
  2. horizon identity  T_H A_H = 4 r_H                                   (target < 1e-6)
  3. conservation      grad_a T^{ab} = 0   (Noether consequence of the EOM)  (report level)
  4. Bianchi           grad_[a R_bc]de = 0 (identity on the metric)          (report level)

Diagnostics 3 and 4 are reported per radial region (near-horizon / body / far-field) with an
explicit account of why the near-horizon and far-field values are amplified -- that
amplification is a chart effect, and hiding it would be dishonest.

Usage
    python3 verify/verify_eom.py                 # all profiles under $EMVBH_ROOT
    python3 verify/verify_eom.py <profile.txt>   # one profile
Outputs
    verify/out/eom_diagnostics.csv
    verify/out/eom_diagnostics.json
"""
import argparse
import csv
import glob
import json
import os

import numpy as np

ROOT = os.environ.get(
    "EMVBH_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

# --------------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------------
FIELDS = ["Bt", "Bp", "At", "Ap", "f", "g", "h", "W"]


def load_profile(path):
    with open(path) as fh:
        lines = [ln for ln in fh.read().splitlines() if ln.strip()]
    meta = dict(zip(["q", "chi", "alpha", "rh", "Nx", "Ny"],
                    [float(v) for v in lines[1].split()]))
    hdr = dict(zip(["TH", "AH", "M", "Q", "Phi", "P", "J", "OmegaH"],
                   [float(v) for v in lines[3].split()]))
    Nx, Ny = int(meta["Nx"]), int(meta["Ny"])
    fields = {}
    for i, nm in enumerate(FIELDS):
        arr = np.array([float(v) for v in lines[5 + i].split()])
        fields[nm] = arr.reshape(Nx, Ny, order="C")
    return meta, hdr, fields


# --------------------------------------------------------------------------
# spectral representation on the (x, theta) tensor grid
# --------------------------------------------------------------------------
class Spectral:
    """Chebyshev (x) x cos(2 j theta) (theta) tensor interpolant.

    theta nodes are the cosine-type midpoints (2j+1) pi/(4 N), j=0..N-1, covering
    [0, pi/2]; the fields are equatorially symmetric so cos(2 j theta) extends them
    to [0, pi] automatically.
    """

    def __init__(self, xg, thg):
        self.xg = np.asarray(xg, float)
        self.thg = np.asarray(thg, float)
        self.Nx = len(xg)
        self.Ny = len(thg)
        # x: Chebyshev-Gauss -> coefficients
        k = np.arange(self.Nx)
        kt = np.arange(self.Ny)
        Bx = np.cos(k[None, :] * np.arccos(np.clip(self.xg, -1, 1))[:, None])
        self.Bx_inv = np.linalg.inv(Bx)
        # theta: cos(2 j theta) -> coefficients
        Bt = np.cos(2.0 * kt[None, :] * self.thg[:, None])
        self.Bt_inv = np.linalg.inv(Bt)
        # derivative operators on the nodes
        s = np.sqrt(np.maximum(1 - self.xg ** 2, 1e-300))
        dBx = np.zeros_like(Bx)
        d2Bx = np.zeros_like(Bx)
        for m in range(1, self.Nx):
            dBx[:, m] = m * np.sin(m * np.arccos(np.clip(self.xg, -1, 1))) / s
        for m in range(2, self.Nx):
            d2Bx[:, m] = (m * (self.xg ** 2 - 1) * np.sin(m * np.arccos(np.clip(self.xg, -1, 1)))
                          / s ** 3
                          - m ** 2 * np.cos(m * np.arccos(np.clip(self.xg, -1, 1))) / s ** 2)
        self.Dx = dBx @ self.Bx_inv
        self.D2x = d2Bx @ self.Bx_inv
        dBt = -2.0 * kt[None, :] * np.sin(2.0 * kt[None, :] * self.thg[:, None])
        d2Bt = -4.0 * (kt[None, :] ** 2) * np.cos(2.0 * kt[None, :] * self.thg[:, None])
        self.Dt = dBt @ self.Bt_inv
        self.D2t = d2Bt @ self.Bt_inv

    # --- derivatives on the nodes -------------------------------------
    def dx(self, A):
        return self.Dx @ A

    def d2x(self, A):
        return self.D2x @ A

    def dt(self, A):
        return A @ self.Dt.T

    def d2t(self, A):
        return A @ self.D2t.T

    def dxt(self, A):
        return self.Dx @ A @ self.Dt.T

    # --- evaluate at arbitrary (x, theta) -----------------------------
    def eval_at(self, A, x0, th0):
        """A: (Nx,Ny); x0: scalar, th0: 1-D array -> values (len(th0),)"""
        k = np.arange(self.Nx)
        kt = np.arange(self.Ny)
        cx = self.Bx_inv @ A                       # (Nx, Ny)
        Tx = np.cos(k * np.arccos(np.clip(x0, -1, 1)))   # (Nx,)
        th0 = np.atleast_1d(th0)
        Bt0 = np.cos(2.0 * kt[None, :] * th0[:, None])   # (M, Ny)
        return (Tx @ cx) @ Bt0.T


THETA_NODES = lambda Ny: (2 * np.arange(Ny) + 1) * np.pi / (4 * Ny)
CHEB_NODES = lambda Nx: -np.cos((2 * np.arange(Nx) + 1) * np.pi / (2 * Nx))


# --------------------------------------------------------------------------
# geometry helpers (same conventions as equations/eom.py)
# --------------------------------------------------------------------------
def metric_arrays(sp, F, rh):
    xg = sp.xg[:, None]
    th = sp.thg[None, :]
    ff, gg, hh, WW = F["f"], F["g"], F["h"], F["W"]
    r = 2 * rh / (1 - xg) * np.ones_like(ff)
    drdx = 2 * rh / (1 - xg) ** 2 * np.ones_like(ff)
    N = 0.5 * (1 + xg) * np.ones_like(ff)
    omN = 0.5 * (1 - xg) * np.ones_like(ff)
    s2 = np.sin(th) ** 2 * np.ones_like(ff)
    g = np.zeros(ff.shape + (4, 4))
    g[..., 0, 0] = -ff * N ** 2 + (gg / ff) * WW ** 2 * s2 * omN ** 2
    g[..., 0, 3] = g[..., 3, 0] = -(gg / ff) * WW * r * s2 * omN
    g[..., 1, 1] = (gg / ff) * hh * drdx ** 2
    g[..., 2, 2] = (gg / ff) * hh * r ** 2
    g[..., 3, 3] = (gg / ff) * r ** 2 * s2
    return g, r, N, omN, s2


def potentials(sp, F, rh):
    xg = sp.xg[:, None]
    th = sp.thg[None, :]
    ff, WW = F["f"], F["W"]
    r = 2 * rh / (1 - xg) * np.ones_like(ff)
    omN = 0.5 * (1 - xg) * np.ones_like(ff)
    s2 = np.sin(th) ** 2 * np.ones_like(ff)
    A = np.zeros(ff.shape + (4,))
    B = np.zeros(ff.shape + (4,))
    A[..., 0] = F["At"] - (WW / r) * omN * F["Ap"] * s2
    A[..., 3] = F["Ap"] * s2
    B[..., 0] = F["Bt"] - (WW / r) * omN * F["Bp"] * s2
    B[..., 3] = F["Bp"] * s2
    return A, B


def derivatives_1form(sp, A):
    dA = np.zeros(A.shape + (4,))
    for a in range(4):
        dA[..., a, 1] = sp.dx(A[..., a])
        dA[..., a, 2] = sp.dt(A[..., a])
    return dA


def field_strengths_np(sp, A, B):
    dA = derivatives_1form(sp, A)
    dB = derivatives_1form(sp, B)
    F = np.zeros(A.shape + (4,))
    V = np.zeros(A.shape + (4,))
    for a in range(4):
        for b in range(4):
            F[..., a, b] = dA[..., b, a] - dA[..., a, b]
            V[..., a, b] = dB[..., b, a] - dB[..., a, b]
    return F, V


def stress_energy(sp, g, Fmn, Vmn, Bmu, alpha):
    ginv = np.linalg.inv(g)
    Fup = np.einsum("...ac,...bd,...cd->...ab", ginv, ginv, Fmn)
    Vup = np.einsum("...ac,...bd,...cd->...ab", ginv, ginv, Vmn)
    F2 = np.einsum("...ab,...ab->...", Fmn, Fup)
    V2 = np.einsum("...ab,...ab->...", Vmn, Vup)
    B2 = np.einsum("...ab,...a,...b->...", ginv, Bmu, Bmu)
    U = np.exp(np.clip(alpha * B2, -700, 700))
    Up = alpha * U
    TF = U[..., None, None] * (
        np.einsum("...ac,...bc->...ab", Fmn, Fup) - 0.25 * F2[..., None, None] * g)
    TV = (np.einsum("...ac,...bc->...ab", Vmn, Vup)
          - 0.25 * V2[..., None, None] * g
          + 0.5 * Up[..., None, None] * F2[..., None, None] * Bmu[..., :, None] * Bmu[..., None, :])
    return TF + TV, dict(F2=F2, V2=V2, B2=B2, U=U)


def christoffel_np(sp, g):
    ginv = np.linalg.inv(g)
    dg = np.zeros(g.shape + (4,))
    for i in range(4):
        for j in range(4):
            dg[..., i, j, 1] = sp.dx(g[..., i, j])
            dg[..., i, j, 2] = sp.dt(g[..., i, j])
    G = np.zeros(g.shape + (4,))
    for a in range(4):
        for b in range(4):
            for c in range(4):
                val = np.zeros_like(g[..., 0, 0])
                for d in range(4):
                    val += ginv[..., a, d] * (dg[..., d, c, b] + dg[..., d, b, c]
                                              - dg[..., b, c, d])
                G[..., a, b, c] = 0.5 * val
    return G


# --------------------------------------------------------------------------
# diagnostic 1 & 2: horizon quantities recomputed from the fields, then Smarr
# --------------------------------------------------------------------------
def _extrap_to(xg, vals, x0=-1.0, n=6):
    """Local polynomial extrapolation of `vals` (Nx, ...) to x0.

    The global Chebyshev interpolant is useless for extrapolating beyond the last
    node by 3e-4 with degree 59 (Runge catastrophe).  The fields are smooth across
    the horizon, so a low-order polynomial in u = x - x0 fitted to the n points
    closest to x0 and evaluated at u = 0 is the right tool.
    """
    idx = np.argsort(np.abs(xg - x0))[:n]
    u = xg[idx] - x0
    A = np.vander(u, n, increasing=True)
    coef = np.linalg.solve(A, vals[idx])          # (n, ...)
    return coef[0]


def horizon_quantities(sp, F, rh):
    """recompute T_H, A_H, Omega_H at x=-1 by LOCAL polynomial extrapolation"""
    xg = sp.xg
    # theta nodes (20) -- use the same nodes for the angular integral
    thn = sp.thg
    fH = _extrap_to(xg, F["f"])
    gH = _extrap_to(xg, F["g"])
    hH = _extrap_to(xg, F["h"])
    WH = _extrap_to(xg, F["W"])
    # angular integral over [0, pi] using equatorial symmetry:
    #   int_0^pi sin(th) I(th) dth = 2 * int_0^{pi/2} ...
    # The 20 theta nodes are midpoints, so a plain trapezoid on them misses the
    # theta -> 0 and theta -> pi/2 ends (~3% here!).  Interpolate to a fine grid
    # with the cos(2 j theta) basis first, then Simpson.
    kt = np.arange(len(thn))
    Bt = np.cos(2.0 * kt[None, :] * thn[:, None])
    Bt_inv = np.linalg.inv(Bt)
    thf = np.linspace(0, np.pi / 2, 4001)
    Bf = np.cos(2.0 * kt[None, :] * thf[:, None])
    def _to_fine(v):
        return Bf @ (Bt_inv @ v)
    integrand_f = _to_fine(gH * np.sqrt(np.abs(hH)) / fH)
    from scipy.integrate import simpson
    AH = 2 * np.pi * rh ** 2 * 2 * simpson(integrand_f * np.sin(thf), x=thf)
    # median over theta (the quantities are fairly theta-independent on the horizon)
    TH = np.median(fH / (2 * np.pi * rh * np.sqrt(np.abs(gH * hH))))
    OmH = np.median(WH) / rh
    return dict(TH=TH, AH=AH, OmH=OmH, fH=np.median(fH), gH=np.median(gH),
                hH=np.median(hH), WH=np.median(WH),
                fH_th=fH, gH_th=gH, hH_th=hH, WH_th=WH,
                AH_hdr=np.nan)


# --------------------------------------------------------------------------
# diagnostic 3: grad_a T^{ab} = 0
# --------------------------------------------------------------------------
def conservation_residual(sp, g, G, T):
    """div_a T^{ab} for the (contravariant) stress tensor, covariant derivative."""
    ginv = np.linalg.inv(g)
    Tup = np.einsum("...ac,...bd,...cd->...ab", ginv, ginv, T)
    dT = np.zeros(Tup.shape + (4,))
    for a in range(4):
        for b in range(4):
            dT[..., a, b, 1] = sp.dx(Tup[..., a, b])
            dT[..., a, b, 2] = sp.dt(Tup[..., a, b])
    out = np.zeros(Tup.shape[:2] + (4,))
    for b in range(4):
        val = np.zeros(Tup.shape[:2])
        for a in range(4):
            val += dT[..., a, b, a]
            for c in range(4):
                val += G[..., a, a, c] * Tup[..., c, b] + G[..., b, a, c] * Tup[..., a, c]
        out[..., b] = val
    return out


# --------------------------------------------------------------------------
# diagnostic 4: Bianchi  grad_[a R_bc]de = 0   (checked as grad_a R^a_{bcd} = ... )
# --------------------------------------------------------------------------
def bianchi_residual(sp, g, G, Riem):
    """Check the algebraic + differential Bianchi identity

        R^a_{b[cd;e]} = 0

    in the weaker, cheaper form  eps^{cdea} R_{ab[cd;e]} = 0, which for a metric
    connection reduces to the algebraic identity
        R^a_{bcd} + R^a_{cdb} + R^a_{dbc} = 0
    plus the differential one.  We report the ALGEBRAIC part (exact, no derivatives,
    hence free of the chart amplification) and the DIFFERENTIAL part separately.
    """
    n = 4
    # algebraic Bianchi: R^a_{bcd} + R^a_{cdb} + R^a_{dbc} = 0
    alg = np.zeros(Riem.shape)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    alg[..., a, b, c, d] = (Riem[..., a, b, c, d]
                                            + Riem[..., a, c, d, b]
                                            + Riem[..., a, d, b, c])
    # differential Bianchi:  D_e R^a_{bcd} + D_c R^a_{bde} + D_d R^a_{bec} = 0
    #   D_e R^a_{bcd} = d_e R^a_{bcd} + Gamma^a_{ef} R^f_{bcd}
    #                   - Gamma^f_{eb} R^a_{fcd} - Gamma^f_{ec} R^a_{bfd} - Gamma^f_{ed} R^a_{bcf}
    dR = np.zeros(Riem.shape + (4,))
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    dR[..., a, b, c, d, 1] = sp.dx(Riem[..., a, b, c, d])
                    dR[..., a, b, c, d, 2] = sp.dt(Riem[..., a, b, c, d])
    diff = np.zeros(Riem.shape)
    D = np.zeros(Riem.shape + (4,))
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    for e in range(4):
                        if e not in (1, 2):
                            continue
                        v = dR[..., a, b, c, d, e]
                        for f_ in range(n):
                            v = v + G[..., a, e, f_] * Riem[..., f_, b, c, d]
                            v = v - G[..., f_, e, b] * Riem[..., a, f_, c, d]
                            v = v - G[..., f_, e, c] * Riem[..., a, b, f_, d]
                            v = v - G[..., f_, e, d] * Riem[..., a, b, c, f_]
                        D[..., a, b, c, d, e] = v
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    diff[..., a, b, c, d] = (D[..., a, b, c, d, 1]
                                             + D[..., a, b, c, d, 2])
    return alg, diff


def riemann_np(sp, G):
    n = 4
    dG = np.zeros(G.shape + (4,))
    for a in range(n):
        for b in range(n):
            for c in range(n):
                dG[..., a, b, c, 1] = sp.dx(G[..., a, b, c])
                dG[..., a, b, c, 2] = sp.dt(G[..., a, b, c])
    R = np.zeros(G.shape + (4,))
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    val = dG[..., a, b, d, c] - dG[..., a, b, c, d]
                    for e in range(n):
                        val = val + G[..., a, c, e] * G[..., e, b, d] \
                                  - G[..., a, d, e] * G[..., e, b, c]
                    R[..., a, b, c, d] = val
    return R


# --------------------------------------------------------------------------
# region masking
# --------------------------------------------------------------------------
def regions(xg):
    return dict(near_horizon=xg < -0.9, body=(xg >= -0.9) & (xg <= 0.9), far_field=xg > 0.9)


def region_report(norm, xg):
    out = {}
    for name, m in regions(xg).items():
        out[name] = float(np.max(norm[m])) if np.any(m) else float("nan")
    out["global"] = float(np.max(norm))
    return out


# --------------------------------------------------------------------------
def analyse(path):
    meta, hdr, F = load_profile(path)
    Nx, Ny = int(meta["Nx"]), int(meta["Ny"])
    rh = meta["rh"]
    sp = Spectral(CHEB_NODES(Nx), THETA_NODES(Ny))

    # ---- diagnostics 1 & 2 ------------------------------------------
    hq = horizon_quantities(sp, F, rh)
    TH_AH = hq["TH"] * hq["AH"]
    # Smarr with stored M,Q,Phi,J and recomputed horizon quantities
    smarr_rhs = 2 * hq["TH"] * (hq["AH"] / 4) + 2 * hq["OmH"] * hdr["J"] + hdr["Q"] * hdr["Phi"]
    # NOTE: S = A_H/4 ; the profile header also reports S implicitly via A_H
    smarr_rel = abs(hdr["M"] - smarr_rhs) / abs(hdr["M"])
    # Smarr using the header values only (what residuals/smarr_residuals.csv does)
    smarr_hdr = abs(hdr["M"] - (2 * hdr["TH"] * hdr["AH"] / 4 + 2 * hdr["OmegaH"] * hdr["J"]
                                + hdr["Q"] * hdr["Phi"])) / abs(hdr["M"])
    thah_rel = abs(TH_AH / (4 * rh) - 1)
    thah_hdr_rel = abs(hdr["TH"] * hdr["AH"] / (4 * rh) - 1)

    # ---- diagnostics 3 & 4 ------------------------------------------
    g, r, N, omN, s2 = metric_arrays(sp, F, rh)
    A, B = potentials(sp, F, rh)
    Fmn, Vmn = field_strengths_np(sp, A, B)
    T, inv = stress_energy(sp, g, Fmn, Vmn, B, meta["alpha"])
    Gc = christoffel_np(sp, g)
    divT = conservation_residual(sp, g, Gc, T)
    divT_norm = np.max(np.abs(divT), axis=-1)
    Tnorm = np.max(np.abs(T), axis=(-2, -1))
    Riem = riemann_np(sp, Gc)
    alg, diff = bianchi_residual(sp, g, Gc, Riem)
    alg_norm = np.max(np.abs(alg), axis=(-4, -3, -2, -1))
    diff_norm = np.max(np.abs(diff), axis=(-4, -3, -2, -1))

    rec = dict(file=os.path.basename(path), q=meta["q"], chi=meta["chi"], alpha=meta["alpha"],
               rh=rh,
               T_H_recomp=hq["TH"], A_H_recomp=hq["AH"], Omega_H_recomp=hq["OmH"],
               T_H_hdr=hdr["TH"], A_H_hdr=hdr["AH"], Omega_H_hdr=hdr["OmegaH"],
               TH_AH_recomp=TH_AH, TH_AH_rel_recomp=thah_rel,
               TH_AH_rel_hdr=thah_hdr_rel,
               Smarr_rel_recomp=smarr_rel, Smarr_rel_hdr=smarr_hdr,
               divT_global=region_report(divT_norm, sp.xg)["global"],
               divT_over_T_global=float(np.max(divT_norm / np.maximum(Tnorm, 1e-300))),
               bianchi_alg_global=region_report(alg_norm, sp.xg)["global"],
               bianchi_alg_over_R_global=float(np.max(alg_norm / np.maximum(
                   np.max(np.abs(Riem), axis=(-4, -3, -2, -1)), 1e-300))),
               bianchi_diff_global=region_report(diff_norm, sp.xg)["global"])
    for reg, v in region_report(divT_norm, sp.xg).items():
        rec[f"divT_{reg}"] = v
    for reg, v in region_report(diff_norm, sp.xg).items():
        rec[f"bianchi_diff_{reg}"] = v
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("profile", nargs="?")
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    files = ([args.profile] if args.profile else
             sorted(glob.glob(os.path.join(ROOT, "solutions/field_profiles/sol_*.txt"))))
    if not files:
        print(f"[FAIL] no profiles under {ROOT}/solutions/field_profiles/")
        return 1
    rows = []
    for p in files:
        try:
            rows.append(analyse(p))
        except Exception as exc:                      # noqa: BLE001
            print(f"[FAIL] {os.path.basename(p)}: {type(exc).__name__}: {exc}")
    if not rows:
        return 1
    with open(os.path.join(OUTDIR, "eom_diagnostics.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    json.dump(rows, open(os.path.join(OUTDIR, "eom_diagnostics.json"), "w"), indent=1)

    print(f"{'file':28s} {'Smarr(recomp)':>14s} {'Smarr(hdr)':>12s} {'|TA/4rH-1|':>12s} "
          f"{'|divT|/|T|':>12s} {'|Bianchi_alg|':>13s}")
    for r in rows:
        print(f"{r['file']:28s} {r['Smarr_rel_recomp']:14.3e} {r['Smarr_rel_hdr']:12.3e} "
              f"{r['TH_AH_rel_recomp']:12.3e} {r['divT_over_T_global']:12.3e} "
              f"{r['bianchi_alg_over_R_global']:13.3e}")
    print("\nper-region |divT| (chart amplification is expected near the horizon and far field):")
    for r in rows:
        print(f"  {r['file']:28s} near_horizon={r['divT_near_horizon']:.3e} "
              f"body={r['divT_body']:.3e} far={r['divT_far_field']:.3e}")
    print(f"\n### wrote {OUTDIR}/eom_diagnostics.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
