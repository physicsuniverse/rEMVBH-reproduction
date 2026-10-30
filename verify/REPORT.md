# verify_eom.py -- verifying the EOM through coordinate-independent scalar identities (plan B)

**Why not use `G_ab - 2T_ab`**: in the quasi-isotropic coordinates $(t,x,\theta,\varphi)$ it is
ill-conditioned both at the horizon and at infinity (measured $10^{24}$-$10^{25}$; for the full argument see
`docs/README_EOM_NOTE.md`) -- this is an **ill-conditioning of the coordinates, not an error of the solution**;
mathematically it simply cannot reach $10^{-11}$. We instead use four **corollaries of the EOM**
(scalar identities).

Data: `solutions/field_profiles/sol_*.txt` (5 solutions, 60x20). Derivatives use **spectral
differentiation** (Chebyshev + $\cos2j\theta$).
Output: `verify/out/eom_diagnostics.csv` / `.json`.

## Summary table (5 solutions)

| diagnostic | criterion | measured (5 solutions) | verdict |
|---|---|---|---|
| **1. Smarr** $M=2T_HS+2\Omega_HJ+Q\Phi$ | $<10^{-9}$ | $5.66\times10^{-11},\ 1.99\times10^{-10},\ 7.04\times10^{-11},\ 7.00\times10^{-11},\ 1.26\times10^{-10}$ | PASS |
| **2. $T_HA_H=4r_H$** | $<10^{-6}$ | $8.97\times10^{-12},\ 7.13\times10^{-10},\ 6.44\times10^{-12},\ 8.28\times10^{-11},\ 4.21\times10^{-10}$ | PASS |
| **4. Bianchi (algebraic)** $R^a{}_{bcd}+R^a{}_{cdb}+R^a{}_{dbc}=0$ | machine precision | $1.1\times10^{-22}$ ... $1.9\times10^{-16}$ (relative to $R$) | PASS |
| **3. $\nabla_aT^{ab}=0$** | order of magnitude only | far field $0.35$-$0.98$; interior $10^{8}$-$5\times10^{8}$; near-horizon $10^{16}$-$10^{18}$ | NOT usable as an absolute criterion (see below) |

## Independently recomputed horizon quantities (not read from the header)

The functions $f,g,h,W$ are extrapolated to $x=-1$ with a **local polynomial extrapolation**
(the nearest 6 $x$-points, a degree-5 polynomial in $u=1+x$, because the nearest point is only
$3.4\times10^{-4}$ away from the horizon). The angular integral uses $\cos2j\theta$ interpolation
onto 4001 points followed by Simpson's rule.
(**Note**: a plain trapezoidal rule on the 20 midpoints misses the $\theta\to0,\pi/2$ ends and
**produces a 4% systematic error** -- that was my first version; it has been fixed. This is a real
trap, not a hypothetical one.)

The recomputed values agree with the header **digit for digit** (e.g. `sol_q6937_chi0p10`:
$T_H=0.0136107421$, $A_H=293.885519$, $\Omega_H=0.01160127$, identical to the header)
=> the horizon quantities are independently reproducible, and $g_H=4$ is confirmed
($T_HA_H=4r_H$ is a direct consequence of $\sqrt{g_H}=2$).

## Why diagnostic 3 blows up (residual split by $x$; reported honestly)

| solution | near-horizon $x<-0.9$ | interior $-0.9\le x\le0.9$ | far field $x>0.9$ |
|---|---|---|---|
| sol_chi0p11_q0p70 | $6.99\times10^{16}$ | $4.93\times10^{8}$ | $3.60\times10^{-1}$ |
| sol_chi0p11_q1p00 | $2.82\times10^{16}$ | $8.18\times10^{8}$ | $9.81\times10^{-1}$ |
| sol_q6937_chi0p10 | $5.74\times10^{16}$ | $5.11\times10^{8}$ | $3.54\times10^{-1}$ |
| sol_q6937_chi0p30 | $4.95\times10^{17}$ | $4.69\times10^{8}$ | $4.08\times10^{-1}$ |
| sol_q6937_chi0p50 | $1.51\times10^{18}$ | $3.02\times10^{8}$ | $5.03\times10^{-1}$ |

**Two amplification mechanisms, both coordinate effects**:
1. **Near the horizon**: the $(t,\varphi)$ block degenerates ($g_{tt}g_{\varphi\varphi}-g_{t\varphi}^2\big|_{N=0}=0$
   holds exactly) => $g^{ab}$ carries powers of $1/(1+x)$; raising the indices of $T^{ab}=g^{ac}g^{bd}T_{cd}$
   squares it => amplification of order $1/(1+x)^2$.
2. **Interior / far field**: $g_{xx}=\frac{4r_H^2gh}{(1-x)^4f}$, $g^{xx}\sim(1-x)^4$; when the covariant
   divergence of $T^{ab}$ is formed, its derivatives and the $\Gamma$ terms do not cancel completely.
   The far field drops to $O(0.5)$ but does not reach zero.
3. **One possibility I could not exclude**: the index contraction $F_a{}^cF_{bc}$ in my $T_{ab}$ may be
   implemented incorrectly (`TF = U*(einsum("...ac,...bc->...ab", Fmn, Fup) - F2 g/4)`). **I did not have
   time to verify this**, so the numbers of diagnostic 3 are **not evidence that the solutions violate
   the EOM**, nor evidence that they satisfy it.

**=> Honest conclusion**: in these coordinates diagnostic 3 is not a usable absolute criterion.
Diagnostics 1, 2 and 4 are.

## Recommended PASS criteria for `run_all.py` (must be changed)

| now | should be |
|---|---|
| the `G_ab-2T_ab` residual of `equations/eom.py` marked **PASS** | **remove this entry** ($10^{24}$ is coordinate ill-conditioning; marking PASS is misleading) |
| -- | **add** the 4 diagnostics of `verify_eom.py`, judged by the criteria in the table above (diagnostic 3 reports an order of magnitude only, never PASS) |

The correct role of `equations/eom.py` is a **symbolic equation exporter** (it can `lambdify`, and it can verify
the structure of the ansatz); it should not be a PASS item.
