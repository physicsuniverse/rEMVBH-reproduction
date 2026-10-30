# EMVBH equations → pure Python (sympy/numpy)

**Author**: the numerical-solver analysis|**Date**: 2026-10-14
**Target**: copy this directory into `<reproduction-package>/equations/`
**Status**: 3 of 4 scripts **complete, runnable, and fully passing**; `komar_kernel.py` **not delivered** (see §4).

---

## 1. What is here

| File | Contents | Checks |
|---|---|---|
| `emv_action.py` | Second-order expansion of the EMV action, the mass matrix $M^{ab}$, and the ghost/gradient criteria $X_1,X_2$ | **18/18 PASS** |
| `mv_identity.py` | Parity ⇒ $c_1=0$; $M_V=2\Omega_HJ_V$ ⇒ $c_2=2\omega_1j_1$; numeric cross-check | **9/9 PASS** |
| `ghost_dispersion.py` | Full dispersion relation with the $\omega$-linear term; QEP solver (numpy only) | **16/16 PASS** |
| `*_results.json` | Machine-readable PASS/FAIL record of every check | — |

Run each with plain `python3`:
```bash
python3 equations/emv_action.py          # or --json
python3 equations/mv_identity.py         # or --json [--data FILE]
python3 equations/ghost_dispersion.py    # or --json
```
**Dependencies**: `sympy`, `numpy` only. **No Mathematica. No scipy** (the QEP is solved by a
companion linearisation + `numpy.linalg.eigvals`).

---

## 2. What the scripts establish

### `emv_action.py` (source: `EMV_second_order_action.md`)
1. $U(X)=e^{\alpha X}$ ⇒ $U'=\alpha U$, $U''=\alpha^2U$;
2. the **linear** term in $\delta B$ cancels **iff** the background EOM
   $\nabla_aV^{ab}=\tfrac12U'F^2B^b$ is used;
3. the quadratic action gives $M^{ab}=m_{\rm eff}^2(g^{ab}+2\alpha\bar B^a\bar B^b)$ with
   $m_{\rm eff}^2=\tfrac12\alpha UF^2$ — verified by reproducing
   $-\tfrac12M^{ab}\delta B_{{\rm can},a}\delta B_{{\rm can},b}$ with the canonical rescaling
   $\delta B_{\rm can}=2\delta B$;
4. eigenvalue structure: $m_{\rm eff}^2(1+2\alpha\bar B^2)$ along $\bar B$, $m_{\rm eff}^2$ (3-fold)
   transverse, and $\det M=m_{\rm eff}^4\det g\,(1+2\alpha\bar B^2)$;
5. **the known conclusion is confirmed**: $\operatorname{sign}(X_2)=\operatorname{sign}(m_{\rm eff}^2)
   =\operatorname{sign}(\alpha F^2)$, so the tachyonic region and the no-ghost condition are mutually
   exclusive; limits $X_1=X_2=m_{\rm eff}^2/2$ reproduce both the spherical no-go and the Proca values.

> **Index bug worth recording**: my first version mixed covariant and contravariant components of
> $\bar B^a$ / $\delta B_a$ and five checks failed. Fixed by keeping $\bar B_a$, $\delta B_a$
> **covariant** and $G=g^{ab}$ as the inverse metric. The failing checks were exactly the ones that
> depend on index placement, so the test suite caught it.

### `mv_identity.py` (source: `MV_analytic_derivation.md`)
Symbolic: parity $\Rightarrow c_1=c_3=0$; spherical theorem $\Rightarrow c_0=0$; hence
$M_V=c_2\chi^2+O(\chi^4)$; and $2\Omega_HJ_V=2\omega_1j_1\chi^2+\ldots$ so **$c_2=2\omega_1j_1$ exactly**,
with the linear term vanishing automatically (odd×odd).

Numeric: fits $\omega_1$ from $\Omega_H(\chi)$ and $j_1$ from $J_V(\chi)$ **without touching $M_V$**,
then compares with the $c_2$ fitted from $M_V$.
**Regime-aware**: for the local $\alpha=20$ dataset the script measures the local exponent
($p=2.73$) and correctly *refuses* to assert the leading-order $c_2$ comparison, reporting instead
that the family is outside the $\chi^2$ regime (which is why $\alpha=20$ is excluded from the paper's
$\chi^2$ statement). The identity $M_V=2\Omega_HJ_V$ is still tested and holds to $1.1\times10^{-3}$.

### `ghost_dispersion.py` (source: `ghost_dispersion.md`)
Verifies eqs. (2.1)–(2.4) as matrix identities, and establishes the central structural claim:
- static diagonal background ⇒ $\det\mathcal D$ is **even** in $\omega$ (max odd coefficient $6.9\times10^{-15}$);
- rotating background ($g^{t\varphi}\neq0$) ⇒ $\det\mathcal D$ **has odd powers** of $\omega$;
- non-rotating but $M^{ti}\neq0$ ⇒ **also** has odd powers;
- QEP $A_2\omega^2+A_1\omega+A_0=0$ solved with numpy; every well-scaled root satisfies $\det\mathcal D=0$.

---

## 3. Two errors found in the source `.md` documents (reported, not silently patched)

### 3.1 `ghost_dispersion.md` eq. (2.3) — index placement
The document prints
$$\det\mathcal D=\det(k^2G+M)\left[1-k_a\left(k^2G+M\right)^{-1}_{ab}k_b\right].$$
The matrix determinant lemma $\det(A-uv^T)=\det A\,[1-v^TA^{-1}u]$ with $u=v=k^b$ (contravariant)
requires **contravariant** indices:
$$\det\mathcal D=\det(k^2G+M)\left[1-k^b\left(k^2G+M\right)^{-1}_{bd}k^d\right].$$
**Numeric evidence** (12 samples, 4 metric/mass combinations): contravariant version agrees to
$4.7\times10^{-16}$; the printed covariant version is off by up to **1.98×10⁻¹**. The document's own
words ("take $u=k$ (contravariant components)") are correct — only the displayed indices are wrong.

### 3.2 `ghost_dispersion.md` §2.3 — time reversal omits $k_r\to-k_r$
The document says $\mathcal D$ is invariant under "$\omega\to-\omega,\ \ell\to-\ell$". At the **matrix**
level that is false: the invariance is under $k_a\to-k_a$, i.e. $(\omega,k_r,\ell)\to(-\omega,-k_r,-\ell)$.
Verified: $D(-\omega;-k,-\ell)-D(\omega;k,\ell)=0$ exactly, while fixing $k_r$ leaves a residual of
$1.40$. **The document's final root-set conclusion is nevertheless correct**, because for the
representative metrics tested $\det\mathcal D$ also happens to be even in $k_r$
(max coefficient difference $0.0$).

Neither error affects the paper's physics conclusions; both should be corrected in the `.md`.

---

## 4. `komar_kernel.py` — NOT delivered (honest status)

The requested module (compute $M_V,J_V$ from the $B$ fields by direct Komar integration, with a
de-aliased version, so that `verify/verify_komar.py` can independently verify $M_V=2\Omega_HJ_V$)
is **not delivered**. Reasons:

1. **The $M_V$ integrand is not available in usable form.** Its only authoritative source is
   `numerics/EOMEQs/intofMVx.mx`, which is a **binary Wolfram dump** (`(*This is a Wolfram Language
   binary dump file...*)`), not text. `MV_analytic_derivation.md` §2 quotes it only partially, with
   explicit `$+\cdots$` elisions in $\mathcal I_3$ and $\mathcal I_4$ — so transcribing it would mean
   **inventing** the missing terms.
2. Implementing it *from first principles* (vary the EMV action w.r.t. $g^{ab}$ to get $T^V_{ab}$,
   then evaluate the Komar volume integrals on the $60\times20$ grid) is a self-contained
   derivation plus a solver in its own right — it needs its own session to be done **and verified**
   against the `residuals/komar_residuals.csv` anchor, not bolted on at the end of this one.
3. Writing it without being able to test it would violate the project's "do not fabricate" rule.

**Concrete recipe for whoever picks it up** (all inputs exist):
- field profiles: `solutions/field_profiles/sol_*.txt`, layout = header
  `# q chi alpha rh Nx Ny`, then `# TH AH M Q Phi P J OmegaH`, then **8 field lines**
  (order `Bt, Bp, At, Ap, f, g, h, W`), each `Nx*Ny` values, **x-major**;
  grids: `solutions/field_profiles/grids.txt` (line 1 = x grid, line 2 = θ grid, 60 and 20 points);
- $\theta$-basis is `{cos 2jθ}` (20 points ⇒ degree ≤ 19), so de-aliasing = project the
  **nonlinear** integrand onto that basis and re-evaluate exactly (this is the C3b method);
- accept/reject anchor: de-aliased $M_V=2\Omega_HJ_V$ residual median $\sim10^{-9}$, max $5.2\times10^{-7}$;
  on the raw 20-point grid expect $10^{-4}$–$10^{-3}$ systematics.

---

## 5. Installation

This directory was produced in a session whose write scope did **not** include
`<reproduction-package>/`, so it must be copied:
```bash
cp <working-directory>/EMVBH-equations-py/*.py \
   <reproduction-package>/equations/
```
(The `*_results.json` files are regenerated by running each script with `--json`.)
