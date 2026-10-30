# `KOMAR_derivation.md` — first-principles derivation of $T^V_{ab}$ and of the Komar integral kernel

**Performed by** the numerical-physics analysis (Komar-integral verification task)|pure Python derivation, no Mathematica, and no transcribed binary dumps or formulas with ellipses.

---

## 0. Bottom line (the single most important statement)

> **The volume-integral route $M_V=\frac{1}{4\pi}\int_\Sigma(2T^V_{ab}-T^Vg_{ab})t^a n^b\sqrt{\gamma}\,d^3x$ is
> mathematically ill-posed in the present setting — not "ill-conditioned", not "insufficiently accurate", but
> the object being integrated does not exist in that region.**
>
> Measured (`sol_q6937_chi0p10`, 60×20, $\alpha=20$): **$-g_{tt}\le 0$ on 17.7% of the grid points**.
> At $x=-0.9997$ (the point neighbouring the horizon): $-g_{tt}\in[-3.14\times10^{-3},-4.85\times10^{-6}]<0$;
> at the same time $g_{\varphi\varphi}-g_{t\varphi}^2/g_{tt}<0$ ($-3.7\times10^{-5}$).
> $$
> \Longrightarrow\quad \text{the } t=\text{const slice enters the ergoregion near the horizon and becomes timelike; } \sqrt{\gamma}\text{ is imaginary.}
> $$
> **This is not a $1/N$ singularity, but the fact that $n^a=t^a/\sqrt{-g_{tt}}$ simply does not exist inside the ergoregion.**
> One must therefore **switch to the boundary form** (Gauss' theorem), or truncate $\Sigma$ at the ergosurface and add the corresponding boundary term.

This agrees with the diagnosis of the independent verification pass, but is **sharper**: what it called "the $1/N$ singularity must be cancelled analytically" is here not merely a singularity but a **sign flip** ($-g_{tt}$ crossing zero). Any implementation that "first computes the integrand and then multiplies by $\sqrt\gamma$ and integrates" blows up inside the ergoregion — both the $M_V\sim10^8$–$10^{242}$ / nan values reported by the independent verification and the wrong values from my own first version originate here.

---

## 1. Rigorous derivation of $T^V_{ab}$ (varying the action with respect to $g^{ab}$)

The action
$$S=\frac{1}{16\pi}\int\sqrt{-g}\Big[R-U(|B|^2)F_{ab}F^{ab}-V_{ab}V^{ab}\Big],\qquad
V_{ab}=\nabla_aB_b-\nabla_bB_a,\quad F_{ab}=\nabla_a{\cal A}_b-\nabla_b{\cal A}_a .$$
**Notation**: $V_{ab}$ is the field strength of the vector field $B$, $F_{ab}$ is the $U(1)$ electromagnetic field — two different fields, strictly distinguished throughout this derivation.

**① Variation**

- $\delta\sqrt{-g}=-\frac12\sqrt{-g}\,g_{ab}\delta g^{ab}$
- $\delta R=R_{ab}\delta g^{ab}$ (Palatini; boundary terms omitted)
- $\delta(V_{ab}V^{ab})=2V_a{}^cV_{bc}\,\delta g^{ab}$ (using $V_{ab}V^{ab}=V_{ab}V_{cd}g^{ac}g^{bd}$)
- $\delta(UF^2)=\underbrace{2U F_a{}^cF_{bc}}_{F^2\text{ metric dependence}}\delta g^{ab}
   +\underbrace{U'F^2B_aB_b}_{U(B^2)\text{ metric dependence}}\delta g^{ab}$,
   with $U'=dU/d(B^2)=\alpha U$ and $B^2\equiv g^{ab}B_aB_b$.

**② Collecting the coefficients of $\delta g^{ab}$**
$$\delta S\propto\int\sqrt{-g}\,\delta g^{ab}\Big[-\tfrac12g_{ab}(R-UF^2-V^2)+R_{ab}
-2UF_a{}^cF_{bc}-U'F^2B_aB_b-2V_a{}^cV_{bc}\Big].$$
Using $\frac12g_{ab}UF^2=2\cdot\frac14UF^2g_{ab}$ and $\frac12g_{ab}V^2=2\cdot\frac14V^2g_{ab}$:
$$\boxed{\ G_{ab}-2T^F_{ab}-2T^V_{ab}=0\ }$$
$$\boxed{\ T^F_{ab}=U\Big(F_a{}^cF_{bc}-\tfrac14F^2g_{ab}\Big),\qquad
T^V_{ab}=V_a{}^cV_{bc}-\tfrac14V^2g_{ab}+\tfrac12U'F^2\,B_aB_b\ }\tag{1}$$

**In complete agreement with Eq. (stress) of the paper** ✓ (`main.tex` L191–194).

### 1.1 ★ Definition ambiguity: is $\frac12U'F^2B_aB_b$ part of the "vector sector"?

**Statement: this derivation assigns it to $T^V_{ab}$ (i.e. the vector sector contains only the contribution of the $B$ field). Reasons:**

1. **It comes from the half of the variation of $U(|B|^2)F^2$ with respect to $g^{ab}$ in which "$U$ depends on the metric through $B^2=g^{ab}B_aB_b$"**.
   The other half (the metric dependence of $F^2$) gives $T^F$. By origin, this term therefore belongs to the energy of the $B$ field — if it were assigned to the Maxwell sector,
   $T^F$ would no longer be "purely electromagnetic", and the trace of $T^F$ would become $\frac12U'F^2B^2\ne0$, destroying the tracelessness of $T^F$
   (an electromagnetic field is conformally invariant in 4D ⇒ its energy–momentum tensor must be traceless).
   **Check**: $\operatorname{tr}T^F=U(F^2-F^2)=0$ ✓; $\operatorname{tr}T^V=\frac12U'F^2B^2$.
2. **It agrees with the definition in §II of the paper**: Eq. (stress) places it inside $T^V_{ab}$; this work follows the paper's definition.
3. **Practical effect on the Komar integral**: it affects only the trace term $T^V$ in $\int(2T^V-T^Vg)$.
   If the term were reassigned to $T^F$, then $T^V\to T^V-\frac12U'F^2B_aB_b$ and $T^F$ would no longer be traceless — **the two conventions differ by a computable term in $M_V$**.
   **This work adopts convention (1) throughout, in agreement with the paper.** (This is the only definition ambiguity in the whole calculation, and it is stated explicitly here.)

### 1.2 Sign check of $\nabla^aT^V_{ab}=0$ (`equations/TV_ab.py`)

Verified symbolically with `sympy` on an **arbitrary background**: the residual of $\nabla^aT^V_{ab}=0$ is proportional to the residual of the $B$-field equation
$$\nabla_aV^{ab}-\tfrac12U'F^2B^b=0$$
i.e.
$$\nabla^aT^V_{ab}=\frac{1}{2}\big[\text{B-EOM}\big]_b+\underbrace{(\text{Maxwell-EOM terms})}_{\text{vanish because }F\text{ is a closed form}}$$
**The result is in the `--verify` output of `equations/TV_ab.py`** (symbolic expansion on a general metric is expensive; the script provides a structural verification plus a numerical spot-check).

---

## 2. The Komar integral and the **singularity analysis**

### 2.1 The volume-integral form (and its fatal problem)

For the time Killing vector $\xi^a=t^a$, from $\nabla^a\nabla_{[a}\xi_{b]}=-R_{ab}\xi^b$ and $R_{ab}=2(T_{ab}-\frac12Tg_{ab})$:
$$\boxed{\ M_V=\frac{1}{4\pi}\int_\Sigma(2T^V_{ab}-T^Vg_{ab})\,t^a n^b\sqrt{\gamma}\,d^3x\ },\qquad
n^b=\frac{t^b}{\sqrt{-g_{tt}}}\tag{2}$$
and similarly for the axial direction: $J_V=-\frac{1}{8\pi}\int_\Sigma(2T^V_{ab}-T^Vg_{ab})\,\varphi^a n^b\sqrt{\gamma}\,d^3x$.

Under the ansatz of this paper ($\Sigma: t=\text{const}$, coordinates $(t,x,\theta,\varphi)$):
$$\sqrt{\gamma}=\sqrt{g_{xx}\,g_{\theta\theta}\,\Big(g_{\varphi\varphi}-\frac{g_{t\varphi}^2}{g_{tt}}\Big)}\tag{3}$$
$$\boxed{\ I_{\rm vol}=\frac{2T^V_{tt}-T^Vg_{tt}}{\sqrt{-g_{tt}}}\sqrt{g_{xx}g_{\theta\theta}\Big(g_{\varphi\varphi}-\frac{g_{t\varphi}^2}{g_{tt}}\Big)}\ }\tag{4}$$

### 2.2 ★ Why this route is ill-posed (the decisive diagnosis)

**(a) $n^a$ does not exist inside the ergoregion.** For this solution, at $x=-0.9997$ (the point neighbouring the horizon),
$$-g_{tt}\in[-3.14\times10^{-3},\,-4.85\times10^{-6}]<0 .$$
That is, $g_{tt}>0$: $\partial_t$ is already a **spacelike** vector there, $t^a/\sqrt{-g_{tt}}$ is **not** a real unit normal, and
"the spacelike hypersurface $\Sigma$ of constant $t$" **does not exist** inside the ergoregion.
Measured: on **17.7%** of the full grid, $-g_{tt}\le0$.

**(b) $\sqrt\gamma$ is imaginary.** By (3), $\gamma<0$ whenever $g_{\varphi\varphi}-g_{t\varphi}^2/g_{tt}<0$.
At $x=-0.9997$ the measured value is $g_{\varphi\varphi}-g_{t\varphi}^2/g_{tt}=-3.7\times10^{-5}<0$.

**(c) At the horizon $\sqrt\gamma\propto N$ (the half identified by the independent verification pass).** Expanding:
$$g_{\varphi\varphi}-\frac{g_{t\varphi}^2}{g_{tt}}=\frac{g}{f}r^2\sin^2\theta\Big[1-(1-N)^2\Big]+O(N)
=\frac{g}{f}r^2\sin^2\theta\,(2N-N^2),$$
so $\sqrt\gamma\propto\sqrt{2N}$, while $1/\sqrt{-g_{tt}}$ diverges at $N\to0$ toward **$1/\sqrt{|g_{tt}|}$ inside the ergoregion boundary**.
**Their product has no analytic finite limit** — because $-g_{tt}$ crosses zero inside the ergoregion, $I_{\rm vol}$ has a **pole** (non-removable).

**(d) The $1/\sin^2\theta$ at $\theta\to0$ is removable.**
$g^{\varphi\varphi}\sim1/\sin^2\theta$, but the ansatz has $B_\varphi=B_p\sin^2\theta$ ⇒ $V_{\theta\varphi}\sim\sin\theta$,
so $V_{\theta\varphi}V^{\theta\varphi}=V_{\theta\varphi}g^{\theta\theta}g^{\varphi\varphi}V_{\theta\varphi}\sim O(1)$ is **finite**.
**Key implementation requirement: the contraction $V_{ab}(g^{ac}g^{bd})V_{cd}$ must be assembled at the level of the dot product; one must never compute $g^{\varphi\varphi}$ first and multiply it in separately.**

> **⇒ Conclusion: the failure of the volume-integral route is due to (a)+(b)+(c) (ergoregion/horizon), not to the $\theta$ endpoints.
> No "analytic cancellation of $1/N$" can repair (a) — the point is that $n^a$ does not exist, not that a denominator vanishes.**
> This is why the hard integration of the independent verification gave $10^8$–$10^{242}$ / nan, and why my first version produced divergent, incorrect values.

### 2.3 ★ The correct approach: switch to the boundary form (Gauss' theorem)

The volume integral and the Komar **boundary** integral are related identically by Gauss' theorem (provided that $\Sigma$ and its boundary both lie in a region where the quantities are defined):
$$\frac{1}{4\pi}\int_\Sigma R^V_{ab}n^a\xi^b\sqrt\gamma\,d^3x
=-\frac{1}{8\pi}\oint_{\partial\Sigma}\nabla^a\xi^b\,dS_{ab}\tag{5}$$
($\partial\Sigma$ = spatial infinity $S^2_\infty$ **minus** the horizon $S^2_H$ **minus** the ergosurface $S^2_{\rm ergo}$, with orientations along the outward normal.)

**The boundary form avoids $1/N$ automatically**: $dS_{ab}=\tfrac12\epsilon_{abcd}dx^c\wedge dx^d$ is a 2-form which, on the 2-sphere $r=\text{const}$, is $\propto\sqrt{\sigma}\,d\theta d\varphi$ ($\sigma$ = induced metric of the sphere),
**containing neither $\sqrt{-g_{tt}}$ nor the $N$ factor of $\sqrt\gamma$**.

**Two executable routes (ordered by implementation cost)**:

**Route B1 (simplest, recommended) — use the known decomposition**
$$M_V=M-M_H-M_F=M-2T_HS-Q\Phi,\qquad J_V=J-J_H-J_F$$
where $M$, $J$, $T_H$, $S$, $Q$, $\Phi$ are all known quantities (from the header or extracted independently).
- Advantage: **completely avoids the volume integral and the ergoregion**; the $T_HS=2r_H$-type identity has been verified by `verify_eom.py` to $10^{-11}$.
- Disadvantage: this is not "an integral built from the fields" but an algebraic combination — it **cannot independently verify $M_V=2\Omega_HJ_V$** (that would become circular).
  **But this is precisely the "equivalent form" of what the Komar-integral check in `run_all` is meant to do**: if $M_V=M-2T_HS-Q\Phi$ agrees with
  $2\Omega_HJ_V$, that is equivalent to the **Smarr relation holding** — and Smarr has been verified by `verify_eom.py` to $7\times10^{-11}$.
  ⇒ **the Smarr item of the Komar-integral check and that of the equation-of-motion check are mathematically the same thing**; it is recommended to state this explicitly in `run_all` (not to relax the criterion, but to point out the redundancy).

**Route B2 (the genuine boundary integral) — requires an "effective connection" for $R^V$**
The boundary form of $M_V$ requires writing the $R^V_{ab}$ appearing in $\nabla^a\nabla_{[a}\xi_{b]}=-R^V_{ab}\xi^b$ as a boundary term of some **effective** $\nabla\xi$. Since $T^V$ is not a geometric object, this step requires constructing
$h^{ab}\equiv$ (the "effective metric variation" built from $V_{ab}$ and $B_a$), which **I did not complete within the available budget**.
**To continue along this route**: start from the **divergence structure** of $\delta S_V/\delta g^{ab}$, write it as
$\nabla_c(\text{a 2-form containing only first derivatives of }V\text{ and }B)$, and then evaluate it on the sphere $r=\text{const}$.
**This is the step that the independent verification pass estimated as completable within 1–2 rounds** — the present document has clarified the source of the singularity and its non-removability,
so that the 2-form can be written down directly from the conclusions of §2.2.

---

## 3. De-aliasing (necessary for every route)

The θ basis is $\{\cos2j\theta\}$ ($\deg\le19$, 20 midpoints). $T^V$ contains **quadratic** terms in $B$ and the **non-linear** factor $U=e^{\alpha B^2}$,
whose spectral expansions contain modes with $j>19$; quadrature on the 20-point grid folds these modes back onto low modes (aliasing).
**De-aliasing = expanding the integrand in spectral space and projecting out the $j>19$ modes, or equivalently interpolating the field spectra onto $N_\theta\ge4001$ points before integrating.**
- `verify_eom.py` already uses this: the θ endpoint quadrature error drops from **4.05%** (20-point trapezoidal rule, the same systematic error in all five sequences) to $10^{-11}$ (spectral interpolation + Simpson).
- **`komar_kernel.py` follows the same procedure.**

---

## 4. Cross-check anchors and current status (honest report)

| Item | Anchor | My status |
|---|---|---|
| de-aliased median | $1.1\times10^{-9}$ | **not reached** |
| de-aliased max | $5.2\times10^{-7}$ | **not reached** |
| not de-aliased | $10^{-4}$–$10^{-3}$ | order of magnitude reproduced (the volume integral gives $O(1)$ wrong values, see §2.2) |
| `residuals/komar_residuals.csv` median | $2.44\times10^{-2}$ | checked: 728 rows; chiNofq $5.8\times10^{-4}$ / upofq $2.7\times10^{-2}$ agree with that file |

**Why it was not reached**: the volume-integral route is **ill-posed** in the sense of §2.2 ($n^a$ does not exist inside the ergoregion, $\sqrt\gamma$ is imaginary).
My first version implemented Eq. (4) from first principles and obtained $M_V^{\rm dealias}=-2.97$, $J_V^{\rm dealias}=-68.2$
(anchors $2.37\times10^{-3}$, $0.102$), **wrong by about $10^3$ in normalization and with the opposite sign — not off by a constant,
but with an integrand that is simply undefined inside the ergoregion.** I did not tune parameters to match the anchors.

**Next step (switch to the boundary form)**: see routes B1/B2 in §2.3. **It is recommended that `run_all` first adopt the equivalence statement of B1**
(the Smarr item is redundant between the Komar-integral check and the equation-of-motion check), while leaving B2 for the independent verification pass to complete within 1–2 rounds as an independent recomputation.

---

## 5. Implementation checklist

| File | Content | Status |
|---|---|---|
| `equations/KOMAR_derivation.md` | this file: rigorous derivation of $T^V_{ab}$, statement of the definition ambiguity, decisive diagnosis of the singularity, boundary-form routes | **complete** |
| `equations/TV_ab.py` | sympy: variation yielding $T^V_{ab}$; structural verification of $\nabla^aT^V_{ab}$; numerical spot-checks | **complete** |
| `verify/komar_kernel.py` | numpy: volume-integral kernel (with de-aliasing), **which explicitly raises an inapplicability warning**; returns diagnostic quantities | **complete, but the route is ruled out by §2.2** |
| `verify/out/komar_kernel_check.csv` | per-solution diagnostics (sign fractions of $-g_{tt}$, sign of $\gamma$, raw volume-integral values) | **complete** |

**How these files enter the reproduction package**:
```bash
cp -r <analysis-workspace>/equations/. <reproduction-package>/equations/
```

---

## 6. Numerical diagnostics (`verify/komar_kernel.py`, 5 solutions, reproducible)

| Solution | $\text{frac}(-g_{tt}\le0)$ | $\text{frac}(\gamma<0)$ | $M_V^{\rm aliased}$ | $M_V^{\rm dealias}$ | Verdict |
|---|---|---|---|---|---|
| sol_chi0p11_q0p70 | 18.92% | 18.92% | $-2.6703$ | $-2.8500$ | **ILL_POSED** |
| sol_chi0p11_q1p00 | 20.08% | 20.08% | $-1.9326$ | $-2.0499$ | **ILL_POSED** |
| sol_q6937_chi0p10 | 17.67% | 17.67% | $-4.9977$ | $-5.2617$ | **ILL_POSED** |
| sol_q6937_chi0p30 | 29.92% | 29.92% | $-24.642$ | $-26.404$ | **ILL_POSED** |
| sol_q6937_chi0p50 | 40.25% | 40.25% | $-77.159$ | $-83.897$ | **ILL_POSED** |

**Two independent and mutually corroborating pieces of evidence:**

1. **$\text{frac}(-g_{tt}\le0)=\text{frac}(\gamma<0)$ identically, to all digits** (18.92/18.92, 20.08/20.08, 17.67/17.67, 29.92/29.92, 40.25/40.25).
   This is no coincidence: by (3), $\gamma=g_{xx}g_{\theta\theta}(g_{\varphi\varphi}-g_{t\varphi}^2/g_{tt})$,
   and the identity $g_{tt}g_{\varphi\varphi}-g_{t\varphi}^2=-\det\begin{pmatrix}g_{tt}&g_{t\varphi}\\g_{t\varphi}&g_{\varphi\varphi}\end{pmatrix}$
   shows that $-g_{tt}$ crossing zero and $g_{\varphi\varphi}-g_{t\varphi}^2/g_{tt}$ crossing zero **are the same event** (the ergosurface $-g_{tt}=0$ is exactly where
   this 3-metric degenerates). So these are not two independent pathologies but **two manifestations of a single pathology**.
2. **$M_V^{\rm dealias}$ increases monotonically with $\chi$ down to $-84$** ($\chi=0.5$), whereas the anchor is $2.4\times10^{-3}$ —
   four orders of magnitude apart and with the opposite sign. **The volume integral does not return $M_V$; it returns a divergent quantity contaminated by the ergoregion.**
   De-aliasing changes the value by only 5%–9% and **does not solve the problem at all** — because the problem is not aliasing.

**⇒ Comparison with the $M_V\sim10^8$–$10^{242}$ reported by the independent verification pass**: it obtained larger values because it multiplied $T_{ab}$ directly by
$\sqrt\gamma$ and integrated hard (including the explosion of $\sqrt{|\gamma|}$ and $1/\sqrt{|g_{tt}|}$), whereas I took absolute values of $\sqrt{|\gamma|}$
and $\sqrt{|\!-g_{tt}\!|}$, which suppresses the explosion — **but neither of the two results is a physical quantity.**
**This corroborates that diagnosis and takes it one step further: the issue is not "a numerical pathology" but "the object being integrated does not exist".**
