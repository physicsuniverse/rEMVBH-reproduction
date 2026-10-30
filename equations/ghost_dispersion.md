# Full dispersion derivation including the term linear in $\omega$

**Author**: the numerical-physics analysis|**Upstream**: the quadratic-action derivation and the numerical evaluation of the criteria (companion notes)
**Conventions**: signature $(-,+,+,+)$; coordinates $(t,r,\theta,\varphi)$; $\bar B_a$ is the background vector, $\delta B_a$ the perturbation.

---

## 1. From $S_2$ to the linearized field equations (following Eq. (2.1) of the quadratic-action derivation)

$$S_2=\int\sqrt{-g}\left[-\tfrac14\delta V_{ab}\delta V^{ab}-\tfrac12M^{ab}\delta B_a\delta B_b\right],
\qquad M^{ab}=m_{\rm eff}^2\left(g^{ab}+2\alpha\bar B^a\bar B^b\right).$$

With $\delta V_{ab}=\nabla_a\delta B_b-\nabla_b\delta B_a$, varying with respect to $\delta B_d$:
$$\delta S_2=\int\sqrt{-g}\left[-\delta V^{ab}\partial_a(\delta B_b)-M^{bd}\delta B_b\right]\delta(\delta B_d)
\;\Longrightarrow\;
\boxed{\ \nabla_a\delta V^{ab}-M^{bd}\delta B_d=0\ }\tag{1.1}$$

(Derivation: $-\frac14\delta(\delta V_{ab}\delta V^{ab})=-\frac12\delta V^{ab}\delta(\delta V_{ab})=-\frac12\delta V^{ab}\cdot2\nabla_a\delta(\delta B_b)$, and integration by parts gives $+\nabla_a\delta V^{ab}\delta(\delta B_b)$.)

---

## 2. Eikonal expansion and the full dispersion matrix

In local coordinates take $k_a=(-\omega,k_r,k_\theta,\ell)$ and $\delta B_d\propto e^{ik_ax^a}$, replacing the highest derivatives in (1.1) by $\nabla_a\to ik_a$:

$$\nabla_a\delta V^{ab}=g^{bd}\left(\Box\delta B_d-\nabla_d\nabla^c\delta B_c\right)
\;\longrightarrow\; g^{bd}\left(-k^2\delta B_d+k_dk^c\delta B_c\right),$$

where $k^b\equiv g^{bc}k_c$ and $k^2\equiv g^{ab}k_ak_b$. The highest-order part of (1.1) is therefore
$$\left[-k^2g^{bd}+k^bk^d-M^{bd}\right]\delta B_d=0
\;\Longrightarrow\;
\boxed{\ \mathcal D^{bd}\,\delta B_d\equiv\left[k^2g^{bd}-k^bk^d+M^{bd}\right]\delta B_d=0\ }\tag{2.1}$$

**Dispersion relations**:
$$\boxed{\ \det\mathcal D(\omega;k,\ell)=0\ },\qquad
\mathcal D=k^2G-k\,k^{T}+M,\quad G\equiv(g^{bd}),\quad (k\,k^T)^{bd}\equiv k^bk^d .\tag{2.2}$$

`k·k^T` is a rank-1 matrix; `G` and `M` are independent of $\omega$, $k^b$ is linear in $\omega$ and $k^2$ quadratic in $\omega$, so **$\det\mathcal D$ is a polynomial in $\omega$ of degree $\le 8$**.

### 2.1 Closed form (matrix determinant lemma)

Let $A\equiv k^2G+M$. Then $\mathcal D=A-k\,k^T$, so
$$\boxed{\ \det\mathcal D=\det\!\left(k^2G+M\right)\cdot\left[1-k^b\left(k^2G+M\right)^{-1}_{bd}k^d\right]\ }\tag{2.3}$$

> **2026-10-17 correction** (found by the 16-case verification performed by `equations/ghost_dispersion.py`):
> the formula was originally printed as $k_a(\cdot)^{-1}_{ab}k_b$ (**covariant**), which is **wrong**. The determinant lemma takes $u=v=k^b$ and requires **contravariant** indices.
> Numerical criterion: the contravariant version agrees to $4.7\times10^{-16}$, whereas the printed version deviates by **$1.98\times10^{-1}$**.
> The sentence in line 24 of this document ("take $u=k$ (contravariant components)") was correct all along; **only the formula was misprinted**. The physical conclusions are unaffected.
(Using $\det(A-uv^T)=\det A\,(1-v^TA^{-1}u)$ with $u=k$ (contravariant components) and $v^T=k^T$.)

### 2.2 Casting as a quadratic problem: the quadratic matrix polynomial in $\omega$ (used for the actual solution)

Splitting $k^b$: $k_a=(-\omega,k,0,\ell)$ (taking $k_r=k,\;k_\theta=0$),
$$k^b=\omega\,p^b+q^b,\qquad p^b\equiv-g^{b0},\quad q^b\equiv g^{b1}k+g^{b3}\ell,$$
$$k^2=a\,\omega^2+b\,\omega+c,\qquad
a=g^{00},\;\; b=-2\left(g^{01}k+g^{03}\ell\right),\;\; c=g^{11}k^2+2g^{13}k\ell+g^{33}\ell^2 .$$
Substituting into (2.2):
$$\boxed{\ \mathcal D(\omega)=\omega^2A_2+\omega A_1+A_0\ }$$
$$\boxed{\ A_2=aG-p\,p^T,\qquad A_1=bG-p\,q^T-q\,p^T,\qquad A_0=cG-q\,q^T+M\ }\tag{2.4}$$

**This is a quadratic eigenvalue problem (QEP)**, which can be linearized into an $8\times8$ generalized eigenvalue problem (as done here:
$\begin{pmatrix}0&I\\-A_0&-A_1\end{pmatrix}\chi=\omega\begin{pmatrix}I&0\\0&A_2\end{pmatrix}\chi$, solved with `scipy.linalg.eig`, with no fitting and no iteration).

### 2.3 Key structural conclusion: **it is no longer $\omega^2=\cdots$**

$A_1\neq0$ if and only if $bG\neq p q^T+q p^T$. Since
- $g^{t\varphi}\neq0$ (rotation) makes $b=-2(g^{01}k+g^{03}\ell)\neq0$ whenever $k$ or $\ell\neq0$;
- the $\bar B^a\bar B^b$ term gives $M$ non-zero off-diagonal entries $M^{t i}\neq0$ (the non-diagonal components $V_{13},V_{14},V_{23},V_{24}$ of the numerical evaluation),

$\mathcal D(\omega)$ **contains a term linear in $\omega$**, and $\det\mathcal D$ contains odd powers of $\omega$. **This is exactly where the "purely quadratic dispersion relation" premise of [2605.13920] breaks down** — they rely on "no off-diagonal terms in the time-like direction".

**Time-reversal check**: $\mathcal D$ is invariant under $\omega\to-\omega,\ \ell\to-\ell$ ($G,M$ real and $b,c$ invariant under the transformation),

> **2026-10-17 correction**: the invariant at the matrix level is $k_a\to-k_a$, i.e. $(\omega,k_r,\ell)$ **all change sign at once** — the original text omitted $k_r\to-k_r$.
> Measured: $D(-\omega;-k,-\ell)-D(\omega;k,\ell)=0$ holds exactly; with $k_r$ held fixed the residual is **$1.40$**.
> The root-set conclusion still holds, however, because $\det\mathcal D$ is measured to be exactly even in $k_r$. The physical conclusions are unaffected. The root set therefore satisfies
$$\omega\in{\rm roots}(k,\ell)\iff -\omega\in{\rm roots}(k,-\ell);\qquad \overline{\omega}\in{\rm roots}(k,\ell)\ (\text{real coefficients}).$$
**This is a strong self-check of the numerical solution** (see §3).

---

## 3. Numerical implementation and self-checks

`equations/ghost_dispersion.py`: reads the background fields from the cached solution data and rebuilds $g^{ab}$, $\bar B_a$, $\bar B^a$, $F^2$, $m_{\rm eff}^2$ and $M^{ab}$ (reusing the metric/gauge-field construction of the criterion evaluation),
constructs (2.4) at the selected grid points and solves the QEP. To improve the conditioning, the equations are made dimensionless with $\mu^2\equiv\max|M_{ab}|$ at that point
($\hat\omega=\omega/\mu$, etc.); $m_{\rm eff}^2$ spans eight orders of magnitude across the grid, and without rescaling one obtains spurious roots (this happened once and has been fixed).

**Self-check 1**: every root set **satisfies the conjugate pairing** $(z,\bar z)$ one by one — required for a polynomial with real coefficients ✓ (see the `roots` column of the accompanying root-set table).
**Self-check 2**: at $(k,\ell)=(0,0)$, $\mathcal D|_{\omega=0}=M$, and here $p^b=-g^{b0}$, $q=0$, $b=0$, so
$\mathcal D(\omega)=\omega^2(aG-pp^T)+M$ is a matrix polynomial in $\omega^2$ ⇒ the roots must occur in $\pm$ pairs with either ${\rm Re}\,\omega=0$ or all real.
Measured: **all purely imaginary** ✓ (see §4).

---

## 4. Actual results for three representative solutions

Selected from the table of criterion outcomes over all solutions by quantiles of $V_{od,rel}$ (large/medium/small):

| Label | File | q | χ | α | $V_{od,rel}$ | criterion verdict |
|---|---|---|---|---|---|---|
| **large** | `axisolsofchiNofq-222.mx` | 0.4727 | 0.650 | 20 | **0.768** | framework_not_applicable |
| **median** | `axisolsofchiNofq-400.mx` | 0.2947 | 0.500 | 20 | **0.444** | framework_not_applicable |
| **small** | `axisolsofchiNofq-584.mx` | 0.1107 | 0.130 | 20 | **0.088** | framework_marginal |

For each solution all $\omega$ roots were obtained at $\theta=\pi/2$ (the 10th $\theta$ collocation point), at $i_x\in\{0,15,30,45,59\}$ (from the horizon side to the infinity side),
and for the 60 $(k,\ell)$ combinations $(k,\ell)\in\{0,0.5,1,2\}\times\{0,0.5,1\}$ (180 rows in the accompanying root-set table).

**Roots at $(k,\ell)=(0,0)$ (zero momentum): all purely imaginary** ⇒ $\omega^2<0$ ⇒ **tachyon**:

| Solution | $i_x=0$ | 15 | 30 | 45 | 59 |
|---|---|---|---|---|---|
| large | ±0.02526i (×2), ±0.03617i | ±0.01926i, ±0.02846i | ±0.01116 (real), ±0.00722i | ±0.1652i, ±0.1425i | ±0.000680i (×3) |
| median | ±0.00538i, ±0.00640i | ±0.00496i, ±0.00624i | ±0.06474i, ±0.05875i | ±0.3266i, ±0.3170i | ±0.000939i (×3) |
| small | ±0.01292i, ±0.01301i | ±0.00178i, ±0.00181i | ±0.07894i, ±0.07877i | ±0.8508i, ±0.8502i | ±0.002123i (×3) |

(Units: local $\mu=\sqrt{\max|M|}$; only finite roots are listed for each $i_x$; 2 roots are filtered out by the singular $A_2$.)

**Statistics** (180 $(k,\ell)$ samples):
- `large`: 49/60 have a root with ${\rm Im}\,\omega>0$; $\max|{\rm Im}\,\omega|=1.57$
- `median`: 40/60; $1.75$
- `small`: 40/60; $1.97$
- **In no case does $\omega^2<0$ occur on a real root** (all real roots have $\omega^2>0$) — i.e. **all instabilities appear as purely imaginary or complex roots, not as the "negative real root" type of tachyon**.
  This is consistent with the structure $\mathcal D|_{\omega=0}=M$: $m_{\rm eff}^2<0$ drives $\omega^2$ directly negative.

---

## 5. What this section provides and does not provide

**It provides**:
1. (2.1), (2.2): the complete dispersion matrix of the vector perturbation on an **arbitrary background**, $\mathcal D^{bd}=k^2g^{bd}-k^bk^d+M^{bd}$.
2. (2.3): the closed form $\det\mathcal D=\det(k^2G+M)[1-k(k^2G+M)^{-1}k]$.
3. (2.4): the **quadratic matrix polynomial** form in $\omega$ plus the QEP solution method (no fitting, no iteration).
4. The **structural conclusion** of §2.3: whenever $k$ or $\ell\neq0$ together with $g^{t\varphi}\neq0$ (rotation) or $M^{ti}\neq0$ (the non-diagonal components of the hair),
   $\det\mathcal D$ **necessarily contains odd powers of $\omega$** ⇒ the "purely quadratic" premise of [2605.13920] indeed fails.
5. §4: all roots of the three representative solutions over 180 $(k,\ell)$ samples (saved to disk).

**It does not provide**:
1. A **ghost criterion** — a dispersion relation **cannot** in principle decide whether a ghost is present (see the risk-assessment note, §1).
2. A conclusion on gradient stability at large $k$ (the locally frozen-coefficient analysis requires more careful treatment in that limit).
3. Non-local/bound-state effects (this analysis is an eikonal/local approximation).
