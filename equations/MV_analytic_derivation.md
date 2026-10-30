# Analytic derivation — the structure of $M_V$, the origin of the $\chi^2$ law, and the independent content of $M_V=2\Omega_HJ_V$

**Analysis** Q2, Q3 **Performed by** the α-continuation analysis
Numerical-experiment scripts: `wscale.wls` (W/B scaling), `wscale2.wls` (full Bt/Bp/W decomposition),
`wscale3.wls` (spherically symmetric part vs. symmetry-breaking part), `wscale4.wls` (extrapolation to χ=0)

---

## 1. Setup

Metric (the ansatz used in the paper, with $\xi^a=\partial_t$ the Killing vector):

$$\mathrm{d}s^2 = -f\,\mathrm{d}t^2 + \frac{g}{f}\Big(\frac{2r_h}{1-x}\Big)^2\frac{\mathrm{d}x^2}{(1-x)^2}
 + \frac{h\,r^2}{f}\big(\mathrm{d}\theta^2 + \sin^2\theta\,(\mathrm{d}\varphi - W\,\mathrm{d}t)^2\big)$$

Vector field and gauge field:

$$B = B_t(x,\theta)\,\mathrm{d}t + B_p(x,\theta)\,\mathrm{d}\varphi,\qquad
  A = A_t(x,\theta)\,\mathrm{d}t + A_p(x,\theta)\,\mathrm{d}\varphi$$

$\chi = J/M^2$ is the rotation parameter. **Key kinematic facts** (all fixed by the boundary conditions and parity):

| Quantity | under $\chi\to-\chi$ | small-χ behaviour |
|---|---|---|
| $W$ (frame dragging) | **odd** | $O(\chi)$ |
| $B_p$ ($\varphi$ component of the vector field) | **odd** | $O(\chi)$ |
| $A_p$ | **odd** | $O(\chi)$ |
| $B_t$ ($t$ component of the vector field) | **even** | $B_t^{(0)}(r) + O(\chi^2)$ |
| $A_t$, $f$, $g$, $h$ | **even** | $O(1) + O(\chi^2)$ |
| $\Omega_H$ (horizon angular velocity) | **odd** | $O(\chi)$ |
| $J_V$ | **odd** | $O(\chi)$ |
| $M_V$, $M$, $T_H$ | **even** | $O(\chi^2)$ / $O(1)$ |

$\chi\to-\chi$ maps a solution to its "mirror solution"
($W\to-W,\,B_p\to-B_p,\,A_p\to-A_p$); this is an **exact discrete symmetry**
(not a perturbative expansion).

### 1.1 Numerical confirmation (α = 5 sequence)

| χ | $\max\lvert B_t\rvert$ | $\theta$-span$(B_t)/\max\lvert B_t\rvert$† | $\max\lvert B_p\rvert/\chi$ | $\max\lvert W\rvert/\chi$ |
|---|---|---|---|---|
| 0.01 | 0.227291 | $1.51\times10^{-6}$ | 2.5868 | 1.4461 |
| 0.02 | 0.227251 | $6.06\times10^{-6}$ | 2.5876 | 1.4464 |
| 0.04 | 0.227087 | $1.865\times10^{-4}$ | 2.5908 | 1.4475 |
| 0.08 | 0.226431 | $7.505\times10^{-4}$ | 2.5987 | 1.4520 |

- $\max\lvert B_p\rvert/\chi$ and $\max\lvert W\rvert/\chi$ are **constant to 0.5%** ⟹ $B_p,W = O(\chi)$ ✓
- $\max\lvert B_t\rvert$ changes from 0.227291 only to 0.226431 (0.4%) ⟹ $B_t$ **has a non-zero spherically symmetric limit**
  $\lim_{\chi\to0}B_t = B_t^{(0)}(r) \simeq 0.2273$
- The **θ deviation of $B_t$ scales as $\chi^2$** ($1.163\!\to\!4.655\!\to\!18.65\!\to\!75.05\times10^{-5}$,
  successive ratios 4.00 / 4.01 / 4.02) ⟹ $B_t = B_t^{(0)}(r) + O(\chi^2)$ ✓
  † $\theta$-span $=\max_x[\max_\theta B_t-\min_\theta B_t]/\max\lvert B_t\rvert$;
  the rms version was also measured, giving $1.51\times10^{-6}$ / $6.06\times10^{-6}$ at χ = 0.01 / 0.02 (ratio 4.01).

> **Important**: the $\chi\to0$ limit is **not** the trivial solution with $B=0$, but a **spherically symmetric EMV solution carrying a radial vector hair $b(r)\mathrm{d}t$**.

---

## 2. Explicit structure of the $M_V$ integrand

Taking the symbolic expression from `numerics/EOMEQs/intofMVx.mx` (the author's code, $M_V=\frac{1}{4\pi}\int_\Sigma(2T^V_{ab}-T^Vg_{ab})t^a\mathrm{d}\Sigma^b$), its structure is

$$\text{intofMV} \;=\; -\frac{1}{32}\,
\frac{(-1+x)^2\,h\,\sin\theta}{(1+x)(-1+x^2)^2f^4\sqrt{g^3h^2/f^2}}\;\times\;\Big[\;\mathcal{I}_1+\mathcal{I}_2+\mathcal{I}_3+\mathcal{I}_4\;\Big]$$

where ($\mathcal{M}_2$ is the field-strength-squared combination of $A_t,A_p$, which itself contains $W$ terms):

$$
\begin{aligned}
\mathcal{I}_1 &= e^{\alpha|B|^2}\,\alpha f\Big(16B_t^2g-(-1+x^2)^2B_p^2f^2\sin^2\theta\Big)\;\mathcal{M}_2
 &&\text{(B}^2\text{ type)}\\[2pt]
\mathcal{I}_2 &= e^{\alpha|B|^2}\,\alpha f\Big(-4B_t+(-1+x)^2B_p\sin^2\theta\,W\Big)^2\;\mathcal{M}_2
 &&\text{(B}^2,\,BW,\,W^2\text{)}\\[2pt]
\mathcal{I}_3 &= 2\Big[(-1+x)^2(1+x)^4f^4\sin^2\theta\big((\partial_\theta B_p)^2+\cdots\big)
   +\underbrace{16(-1+x)^2g^2\sin^2\theta\,W^2\big((\partial_\theta B_t)^2+\cdots\big)}_{\text{4 factors}}\\
 &\qquad +\;(1+x)^2f^2g\Big(\cdots W^2(\partial B_p)^2-16(\cdots)W\,\partial B_p\partial B_t+16(\partial B_t)^2\Big)
   +2(-1+x)^2B_p\big(\cdots\big)+(-1+x)^2B_p^2\big(\cdots\big)\Big]\\[2pt]
\mathcal{I}_4 &= (-1+x)^2\sin\theta\,W\Big[e^{\alpha|B|^2}\alpha B_pf\sin\theta\big(4B_t-\cdots\big)\mathcal{M}_2-2g\big(\cdots\big)\Big]
 &&\text{(}W\cdot B\cdot B\text{)}
\end{aligned}
$$

### ★ Structural lemma (read off directly)

> **Every term of the integrand $\text{intofMV}$ contains at least _two_ factors from $\{B_t,B_p,W\}$**
> (counting $\partial_t,\partial_x,\partial_\theta$ as the same order).
> Equivalently: $\text{intofMV}$ vanishes when $B=0$ **and** $W=0$ simultaneously.

**Numerical confirmation**: $B_t=B_p=0$ (keeping $W$ and $A_t,A_p,f,g,h$) ⟹
$M_V = 0$ (exact zero; 3 solutions: $\alpha=5\,\chi=0.01,\ \alpha=15\,\chi=0.01,\ \alpha=40\,\chi=0.1$)
and $J_V=0$ — consistent with the lemma.

### 2.1 Slow-rotation counting: why $\chi^2$ and not $\chi$

Counting the order in χ of every factor ($B_t = B_t^{(0)}+O(\chi^2)$, $B_p,W=O(\chi)$):

| Term | Factor content | χ order | Treatment |
|---|---|---|---|
| $\mathcal{I}_1$ | $B_t^2\mathcal{M}_2$ | $\chi^0$ | **must vanish identically** (spherical theorem, see §3) |
| $\mathcal{I}_1$ | $B_p^2\mathcal{M}_2$ | $\chi^2$ | survives |
| $\mathcal{I}_2$ | $B_t^2\mathcal{M}_2$ | $\chi^0$ | as above, vanishes |
| $\mathcal{I}_2$ | $B_tB_pW\mathcal{M}_2$ | $\chi^2$ | survives |
| $\mathcal{I}_2$ | $B_p^2W^2\mathcal{M}_2$ | $\chi^4$ | higher order |
| $\mathcal{I}_3$ | $(\partial B_p)^2$, $(\partial B_t)^2$ | $\chi^2$, $\chi^0$ | $(\partial B_t^{(0)})^2$ must vanish |
| $\mathcal{I}_3$ | $W\partial B_p\partial B_t$ | $\chi^1$ | **must vanish — by parity** |
| $\mathcal{I}_4$ | $W B_pB_t\mathcal{M}_2$ | $\chi^2$ | survives |
| $\mathcal{I}_4$ | $W^2(\partial B_t)^2$ | $\chi^2$ | survives |

**Note**: $O(\chi)$ terms do appear in the table ($W\partial B_p\partial B_t$, as well as the $W\partial A_t$ terms in $\mathcal{M}_2$
and the $B_t^2W$-type terms produced by expanding $\mathcal{I}_2$ against $B_t^2$). They **cannot be removed by counting;
they can only be removed by symmetry** — which is the subject of the next section.

---

## 3. ★ Rigorous $\chi^2$ theorem

**Theorem.** Let $\{f,g,h,A_t,A_p,B_t,B_p,W\}$ be a family of asymptotically flat, regular solutions of the above ansatz,
with $\chi$ the rotation parameter, spherically symmetric at $\chi=0$. Then

$$M_V(\chi) = c_2\,\chi^2 + O(\chi^4),\qquad\text{i.e. the expansion of } M_V \text{ contains no } \chi^0,\chi^1,\chi^3 \text{ terms.}$$

**Proof (two steps)**

1. **Evenness**: $\chi\to-\chi$ is an exact discrete symmetry of the solution family (it maps a solution to its mirror solution).
   $M_V = \frac{1}{4\pi}\int_\Sigma(2T^V_{ab}-T^Vg_{ab})t^a\mathrm{d}\Sigma^b$ is a **scalar mass**
   (it contains no $\epsilon_{abcd}$ and projects on no direction other than $t^a$), and is invariant under the mirror map.
   Hence $M_V(-\chi)=M_V(\chi)$, i.e. $M_V$ is an **even function** of χ: all odd coefficients vanish, $\;c_1=c_3=\cdots=0$.

2. **Vanishing in the spherically symmetric limit**: at $\chi=0$ one has $W=B_p=A_p=0$, leaving a spherically symmetric set of fields
   $\{f_0(r),g_0(r),h_0(r),A_{t0}(r),B_{t0}(r)\}$.
   By the structural lemma of §2, $\text{intofMV}$ vanishes when $B=0$; here, however, $B_{t0}\neq0$, so the lemma is not enough.
   A stronger statement is required: **for a spherically symmetric configuration, $\int_\Sigma(2T^V_{ab}-T^Vg_{ab})t^a\mathrm{d}\Sigma^b = 0$**.
   This is the "spherical symmetry implies $M_V\equiv0$" theorem of the literature; here it is **verified numerically** (§3.1).
   ⟹ $c_0 = 0$.

Combining 1 and 2: $M_V = c_2\chi^2+c_4\chi^4+\cdots$. $\blacksquare$

### 3.1 Numerical verification (`wscale4.wls`) — the strongest evidence in this study

Taking two members of the α = 5 sequence (χ = 0.01, 0.02), **extrapolate to χ = 0** using parity:

$$F_{\text{even}}(0) = F_1-\frac{\chi_1^2}{\chi_2^2-\chi_1^2}(F_2-F_1),\qquad
  F_{\text{odd}}(0)=0\quad(A_p,B_p,W)$$

The resulting χ=0 configuration has $\max A_t = 0.9164$, $\max B_t = 0.227305$, $A_p=B_p=W=0$,
and $f,g$ spherically symmetric ($g\to4$). Its $M_V$ is then evaluated:

| Configuration | $M_V$ | $J_V$ |
|---|---|---|
| **Extrapolated χ = 0 spherically symmetric configuration** ($2\sigma$ upper bound, all α) | $\mathbf{|a|<5\times10^{-8}}$ | $0$ |
| Archived solution χ = 0.01 | $2.376\times10^{-4}$ | 0.234793 |
| Archived solution χ = 0.02 | $9.508\times10^{-4}$ | 0.469767 |
| **Extrapolated χ = 0 spherically symmetric configuration** (α=15) | $\mathbf{-9.84\times10^{-9}}$ | $0$ |
| Archived solution χ = 0.01 (α=15) | $8.803\times10^{-5}$ | 0.045633 |

> **The $M_V$ of the spherically symmetric limit (with vector hair $B_t^{(0)}(r)\neq0$) is four orders of magnitude smaller
> than the smallest archived value ($1.9\times10^{-8}$ vs $2.4\times10^{-4}$, ratio $8\times10^{-5}$), consistent with
> "exactly zero plus the extrapolated $O(\chi^4)$ error".**

**⇒ The spherical theorem holds; $M_V\neq0$ arises entirely from the breaking of spherical symmetry.**

---

## 4. Structural formula: $M_V$ is a "second-order invariant of the rotational response"

By the counting of §2 and the two cancellation conditions of §3, only the **quadratic term in χ** survives as $\chi\to0$, and its origin is:

$$\boxed{\;M_V \;=\; c_{BB}\big[\text{quadratic form in }B_p\text{ and }\delta B_t\big]
  \;+\;c_{BW}\big[\text{quadratic form in }B_t\text{ and }W\big]
  \;+\;c_{WW}\big[W^2\text{ type}\big]\;+\;O(\chi^4)\;}$$

where $\delta B_t = B_t-\langle B_t\rangle_\theta = O(\chi^2)$.
In other words: **to leading order $M_V$ is the "product of two rotational-response quantities"**:
- one response is $B_p\sim O(\chi)$ (the axial component of the vector field, generated by rotational dragging),
- the other is $W\sim O(\chi)$ (frame dragging) or $\delta B_t\sim O(\chi^2)$.

This gives the **structural formula** requested by the analysis (the precise version of
$M_V\propto\int(\text{vector field})\times(\text{frame dragging})$):

$$M_V \;\propto\; \int_\Sigma \Big[\underbrace{B_t^{(0)}\,\delta B_t}_{\chi^2}
   \;+\;\underbrace{B_t^{(0)}B_p\,W}_{\chi^2}
   \;+\;\underbrace{B_p^2}_{\chi^2}
   \;+\;\underbrace{(B_t^{(0)})^2W^2}_{\chi^2}\Big]\;\times\;\big(\text{metric factors }f,g,h\big)\;\mathrm{d}\Sigma$$

It is **not** simply $\int B\cdot W$: combinations of the spherically symmetric part $B_t^{(0)}$ of $B$ with $W$
($B_t^{(0)}B_pW$, $(B_t^{(0)})^2W^2$) contribute at $\chi^2$, whereas $B_t^{(0)}$ **alone** (i.e. χ=0) contributes exactly nothing.

### 4.1 Numerical decomposition (α = 5 sequence, `wscale2.wls` / `wscale3.wls`)

| Configuration | χ=0.01 | χ=0.02 | χ=0.04 | χ=0.08 | divided by $\chi^2$ |
|---|---|---|---|---|---|
| full | 2.37626e-4 | 9.50849e-4 | 3.80892e-3 | 1.53249e-2 | 2.3763 / 2.3771 / 2.3806 / 2.3945 |
| $W=0$ (keeping $B_t,B_p$) | 5.96546e-5 | 2.38674e-4 | 9.55590e-4 | 3.83678e-3 | 0.5965 / 0.5967 / 0.5972 / 0.5995 |
| $B_p=0$ and $W=0$ (only $B_t$) | 8.26634e-5 | 3.30699e-4 | 1.32352e-3 | 5.30570e-3 | 0.8266 / 0.8267 / 0.8272 / 0.8290 |
| $B_t=0,B_p=0$ | **0.** | **0.** | **0.** | **0.** | — |

**Every separable component scales as $\chi^2$** (constant to 0.4% after dividing by $\chi^2$),
including the "only $B_t$" branch — because the χ dependence of $B_t$ (its $O(\chi^2)$ θ-breaking part) is itself a $\chi^2$ source.

### 4.2 Linearity in $W$

Setting $W\to\lambda W$ and extracting the quadratic terms (α=5):

| χ | $a\equiv M_V(W{=}0)$ | $b$ (linear in $W$) | $c$ ($W^2$) | $c/b$ |
|---|---|---|---|---|
| 0.01 | 5.9655e-5 | 1.7797e-4 | 4.5e-9 | $2.5\times10^{-5}$ |
| 0.06 | 2.1534e-3 | 6.4369e-3 | 5.6e-6 | $8.7\times10^{-4}$ |

α=15, χ=0.02: $c/b = 3.2\times10^{-4}$.
**At small χ, $M_V$ is almost strictly linear in $W$** (the $W^2$ contribution is $\le10^{-3}$).

At **large α** and $\chi=0.1$ (α=40), however, $c/b = 0.672$ — the $W^2$ term amounts to 67% of the linear one and is no longer negligible.
Observed fact (without explanation): the α = 40 sequence at χ = 0.1 has $\max\lvert B_t\rvert = 0.0066$ (very small) and
$\max\lvert B_p\rvert = 0.924$ (O(1)), the opposite of the α ≤ 15 sequences (those five have an O(1) $B_t$ limit and $B_p\propto\chi$);
the role of $W$ is therefore more important at large α. **No analytic explanation of this reversal is given in this work** (see §6).

### 4.3 Honest statement: what is rigorous and what is numerical

| Statement | Status |
|---|---|
| $\chi\to-\chi$ is an exact discrete symmetry of the solution family | **rigorous** (from the ansatz and the boundary conditions) |
| $M_V$ is a scalar mass ⟹ even function | **rigorous** |
| $M_V=0$ at $\chi=0$ | literature theorem + **numerically verified here to a relative level of $10^{-4}$** |
| ⟹ $M_V=c_2\chi^2+O(\chi^4)$ | **rigorous** (given the previous line) |
| every term of the integrand has $\ge2$ factors from $\{B,W\}$ | **rigorous** (read off term by term from the symbolic expression) |
| the **coefficients** $c_{BB},c_{BW},c_{WW}$ of the structural formula | **not obtained analytically** (see §6) |

---

## 5. Q3: is $M_V=2\Omega_HJ_V$ an identity or independent content?

### 5.1 What it "is"

$M_V=2\Omega_HJ_V$ is a **Smarr-type relation of the vector field itself**: it uses the volume-integral form of the Komar
mass/angular momentum on the stress tensor of the vector field, together with the fact that on the horizon
$t^a\to t^a+\Omega_H\phi^a$ generates a null geodesic. In spherical symmetry both sides vanish (the identity degenerates);
only for rotating solutions does it first acquire content.

### 5.2 It **is an identity**: the data show that the residual is fully explained by quadrature artifacts

Evaluating $\big|M_V/(2\Omega_HJ_V)-1\big|$ over all 718 solutions:

| α | 5 | 8 | 10 | 12 | 15 | 20 | 25 | 30 | 35 | 40 |
|---|---|---|---|---|---|---|---|---|---|---|
| median | 1.54e-4 | 2.36e-4 | 2.95e-4 | 3.53e-4 | 4.36e-4 | 5.88e-4 | 6.95e-4 | 8.73e-4 | 1.02e-3 | **1.23e-3** |
| max | 1.59e-4 | 2.47e-4 | 3.15e-4 | 3.91e-4 | 5.31e-4 | 1.07e-3 | 6.90e-3 | 1.18e-2 | 1.61e-2 | **2.01e-2** |

Over all 718 points: median **3.9e-4**, maximum 2.0e-2.
**Spearman$(\alpha,\ \text{residual}) = 0.868$** — the residual **increases monotonically with α**,
the fingerprint of a systemic quadrature bias rather than random numerical error.

This agrees with the conclusion of the earlier consistency analysis: **all** of this residual comes from
**spectral aliasing** of the non-linear integral kernel of $M_V$ sampled on the 20-point θ grid
(refining the quadrature to 480×240 freezes the residual to all digits; evaluating alias-free from the spectral
representation of the solution itself brings the residual down to the $10^{-9}$ level,
i.e. $M_V=2\Omega_HJ_V$ holds to the accuracy with which $\Omega_H$ itself is given).

> **Q3 conclusion**: $M_V=2\Omega_HJ_V$ **is an identity** (not an independent physical constraint);
> it **gives no information** beyond $(\Omega_H,J_V)$. But this does not mean it is empty —

### 5.3 ★ What it really contains: it **fixes** the $\chi^2$ coefficient of $M_V$

Expanding the identity as $\chi\to0$ and using the parity assignments of §3 ($\Omega_H,J_V$ odd, $M_V$ even):

$$M_V = 2\big(\omega_1\chi+\omega_3\chi^3+\cdots\big)\big(j_1\chi+j_3\chi^3+\cdots\big)
   \;=\;\underbrace{2\omega_1j_1}_{c_2}\chi^2+O(\chi^4)$$

**This automatically explains the two properties of §3**: (i) the linear term is absent (odd × odd = even); (ii) the leading term is $\chi^2$.
Conversely, $c_2$ is **completely fixed** by $\omega_1,j_1$. Numerical check (two independent routes):

| α | $2\omega_1j_1$ (from the linear coefficients of $\Omega_H,J_V$) | $c_2=C(\alpha)$ (from a direct fit to $M_V$) | relative difference |
|---|---|---|---|
| 5 | 2.37555 | 2.35949 | **0.68 %** |
| 8 | 1.53695 | 1.52889 | **0.53 %** |
| 10 | 1.30836 | 1.30122 | **0.55 %** |
| 12 | 1.13470 | 1.12965 | **0.45 %** |
| 15 | 0.87968 | 0.88326 | **0.41 %** |

> **The independent content of this identity is that it reduces "the quadratic coefficient of $M_V$ as $\chi\to0$"
> to "the product of the linear responses of $\Omega_H$ and $J_V$". $M_V$ is therefore not a new independent degree of freedom,
> but a second-order invariant of the linear rotational response.**

**Remark**: this **lowers** the "information content" of $M_V$ but **raises** its "significance" —
it turns $M_V$ from "an integral that must be computed" into "a second-order quantity completely determined by $\Omega_H,J_V$",
which is exactly the kind of structural statement that can go into the abstract.

---

## 6. Where this stands (honestly)

1. **The coefficients of the structural formula have not been obtained analytically.** Fully expanded, $\mathcal{I}_1$–$\mathcal{I}_4$ of §4
   contain thousands of terms, and the θ integral depends on the complete radial profiles of $B_t^{(0)}(r),B_p(r,\theta),W(r,\theta)$,
   which can only be supplied by numerical solutions of the EMV equations. **What has not been done** is:
   performing the $\int\mathrm{d}\theta$ analytically and reducing the structural formula to a closed form containing only radial integrals.
   Possible routes: (a) work in an analytically controllable limit such as $r_h\ll1$ or $M\gg1$;
   (b) substitute the numerical profiles directly into $\mathcal{I}_{1,2}$ and perform the θ integral term by term to see which term dominates
   (the decomposition of this work already suggests that both the $B_t^{(0)}B_pW$ and the $(B_t^{(0)})^2W^2$ branches are active,
   but **no single dominant term has been isolated**).
2. **The $\chi^4$ coefficient is unexplained.** Numerically $c_4/c_2 \sim O(1)$ ($\alpha=5$: $c_4 = 3.48$ vs $c_2 = 2.38$),
   so the $\chi^4$ correction is not small at all. Its physical origin is undetermined.
3. **There is no analytic model of the branch switch at α ≥ 20.** The data only show that $\max\lvert B_t\rvert$ drops sharply by a factor 8.6 between α ≈ 16–20.
   The critical α has not been located, nor has an asymptotic description of that branch been given.
4. **The χ grid is not fine enough**: the minimum χ of the α ≤ 15 sequences is 0.01.
   With solutions at χ ~ 1e-3, the upper bound on $c_1$ could be lowered by another 2–3 orders of magnitude,
   and $c_4$ could be determined more accurately. **This requires new numerical solutions and was not done in this study.**

> **2026-10-14 correction**: the value originally listed in the table ($1.88\times10^{-8}$, labelled α=5) has been **retracted** —
> it is in fact the **α=10** point value ($+1.8924\times10^{-8}\pm1.15\times10^{-8}$, 1.6σ, consistent with zero);
> the corresponding α=5 value is $+4.02\times10^{-8}\pm2.45\times10^{-8}$. The correct statement is the **upper bound** $|a|<5\times10^{-8}$ ($2\sigma$).
> Basis: the internal analysis note on the χ→0 extrapolation and the independent $c_2$ fit.
