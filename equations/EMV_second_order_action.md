# The EMV quadratic action and the ghost/gradient criteria (step-by-step derivation)

**Author**: the paper-reading analysis|**Companion notes** (separate files): the literature-criteria excerpt, the numerical-criteria interface, the worst-case assessment
**Conventions**: signature $(-,+,+,+)$; $S=\frac{1}{16\pi}\int\sqrt{-g}\left[R-U(|B|^2)F_{ab}F^{ab}-V_{ab}V^{ab}\right]$, $U=e^{\alpha|B|^2}$, $V_{ab}=\partial_aB_b-\partial_bB_a$, $\alpha>0$, $F^2\equiv F_{ab}F^{ab}$ (electrically dominated: $F^2<0$).

---

## 1. Background field equations and the definition of the "effective mass"

Varying with respect to $B_b$:
$$\delta S_B=\int\sqrt{-g}\left[-2V^{ab}\delta V_{ab}-U'(|B|^2)\,F^2\,2B_b\delta B^b\right],\qquad U'=\frac{dU}{d(|B|^2)}=\alpha U .$$
Integrating by parts, $\int\sqrt{-g}\,(-2V^{ab}\delta V_{ab})=\int\sqrt{-g}\,(-4V^{ab}\partial_a\delta B_b)\to 4\int\sqrt{-g}\,\nabla_aV^{ab}\,\delta B_b$. Setting the total variation to zero:
$$\boxed{\ \nabla_aV^{ab}=\tfrac12U'(|B|^2)F^2B^b=\alpha U F^2B^b/2\ }\tag{1.1}$$
**Define the effective mass squared** (same dimensions and same sign as $m_{\rm eff}^2=\omega^2|_{k,\ell=0}$ of [2605.13920](https://arxiv.org/abs/2605.13920)):
$$\boxed{\ m_{\rm eff}^2\equiv\tfrac12U'F^2=\tfrac12\alpha\,U(|B|^2)\,F_{ab}F^{ab}\ }\tag{1.2}$$
- Pure-Proca comparison: $\mathcal L=-\frac14F^2-\frac12\mu^2A^2\Rightarrow\nabla_aF^{ab}=\mu^2A^b$, so the coefficient on the right-hand side of (1.1) is precisely $m_{\rm eff}^2$ ✓.
- **Electrically dominated background ($F^2<0$) with $\alpha>0$ ⇒ $m_{\rm eff}^2<0$**: this is exactly the vectorization driver (a negative effective mass squared).

---

## 2. The quadratic action ($B=\bar B+\delta B$, arbitrary background $\bar B$)

Expanding term by term:

**(a) Kinetic term** $ -V_{ab}V^{ab}=-(\bar V+\delta V)^2=-\bar V^2-2\bar V^{ab}\delta V_{ab}-\delta V_{ab}\delta V^{ab}$,
and $2\bar V^{ab}\delta V_{ab}=4\bar V^{ab}\partial_a\delta B_b$.

**(b) Coupling term** Write $B^2\equiv B_aB^a$. For $U(B^2)=e^{\alpha B^2}$ at $\bar B$:
$$U=\bar U+\bar U'\left(2\bar B\!\cdot\!\delta B+\delta B^2\right)+\tfrac12\bar U''\left(2\bar B\!\cdot\!\delta B\right)^2+\cdots,\qquad \bar U'=\alpha\bar U,\quad \bar U''=\alpha^2\bar U .$$
$$\Rightarrow\ -UF^2=-\bar U\bar F^2-\bar F^2\left[2\bar U'(\bar B\!\cdot\!\delta B)+\bar U'\delta B^2+2\bar U''(\bar B\!\cdot\!\delta B)^2\right]+\cdots$$

**(c) The linear terms cancel each other** (a necessary consequence of the background field equation (1.1)):
$$-4\bar V^{ab}\partial_a\delta B_b\to 4\nabla_a\bar V^{ab}\delta B_b=4\cdot\tfrac12\bar U'\bar F^2(\bar B\!\cdot\!\delta B)=2\bar U'\bar F^2(\bar B\!\cdot\!\delta B),$$
$$-\bar F^2\cdot2\bar U'(\bar B\!\cdot\!\delta B)=-2\bar U'\bar F^2(\bar B\!\cdot\!\delta B)\ \Rightarrow\ \text{their sum vanishes}\ \checkmark$$

**(d) The quadratic action** (writing the kinetic term in the standard normalization $-\frac14\delta V^2$, i.e. setting $\delta B_{\rm can}=2\,\delta B$):
$$S_2=\int\sqrt{-g}\left[-\tfrac14\,\delta V_{ab}\delta V^{ab}-\tfrac12\,M^{ab}\,\delta B_a\delta B_b\right],$$
$$\boxed{\ M^{ab}=m_{\rm eff}^2\left(g^{ab}+2\alpha\,\bar B^a\bar B^b\right)\ }\tag{2.1}$$
(Derivation: the mass term is $=-\bar F^2\bar U'\left[\delta B^2+2\alpha(\bar B\!\cdot\!\delta B)^2\right]=-\tfrac14\bar F^2\bar U'\left[\delta B_{\rm can}^2+2\alpha(\bar B\!\cdot\!\delta B_{\rm can})^2\right]$; comparing with $-\frac12M^{ab}\delta B_a\delta B_b$ gives (2.1).)

**Self-check against the linearized field equations**: linearizing (1.1) directly, the linear part of $U'(B^2)B^b$ is $=\bar U'\bar B^b+\alpha\bar U\left[\delta B^b+2(1{+}\alpha)(\bar B\!\cdot\!\delta B)\bar B^b\right]$; subtracting the background equation gives
$$\nabla_a\delta V^{ab}=m_{\rm eff}^2\left[\delta B^b+2\alpha(\bar B\!\cdot\!\delta B)\bar B^b\right],$$
in complete agreement with the $M^{ab}$ of (2.1) ✓ (**note**: the term $\bar U'\bar B^b$ must be subtracted, otherwise a spurious $(1+\alpha)$ factor appears — this is the easiest place to make a mistake in this derivation).

### Key structure
$$M^{ab}=m_{\rm eff}^2\left(g^{ab}+2\alpha\bar B^a\bar B^b\right)$$
- Eigenvalue along $\bar B$: $m_{\rm eff}^2(1+2\alpha\bar B^2)$; transverse directions (3-fold): $m_{\rm eff}^2$.
- When $\bar B$ is timelike ($\bar B^2<0$) and $2\alpha|\bar B^2|>1$, **the eigenvalue along $\bar B$ flips to positive** — this is the only analytic difference between the rotating and the spherically symmetric solutions.
- **The term $\bar B^a\bar B^b$ is absent from Eq. (1) of [2605.13920](https://arxiv.org/abs/2605.13920)** (their $\mathcal V^{\mu\nu}=\frac12\mu^2g^{\mu\nu}-\beta G^{\mu\nu}$ allows only $g$ and $G$). Because they linearize only a "hairless background" ($A_\mu=0$), no background vector is available to span $\bar B^\mu\bar B^\nu$. **This is our structural loophole.**

---

## 3. Mapping onto the criteria of 2605.13920

Their (1) and (3): $\mathcal L=F\mathcal G F-A\mathcal V A$, $X_1\equiv-\mathcal V^{\mathbf{12}}$, $X_2\equiv\mathcal V^{\mathbf{34}}$, $F_1\equiv4\mathcal G^{\mathbf{1212}}$, etc., in the Newman–Penrose convention $g_{\mathbf{34}}=-g_{\mathbf{12}}=1$.

Our kinetic term is precisely $\mathcal G^{\mu\nu\rho\sigma}=-\frac14g^{\mu\rho}g^{\nu\sigma}$ (their $\alpha_{\rm curv}=\beta=0$) ⇒
$$F_1=F_2=1>0,\qquad F_3=F_4=0 .$$
Our mass term: $\mathcal L\supset-\frac12M^{ab}\delta B_a\delta B_b=-A_\mu\mathcal V^{\mu\nu}A_\nu$ ⇒ $\mathcal V^{\mu\nu}=\frac12M^{\mu\nu}$. Substituting (2.1):
$$\boxed{\ X_1=-\mathcal V^{\mathbf{12}}=\tfrac12m_{\rm eff}^2\left(1-2\alpha\bar B^{\mathbf 1}\bar B^{\mathbf 2}\right),\qquad X_2=\mathcal V^{\mathbf{34}}=\tfrac12m_{\rm eff}^2\left(1+2\alpha\bar B^{\mathbf 3}\bar B^{\mathbf 4}\right)\ }\tag{3.1}$$
**Sign self-check (Proca limit)**: for $\bar B=0$ and $m_{\rm eff}^2\to\mu^2>0$, $X_1=X_2=\frac12\mu^2>0$ ✓, exactly their positive value in the Proca case.

Their **Eq. (7)**: $F_1>0,\ F_2>0,\ F_3\ge0,\ X_1>0,\ X_2\ge0$ ⇒ tachyonically stable (i.e. **no** vectorization). Here the $F$ conditions are automatic, so the criterion reduces to
$$\boxed{\ \text{no ghost}+\text{no gradient instability}\iff X_1>0\ \text{and}\ X_2\ge0\ }\tag{3.2}$$

### 3.1 Immediate analytic conclusion (**hairless background, $\bar B=0$**)
$$X_1=X_2=\tfrac12 m_{\rm eff}^2=\tfrac14\alpha F^2 .$$
With $\alpha>0$ and electrical dominance $F^2<0$ ⇒ $X_1=X_2<0$ ⇒ **(3.2) is violated**.
⇒ **"hairless KN background + EMV coupling" falls exactly within the scope of their theorem, and the conclusion is: the tachyon (the vectorization driver) necessarily comes together with a ghost or a gradient instability.** This agrees completely with the spherically symmetric conclusion of [2310.18399](https://arxiv.org/abs/2310.18399) (our mapping passes this consistency check).

### 3.2 Hairy, rotating background ($\bar B\neq0$)
$$\operatorname{sign}(X_2)=\operatorname{sign}\!\left[m_{\rm eff}^2\left(1+2\alpha\bar B^{\mathbf3}\bar B^{\mathbf4}\right)\right].$$
In the Newman–Penrose convention $\bar B^{\mathbf3}\bar B^{\mathbf4}=(\bar B\!\cdot\! m)(\bar B\!\cdot\!\bar m)=|\bar B\!\cdot\! m|^2\ge0$ (real vector + complex null tetrad). Hence
$$\boxed{\ \operatorname{sign}(X_2)=\operatorname{sign}(m_{\rm eff}^2)=\operatorname{sign}(\alpha F^2)\ }\tag{3.3}$$
⇒ **whenever the tachyon region exists ($\alpha F^2<0$), $X_2<0$ necessarily holds**, independently of rotation and of the magnitude of $\alpha$. **In the sense of (3.1), the "allowed region with $X_2\ge0$" is empty.**

$X_1$ can be rescued by rotation: $X_1>0\iff 2\alpha\bar B^{\mathbf1}\bar B^{\mathbf2}>1$ (since $m_{\rm eff}^2<0$). This is the **only** quantity that can be flipped. But note the structure of their theorem: **if $X_1>0$ while $X_2<0$, the longitudinal gradient condition $X_2/X_1\ge0$ necessarily fails** ⇒ at large $\ell$, $\omega_L^2=k^2+(X_2/X_1)\ell^2+X_2<0$ ⇒ **gradient instability**.
⇒ Therefore: **within the framework of (3.1)/(3.2), "tachyon + no ghost + no gradient instability" are mutually exclusive (the allowed region is empty): one gets either a ghost or a gradient instability.**

### 3.3 Explicit dispersion relations (the specialization $F_3=F_4=0$, convenient for cross-checks on the numerical side)
Substituting $F_1=F_2=1,\ F_3=F_4=0$ into their (5) and (6):
$$\omega_L^2=k^2+\frac{X_2}{X_1}\ell^2+X_2,\qquad \mathcal R_L=\frac{\ell^2}{X_1}+1,$$
$$\omega_-^2=k^2+X_2,\qquad \omega_+^2=k^2+X_1+\ell^2 .$$
(The third relation holds in the branch $\ell^2>X_2-X_1$; that branch is precisely the large-$\ell$ branch when $X_2<0<X_1$.)
**How to read these**: the constant term of $\omega_-^2$ is $X_2$ — $X_2<0$ gives a tachyon; the $\ell^2$ coefficient of $\mathcal R_L$ is $1/X_1$ — $X_1<0$ gives a ghost; the $\ell^2$ coefficient of $\omega_L^2$ is $X_2/X_1$ — opposite signs of the two give a gradient instability.

---

## 4. Three escape routes (listed honestly; only the first can be tested by us)

| # | Escape route | Basis | What we can do |
|---|---|---|---|
| **1** | **The $\bar B^a\bar B^b$ term of $M^{ab}$ is not contained in their Eq. (1)**; moreover this term generates **non-diagonal "timelike" components** in $\mathcal V^{\mu\nu}$ ($\mathcal V^{\mathbf{13}},\mathcal V^{\mathbf{14}},\mathcal V^{\mathbf{23}},\mathcal V^{\mathbf{24}}$), destroying the premise "no off-diagonal terms in the time-like direction" on which their dispersion relations rest | Their text states explicitly: "This structure guarantees that $\mathcal G^{\mu\nu\rho\sigma}$ has no off-diagonal terms in the time-like direction, thus leading to a purely quadratic dispersion relation"; and "Whether they are also necessary conditions **remains an open question**" | ✅ **This is precisely what the numerical side must compute** (see the numerical-criteria interface, §1.4): if these components are significantly non-zero ⇒ criterion (3.2) **does not apply**, and a full dispersion analysis including the term linear in $\omega$ is required |
| 2 | Their assumption that "the theory must not contain extra degrees of freedom beyond spin 1" concerns the most general linear theory after maximizing over a **hairless background**; the expansion of $U(|B|^2)F^2$ in EMV contains infinitely many operators such as $\alpha^2|B|^4F^2$ | Our derivation of (2.1) shows that on a **hairy background** the perturbations do acquire extra structure | This is the same thing as #1 seen from the other side; one cannot escape by "having more operators" per se (the linear order is unaffected), but one **can** escape through the $\bar B\bar B$ term |
| 3 | A "universal time direction" (a direction in which the Hamiltonian can be defined) may not exist for rotating solutions | They themselves admit that their conditions are sufficient but not necessary | ❌ We cannot test this; it can only be recorded as a caveat |

**It must be emphasised**: escape route 1 is not "an escape we can take at will". The sign of $X_2$ **does not depend** on any component of $\bar B$ ((3.3)), so even if the $\bar B\bar B$ term is large, $X_2<0$ still holds; the only change is that criterion (3.2) itself becomes invalid and a complete dispersion analysis (including the term linear in $\omega$) must be redone — a **new and harder** calculation, about 1–2 weeks of work, with an unknown outcome.

---

## 5. Interface with our paper and data (an executable form of the numerical evaluation)

See the numerical-criteria interface note. The three core quantities are
$$\texttt{F2}=F_{ab}F^{ab},\quad \texttt{B2}=B_aB^a,\quad \texttt{meff2}=\tfrac12\alpha e^{\alpha\,\texttt{B2}}\,\texttt{F2};$$
$$\texttt{X1}=\tfrac12\texttt{meff2}\left(1-2\alpha\bar B^{\mathbf1}\bar B^{\mathbf2}\right),\qquad \texttt{X2}=\tfrac12\texttt{meff2}\left(1+2\alpha\bar B^{\mathbf3}\bar B^{\mathbf4}\right).$$
**Unit test** ($\chi\to0$): one must obtain $X_1<0$ and $X_2<0$ — this is exactly 2310.18399; if $X_1>0$ is obtained, first re-check the sign conventions.

---

## 6. Conclusions (what this derivation does and does not provide)

**It provides (directly usable in the paper or in follow-up work)**:
1. (2.1): the quadratic action of the vector perturbation in EMV around an **arbitrary** background, $M^{ab}=m_{\rm eff}^2(g^{ab}+2\alpha\bar B^a\bar B^b)$.
2. (3.1), (3.3): the **analytic conclusions** after mapping onto the literature criteria — the sign of $X_2$ equals $\operatorname{sign}(\alpha F^2)$, so the **tachyon region and the "no-ghost criterion" are mutually exclusive**.
3. Agreement with the spherically symmetric conclusion of [2310.18399](https://arxiv.org/abs/2310.18399) (an independent check that the mapping is correct).
4. **An executable form of the criteria**, delivered to the numerical evaluation.

**It does not provide (stated honestly)**:
1. whether escape route 1 actually holds (this requires the numerical non-diagonal components of $\mathcal V$ plus a complete dispersion analysis);
2. whether the ghost region can be "crossed" at the non-linear level ([2204.10868](https://arxiv.org/abs/2204.10868) concludes that it **cannot**; they performed a time evolution on Kerr);
3. we have **not** independently reproduced Eqs. (5) and (6) of 2605.13920 (that would require the derivation of their Appendix A); this document merely **transplants their conclusions onto EMV** and uses the conditions they state; if their appendix contains additional assumptions, our mapping must be tightened accordingly.
