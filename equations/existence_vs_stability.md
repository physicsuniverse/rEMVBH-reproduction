# Stationary existence ≠ dynamical stability — a three-layer argument and self-attack

**Author**: the paper-reading analysis|**Prerequisite**: the conclusions of the quadratic-action derivation (companion note)
**Project judgement to be tested**: "It is a static solution (a rotating static solution); we have not yet gone as far as the dynamics; the dynamical problem requires a response, but at least it shows that the static solution exists."
**Purpose of this document**: to confirm or refute that judgement. **Conclusions first**:

| Layer | Verdict |
|---|---|
| ① existence and stability are logically independent | ✅ **defensible** (but the wording must be revised, see §1.3) |
| ② the applicability conditions of the no-go do not cover our hairy solutions | ⚠️ **wording must be weakened** (not "strictly inapplicable", but "the premise does not hold ⇒ the conclusion is not established", see §2.4) |
| ③ the tension of "existing but not dynamically reachable" | ✅ **defensible, and worth stating** (supported by precedents, see §3.4) |
| **Overall** | ✅ **the separation argument may be used in the paper, but one must not say "we evade the no-go"** |

---

## 1. First layer: existence and stability are two logically independent matters

### 1.1 The argument

1. **Our equations form an elliptic boundary-value problem.** Under the stationary axisymmetric ansatz ($\partial_t,\partial_\phi$ Killing, discrete symmetry $t\to-t,\ \phi\to-\phi$), the Einstein–Maxwell–vector equations reduce to an **elliptic** PDE system in the two variables $(x,\theta)$, together with boundary conditions (regularity at the horizon, asymptotic behaviour at infinity, regularity on the axis, symmetry about the equator).
2. **An elliptic BVP is well posed**: given the BCs, existence of a solution is demonstrated by the convergence of Newton–Raphson plus residuals $\sim10^{-9}$ (measured); this has nothing to do with "whether the solution is stable under time evolution" — **there is simply no time in a BVP.**
3. **A stationary solution = a critical point of the field equations** ($\delta S=0$ restricted to the space of stationary axisymmetric configurations).
4. **A ghost/gradient instability is a property of the "Hessian/fluctuation operator at that critical point"**: it lives in the **time-dependent (hyperbolic) sector** and is determined by the sign of the kinetic term of the quadratic action $\delta^2S$ ($X_1,X_2$ of the quadratic-action derivation).
5. ⇒ **"$u_0$ satisfies $E[u_0]=0$" and "$\delta^2S|_{u_0}$ has no ghost" are two propositions about different objects**: the former concerns a **point**, the latter the **spectrum at that point**.

**Analogy (useful for the writing)**: a cubic potential $V(\phi)=\phi^3$ has a Hamiltonian that is unbounded below, yet its critical point (the saddle at $\phi=0$) is still a legitimate solution of the equation $\phi^2=0$; "the solution exists" never requires the action to be bounded below.

### 1.2 Self-attack (the place where this layer is most likely to be broken)

**Attack A: ghost ⇒ the action is unbounded below ⇒ the configuration is not a legitimate state ⇒ the "solution" is meaningless.**
- **Partly valid.** The valid part is: unbounded below ⇒ it cannot be quantized, nor treated as a stable final state. The invalid part is the logical jump: **"unbounded below" and "a critical point exists" can coexist** (see above). So what a referee can destroy is the **physical interpretation**, not the **existence**.
- **Reinforcement**: we **do not claim** that these solutions are "final states/astrophysically realisable black holes". What we claim is: **they are stationary solutions of this EFT**, and their dynamical viability is an independent question.

**Attack B: your ghost is "background-induced", not a genus defect of the fundamental Lagrangian — doesn't that amount to saying the problem is not serious?**
- **This is a temptation and must be refused.** Our fundamental kinetic term is indeed canonical ($-\frac14\delta V^2$), and the ghost appears only on backgrounds with $F^2\neq0$ ($m_{\rm eff}^2=\frac12\alpha UF^2<0$). But [2204.10868](https://arxiv.org/abs/2204.10868) has already shown, using a **time evolution on Kerr**, that **the instability is indeed triggered** when the self-interaction becomes strong ("instabilities are triggered in a finite time once the self-interaction becomes significant"), and states explicitly: "**It has been suggested that nonperturbative dynamics may drive systems away from such instabilities. We demonstrate that this is not the case**".
- ⇒ **Any phrasing such as "background-induced, hence harmless" or "non-linearities can cure it" is wrong** and will be caught immediately.

**Attack C: could your Newton solutions be mere numerical artifacts?**
- Our solutions agree to all digits with the independent archive (relative difference in $M$ of $2.4\times10^{-11}$), satisfy $\Delta_{\rm Smarr}\le2.6\times10^{-10}$, $J=\chi M^2$ and $Q=qM$ to $10^{-10}$, and the 100×50 and 60×20 grids agree. ⇒ The attack fails.

### 1.3 Verdict and recommended wording

**Defensible**, but it must be phrased as "**we restrict the problem to the stationary sector**", not as "the ghost does not matter":

> ✅ Acceptable: "We construct stationary solutions; we do not analyse their dynamical stability. Existence of a stationary solution and the absence of ghosts in its fluctuation operator are logically independent statements: the former is a property of the elliptic boundary-value problem we solve, the latter of the hyperbolic sector."
> ❌ Not acceptable: "the ghost is harmless", "non-linearities may cure it", "the no-go only concerns dynamics so our solutions are unaffected".

---

## 2. Second layer: the applicability conditions of the no-go do not cover our **hairy** background

### 2.1 Facts (from the quadratic-action derivation; directly citable)

| Fact | Source |
|---|---|
| [2605.13920](https://arxiv.org/abs/2605.13920) studies "linear stability of a massive vector field $A_\mu$ **about a GR background … with vanishing $A_\mu$**" | their §II |
| In their most general Lagrangian (1), $\mathcal V^{\mu\nu}=\frac12\mu^2g^{\mu\nu}-\beta G^{\mu\nu}$ contains **only $g$ and $G$** | their (1), (2) |
| Their dispersion relations rely on "$\mathcal G^{\mu\nu\rho\sigma}$ has **no off-diagonal terms in the time-like direction**, thus leading to a purely quadratic dispersion relation" | their §II |
| Our quadratic action for an **arbitrary (hairy) background** gives $M^{ab}=m_{\rm eff}^2(g^{ab}+2\alpha\bar B^a\bar B^b)$, **which contains $\bar B^a\bar B^b$ and is not covered by their Eq. (1)** | quadratic-action derivation, Eq. (2.1) |
| The content of their theorem: Eq. (7) holding everywhere ⇒ tachyonically stable; **dynamical reachability of perturbations linearized about the hairless state** | their §IV |

### 2.2 How strong is "outside their framework" really? Two sub-questions

**(a) Might the physical conclusion be the same anyway?**
- **Our calculation says: very likely yes.** Because we find $\operatorname{sign}(X_2)=\operatorname{sign}(m_{\rm eff}^2)<0$ (quadratic-action derivation, Eq. (3.3)), i.e. "naively transplanting their criterion" still fails.
- ⇒ **One cannot claim "because we are outside their framework, we are safe".** This is the easiest step at which to fool oneself.

**(b) What, then, is left of "outside their framework"?**
- What is left is this: **the premises of their derivation do not hold on our hairy background ⇒ their conclusion is "not established" in our case**, not "refuted".
- Concretely, $\bar B^a\bar B^b$ generates non-diagonal "timelike" components of $\mathcal V$ ($\mathcal V^{\mathbf{13}},\mathcal V^{\mathbf{14}},\mathcal V^{\mathbf{23}},\mathcal V^{\mathbf{24}}$). When these components are non-zero, the dispersion relation acquires a term **linear in $\omega$**, the definition of a "ghost" itself must be replaced by a more general criterion (the framework of Babichev et al. 2018 cited by them), and the conclusion is **unknown**.
- **The numerical evaluation is currently computing these four components** (see the numerical-criteria interface of the quadratic-action analysis, §1.4). This is the only empirical support for this line of defence.

### 2.3 Self-attack

**Attack D: "Your solution has a hairy background, but the whole motivation of the no-go is 'how does the hair grow'; you want to dodge by saying 'I only have the solution and say nothing about hair growth' — that is evasive."**
- **Partly valid.** Response: we **indeed say nothing about hair growth**, and we state it as a limitation rather than hiding it in a footnote. At the same time we **state what should be checked next** (the $M^{ab}$ criteria + the non-diagonal components), which is more honest than evasion.
- There is also a genuine theoretical distinction worth stating: **the no-go conclusion concerns the "linearized tachyon channel"**, whereas our solutions are **non-linear stationary states** obtained in one shot. "A non-linear stationary state exists" and "the linear channel is clean" are indeed two different things — this is still layer ①, but it is appropriate to use it here.

**Attack E: "Your $\bar B\bar B$ term is non-zero on a hairy background, but in the spherically symmetric limit it yields exactly the ghost ($X_1=X_2=\frac12m_{\rm eff}^2<0$, in agreement with [2310.18399](https://arxiv.org/abs/2310.18399)); what makes rotation change the conclusion?"**
- **Valid, and this is the core of the debate.** We **cannot** guarantee that rotation changes the conclusion; all we can say is that "the premises of the criterion may fail". **So the wording must be "not established", not "refuted".**

### 2.4 Verdict and recommended wording

**The wording must be weakened.** Recommended:

> ✅ "The no-go theorem of Ref. [2605.13920] is derived for linearized perturbations of a **hairless** background, for which the most general fluctuation operator has the form of their Eq. (1). Our hairy backgrounds generate an additional term $\propto\bar B^\mu\bar B^\nu$ that lies outside that class; consequently the theorem does not establish a ghost instability for these backgrounds, and the question must be settled by the generalized analysis we outline."
> ❌ "Our solutions evade the no-go." (the inference does not follow)
> ❌ "The no-go applies only to the dynamical formation, not to stationary solutions." (the first half is right, the second wrong: a stationary solution is itself the background of the fluctuation analysis)

---

## 3. Third layer: the tension of "existing but not dynamically reachable"

### 3.1 Physical content
- The no-go line of argument says: a "clean" tachyon growth starting from a hairless state is impossible; our solutions are **non-linear stationary states**, so **one cannot claim that they are dynamically reachable through the linear instability**.
- This is not at all unusual in physics: **non-attractor stationary states** (unstable equilibria, negative-mode black holes, small-AdS-black-hole branches, thermodynamically unstable branches) are ubiquitous in GR, and their stationary solutions are nevertheless constructed and studied.

### 3.2 Precedents (see the literature-precedent survey, companion note)
- **[1403.2757](https://arxiv.org/abs/1403.2757)** (Herdeiro & Radu, PRL 2014): only existence plus physical properties ($J/M^2>1$, quadrupole moment, ISCO), **with no treatment of dynamical stability** — this paper established the paradigm of "construct the stationary state first, discuss stability separately".
- **[2004.09536](https://arxiv.org/abs/2004.09536)**: the domain of stationary solutions and phase space; its abstract explicitly lists "**fully non-linear dynamical evolutions** … have shown the formation of a new equilibrium state" as **someone else's** work (i.e. the dynamical verification is subsequent, independent work).
- **[2604.04568](https://arxiv.org/abs/2604.04568)** (Kleihaus & Kunz, PLB 2026): rotating vectorized BHs, treating only the domain of existence and thermodynamics, **published after all three no-go papers**, and still without any treatment of the dynamics.

### 3.3 Self-attack
**Attack F: "The referee will ask: how does this solution form?"**
- **Valid, but answerable**: (i) we do not claim a formation mechanism; (ii) the existence of a stationary solution is an independent mathematical fact and a candidate endpoint of **any** formation scenario; (iii) we state it as an explicit limitation + agenda.
**Attack G: "If your solution can never be reached, what physical significance does it have?"**
- Three answers (in decreasing order of strength): ① it is **a solution that the EFT must acknowledge** (every UV completion contains it), and its domain of existence/phase diagram is **the target region for dynamical studies**; ② the two headline results (rotation reversing the undercharging, $M_V=2\Omega_HJ_V$) **do not depend on stability at all**; ③ **"can the mechanism survive" is itself the hottest question in the field at present** (all three no-go papers contest exactly this), and our stationary solutions + criteria are precisely an input to that debate.

### 3.4 Verdict
**Defensible, and worth stating.** But three conditions must be met: state the limitation explicitly; do not imply dynamical reachability; provide the agenda.

---

## 4. One-sentence summary (for the project record)

**The separation argument is usable; its core logic is**: *we solve an elliptic boundary-value problem (a stationary critical point), whereas the ghost criterion concerns the fluctuation spectrum of the hyperbolic sector; the objects are different, hence the two are logically independent.* **But it protects only "existence" and the two results that are stability-independent; it does not protect the proposition of "clean formation", which we have never claimed.**
**Three sentences must be avoided**: ① "the ghost is harmless"; ② "non-linearities can cure it"; ③ "we evade the no-go".
