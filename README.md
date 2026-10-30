# rEMVBH-reproduction

[![reproducibility](https://img.shields.io/badge/run__all.py-PASS%3D46%20%7C%20FAIL%3D0-brightgreen)](#quick-start)
[![python](https://img.shields.io/badge/python-%E2%89%A53.9-blue)](https://www.python.org/)
[![license](https://img.shields.io/badge/license-MIT-blue)](LICENSE)
[![Mathematica](https://img.shields.io/badge/Mathematica-not%20required-lightgrey)](#quick-start)
[![arXiv](https://img.shields.io/badge/arXiv-2610.40179-b31b1b)](https://arxiv.org/abs/2610.40179)

Data, derived equations, verification scripts, unit tests and plotting code for

> **Rotation overcharges vectorized black holes**
> Shengjun Guo, [Guangzai Ye](https://orcid.org/), Peng Liu
> [arXiv:2610.40179](https://arxiv.org/abs/2610.40179) *(to appear in Physical Review D)*

**One command reproduces the entire verification pipeline:**

```bash
git clone https://github.com/physicsuniverse/rEMVBH-reproduction.git
cd rEMVBH-reproduction
pip install numpy scipy sympy
python3 run_all.py          # -> PASS=46  FAIL=0  SKIPPED=0  ERROR=0  (~8 min)
```

---

## What this package is

A **self-contained reproduction package** for the core numerical results of the paper.
It contains:

| Directory | Contents |
|---|---|
| `equations/` | The **derived equations** used by the solver (EOM, second-order action, dispersion relation) |
| `solutions/` | **Solution files**: metric tensors and field profiles, plus the physical-quantity table |
| `data/` | Numerical data behind Figs. 1–4 (mass, existence domain, shadows) |
| `residuals/` | **Solver residual diagnostics**: Smarr/Komar residuals, spectral smoothness |
| `plotting/` | Python plotting code for all figures. `fig4_shadows.py` and `fig10_spectral.py` run self-contained from the shipped data; `fig1`--`fig3` additionally need the raw parameter-scan archive, which is too large to ship (set `EMVBH_ARCHIVE` to its location) |
| `docs/` | Verification notes: why direct substitution into the field equations is ill-conditioned in these coordinates, and how the EOM is verified instead |
| **`run_all.py`** | **One-click driver for the pure-Python verification pipeline** (see below) |
| **`equations/*.py`** | **sympy derivations** that export the equations as lambdified numpy callables |
| **`verify/*.py`** | **numpy/scipy verification** of the archived solutions (Smarr, Komar, `T_H A_H`, EOM residuals) |
| **`tests/`** | **Known-limit unit tests**: flat, Schwarzschild, Kerr, $\chi\to0$, exact Smarr, circle $\Rightarrow\varepsilon=0$ |

> **This package deliberately does NOT contain the spectral solver itself.**
> The solver, its Mathematica implementation, and the solution-generation pipeline are not
> distributed here. What is provided are the **solutions** (as data), the **equations**
> they satisfy, the **residuals** that certify them, and the **code that turns the data
> into the published figures** — i.e. everything needed to *verify and replot* the results,
> not to *regenerate* the solutions from scratch.

---

## Core results reproduced here

### 1. Rotation reverses the spherical undercharging

In the spherical Einstein–Maxwell–Vector (EMV) model, a vectorized black hole can **never**
reach the extremal charge-to-mass ratio. Rotation changes this: the rotating branch crosses
the Kerr–Newman extremal bound.

| Quantity | Value |
|---|---|
| Largest violation found | $q^2+\chi^2 = \mathbf{1.8425}$ |
| Location | $(q,\chi) = (1.3537,\ 0.10)$, `highq_grid_0.csv` line 661 |
| Fraction of solutions beyond the bound | $18.1\%$ of all 72,480 archived ($\mathbf{17.6\%}$ restricting to accuracy-qualified solutions) |
| Extremal KN at $\chi=0.10$ | $\sqrt{1-0.1^2} = 0.99499$ |

⚠️ **The violating region is not a naked-singularity region**: all 13,085 beyond-bound
solutions have $A_H>0$, $T_H>0$ (non-extremal) and $S>0$.

### 2. The vector mass is a *second-order* response to rotation

$$M_V = c_1\chi + c_2\chi^2 + O(\chi^3),\qquad c_1 = 0\ \text{(exact, by parity)}$$

| Quantity | Value | Note |
|---|---|---|
| Exponent (narrow window) | $p = 2.0007 \pm 0.0010$ | $\alpha=5$–$15$, $0.01\le\chi\le0.05$, fitting $M_V/M$ |
| Exponent (wide window) | $p \approx 2.003$ | $0.01\le\chi\le0.08$ |
| $\lvert c_1\rvert/c_2$ bound | $1.3\times10^{-7}$ | two orders below a naive estimate |
| $c_2$ vs. the analytic value $2\omega_1 j_1$ | agrees to $0.112$–$0.192\%$ | $\omega_1,j_1$ fitted independently, **never using $M_V$** |

**The two windows are different quantities and must not be conflated** — see the note in
the paper §III A.

### 3. Rotating hairy solutions exist over a two-dimensional domain

The existence domain in the $(q,\chi)$ plane, and its boundary against the extremal KN curve,
is shown in Figs. 2–3. The data are `data/fig2_lowq.csv` and `data/fig3_highq.csv`.

### 4. Shadow ellipticity carries the hair signature

$$\varepsilon \equiv H/W - 1 = 0.182\,\chi^{1.954}$$

| Contribution to $\varepsilon$ | At $(q,\chi)=(0.9947,0.10)$ |
|---|---|
| Vector hair | $8.30\%$ |
| Electric charge | $0.44\%$ |
| Rotation | $0.06\%$ |
| **Ratio hair/charge** | **19** |

⚠️ The alternative two-digit numbers $8.8\%$ / factor 20 that appear in some earlier
drafts were a **convention mix-up** (total vs. increment). The numbers above are the
self-consistent increment convention and are what the paper uses.

---

## Solutions and their accuracy

The table `solutions/physical_quantities.csv` lists, for every archived solution, the
independent extraction of $(M, J, \chi, \Omega_H, T_H, A_H, S, \tilde T_H, \tilde S)$ from
the metric file, together with the deviation from the header values.

**Accuracy certification** (`residuals/`):

| Diagnostic | Result |
|---|---|
| Smarr relation $M = 2T_HS + 2\Omega_HJ + Q\Phi$ | median relative deviation $1.8\times10^{-11}$ |
| Komar identity $M_V = 2\Omega_H J_V$ | **un-de-aliased**, the stored column gives median $2.4\times10^{-2}$ (chiNofq family $5.8\times10^{-4}$, upofq family $2.7\times10^{-2}$); de-aliasing is a **separate step, verified at 2 points only**, and reaches $\sim10^{-9}$ there |
| Exact identity $T_HA_H = 4r_H$ | median $7.4\times10^{-9}$ for accuracy-qualified solutions; 45% violate it by >1% among the gray (unqualified) ones |
| Spectral smoothness | Chebyshev coefficients decay by **15.6 / 15.2 orders of magnitude** (`residuals/spectral_smoothness.pdf`) |

The last two are used to **justify the colour/gray partition** in Figs. 2–3: the gray
region is not a plotting choice, it is the set where these identities fail.

---

## Quick start

Everything is plain Python 3 + numpy + matplotlib (LaTeX rendering requires a TeX install).

```bash
# 1. Inspect the data
head -3 data/fig1_massofV.csv
head -3 solutions/physical_quantities.csv

# 2. Replot every figure
cd plotting
python3 fig1_massofV.py        # Fig. 1  (M_V landscape + quadratic scaling)
python3 fig2_lowq.py           # Fig. 2  (low-charge existence domain)
python3 fig3_highq.py          # Fig. 3  (high-charge existence domain)
python3 fig4_shadows.py        # Fig. 4  (shadows + ellipticity)
python3 fig10_spectral.py      # spectral smoothness
```

`plotting/figstyle.py` holds the shared style (figure sizes, colour maps, fonts) so that
the output matches the published figures.

## What runs out of the box

`python3 run_all.py` is fully self-contained. The figure scripts differ:

| script | reproducible from the package alone? |
|---|---|
| `plotting/fig4_shadows.py` | **yes** — reads `data/shadow_boundaries/` (16 CSV files shipped) |
| `plotting/fig10_spectral.py` | **yes** — reads `solutions/field_profiles/` |
| `plotting/fig1_massofV.py` | no — needs the raw alpha scan; set `EMVBH_ARCHIVE=/path/to/archive` |
| `plotting/fig2_lowq.py` | no — needs the raw low-$q$ scan; set `EMVBH_ARCHIVE` |
| `plotting/fig3_highq.py` | no — needs the raw high-$q$ scan; set `EMVBH_ARCHIVE` |
| `plotting/fig_solutions.py` | no — needs the raw `.mx` exports; set `EMVBH_ARCHIVE` |

The **derived data** behind every figure *is* shipped (`data/*.csv`), so the figures can be
inspected and re-checked even where the raw archive is not distributed. This is a deliberate
size trade-off, and it is the only respect in which the package is not self-contained.

---

## Conventions

* Units: $M_{\rm ADM}$ is set by the asymptotic falloff; $q\equiv Q/M$, $\chi\equiv J/M^2$.
* The solution files use quasi-isotropic coordinates $x = 1 - 2r_H/r$, with the metric
  functions stored as $f,g,h,W$ and the Maxwell/vector fields as $B_t,B_\varphi,A_t,A_\varphi$.
* $\theta\in[0,\pi/2]$ is the fundamental domain; the solution is reflection-symmetric
  about the equatorial plane ($\theta\to\pi-\theta$).
* $M_V$ is quoted to **at most 3 significant figures**, limited by the $\theta$-grid
  aliasing of the Komar integral (de-aliased values support more digits; see
  `residuals/`).

---

## Python verification pipeline

Everything above can be checked **without Mathematica**, from a single entry point:

```bash
python3 run_all.py            # run everything, print a PASS/FAIL table
python3 run_all.py --only tests   # only the known-limit unit tests
python3 run_all.py --list         # show what is registered / on disk
```

Requirements: `numpy`, `scipy`, `sympy` (Python >= 3.9). No Mathematica anywhere in this path.

### What it runs

```
equations/*.py     symbolic derivation in sympy (metric -> Christoffel -> Riemann -> Ricci
                   -> Einstein; EMV action -> second-order expansion -> mass matrix M^ab;
                   M_V = 2 Omega_H J_V; dispersion relation), exported via sympy.lambdify
        |
        v  the exported equations
verify/*.py        numerical verification in numpy/scipy: independent extraction of
                   M, J, chi, Omega_H, T_H, A_H, S, Phi, Q from the raw field files, then
                   residual tests of Smarr, Komar, T_H A_H = 4 r_H and the EOM
        |
        v
tests/test_*.py    known-limit unit tests (see tests/README.md)
        |
        v
run_all_report.md  PASS/FAIL table: id | script | item | value | criterion | status
run_all_results.json   the same, machine readable
```

### The pass criteria

| test | criterion |
|---|---|
| Smarr relation | median relative residual $<10^{-9}$ |
| Komar identity $M_V=2\Omega_HJ_V$ | median relative residual $<10^{-6}$ |
| $T_HA_H=4r_H$ (exact for this ansatz) | qualified median $<10^{-6}$ |
| EOM residual | median $<10^{-3}$ |
| known-limit unit tests | see `tests/README.md` |

A script that has not been written yet is reported as **SKIPPED** rather than failing the
run: the driver is designed to be usable while the pipeline is still being assembled.
Pass `--strict` if you want a non-zero exit code whenever anything fails.

> **Why known-limit tests matter.** A closed curve is not automatically the right curve and
> a small number is not automatically the right number. Two concrete failures caught this
> way: (i) a hand-rolled "circle" shadow returned $\varepsilon=0.001\%$ instead of $0$, which
> masked a $0.056\%$ Kerr signal; (ii) plots that looked like smooth grey polygons were an
> artefact of convex-hull interpolation rather than data. `tests/test_ellipse_circle.py`
> now pins the first one: *a perfect circle must give $\varepsilon=0$ to machine precision*.

---

## Verification notes

`docs/README_EOM_NOTE.md` explains why substituting an archived solution directly into
$G_{ab}-2T_{ab}$ in the quasi-isotropic coordinates $(t,x,\theta,\varphi)$ returns residuals of
$10^{24}$–$10^{25}$: this is an **ill-conditioning of the coordinate system, not an error of the
solution** (the horizon is a degenerate point of the 2-metric, so $g^{ab}$ diverges, and at
infinity $g_{xx}\sim(1-x)^{-4}$). No finite-difference or spectral method at any resolution can
push that route to $10^{-11}$.

The equations of motion are therefore verified through **coordinate-independent scalar
identities** — the four diagnostics in `verify/verify_eom.py`: the Smarr relation,
$T_HA_H=4r_H$, $\nabla_aT^{ab}=0$ (order of magnitude only), and the algebraic Bianchi identity.
These reach $10^{-10}$, $10^{-11}$ and $10^{-16}$ respectively; see `run_all_report.md`.

---

## Licence / citation

Code: **MIT** (see [LICENSE](LICENSE)). Data and figures: freely usable with attribution.

If you use this package, please cite both the paper and the repository:

```bibtex
@article{Guo:2026overcharge,
  author        = {Guo, Shengjun and Ye, Guangzai and Liu, Peng},
  title         = {Rotation overcharges vectorized black holes},
  year          = {2026},
  eprint        = {2610.40179},
  archivePrefix = {arXiv},
  primaryClass  = {gr-qc}
}

@misc{rEMVBH-reproduction,
  author       = {Guo, Shengjun and Ye, Guangzai and Liu, Peng},
  title        = {{rEMVBH-reproduction}: reproduction package for
                  ``Rotation overcharges vectorized black holes''},
  year         = {2026},
  howpublished = {\url{https://github.com/physicsuniverse/rEMVBH-reproduction}}
}
```

---

## Contact

Peng Liu (corresponding author) — phylp@email.jnu.edu.cn

Questions, bug reports and pull requests are welcome via the GitHub issue tracker.
