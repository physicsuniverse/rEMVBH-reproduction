# Known-limit unit tests

Six tests that pin the pipeline to limits where the answer is known **analytically**.
They are the credibility backbone of the Python path: a derivation that cannot reproduce
flat space, Schwarzschild, Kerr and an exact circle is not trustworthy anywhere else.

Run all of them, or one at a time:

```bash
python3 run_all.py --only tests          # via the driver, PASS/FAIL table
python3 tests/test_ellipse_circle.py     # or individually (exit code 0 = pass)
```

Requirements: `numpy`, `scipy`, `sympy`. `pytest` is optional — every file is a plain
script with a `__main__` runner, and the functions are named `test_*` so `pytest tests/`
works too if pytest is installed.

| file | limit tested | criterion |
|---|---|---|
| `test_flat_limit.py` | all fields off, $r_H\to0$ | $R_{ab}=0$ identically (symbolic) **and** $<10^{-6}$ (independent finite-difference Ricci) |
| `test_schwarzschild_limit.py` | $Q=B=\chi=0$ | vacuum ansatz is Ricci-flat; $g_{tt}g_{RR}=-1$; $M=2r_H$; $T_H=1/(16\pi r_H)$; $A_H=64\pi r_H^2$; $T_HA_H=4r_H$; Smarr $M=2T_HS$ |
| `test_kerr_limit.py` | $Q=B=0$, spin $\chi$ | Kerr is Ricci-flat; horizon $\Delta(r_+)=0$; extremal at $a=M$ with $T_H\to0$; $g_{t\varphi}\to0$ as $a\to0$; $a\to0$ reproduces Schwarzschild |
| `test_chi_to_zero.py` | $\chi\to0$ | $\mathrm{d}\ln M_V/\mathrm{d}\ln\chi=2$ within 2% for $\alpha=5$–$15$; no linear term ($|c_1|/c_2<10^{-3}$) |
| `test_smarr_exact.py` | exact analytic metrics | $M=2T_HS+2\Omega_HJ+Q\Phi$ to $10^{-13}$ for Schwarzschild, RN, Kerr, including extremal cases |
| `test_ellipse_circle.py` | circle / Kerr shadow | perfect circle $\Rightarrow\varepsilon=0$ to $10^{-14}$; Kerr $a=0.10$ $\Rightarrow\varepsilon=0.0558\%$; $\varepsilon\propto a^2$ |

## The vacuum solution in the paper's ansatz

Several tests rely on the fact that the paper's ansatz

$$ds^2=-fN^2dt^2+\frac{g}{f}\Big[h\big(dr^2+r^2d\theta^2\big)+r^2\sin^2\theta\big(d\varphi-\tfrac{W}{r}(1-N)dt\big)^2\Big],\qquad N=1-\frac{r_H}{r}$$

has the **exact vacuum solution**

$$g=\Big(1+\frac{r_H}{r}\Big)^2,\qquad f=\frac1g,\qquad h=1,\qquad W=0,$$

whose areal radius is $R=r\,g$ and for which $g_{tt}g_{RR}=-1$ exactly. This is *derived*
symbolically in `test_schwarzschild_limit.py` (the Ricci tensor of the ansatz is computed
and required to vanish identically), not assumed. It also fixes the analytic horizon
relations used elsewhere: $g_H=4$, $M=2r_H$, $T_H=1/(16\pi r_H)$, $A_H=64\pi r_H^2$, hence
$T_HA_H=4r_H$ — the identity the paper uses as an accuracy diagnostic.

## Canaries: each test must be able to fail

Every file contains at least one **canary** that deliberately feeds a plausible but wrong
object through the same machinery and asserts that it is rejected:

| canary | what it shows |
|---|---|
| wrong vacuum embedding ($f=h=1-r_H/r$, $g=1$) | fails $g_{tt}g_{RR}=-1$ — the identity has teeth |
| Kerr with $g_{t\varphi}$ normalised 25% too large | fails the Ricci test (note: flipping the *sign* of $g_{t\varphi}$ is just $a\to-a$, i.e. still Kerr — that is why the canary perturbs the magnitude) |
| all-$\alpha$ fit including $\alpha=20$ | gives $p\simeq2.4$–$2.7$ and is rejected by the $\chi^2$ test |
| Smarr formula with the $2\Omega_HJ$ term dropped | fails badly for Kerr $a=0.9$ |
| naive bounding-box ellipticity on a coarse pixel grid | returns $\varepsilon\neq0$ for a *perfect circle* — exactly the historical bug |

## Conventions

* signature $(-,+,+,+)$; quasi-isotropic coordinates $x=1-2r_H/r$;
* metric functions $\{f,g,h,W\}$, Maxwell/vector fields $\{B_t,B_\varphi,A_t,A_\varphi\}$;
* $\theta\in[0,\pi/2]$ is the fundamental domain (reflection symmetry about the equator);
* the ellipticity is the paper's $\varepsilon=H/W-1$ with $W,H$ the equatorial and polar
  diameters of the shadow.

## Data used

Only files already in this package: `solutions/physical_quantities.csv` is *not* used by
the limit tests (they are analytic), except `data/fig1_massofV.csv` for
`test_chi_to_zero.py`. Symbolic Ricci results are cached under `tests/_cache/`
(delete the cache to force a recomputation).
