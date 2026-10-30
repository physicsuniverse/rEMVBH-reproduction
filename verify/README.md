# Numerical verification (numpy / scipy)

Independent checks that the archived solutions really solve the equations, and that the
physical quantities quoted in the paper can be re-extracted from the raw fields rather
than read off a stored header.

```bash
python3 run_all.py --only verify     # all of them, one PASS/FAIL table
python3 verify/verify_smarr.py       # or individually
```

All scripts write their output under `verify/out/`.

| script | what it does | criterion |
|---|---|---|
| `extract_quantities.py` | extracts $M,J,\chi,\Omega_H,T_H,A_H,S,\Phi,Q$ from the raw field files (`solutions/field_profiles/`, `solutions/metric_files/`) by asymptotic fitting and horizon extrapolation — **never** from the stored header values | writes `out/extracted.csv`; compare against `solutions/physical_quantities.csv` |
| `verify_smarr.py` | checks $M=2T_HS+2\Omega_HJ+Q\Phi$ for every solution with gauge fields | median relative residual $<10^{-9}$ |
| `verify_komar.py` | checks $M_V=2\Omega_HJ_V$ | median relative residual $<10^{-6}$ |
| `verify_identity.py` | checks $T_HA_H=4r_H$, which holds **exactly** for the metric ansatz and therefore measures the solution error directly | qualified median $<10^{-6}$ |
| `verify_eom.py` | substitutes the fields into the equations of motion exported by `equations/eom.py` and reports the pointwise residual | median $<10^{-3}$ |

`_common.py` holds the shared readers and formulas (grid reading, horizon extrapolation,
spectral $\theta$-quadrature for $A_H$, asymptotic fitting for $M,J,Q,\Phi$).

## Cross-checks against the packaged residuals

`residuals/smarr_residuals.csv`, `residuals/komar_residuals.csv` and the `TH_AH_rel_dev`
column of `solutions/physical_quantities.csv` are **independently produced** references.
Each verification script prints its median next to the corresponding reference value, so a
mismatch is visible immediately. Typical values for accuracy-qualified solutions: Smarr
median $1.8\times10^{-11}$, Komar median $1.1\times10^{-9}$ (after de-aliasing), $T_HA_H$
median $7.4\times10^{-9}$.

## Conventions

$(-,+,+,+)$; quasi-isotropic $x=1-2r_H/r$ with $r=2r_H/(1-x)$; $N=1-r_H/r$;
$\theta\in[0,\pi/2]$ fundamental domain; horizon quantities obtained by polynomial
extrapolation in $\delta=1+x$; $A_H=2\pi r_H^2\int_0^\pi\sin\theta\,(g\sqrt h/f)\,d\theta$;
$T_H=f/(2\pi r_H\sqrt{gh})|_H$.
