# Residual diagnostics

These files certify that the solutions in `solutions/` actually solve the equations, and
that the numerical resolution is sufficient for the quoted precision.

| File | Diagnostic |
|---|---|
| `smarr_residuals.csv` | Relative deviation from the Smarr relation $M = 2T_HS + 2\Omega_HJ + Q\Phi$, per solution |
| `komar_residuals.csv` | Relative deviation from the Komar identity $M_V = 2\Omega_H J_V$ — **un-de-aliased** (median $2.4\times10^{-2}$); de-aliased values were obtained for 2 points only ($\sim10^{-9}$) |
| `spectral_smoothness.pdf` | Chebyshev coefficient decay for the spectral expansion of the fields — evidence that the solutions are analytically smooth, and that the plotted surfaces are not interpolated |

## Why these matter

Two of these are used in the paper as **accuracy criteria**, not just as checks:

* $T_HA_H = 4r_H$ holds **exactly** for the metric ansatz used. Its violation therefore
  measures the solution's error directly. Solutions violating it by more than 1% are shown
  in gray in Figs. 2–3; the colour/gray partition is derived from this, not chosen.
* $\Delta_{\rm Smarr} < 10^{-4}$ is the acceptance criterion for a solution to be quoted.

Typical values for accuracy-qualified solutions: Smarr median $1.8\times10^{-11}$,
**Komar (important):** the values in `komar_residuals.csv` are **NOT de-aliased**. Over its 728 rows the median relative deviation is $2.4\times10^{-2}$ (chiNofq family $5.8\times10^{-4}$, upofq family $2.7\times10^{-2}$). De-aliasing is a **separate step that was verified at 2 points only** and reaches $\sim10^{-9}$ there. The stored `rel_residual` is normalised by $|2\Omega_HJ_V|$, **not** by $|M_V|$. $T_HA_H$ median (qualified solutions) is $7.4\times10^{-9}$.
