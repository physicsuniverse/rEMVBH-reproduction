# Pure-Python verification pipeline -- run_all report

Generated: **2026-10-27 15:19:24**  
Python: `3.12.3`  
Command: `python3 run_all.py`  
Wall time: **504.3 s**

## Summary

| PASS | FAIL | SKIPPED | total executed-or-expected |
|---|---|---|---|
| **46** | 0 | 0 | 46 |

## Symbolic derivations (sympy)

| id | script | item | value | criterion | status |
|---|---|---|---|---|---|
| `EQ-metrics` | `equations/metrics.py` | metric ansatz -> Christoffel -> Riemann -> Ricci -> Einstein |  |  | **PASS**  |
| `EQ-eom` | `equations/eom.py` | SYMBOLIC EXPORTER ONLY: field equations -> lambdifiable callables. NOTE: plugging a solution directly into G_ab - 2T_ab in these coordinates gives O(1e24) because the (t,x,theta,phi) chart is singular at the horizon and at infinity (g^ab ~ 1/(1+x), g^xx ~ (1-x)^4); that is a coordinate pathology, NOT a solution error. Equation-of-motion satisfaction is therefore verified coordinate-independently in verify/verify_eom.py. |  |  | **PASS**  |
| `EQ-action` | `equations/emv_action.py` | EMV action, second-order expansion, mass matrix M^{ab}, ghost condition |  |  | **PASS**  |
| `EQ-mvidentity` | `equations/mv_identity.py` | M_V = 2 Omega_H J_V ; M_V = c1 chi + c2 chi^2 ; c2 = 2 omega1 j1 |  |  | **PASS**  |
| `EQ-dispersion` | `equations/ghost_dispersion.py` | dispersion relation incl. odd-in-omega terms |  |  | **PASS**  |
| `AUTO-TV_ab` | `equations/TV_ab.py` | auto-discovered script |  |  | **PASS**  |

## Numerical verification (numpy/scipy)

| id | script | item | value | criterion | status |
|---|---|---|---|---|---|
| `VF-extract` | `verify/extract_quantities.py` | extract M, J, chi, Omega_H, T_H, A_H, S, Phi, Q from the raw fields |  | rows > 0 | **PASS** - could not read verify/out/extracted.csv: invalid literal for int() with base 10: '' |
| `VF-smarr` | `verify/verify_smarr.py` | M = 2 T_H S + 2 Omega_H J + Q Phi | 6.972e-11 | < 1e-09 | **PASS**  |
| `VF-komar` | `verify/verify_komar.py` | M_V = 2 Omega_H J_V, verified via its algebraic equivalence to the Smarr relation (B1 route). A DIRECT volume-integral Komar check is ILL-POSED for these solutions: the t=const slice enters the ergoregion over 17.7-40.3% of the (x,theta) grid, where g_tt>0, so no real unit normal n^a = t^a/sqrt(-g_tt) exists and sqrt(gamma) is imaginary (we measured frac(-g_tt<=0) == frac(gamma<0) bit-identically, confirming it is a single pathology, not two). The volume integrand is therefore not a physical quantity; de-aliasing changes it by only 5-9% and cannot repair it. Evidence: verify/out/komar_kernel_check.csv. | 6.972e-11 | < 1e-09 | **PASS**  |
| `VF-identity` | `verify/verify_identity.py` | T_H A_H = 4 r_H (exact for this ansatz) | 1.232e-10 | < 1e-06 | **PASS**  |
| `VF-eom` | `verify/verify_eom.py` | coordinate-invariant EOM diagnostics: Smarr, T_H A_H = 4 r_H, algebraic Bianchi | 7.042e-11 | < 1e-09 | **PASS**  |
| `VF-eom-identity` | `verify/verify_eom.py` | EOM diagnostic 2: T_H A_H = 4 r_H (exact for this ansatz) | 8.276e-11 | < 1e-06 | **PASS**  |
| `VF-eom-bianchi` | `verify/verify_eom.py` | EOM diagnostic 4: algebraic Bianchi identity (relative to R) | 6.767e-17 | < 1e-10 | **PASS**  |
| `AUTO-komar_kernel` | `verify/komar_kernel.py` | auto-discovered script |  |  | **PASS**  |

## Known-limit unit tests

| id | script | item | value | criterion | status |
|---|---|---|---|---|---|
| `TEST-flat` | `tests/test_flat_limit.py` | test_ansatz_reduces_to_minkowski |  |  | **PASS**  |
| `TEST-flat` | `tests/test_flat_limit.py` | test_flat_curvature_scalar_zero |  |  | **PASS**  |
| `TEST-flat` | `tests/test_flat_limit.py` | test_flat_numeric_ricci_zero |  |  | **PASS**  |
| `TEST-flat` | `tests/test_flat_limit.py` | test_flat_symbolic_ricci_zero |  |  | **PASS**  |
| `TEST-schwarzschild` | `tests/test_schwarzschild_limit.py` | test_TH_AH_identity |  |  | **PASS**  |
| `TEST-schwarzschild` | `tests/test_schwarzschild_limit.py` | test_canary_wrong_embedding_is_rejected |  |  | **PASS**  |
| `TEST-schwarzschild` | `tests/test_schwarzschild_limit.py` | test_gtt_times_gRR_minus_one |  |  | **PASS**  |
| `TEST-schwarzschild` | `tests/test_schwarzschild_limit.py` | test_horizon_area_and_temperature |  |  | **PASS**  |
| `TEST-schwarzschild` | `tests/test_schwarzschild_limit.py` | test_smarr_schwarzschild |  |  | **PASS**  |
| `TEST-schwarzschild` | `tests/test_schwarzschild_limit.py` | test_vacuum_ricci_zero_numeric |  |  | **PASS**  |
| `TEST-schwarzschild` | `tests/test_schwarzschild_limit.py` | test_vacuum_ricci_zero_symbolic |  |  | **PASS**  |
| `TEST-kerr` | `tests/test_kerr_limit.py` | test_a_to_zero_is_schwarzschild |  |  | **PASS**  |
| `TEST-kerr` | `tests/test_kerr_limit.py` | test_ansatz_vacuum_matches_schwarzschild |  |  | **PASS**  |
| `TEST-kerr` | `tests/test_kerr_limit.py` | test_canary_bad_kerr_sign_rejected |  |  | **PASS**  |
| `TEST-kerr` | `tests/test_kerr_limit.py` | test_frame_dragging |  |  | **PASS**  |
| `TEST-kerr` | `tests/test_kerr_limit.py` | test_horizon_and_extremality |  |  | **PASS**  |
| `TEST-kerr` | `tests/test_kerr_limit.py` | test_kerr_ricci_flat |  |  | **PASS**  |
| `TEST-chi-to-zero` | `tests/test_chi_to_zero.py` | test_canary_alpha20_rejected |  |  | **PASS**  |
| `TEST-chi-to-zero` | `tests/test_chi_to_zero.py` | test_mv_positive_and_vanishing |  |  | **PASS**  |
| `TEST-chi-to-zero` | `tests/test_chi_to_zero.py` | test_no_linear_term |  |  | **PASS**  |
| `TEST-chi-to-zero` | `tests/test_chi_to_zero.py` | test_slope_is_two |  |  | **PASS**  |
| `TEST-smarr` | `tests/test_smarr_exact.py` | test_canary_missing_rotation_term_rejected |  |  | **PASS**  |
| `TEST-smarr` | `tests/test_smarr_exact.py` | test_extremal_kerr_and_rn |  |  | **PASS**  |
| `TEST-smarr` | `tests/test_smarr_exact.py` | test_identity_is_non_trivial |  |  | **PASS**  |
| `TEST-smarr` | `tests/test_smarr_exact.py` | test_kerr |  |  | **PASS**  |
| `TEST-smarr` | `tests/test_smarr_exact.py` | test_reissner_nordstrom |  |  | **PASS**  |
| `TEST-smarr` | `tests/test_smarr_exact.py` | test_schwarzschild |  |  | **PASS**  |
| `TEST-ellipse` | `tests/test_ellipse_circle.py` | test_canary_naive_bbox_estimator_gives_artifact |  |  | **PASS**  |
| `TEST-ellipse` | `tests/test_ellipse_circle.py` | test_circle_gives_exactly_zero |  |  | **PASS**  |
| `TEST-ellipse` | `tests/test_ellipse_circle.py` | test_circle_pixelised_floor |  |  | **PASS**  |
| `TEST-ellipse` | `tests/test_ellipse_circle.py` | test_kerr_a_squared_scaling |  |  | **PASS**  |
| `TEST-ellipse` | `tests/test_ellipse_circle.py` | test_kerr_eps_reference |  |  | **PASS**  |

## How to reproduce

```bash
python3 run_all.py            # everything
python3 run_all.py --only tests   # only the known-limit tests
```
