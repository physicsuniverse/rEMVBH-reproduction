# Unqualified / gray-zone profiles (kept deliberately)

These two profiles are kept deliberately: they are the two that **fail** the accuracy criteria,
and they demonstrate that the criteria are not vacuous.

| profile | Smarr residual (rel.) | $T_H A_H = 4 r_H$ rel. dev. | $g_H$ | status |
|---|---|---|---|---|
| `sol_chi0p11_q1p34.txt` | $1.34\times10^{-2}$ | $5.27\%$ | 3.9948606495 (≠ 4, off by $1.3\times10^{-3}$) | **unqualified** (fails both) |
| `sol_q6937_chi0p73.txt` | $2.17\times10^{-5}$ | $2.06\times10^{-5}$ | 3.9999986237 | **gray zone** |

`sol_chi0p11_q1p34` violates both the Smarr relation ($1.3\times10^{-2}$) and the exact identity
$T_HA_H = 4r_H$ (5.27%, with $g_H = 3.99486$ instead of 4);
`sol_q6937_chi0p73` sits in the gray region (Smarr $2.2\times10^{-5}$).
**Neither should be used for quantitative statements.**
