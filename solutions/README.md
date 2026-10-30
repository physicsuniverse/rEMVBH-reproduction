# Solution files

## `physical_quantities.csv`
One row per archived solution. Physical quantities were extracted **independently** from
the metric file (not read from the header), so that the header values can be checked
against them. Columns include `M`, `J`, `chi_J`, `Omega_H`, `T_H`, `A_H`, `S`,
`T_H_tilde`, `S_tilde`, `TH_AH`, `TH_AH_rel_dev`, and the header-vs-extracted deviations.

## `metric_files/`
Raw metric tensors for the solutions used in the shadow figures. Format: the
quasi-isotropic coordinate grid $(x,\theta)$ and the eight functions
$\{B_t, B_\varphi, A_t, A_\varphi, f, g, h, W\}$.

## `field_profiles/`
Field profiles exported on a uniform grid for a representative set of solutions
(`sol_chi*_q*.txt`), used to produce the appendix figures. `grids.txt` gives the
$(x,\theta)$ grid the profiles are sampled on.

> **Not included:** the solver that generates these files. Only the solutions themselves.
