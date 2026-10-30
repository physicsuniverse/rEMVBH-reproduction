# Plotting code

Each script regenerates one figure from the corresponding CSV in `../data/`.
`figstyle.py` is imported by all of them and holds the shared style, so the output matches
the published figures (sizes are built at the final printed width, so 9 pt type stays 9 pt).

| Script | Figure |
|---|---|
| `fig1_massofV.py` | Fig. 1 — $M_V$ landscape and the quadratic-scaling panel |
| `fig2_lowq.py` | Fig. 2 — low-charge existence domain |
| `fig3_highq.py` | Fig. 3 — high-charge existence domain |
| `fig4_shadows.py` | Fig. 4 — shadows and ellipticity |
| `fig_solutions.py` | Appendix — 3D field profiles |
| `fig10_spectral.py` | Spectral smoothness |

Requirements: `numpy`, `matplotlib` (with a working LaTeX install for `text.usetex`).
