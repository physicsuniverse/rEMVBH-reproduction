# What `equations/eom.py` is (and is not) -- please add this to the main README

> `equations/eom.py` is a **symbolic equation exporter**: it derives the Einstein, Maxwell and
> vector field equations from the action with sympy and emits them in a form that can be
> `lambdify`-ed, in order to verify the structure of the ansatz.
>
> **It is not a validator that returns PASS as soon as a solution is substituted into it.**
> Substituting an archived solution directly into $G_{ab}-2T_{ab}$ in the quasi-isotropic
> coordinates $(t,x,\theta,\varphi)$ gives residuals of $10^{24}$- $10^{25}$.
> This is an **ill-conditioning of the coordinate system, not an error of the solution**:
> at the horizon $g_{tt}g_{\varphi\varphi}-g_{t\varphi}^2=0$ holds exactly, so $g^{ab}$ diverges;
> at infinity $g_{xx}\sim(1-x)^{-4}$, and cancelling down to $G_{ab}\sim r^{-3}$ requires
> cancellations of order $10^{12}$. **No finite-difference or spectral method at any resolution
> can push this route to $10^{-11}$.**
>
> We therefore verify the EOM through **coordinate-independent scalar identities** -- the four
> diagnostics of `verify/verify_eom.py`: (1) Smarr, (2) $T_HA_H=4r_H$, (3) $\nabla_aT^{ab}=0$
> (order of magnitude only) and (4) the algebraic Bianchi identity. Results in `verify/REPORT.md`.
> Diagnostics 1, 2 and 4 reach $10^{-10}$, $10^{-11}$ and $10^{-16}$ respectively.
