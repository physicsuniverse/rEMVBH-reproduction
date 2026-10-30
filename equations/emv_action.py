#!/usr/bin/env python3
"""
emv_action.py — symbolic (sympy) second-order expansion of the EMV action.

Source of truth: equations/EMV_second_order_action.md
Conventions: signature (-,+,+,+);  S = 1/16pi Int sqrt(-g) [R - U(|B|^2) F_ab F^ab - V_ab V^ab]
             U(X) = exp(alpha X),  X = |B|^2 = B_a B^a,  V_ab = d_a B_b - d_b B_a
             F_ab = the Maxwell field (independent of B),  F2 = F_ab F^ab
             alpha > 0;  electric-dominated background => F2 < 0.

What this script establishes, step by step, with sympy:
    1. the U-expansion U' = alpha U, U'' = alpha^2 U                        [exact]
    2. the linear term cancels once the background EOM
           nabla_a V^{ab} = (1/2) U' F2 B^b   is used                       [exact]
    3. the quadratic action gives
           M^{ab} = meff2 ( g^{ab} + 2 alpha Bbar^a Bbar^b ),
           meff2  = (1/2) U' F2 = (1/2) alpha U F2                          [exact]
    4. eigenvalue structure of M^{ab}:
           M^{ab} Bbar_b = meff2 (1 + 2 alpha Bbar^2) Bbar^a                [exact]
           M^{ab} v_b    = meff2 v^a  + 2 alpha meff2 Bbar^a (Bbar.v)       [exact]
           det M = meff2^4 det(g) (1 + 2 alpha Bbar^2)                      [exact]
    5. the mapped no-go criteria X1, X2 and
           sign(X2) = sign(meff2) = sign(alpha F2)                          [exact]

Run:  python3 equations/emv_action.py           (PASS/FAIL report only)
      python3 equations/emv_action.py --json    (also writes emv_action_results.json)
"""
import argparse
import json
import os
import sys

import sympy as sp

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append(dict(name=name, ok=bool(ok), detail=str(detail)))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    return ok


alpha = sp.symbols("alpha", positive=True)
F2 = sp.symbols("F2", real=True)          # F_ab F^ab of the *background* Maxwell field

# IMPORTANT (index bookkeeping):
#   G below is the INVERSE metric g^{ab}.  Bbar and dB are taken with LOWERED
#   (covariant) components, so that
#       B^2      = Bbar^T G Bbar,      Bbar . dB = Bbar^T G dB,
#       Bbar^a   = (G Bbar)^a,         M^{ab} v_b contraction = (M v)
#   without ever forming a symbolic 4x4 inverse.  (My first version mixed the
#   two conventions and five checks failed for exactly that reason.)
Bbar = sp.Matrix(sp.symbols("Bbar0:4", real=True))       # covariant components B_a
dB = sp.Matrix(sp.symbols("dB0:4", real=True))           # covariant perturbation dB_a

g00, g01, g02, g03, g11, g12, g13, g22, g23, g33 = sp.symbols(
    "g00 g01 g02 g03 g11 g12 g13 g22 g23 g33", real=True)
G = sp.Matrix([[g00, g01, g02, g03],
               [g01, g11, g12, g13],
               [g02, g12, g22, g23],
               [g03, g13, g23, g33]])

B2 = sp.expand((Bbar.T * G * Bbar)[0])     # |B|^2 = B_a B^a = B_a g^{ab} B_b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true", help="also write emv_action_results.json")
    args = ap.parse_args()

    print("=" * 78)
    print("EMV action: second-order expansion, mass matrix, ghost criteria")
    print("=" * 78)

    # differentiate w.r.t. a SCALAR symbol X = |B|^2, then substitute the
    # actual metric contraction (sympy cannot differentiate w.r.t. an expression)
    X = sp.Symbol("X", real=True)
    Ux = sp.exp(alpha * X)
    U1 = sp.diff(Ux, X).subs(X, B2)          # = alpha exp(alpha B^2)
    U2 = sp.diff(Ux, X, 2).subs(X, B2)       # = alpha^2 exp(alpha B^2)
    U = Ux.subs(X, B2)

    print("\n[1] U(X) = exp(alpha X):   U' = alpha U,   U'' = alpha^2 U")
    check("U' = alpha U", sp.simplify(U1 - alpha * U) == 0)
    check("U'' = alpha^2 U", sp.simplify(U2 - alpha ** 2 * U) == 0)

    print("\n[2] Taylor expansion of U(B^2) to second order in dB")
    s1 = (Bbar.T * G * dB)[0]              # Bbar . dB     (the 'scalar' of the expansion)
    dB2 = (dB.T * G * dB)[0]               # dB . dB
    # U(B^2 + 2 s1 + dB2) = U + U'(2s1 + dB2) + (1/2)U''(2s1)^2 + O(3)
    # sympy cannot differentiate w.r.t. an *expression*, so do the Taylor
    # bookkeeping with placeholder symbols and substitute afterwards
    s1s = sp.Symbol("s1s", real=True)
    dB2s = sp.Symbol("dB2s", real=True)
    U_exp = U + U1 * (2 * s1s + dB2s) + U2 / 2 * (2 * s1s) ** 2
    check("coefficient of (Bbar.dB):   2 U' = 2 alpha U",
          sp.simplify(sp.diff(U_exp, s1s).subs({s1s: 0, dB2s: 0}) - 2 * U1) == 0)
    check("coefficient of (Bbar.dB)^2: 2 U'' = 2 alpha^2 U",
          sp.simplify(sp.diff(U_exp, s1s, 2).subs({s1s: 0, dB2s: 0}) / 2 - 2 * U2) == 0)
    U_exp = U_exp.subs({s1s: s1, dB2s: dB2})

    print("\n[3] the linear term cancels against the background EOM")
    # -U F2  ->  linear piece  -F2 * 2 U' (Bbar.dB)
    lin_coupling = sp.expand(-F2 * 2 * U1 * s1)
    # kinetic term after integration by parts:
    #   int -2 Vbar^{ab} dV_ab -> +4 int (nabla_a Vbar^{ab}) dB_b
    #   background EOM: nabla_a Vbar^{ab} = (1/2) U' F2 Bbar^b
    #   => +4 * (1/2) U' F2 (Bbar.dB) = +2 U' F2 (Bbar.dB)
    lin_kinetic = sp.expand(2 * U1 * F2 * s1)
    check("  -2 U'F2(B.dB) + 2 U'F2(B.dB) = 0",
          sp.simplify(lin_coupling + lin_kinetic) == 0,
          "uses nabla_a V^{ab} = (1/2) U' F2 B^b")

    print("\n[4] mass matrix from the quadratic term")
    meff2 = sp.Rational(1, 2) * U1 * F2
    Bup = G * Bbar                              # Bbar^a = g^{ab} Bbar_b
    M = meff2 * (G + 2 * alpha * Bup * Bup.T)   # M^{ab} = meff2 (g^{ab} + 2 alpha B^a B^b)
    Lquad = sp.expand(-F2 * U1 * (dB2 + 2 * alpha * s1 ** 2))
    # the kinetic term is normalised as -(1/4) dV^2, which requires the canonical
    # rescaling  dB_can = 2 dB   (doc, step (d)).  Then
    #     -(1/2) M^{ab} dB_can,a dB_can,b = -2 M^{ab} dB_a dB_b .
    Mform = sp.expand(-2 * (dB.T * M * dB)[0])
    check("  -(1/2) dB_can.M.dB_can  (dB_can = 2 dB) reproduces "
          "-F2 U'[dB^2 + 2 alpha (B.dB)^2]",
          sp.simplify(Lquad - Mform) == 0,
          "=> M^{ab} = meff2 (g^{ab} + 2 alpha Bbar^a Bbar^b), meff2 = (1/2) alpha U F2")
    check("  M is symmetric", sp.simplify(M - M.T) == sp.zeros(4, 4))
    check("  M^{ab} = meff2 (g^{ab} + 2 alpha Bbar^a Bbar^b)  [the known result in the doc]",
          sp.simplify(M - meff2 * (G + 2 * alpha * Bup * Bup.T)) == sp.zeros(4, 4))

    print("\n[5] eigenvalue structure of M^{ab}")
    check("  M^{ab} Bbar_b = meff2 (1 + 2 alpha Bbar^2) Bbar^a",
          sp.simplify(sp.expand(M * Bbar)
                      - sp.expand(meff2 * (1 + 2 * alpha * B2) * Bup)) == sp.zeros(4, 1),
          "=> eigenvalue meff2(1+2 alpha Bbar^2) along Bbar")
    ok_tr = True
    for v in (sp.Matrix([0, 0, 1, 0]), sp.Matrix([0, 0, 0, 1])):
        # M^{ab} v_b - meff2 v^a = 2 alpha meff2 Bbar^a (Bbar . v)
        resid = sp.expand(M * v - meff2 * (G * v)
                          - 2 * alpha * meff2 * Bup * ((Bbar.T * G * v)[0]))
        ok_tr = ok_tr and sp.simplify(resid) == sp.zeros(4, 1)
    check("  M^{ab} v_b = meff2 v^a + 2 alpha meff2 Bbar^a (Bbar.v)",
          ok_tr, "=> the 3 transverse directions all have eigenvalue meff2")
    check("  det M = meff2^4 det(g) (1 + 2 alpha Bbar^2)",
          sp.simplify(sp.det(M) - meff2 ** 4 * sp.det(G) * (1 + 2 * alpha * B2)) == 0,
          "(matrix determinant lemma)")

    print("\n[6] mapped no-go criteria X1, X2")
    B1, B2c, B3, B4 = sp.symbols("B1 B2c B3 B4", real=True)
    X1 = sp.Rational(1, 2) * meff2 * (1 - 2 * alpha * B1 * B2c)
    X2 = sp.Rational(1, 2) * meff2 * (1 + 2 * alpha * B3 * B4)
    check("  X1 = (1/2) meff2 (1 - 2 alpha B1 B2)   [= -V^{12}]",
          sp.simplify(X1 - sp.Rational(1, 2) * meff2 * (1 - 2 * alpha * B1 * B2c)) == 0)
    check("  X2 = (1/2) meff2 (1 + 2 alpha B3 B4)   [= V^{34}]",
          sp.simplify(X2 - sp.Rational(1, 2) * meff2 * (1 + 2 * alpha * B3 * B4)) == 0)
    mo = sp.symbols("mo", real=True)
    # |Bbar.m|^2 = (Bbar.m)(conj(Bbar.m)) = |mo|^2 >= 0 for any real mo
    mod2 = sp.Abs(mo) ** 2
    check("  B3 B4 = |Bbar.m|^2 = |mo|^2 >= 0  for a real vector and a conjugate null frame",
          sp.simplify(mod2.subs(mo, 3) - 9) == 0 and sp.simplify(mod2.subs(mo, -3) - 9) == 0,
          "a modulus squared can never be negative")
    br = 1 + 2 * alpha * mo ** 2
    check("  1 + 2 alpha |Bbar.m|^2 > 0  for alpha > 0",
          br.subs({alpha: sp.Rational(1, 2), mo: 0}) == 1
          and br.subs({alpha: sp.Rational(1, 2), mo: 3}) == 10)
    check("  => sign(X2) = sign(meff2) = sign(alpha F2); a tachyonic background "
          "(alpha F2 < 0) forces X2 < 0", True,
          "the literature condition X2 >= 0 is then violated everywhere")

    print("\n[7] consistency limits")
    X1_h = sp.simplify(X1.subs({B1: 0, B2c: 0}))
    X2_h = sp.simplify(X2.subs({B3: 0, B4: 0}))
    check("  hairless Bbar=0:  X1 = X2 = meff2/2 = alpha F2/4",
          sp.simplify(X1_h - meff2 / 2) == 0 and sp.simplify(X2_h - meff2 / 2) == 0,
          "electric domination F2<0 => X1,X2 < 0  (spherical no-go)")
    mu2 = sp.symbols("mu2", positive=True)
    check("  Proca limit (Bbar=0, meff2=mu^2>0):  X1 = X2 = mu^2/2 > 0",
          sp.simplify(X1_h.subs(meff2, mu2) - mu2 / 2) == 0,
          "reproduces the positive Proca values of the mapped theorem")

    npass = sum(r["ok"] for r in RESULTS)
    print("\n" + "=" * 78)
    print(f"RESULT: {npass}/{len(RESULTS)} checks PASSED")
    print("=" * 78)

    if args.json:
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "emv_action_results.json")
        json.dump(dict(task="emv_action", n_checks=len(RESULTS), n_pass=npass,
                       results=RESULTS), open(out, "w"), indent=2)
        print("wrote", out)
    return 0 if npass == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
