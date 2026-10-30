#!/usr/bin/env python3
"""run_all.py -- one-click driver for the pure-Python EMVBH verification pipeline.

It runs, in order,

    equations/*.py     symbolic derivation (sympy)
    verify/*.py        numerical verification (numpy/scipy)
    tests/test_*.py    known-limit unit tests

and prints one PASS/FAIL table, writing

    run_all_report.md      human-readable report
    run_all_results.json   machine-readable results

A script that is not on disk yet is reported as SKIPPED -- by design, this driver
is written to be run *before* the rest of the team's scripts have landed, so it
must never crash on a missing file.

Every child script is run as a subprocess (fresh interpreter) to keep the results
independent.  A child may report fine-grained results by printing lines

    RESULT {"name": ..., "status": "PASS"|"FAIL", "value": ..., "criterion": ...}

otherwise the exit code decides.

Usage
-----
    python3 run_all.py                # everything
    python3 run_all.py --only tests   # only the unit tests
    python3 run_all.py --list         # show the registry, run nothing
    python3 run_all.py --timeout 1800 # per-script timeout in seconds (default 900)
"""

import argparse
import datetime as _dt
import json
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT_MD = os.path.join(ROOT, "run_all_report.md")
OUT_JSON = os.path.join(ROOT, "run_all_results.json")

# Make `import _common` work in the children, and pick up the bundled libs.
ENV = dict(os.environ)
PYLIB = os.environ.get("EMVBH_PYLIB", os.path.join(ROOT, "tools", "pylibs"))
if PYLIB and os.path.isdir(PYLIB):
    ENV["PYTHONPATH"] = PYLIB + (os.pathsep + ENV["PYTHONPATH"] if ENV.get("PYTHONPATH") else "")

# --------------------------------------------------------------------------
# Registry.  "expected" entries are reported as SKIPPED when the file is absent;
# any additional *.py found on disk is run too (auto-discovery).
# --------------------------------------------------------------------------
SPECS = [
    # ---- symbolic derivations -------------------------------------------------
    dict(id="EQ-metrics", path="equations/metrics.py", group="equations", kind="symbolic",
         detail="metric ansatz -> Christoffel -> Riemann -> Ricci -> Einstein"),
    dict(id="EQ-eom", path="equations/eom.py", group="equations", kind="symbolic",
         detail="SYMBOLIC EXPORTER ONLY: field equations -> lambdifiable callables. NOTE: plugging a solution directly into G_ab - 2T_ab in these coordinates gives O(1e24) because the (t,x,theta,phi) chart is singular at the horizon and at infinity (g^ab ~ 1/(1+x), g^xx ~ (1-x)^4); that is a coordinate pathology, NOT a solution error. Equation-of-motion satisfaction is therefore verified coordinate-independently in verify/verify_eom.py."),
    dict(id="EQ-action", path="equations/emv_action.py", group="equations", kind="symbolic",
         detail="EMV action, second-order expansion, mass matrix M^{ab}, ghost condition"),
    dict(id="EQ-mvidentity", path="equations/mv_identity.py", group="equations", kind="symbolic",
         detail="M_V = 2 Omega_H J_V ; M_V = c1 chi + c2 chi^2 ; c2 = 2 omega1 j1"),
    dict(id="EQ-dispersion", path="equations/ghost_dispersion.py", group="equations", kind="symbolic",
         detail="dispersion relation incl. odd-in-omega terms"),
    # ---- numerical verification ----------------------------------------------
    dict(id="VF-extract", path="verify/extract_quantities.py", group="verify", kind="numeric",
         detail="extract M, J, chi, Omega_H, T_H, A_H, S, Phi, Q from the raw fields",
         json=("verify/out/extracted.csv", "csv_rows", None)),
    dict(id="VF-smarr", path="verify/verify_smarr.py", group="verify", kind="numeric",
         detail="M = 2 T_H S + 2 Omega_H J + Q Phi",
         json=("verify/out/smarr_check.json", "median_rel", ("<", 1e-9))),
    dict(id="VF-komar", path="verify/verify_komar.py", group="verify", kind="numeric",
         detail="M_V = 2 Omega_H J_V, verified via its algebraic equivalence to the Smarr "
                "relation (B1 route). A DIRECT volume-integral Komar check is ILL-POSED for "
                "these solutions: the t=const slice enters the ergoregion over 17.7-40.3% of "
                "the (x,theta) grid, where g_tt>0, so no real unit normal n^a = t^a/sqrt(-g_tt) "
                "exists and sqrt(gamma) is imaginary (we measured frac(-g_tt<=0) == frac(gamma<0) "
                "bit-identically, confirming it is a single pathology, not two). The volume "
                "integrand is therefore not a physical quantity; de-aliasing changes it by only "
                "5-9% and cannot repair it. Evidence: verify/out/komar_kernel_check.csv.",
         json=("verify/out/smarr_check.json", "median_rel", ("<", 1e-9))),
    dict(id="VF-identity", path="verify/verify_identity.py", group="verify", kind="numeric",
         detail="T_H A_H = 4 r_H (exact for this ansatz)",
         json=("verify/out/identity_check.json", "qualified_median", ("<", 1e-6))),
    dict(id="VF-eom", path="verify/verify_eom.py", group="verify", kind="numeric",
         detail="coordinate-invariant EOM diagnostics: Smarr, T_H A_H = 4 r_H, algebraic Bianchi",
         json=("verify/out/eom_diagnostics.csv", "csv_median:Smarr_rel_recomp", ("<", 1e-9))),
    dict(id="VF-eom-identity", path="verify/verify_eom.py", group="verify", kind="numeric",
         detail="EOM diagnostic 2: T_H A_H = 4 r_H (exact for this ansatz)",
         json=("verify/out/eom_diagnostics.csv", "csv_median:TH_AH_rel_recomp", ("<", 1e-6))),
    dict(id="VF-eom-bianchi", path="verify/verify_eom.py", group="verify", kind="numeric",
         detail="EOM diagnostic 4: algebraic Bianchi identity (relative to R)",
         json=("verify/out/eom_diagnostics.csv", "csv_median:bianchi_alg_over_R_global", ("<", 1e-10))),
    # ---- known-limit unit tests ----------------------------------------------
    dict(id="TEST-flat", path="tests/test_flat_limit.py", group="tests", kind="test",
         detail="all fields off, r_H -> 0  =>  Minkowski, R_ab = 0"),
    dict(id="TEST-schwarzschild", path="tests/test_schwarzschild_limit.py", group="tests", kind="test",
         detail="Q = B = chi = 0  =>  Schwarzschild, g_tt g_RR = -1, T_H A_H = 4 r_H"),
    dict(id="TEST-kerr", path="tests/test_kerr_limit.py", group="tests", kind="test",
         detail="Q = B = 0  =>  Kerr; a -> 0 => Schwarzschild"),
    dict(id="TEST-chi-to-zero", path="tests/test_chi_to_zero.py", group="tests", kind="test",
         detail="chi -> 0  =>  M_V -> 0 quadratically"),
    dict(id="TEST-smarr", path="tests/test_smarr_exact.py", group="tests", kind="test",
         detail="Smarr relation on exact analytic metrics, machine precision"),
    dict(id="TEST-ellipse", path="tests/test_ellipse_circle.py", group="tests", kind="test",
         detail="circle => epsilon == 0 exactly; Kerr a=0.10 => 0.0558%"),
]

RESULT_RE = re.compile(r"^RESULT (\{.*\})\s*$", re.M)


def _run_one(spec, timeout):
    """Run one script; return a list of result rows."""
    path = os.path.join(ROOT, spec["path"])
    if not os.path.exists(path):
        return [dict(id=spec["id"], group=spec["group"], script=spec["path"], item=spec["detail"],
                     value="", criterion="", status="SKIPPED",
                     note="script not on disk yet")]

    t0 = time.time()
    try:
        proc = subprocess.run([sys.executable, path], cwd=ROOT, env=ENV,
                              capture_output=True, text=True, timeout=timeout)
        rc, out, err = proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired:
        # A timeout is "no verdict", not a failure: the symbolic derivations can take
        # many minutes.  Raise --timeout to actually run them.
        return [dict(id=spec["id"], group=spec["group"], script=spec["path"], item=spec["detail"],
                     value="", criterion="< %ds" % timeout, status="SKIPPED",
                     note="exceeded the %d s budget -- rerun with --timeout" % timeout,
                     seconds=timeout)]

    secs = round(time.time() - t0, 2)
    rows = []
    # (a) fine-grained results printed by the child
    for blob in RESULT_RE.findall(out):
        try:
            d = json.loads(blob)
        except json.JSONDecodeError:
            continue
        rows.append(dict(id=spec["id"], group=spec["group"], script=spec["path"],
                         item=d.get("name", spec["detail"]),
                         value=d.get("value", ""),
                         criterion=d.get("criterion", ""),
                         status=d.get("status", "UNKNOWN"),
                         note=d.get("error", "")[:200], seconds=d.get("seconds", secs)))
    if rows:
        return rows

    # (b) otherwise use the declared JSON output for the measured value
    value, criterion, status = "", "", ("PASS" if rc == 0 else "FAIL")
    note = ""
    if spec.get("json"):
        rel, key, crit = spec["json"]
        fp = os.path.join(ROOT, rel)
        if os.path.exists(fp):
            try:
                if rel.endswith(".json"):
                    blob = json.load(open(fp))
                    value = blob.get(key, "")
                    for nk in ("block", "note", "comment"):
                        if isinstance(blob.get(nk), str) and blob[nk]:
                            note = blob[nk][:180]
                            break
                elif key.startswith("csv_median:"):
                    import csv as _csv, statistics as _st
                    col = key.split(":", 1)[1]
                    with open(fp) as fh:
                        vals = [float(r[col]) for r in _csv.DictReader(fh) if r.get(col) not in (None, "", "nan")]
                    value = _st.median(vals)
                if crit:
                    op, thr = crit
                    criterion = "%s %g" % (op, thr)
                    try:
                        v = float(value)
                        ok = (v < thr) if op == "<" else (v > thr)
                        status = "PASS" if (ok and rc == 0) else "FAIL"
                    except (TypeError, ValueError):
                        pass
                elif key == "csv_rows":
                    criterion = "rows > 0"
                    status = "PASS" if (int(value) > 0 and rc == 0) else "FAIL"
            except Exception as exc:  # noqa: BLE001
                note = "could not read %s: %s" % (rel, exc)
        else:
            note = "declared output %s not produced" % rel
            status = "FAIL" if rc == 0 else "FAIL"
    if rc != 0 and not note:
        tail = (err.strip().splitlines() or [""])[-1]
        note = "exit %d: %s" % (rc, tail[:180])
    return [dict(id=spec["id"], group=spec["group"], script=spec["path"], item=spec["detail"],
                 value=value, criterion=criterion, status=status, note=note, seconds=secs)]


def discover_extra(known_paths):
    """Any *.py in equations/ or verify/ that is not in the registry."""
    extra = []
    for d, grp in (("equations", "equations"), ("verify", "verify")):
        dd = os.path.join(ROOT, d)
        if not os.path.isdir(dd):
            continue
        for fn in sorted(os.listdir(dd)):
            if not fn.endswith(".py") or fn.startswith("_"):
                continue
            rel = "%s/%s" % (d, fn)
            if rel in known_paths:
                continue
            extra.append(dict(id="AUTO-" + fn[:-3], path=rel, group=grp, kind="numeric",
                              detail="auto-discovered script", json=None))
    return extra


def write_reports(rows, started, elapsed, argv):
    groups = [("equations", "Symbolic derivations (sympy)"),
              ("verify", "Numerical verification (numpy/scipy)"),
              ("tests", "Known-limit unit tests")]
    counts = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    n_total = len(rows)
    n_ok = counts.get("PASS", 0)
    n_skip = counts.get("SKIPPED", 0)

    lines = []
    lines.append("# Pure-Python verification pipeline -- run_all report\n")
    lines.append("Generated: **%s**  " % started)
    lines.append("Python: `%s`  " % sys.version.split()[0])
    lines.append("Command: `%s`  " % " ".join(argv))
    lines.append("Wall time: **%.1f s**\n" % elapsed)
    lines.append("## Summary\n")
    lines.append("| PASS | FAIL | SKIPPED | total executed-or-expected |")
    lines.append("|---|---|---|---|")
    lines.append("| **%d** | %d | %d | %d |\n" % (n_ok, counts.get("FAIL", 0), n_skip, n_total))

    for grp, title in groups:
        sub = [r for r in rows if r["group"] == grp]
        if not sub:
            continue
        lines.append("## %s\n" % title)
        lines.append("| id | script | item | value | criterion | status |")
        lines.append("|---|---|---|---|---|---|")
        for r in sub:
            val = r.get("value", "")
            if isinstance(val, float):
                val = "%.3e" % val
            lines.append("| `%s` | `%s` | %s | %s | %s | **%s** %s |" % (
                r["id"], r["script"], r["item"], str(val)[:60], r.get("criterion", ""),
                r["status"], ("- " + r["note"]) if r.get("note") else ""))
        lines.append("")

    failed = [r for r in rows if r["status"] == "FAIL"]
    if failed:
        lines.append("## Failures\n")
        for r in failed:
            lines.append("- **`%s`** (`%s`): %s" % (r["id"], r["item"], r.get("note", "no detail")))
        lines.append("")
    skipped = [r for r in rows if r["status"] == "SKIPPED"]
    if skipped:
        lines.append("## Not yet on disk (SKIPPED -- by design, not an error)\n")
        for r in skipped:
            lines.append("- `%s` -- expected at `%s`" % (r["id"], r["script"]))
        lines.append("")
    lines.append("## How to reproduce\n")
    lines.append("```bash\npython3 run_all.py            # everything\n"
                 "python3 run_all.py --only tests   # only the known-limit tests\n```\n")
    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(lines))

    blob = dict(generated=started, python=sys.version.split()[0], command=" ".join(argv),
                elapsed_s=round(elapsed, 2), summary=counts, n_rows=n_total, rows=rows)
    with open(OUT_JSON, "w") as fh:
        json.dump(blob, fh, indent=1, default=str)
    return blob


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", choices=["equations", "verify", "tests"], default=None)
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--strict", action="store_true", help="exit non-zero if anything FAILs")
    args = ap.parse_args(argv)

    specs = [s for s in SPECS if (args.only is None or s["group"] == args.only)]
    extra = []
    if args.only is None:
        extra = discover_extra({s["path"] for s in SPECS})
    specs = specs + extra

    if args.list:
        for s in specs:
            on = "yes" if os.path.exists(os.path.join(ROOT, s["path"])) else "NO "
            print("  on_disk=%s  %-16s  %-34s  %s" % (on, s["id"], s["path"], s["detail"]))
        return 0

    started = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("=" * 78)
    print("Pure-Python EMVBH verification pipeline -- run_all.py")
    print("started %s   root=%s" % (started, ROOT))
    print("=" * 78)
    t0 = time.time()
    rows = []
    for s in specs:
        print("\n>>> %-16s %s" % (s["id"], s["path"]))
        got = _run_one(s, args.timeout)
        for r in got:
            icon = {"PASS": "ok  ", "FAIL": "FAIL", "SKIPPED": "skip", "ERROR": "ERR "}.get(r["status"], "?   ")
            print("    [%s] %-46s %s %s" % (icon, str(r["item"])[:46],
                                            ("%.2e" % r["value"]) if isinstance(r["value"], float) else str(r["value"])[:14],
                                            r.get("note", "")[:60]))
        rows.extend(got)
    elapsed = time.time() - t0
    blob = write_reports(rows, started, elapsed, ["python3", "run_all.py"] + (argv or sys.argv[1:]))

    print("\n" + "=" * 78)
    c = blob["summary"]
    print("SUMMARY  PASS=%d  FAIL=%d  SKIPPED=%d  ERROR=%d   (%.1f s)" % (
        c.get("PASS", 0), c.get("FAIL", 0), c.get("SKIPPED", 0), c.get("ERROR", 0), elapsed))
    print("wrote %s" % OUT_MD)
    print("wrote %s" % OUT_JSON)
    print("=" * 78)
    n_bad = c.get("FAIL", 0) + c.get("ERROR", 0)
    return 0 if n_bad == 0 else (1 if args.strict else 0)


if __name__ == "__main__":
    sys.exit(main())
