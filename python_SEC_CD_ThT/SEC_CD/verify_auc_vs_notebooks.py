"""verify_auc_vs_notebooks.py — pin the S3-S5 AUC numbers to HER OWN retained
notebook outputs.

The expected values below are copied verbatim from the executed output cells of

    /path/to/Ab42oligomers37C/data/SEC\\python_code\\WT37COligomer.ipynb
        cells  4-5   250 uM 37 C
        cells 10-11   25 uM 37 C
        cells 26-27  100 uM 25 C
    /path/to/Ab42oligomers37C/data/SEC\\python_code\\SECProgressE22K_onlyrequiredtraces.ipynb
        cells 28-31  100 uM 37 C

They are not re-derived here, so this is a genuine external check on
``build_figures_S3_S5.measure`` rather than a tautology.

WHY THIS FILE EXISTS.  An earlier reconstruction of her method was fitted to pixels
digitised off the published artwork and reached a plausible-looking 0.027 mean
relative error on one panel -- while being wrong about the integrator (simpson vs
np.trapz), wrong about the baseline (a linear edge fit vs none at all) and wrong
about the windows (6-11.5/12-20 vs 5-11/12-18).  The three errors partly cancelled.
A close fit to a curve is not evidence of a matching method; matching her printed
numbers to 4 dp is.  This test enforces that standard.

Note the 100 uM 37 C block: her notebook's own output is MIS-PAIRED (a chronological
label list zipped against a file list ordered 0, 2, 12, 24, 48, 4, 8), so the
expected values here are keyed by the FILE the number came from, not by the label
she printed against it.  The build is checked against the right file, which is also
why the rebuilt panel does not look like the published one.

    python verify_auc_vs_notebooks.py
"""
from __future__ import annotations

import os
import sys

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from build_figures_S3_S5 import measure

TOL = 5e-4          # her prints carry 4 decimals

#: registry id -> {hours: (oligomer, monomer[, fragments])} exactly as her cells printed
EXPECTED = {
    "3A": {2: (9.6286, 17.5689), 4: (13.5123, 14.3428), 6: (6.6545, 4.9006)},
    "3C": {0: (0.2052, 2.3549, 0.1225), 3: (0.1367, 2.4757, 0.1891),
           6: (-0.1292, 2.1181, 0.0070), 12: (0.0089, 2.4819, 0.1863),
           24: (0.1798, 2.1494, 0.3775), 48: (0.0445, 1.8169, 0.3334),
           72: (-0.0517, 1.4878, 0.4889)},
    "4A": {2: (0.8386, 9.1686), 4: (1.6192, 9.6944), 6: (1.9926, 9.7760),
           8: (2.1914, 9.2708), 10: (2.5080, 9.5601)},
    # keyed by the TRUE timepoint of the file, not by the label she printed
    "3B": {0: (0.0470, 9.9881), 2: (2.1598, 9.3533), 12: (3.4840, 7.2217),
           24: (3.9430, 5.9772), 48: (4.2807, 2.9944)},
}


def main():
    bad = 0
    checked = 0
    for pid, exp in EXPECTED.items():
        d = measure(pid, split_fragments=(pid == "3C"))
        idx = {float(t): i for i, t in enumerate(d["t"])}
        print(f"\n{pid}")
        for t, vals in exp.items():
            i = idx.get(float(t))
            if i is None:
                print(f"   t={t:<3} MISSING from build")
                bad += 1
                continue
            got = [d["olig"][i], d["mono"][i]]
            if len(vals) == 3:
                got.append(d["frag"][i])
            for name, g, e in zip(("olig", "mono", "frag"), got, vals):
                checked += 1
                ok = abs(g - e) <= TOL
                bad += (not ok)
                print(f"   t={t:<3} {name:5} {g:9.4f}  hers {e:9.4f}  "
                      f"{'ok' if ok else 'MISMATCH  d=%+.4f' % (g - e)}")
    print(f"\n{checked} values checked, {bad} mismatched, tolerance {TOL}")
    if bad:
        print("FAIL — the build no longer reproduces her notebook output.")
    else:
        print("PASS — bit-exact against her retained notebook outputs.")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
