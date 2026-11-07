"""build_fig2_sec_24h.py — the SEC companion to Figure 2: 100 uM, 24 h, three
temperatures, NON-NORMALISED.

Figure 2 claims "Aβ42 exhibits temperature-dependent conversion from monomers to
oligomers" and currently supports it with CD alone (panel 2A: 6 / 25 / 37 °C after
24 h, raw mdeg).  SEC shows the same conversion directly -- monomer peak falling,
oligomer peak rising -- at the same three temperatures and the same timepoint, so
it belongs beside it.

RAW mAU, NOT NORMALISED, and deliberately so.  Min-max normalisation rescales each
trace to its own maximum, which is exactly the quantity that distinguishes the
three temperatures; it would make them look alike.  What the raw traces show
(monomer 12-18 mL, oligomer 7-10 mL, AUC in mAU*mL):

        24 h monomer     24 h oligomer     monomer relative to that series' t = 0
   6 °C   11.381*             0.538*                   0.90
  25 °C    7.251              2.889                    0.76
  37 °C    6.616              4.444                    0.63
        * BOTH cells x5 loop-corrected. An earlier version corrected the
          monomer column only, printing 0.108 for the oligomer -- the two
          columns of one table on different scales, understating the
          oligomer/monomer ratio 0.0472 as 0.0095. The code was always
          right; only this table was wrong.

The oligomer AUC rises monotonically 0.108 -> 2.889 -> 4.444 with temperature.
That is the figure's claim, measured.

TWO THINGS THE READER MUST BE TOLD, both stated in the caption text this script
prints, neither of them correctable by any honest arithmetic:

 1. LOOP.  The 6 °C 24 h reinjection was run on a 100 uL loop; every other run
    here used 500 uL.  Injected volume is not concentration, so that trace is
    multiplied by the NOMINAL LOOP RATIO, 5.  This is a known instrument setting,
    not a fitted scale factor.  Without it the panel would show a false 5x
    monomer collapse at 6 °C -- the opposite of the truth.
 2. PREP.  The three series are different isolations, and their own starting
    monomer AUCs are 12.604 (6 °C), 9.558 (25 °C, the 2 h point -- no zero-time
    reinjection was collected for this series) and 10.495 (37 °C): a spread of
    about 25 %.  No correction is applied for this, because dividing by t = 0
    would BE a normalisation and this panel exists not to be one.  It means
    absolute heights are comparable to about 25 %, while the oligomer/monomer
    ratio WITHIN each trace is exact.  The conversion trend is far larger than
    25 %, so it survives; a reader should still be told.

Colours are Figure 2's own TEMP_INKS, so a temperature is the same colour in the
CD panel and here.

    python build_fig2_sec_24h.py            # bold, matching the rest of Figure 2
    python build_fig2_sec_24h.py --regular  # the measured Griffin reference weights
    python build_fig2_sec_24h.py --letter b
"""
from __future__ import annotations

import argparse
import os
import sys

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from secfigs import style as gs
from secfigs.figure_fig12 import TEMP_INKS          # 6 / 25 / 37 °C, Figure 2's key
from build_fig1_panels import load_trace            # the manuscript's own pipeline

RES = r"/path/to/Ab42oligomers37C/data/SEC\Run_Results"     # READ ONLY
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

MONOMER = (12.0, 18.0)
OLIGOMER = (7.0, 10.0)
trapz = getattr(np, "trapezoid", np.trapz)

#: (label, 24 h run, loop uL, that series' t = 0 run, its loop uL)
TEMPS = [
    ("6 °C",  f"{RES}/E22K/WT 0213 monomer reinj 021425 100ul Loop no fractions  001.csv", 100,
              f"{RES}/E22K/WT immediate reinj 021325 500ul Loop no fractions 001.csv", 500),
    ("25 °C", f"{RES}/RT runs/24h RT WT reinjection 001.csv", 500,
              f"{RES}/RT runs/WT 100uM RT 2hr 122025  001.csv", 500),
    ("37 °C", f"{RES}/WT/WT 24hr reinj 37C 022025 500ul Loop no fractions 001.csv", 500,
              f"{RES}/WT/WT 0 hr reinj 37C 021925 500ul Loop no fractions.csv", 500),
]
REF_LOOP = 500


def area(v, a, window):
    s = (v >= window[0]) & (v <= window[1])
    return float(trapz(a[s], v[s])) if s.sum() > 1 else float("nan")


def build(stem="Fig2_SEC_24h_three_temperatures_nonnormalized", xlim=(6.0, 20.0),
          bold=True, letter=None, formats=("svg", "pdf", "png")):
    recs = []
    for label, p24, loop, p0, loop0 in TEMPS:
        for p in (p24, p0):
            if not os.path.exists(p):
                raise FileNotFoundError(p)
        # normalized=False -> the manuscript pipeline WITHOUT the min-max step:
        # utf-16/tab/skiprows=3, first two numeric columns, sort, 0 <= V <= 25,
        # subtract the minimum only if negative.  Same reader, raw output.
        v, a, _ = load_trace(p24, normalized=False)
        v0, a0, _ = load_trace(p0, normalized=False)
        k = REF_LOOP / loop                       # 5.0 at 6 °C, 1.0 elsewhere
        y = a * k
        recs.append(dict(label=label, v=v, y=y, loop=loop, k=k,
                         mono=area(v, y, MONOMER), olig=area(v, y, OLIGOMER),
                         mono0=area(v0, a0 * (REF_LOOP / loop0), MONOMER)))

    style = gs.bold_variant() if bold else gs.STYLE
    with gs.apply(style):
        fig, ax = gs.figure(style)
        for rec, c in zip(recs, TEMP_INKS):
            ax.plot(rec["v"], rec["y"], color=c, lw=style.trace_linewidth,
                    label=rec["label"], solid_capstyle="round")

        ymax = max(r["y"].max() for r in recs)
        gs.axes(ax, style, normalized=False, x_major=2.0, x_range=xlim,
                y_range=(0.0, ymax * 1.08))
        gs.set_labels(ax, "Elution Volume (ml)", "Absorbance (mAU)", style)
        gs.legend(ax, style)
        if letter:
            gs.panel_letter(ax, letter, style)

        paths = []
        for ext in formats:
            p = os.path.join(OUT, f"{stem}.{ext}")
            fig.savefig(p, dpi=600 if ext == "png" else None, facecolor="white")
            paths.append(p)
    plt.close(fig)

    print("=== Figure 2 SEC companion — 100 µM, 24 h, raw mAU "
          f"({'bold' if bold else 'reference'} weights)")
    for r in recs:
        print(f"   {r['label']:6} loop {r['loop']:3d} µL ×{r['k']:.0f}   "
              f"monomer {r['mono']:7.3f}   oligomer {r['olig']:6.3f}   "
              f"monomer/t0 {r['mono']/r['mono0']:.2f}")
    print("\n   caption facts:")
    print("     oligomer AUC rises %s with temperature"
          % " -> ".join(f"{r['olig']:.3f}" for r in recs))
    print("     6 °C trace ×5, the nominal loop ratio (100 µL run against 500 µL elsewhere)")
    print("     starting monomer AUCs differ by ~25 %% across the three preps "
          "(%s); no correction applied"
          % ", ".join(f"{r['label']} {r['mono0']:.2f}" for r in recs))
    for p in paths:
        print("   wrote", p)
    return paths


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regular", action="store_true",
                    help="use the measured Griffin reference weights instead of bold")
    ap.add_argument("--letter", default=None, help="panel letter, e.g. b")
    ap.add_argument("--stem", default="Fig2_SEC_24h_three_temperatures_nonnormalized")
    a = ap.parse_args(argv)
    build(stem=a.stem, bold=not a.regular, letter=a.letter)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
