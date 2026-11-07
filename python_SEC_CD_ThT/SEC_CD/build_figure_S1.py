"""build_figure_S1.py — Supplementary Figure S1, remade in Nature format from the
artwork that already exists in 10_Collaborator_Review_Clean.docx.

IT ALREADY EXISTS: word/media/image15.png, captioned at paragraph 201:

    "Figure S1. SEC profile measured for injection of synthetic Abeta42 prepped in
     50 mM NaOH at a peptide concentration of 2 mg/mL. The shaded region represents
     monomeric species isolated for all subsequent experiments."

and the surrounding text adds the method: Superdex 75 Increase, 20 mM NaPi mobile
buffer, monomeric Abeta42 (4.5 kDa) eluting near 15 mL.

SOURCE DATA IDENTIFIED BY MATCHING, NOT BY GUESSING.  Her artwork was digitised (the
black trace isolated as the widest connected component -- a naive "topmost dark pixel
per column" locks onto the in-plot title instead) and every WT primary/crude run in
Run_Results was scored against it:

    rho    rms   apex   file
    1.000  0.008 15.74  WT/06-04-25 Results/wt primary crude 060325 5ml loop ...  <-- THIS
    0.991  0.048 15.68  WT/WT crude primary 6mg 021925 5ml loop ...
    0.965  0.085 15.68  RT runs/WT prim run 122025 001.csv
    0.830  0.199 15.20  WT/WT 4mg primary crude 022025 2ml loop ...

The top match is exact.  Note the filename says a 5 mL loop and gives no concentration,
so the caption's "2 mg/mL" is NOT verifiable from the file -- it is taken on trust from
the manuscript.  Digitised reference apex 15.84 mL against 15.74 measured, a 0.10 mL
difference, which is about five pixels of the source image.

SHADED REGION: 14.50 - 16.03 mL, measured from the fill colour (207, 221, 229) in her
artwork, redrawn here with the same colour.

ONE DELIBERATE DEPARTURE.  Her y axis carries tick marks but NO numbers, so the
absorbance scale is unreadable.  This version labels it in mAU.  That adds information
and misrepresents nothing; if the arbitrary-units look is wanted back, pass --no-ynum.

    python build_figure_S1.py
    python build_figure_S1.py --no-ynum
"""
from __future__ import annotations

import argparse
import os
import sys

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True

from dataclasses import replace

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from secfigs import style as gs                                             # noqa: E402
from secfigs.figure_fig12 import _NAT_STYLE, MARGIN, TRIM, MM, FIG_W, RC_EXTRA  # noqa: E402
from build_fig1_panels import load_trace                                    # noqa: E402

OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

CSV = (r"/path/to/Ab42oligomers37C/data/SEC/Run_Results/WT/06-04-25 Results/"
       r"wt primary crude 060325 5ml loop all fractions 001 001.csv")
SHADE = (14.50, 16.03)          # measured off her fill colour (207, 221, 229)
SHADE_RGB = "#CFDDE5"
XLIM = (6.5, 19.0)
ASPECT = 1.70                   # between the Griffin chromatogram 1.270 and a flat band


def build(ynum=True, formats=("svg", "pdf", "png")):
    if not os.path.exists(CSV):
        raise SystemExit("source not found: " + CSV)
    V, A, _ = load_trace(CSV, normalized=False)
    s = (V >= XLIM[0] - 0.5) & (V <= XLIM[1] + 0.5)
    V, A = V[s], A[s]

    avail = FIG_W - MARGIN["left"] - MARGIN["right"]
    ph = avail / ASPECT
    fh = TRIM + MARGIN["top"] + ph + MARGIN["bottom"]
    style = replace(_NAT_STYLE, figsize=(FIG_W, fh))

    with gs.apply(style, **RC_EXTRA):
        fig = plt.figure(figsize=(FIG_W, fh), dpi=style.dpi)
        fx, fy = (lambda v: v / FIG_W), (lambda v: v / fh)
        ax = fig.add_axes([fx(MARGIN["left"]), 1.0 - fy(TRIM + MARGIN["top"] + ph),
                           fx(avail), fy(ph)])

        m = (V >= SHADE[0]) & (V <= SHADE[1])
        ax.fill_between(V[m], 0.0, A[m], color=SHADE_RGB, lw=0, zorder=1)
        ax.plot(V, A, color="black", lw=style.trace_linewidth,
                solid_capstyle="round", zorder=2)

        hi = float(A[(V >= XLIM[0]) & (V <= XLIM[1])].max())
        gs.axes(ax, style, normalized=False, x_major=2.0, x_range=XLIM,
                y_range=(0.0, hi * 1.08))
        gs.set_labels(ax, "Elution Volume (ml)", "Absorbance (mAU)", style)
        gs.set_title(ax, u"Monomer isolation from synthetic Aβ42", style)
        if not ynum:
            ax.set_yticklabels([])
            ax.set_ylabel("Absorbance (a.u.)", fontsize=style.axis_label_size,
                          fontweight=style.axis_label_weight)

        stem = "Figure_S1" + ("" if ynum else "_noynum")
        paths = []
        for ext in formats:
            p = os.path.join(OUT, f"{stem}.{ext}")
            fig.savefig(p, dpi=600 if ext == "png" else None, facecolor="white")
            paths.append(p)
    plt.close(fig)

    apex = float(V[(V >= 14) & (V <= 17)][np.argmax(A[(V >= 14) & (V <= 17)])])
    print(f"=== {stem}   {FIG_W/MM:.4f} x {fh/MM:.4f} mm")
    print(f"   source   {os.path.basename(CSV)}")
    print(f"   monomer apex {apex:.2f} mL   (her artwork digitises to 15.84)")
    print(f"   shaded   {SHADE[0]:.2f}-{SHADE[1]:.2f} mL, colour {SHADE_RGB}")
    print(f"   y axis   {'labelled in mAU' if ynum else 'unlabelled (a.u.), as in her artwork'}")
    print("\nCAPTION (from the manuscript, verbatim):")
    print("   Figure S1. SEC profile measured for injection of synthetic Aβ42 prepped in")
    print("   50 mM NaOH at a peptide concentration of 2 mg/mL. The shaded region")
    print("   represents monomeric species isolated for all subsequent experiments.")
    for p in paths:
        print("   wrote", p)
    return paths


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-ynum", action="store_true",
                    help="drop the y numerals, as in her original artwork")
    a = ap.parse_args(argv)
    build(ynum=not a.no_ynum)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
