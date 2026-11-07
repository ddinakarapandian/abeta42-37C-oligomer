"""build_figure_S3_combined.py — Supplementary Figure S3: SEC peak area vs reinjection
time, all four unique conditions in one figure.

HER LAYOUT, CORRECTED DATA.  The collaborator's consolidated AUC figure puts every
condition in one display item, which is the right call -- but the version supplied has
six panels of which only four are distinct (its c/d and b/e are the same data twice),
and its 100 uM / 37 C panels carry a time-ordering fault.  This build keeps the
one-figure arrangement, drops the duplicates, and uses the corrected values.

THE ORDERING FAULT, for the record.  ``SECProgressE22K_onlyrequiredtraces.ipynb`` cell 28
declares the file list as 0, 2, 12, 24, 48, 4, 8 h while the label list next to it is
chronological, and cell 30 zips the two together -- so from the third point on every area
is attached to the wrong time.  Plotted against ``times_hr_WT37C = [0,2,4,8,12,24,48]``
this makes the monomer appear to collapse to 2.99 at "12 h" and rebound to 8.42 at "24 h",
while the real monomer peak apex falls monotonically 9.74 -> 2.47 mAU.  Read in the right
order the series is monotonic.  ``verify_auc_vs_notebooks.py`` pins the corrected numbers
to her own retained cell outputs, 47 values bit-exact to 4 dp.

Integration follows her code exactly -- raw absorbance, no baseline, ``np.trapz``, windows
5-11 and 12-18 mL (12-16.1 + 16.1-20 at 25 uM), Total = the sum.  See
``build_figures_S3_S5.py`` for the provenance of every one of those choices.

    python build_figure_S3_combined.py
"""
from __future__ import annotations

import os
import sys

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True

from dataclasses import replace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from secfigs import style as gs                                              # noqa: E402
from secfigs.figure_fig12 import (_NAT_STYLE, MARGIN, GAP_W, GAP_H, TRIM,    # noqa: E402
                                  MM, FIG_W, RC_EXTRA)
from build_figures_S3_S5 import measure, draw                                # noqa: E402
from build_figure1_v2 import LETTER_DX, LETTER_DY                            # noqa: E402

OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

#: letter, registry id, title, split-fragments
PANELS = [("a", "3A", u"250 µM, 37 °C", False),
          ("b", "3B", u"100 µM, 37 °C", False),
          ("c", "3C", u"25 µM, 37 °C",  True),
          ("d", "4A", u"100 µM, 25 °C", False)]
NCOL = 2


def build(stem="Figure_S3", formats=("pdf", "svg", "png")):
    avail = FIG_W - MARGIN["left"] - MARGIN["right"]
    pw = (avail - (NCOL - 1) * GAP_W) / NCOL
    ph = pw / 1.30                                  # a touch wider than tall
    nrow = (len(PANELS) + NCOL - 1) // NCOL
    fh = TRIM + MARGIN["top"] + nrow * ph + (nrow - 1) * GAP_H + MARGIN["bottom"]
    style = replace(_NAT_STYLE, figsize=(FIG_W, fh))

    rows = []
    with gs.apply(style, **RC_EXTRA):
        fig = plt.figure(figsize=(FIG_W, fh), dpi=style.dpi)
        fx, fy = (lambda v: v / FIG_W), (lambda v: v / fh)
        for i, (ch, pid, title, split) in enumerate(PANELS):
            r, c = divmod(i, NCOL)
            x = MARGIN["left"] + c * (pw + GAP_W)
            y = TRIM + MARGIN["top"] + r * (ph + GAP_H)
            ax = fig.add_axes([fx(x), 1.0 - fy(y + ph), fx(pw), fy(ph)])
            fig.text(fx(max(x - LETTER_DX, TRIM)), 1.0 - fy(y - LETTER_DY), ch,
                     fontsize=style.panel_letter_size,
                     fontweight=style.panel_letter_weight, ha="left", va="baseline")
            d = measure(pid, split)
            rows.append((ch, title, d))
            draw(ax, d, style, title)

        paths = []
        for ext in formats:
            p = os.path.join(OUT, f"{stem}.{ext}")
            fig.savefig(p, dpi=600 if ext == "png" else None, facecolor="white")
            paths.append(p)
    plt.close(fig)

    print(f"=== {stem}   {FIG_W/MM:.4f} x {fh/MM:.4f} mm   {len(PANELS)} panels")
    for ch, title, d in rows:
        print(f"\n  {ch}  {title}")
        for j, tt in enumerate(d["t"]):
            frag = f"  frag {d['frag'][j]:7.3f}" if d["split"] else ""
            print(f"     t = {tt:5.1f} h   olig {d['olig'][j]:7.3f}   "
                  f"mono {d['mono'][j]:7.3f}{frag}   total {d['total'][j]:7.3f}")
    print("\n  monotonic check, 100 µM 37 °C (the panel that was mis-ordered):")
    d = dict(rows)[u"100 µM, 37 °C"] if False else rows[1][2]
    print("     oligomer " + " -> ".join(f"{v:.2f}" for v in d["olig"]))
    print("     monomer  " + " -> ".join(f"{v:.2f}" for v in d["mono"]))
    for p in paths:
        print("   wrote", p)
    return paths


if __name__ == "__main__":
    raise SystemExit(build() and 0)
