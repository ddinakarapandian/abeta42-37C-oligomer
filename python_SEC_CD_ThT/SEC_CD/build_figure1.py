"""build_figure1.py — Figure 1.

    a   workflow schematic, full width  (slide 2 of the BioRender deck)
        crude peptide -> SEC purification -> incubation at 37 °C -> reisolation
        (SEC) -> monomers & oligomers -> downstream (CD, ThT, NMR, microglia)
    b   SEC, 100 µM, 24 h reinjection at 6 / 25 / 37 °C          raw mAU
    c   CD after 24 h at 6 / 25 / 37 °C   (manuscript Fig 2, registry 2A)
    d   CD, SEC-isolated oligomers, 37 °C                        (registry 1D)
    e   ThT, fibril against oligomer, each read where it is defined

The story runs: how the material is made (a), that it converts with temperature
by size (b) and by secondary structure (c), what the isolated oligomer looks like
(d), and that the oligomer is not a fibril (e).

STYLE.  Everything uses ``figure_fig12._NAT_STYLE`` -- 7 pt type, 8 pt panel
letters, 0.85 pt traces -- with margins and gaps taken from that module's MEASURED
text-block extents.  Not ``gs.STYLE``: that one is calibrated against her
12 x 9.53 in reference figure (21.5 pt ticks, 54 pt letters) and renders about
2.5x oversized on a 180 mm page, which is exactly how the first draft of this file
came out with its titles colliding through each other.

WHERE EACH PANEL'S DATA COMES FROM.  c and d are drawn by
``build_fig12_panels.draw_panel``, the same function the combined Fig 1+2 page
uses.  b reuses ``build_fig1_panels.load_trace``, i.e. sec3.py's formula, with the
min-max step not taken.  e takes its numbers from ``tht_plateau`` -- the analysis
is shared, only the drawing is local, so this page's bar and the standalone ThT
figure can never disagree about a value while still being styled for their own
pages.

PANEL b CARRIES TWO CORRECTIONS THAT MUST REACH THE CAPTION.  The 6 °C 24 h run
used a 100 uL loop against 500 uL everywhere else, so that trace is multiplied by
the nominal loop ratio, 5 -- an instrument setting, not a fit.  And the three
series are different preps whose own t = 0 monomer AUCs span about 25 %; no
correction is applied for that, because dividing by t = 0 would be a
normalisation and this panel exists not to be one.

    python build_figure1.py
    python build_figure1.py --no-schematic
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
import matplotlib.image as mpimg

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from secfigs import style as gs
from secfigs.figure_fig12 import (_NAT_STYLE, MARGIN, GAP_W, GAP_H, TRIM, MM,
                                  ROW_ASPECT, RC_EXTRA, FIG_W, TEMP_INKS)
import build_fig12_panels as p12
from build_fig1_panels import load_trace

OUT = os.path.join(HERE, "figures")
ASSETS = os.path.join(HERE, "assets")
os.makedirs(OUT, exist_ok=True)

SCHEMATIC = os.path.join(ASSETS, "fig1a_schematic.png")
SCHEM_GAP = 0.150
RES = r"/path/to/Ab42oligomers37C/data/SEC\Run_Results"      # READ ONLY
XLIM_SEC = (6.0, 20.0)
MONOMER_WIN = (12.0, 18.0)
trapz = getattr(np, "trapezoid", np.trapz)

#: panel b — the 24 h reinjection at each temperature, with its loop size
B_TEMPS = [
    ("6 °C",  f"{RES}/E22K/WT 0213 monomer reinj 021425 100ul Loop no fractions  001.csv", 100),
    ("25 °C", f"{RES}/RT runs/24h RT WT reinjection 001.csv", 500),
    ("37 °C", f"{RES}/WT/WT 24hr reinj 37C 022025 500ul Loop no fractions 001.csv", 500),
]
REF_LOOP = 500


# --------------------------------------------------------------- panel drawers

def draw_sec_24h(ax, style):
    """b — 100 µM, 24 h, three temperatures, raw mAU."""
    out = []
    for (label, path, loop), colour in zip(B_TEMPS, TEMP_INKS):
        if not os.path.exists(path):
            raise SystemExit("missing SEC run: %s" % path)
        v, a, _ = load_trace(path, normalized=False)   # sec3's formula, no min-max
        k = REF_LOOP / loop                            # 5.0 at 6 °C, else 1.0
        y = a * k
        s = (v >= MONOMER_WIN[0]) & (v <= MONOMER_WIN[1])
        ax.plot(v, y, color=colour, lw=style.trace_linewidth, label=label,
                solid_capstyle="round")
        out.append(dict(label=label, loop=loop, k=k, peak=float(y.max()),
                        mono=float(trapz(y[s], v[s]))))
    ymax = max(o["peak"] for o in out)
    gs.axes(ax, style, normalized=False, x_major=2.0, x_range=XLIM_SEC,
            y_range=(0.0, ymax * 1.08))
    gs.set_labels(ax, "Elution Volume (ml)", "Absorbance (mAU)", style)
    gs.set_title(ax, u"100 µM, 24 h", style)
    gs.legend(ax, style)
    return out


def draw_tht_bar(ax, style):
    """e — ThT of the fibril, the oligomer and the monomer, each read where it is
    defined: fibril at the monomer well's plateau, oligomer at t = 0, monomer at t = 0.

    The monomer bar was added at A.K.P.'s request. It uses THE SAME WELLS as the
    fibril bar, read before conversion, so the fibril/monomer difference is an
    internal comparison rather than a between-well one.

    Numbers come from tht_plateau/tht_bars so they cannot drift from the standalone
    ThT figure; only the drawing is local, because that module styles for its own page.
    """
    sys.path.insert(0, os.path.join(os.path.dirname(HERE), "37C_ThT_Python"))
    import tht_plateau as tp
    import tht_bars as tb
    df, t = tp.load()
    groups = tb.bar_groups(tp.analyse_all(df, t))
    conds = tb.BAR_CONDITIONS

    xs = np.arange(len(groups), dtype=float)
    rng = np.random.default_rng(0)
    top = 0.0
    names = []
    for x, g, c in zip(xs, groups, conds):
        key = tb.read_key(g)
        vals = g._arr(key)
        m, sd = g.mean_sd(key)
        top = max(top, float(vals.max()), m + sd)
        names.append(c.measured_label or c.label)
        ax.bar(x, m, width=0.56, color=c.light, edgecolor=c.dark,
               lw=style.spine_linewidth, zorder=2)
        ax.errorbar(x, m, yerr=sd, fmt="none", ecolor=c.dark,
                    elinewidth=style.spine_linewidth * 1.6,
                    capsize=2.0, capthick=style.spine_linewidth * 1.6, zorder=3)
        ax.scatter(x + (rng.random(vals.size) - 0.5) * 0.20, vals, s=6.0,
                   facecolor="white", edgecolor=c.dark,
                   linewidth=style.spine_linewidth, zorder=4, clip_on=False)

    ax.set_xticks(xs)
    ax.set_xticklabels(names, fontsize=style.tick_labelsize,
                       fontweight=style.tick_label_weight)
    ax.set_xlim(-0.62, len(groups) - 0.38)
    gs.axes(ax, style, normalized=False, x_major=None, x_minor=None,
            y_range=(0.0, top * 1.30))
    ax.xaxis.set_ticks_position("bottom")
    ax.set_xticks(xs)
    ax.set_xticklabels(names, fontsize=style.tick_labelsize,
                       fontweight=style.tick_label_weight)
    gs.set_labels(ax, "", "ThT fluorescence (mAU)", style)
    gs.set_title(ax, u"25 µM, ThT", style)

    # Two brackets, stacked innermost first so they do not cross: oligomer/monomer
    # sits low over those two bars, fibril/oligomer spans above it.
    out = {}
    for k, (i, j) in enumerate(tb.BRACKETS):
        r = tb.pairwise(groups, i, j)
        out["%s_over_%s" % (conds[i].measured_label, conds[j].measured_label)] = r
        y = top * (1.03 + 0.13 * k)
        lo, hi = sorted((xs[i], xs[j]))
        ax.plot([lo, lo, hi, hi], [y * 0.975, y, y, y * 0.975],
                color="black", lw=style.spine_linewidth, zorder=5)
        ptxt = f"p = {r['p']:.3f}" if r["p"] >= 0.001 else "p < 0.001"
        ax.text((lo + hi) / 2.0, y * 1.01, f"{r['fold']:.1f}×  {ptxt}",
                ha="center", va="bottom", fontsize=style.annotation_size,
                fontweight=style.annotation_weight)
    ax.set_ylim(0.0, top * 1.38)          # headroom for the upper bracket
    return out


#: (letter, how to draw it, row, column)
DATA_PANELS = [
    ("b", ("sec24", None), 0, 0),
    ("c", ("registry", "2A"), 0, 1),
    ("d", ("registry", "1D"), 1, 0),
    ("e", ("tht", None), 1, 1),
]


def schematic_ink(path=SCHEMATIC):
    img = mpimg.imread(path)
    ink = (img[..., :3].min(2) < 0.96) if img.ndim == 3 else (img < 0.96)
    if img.ndim == 3 and img.shape[2] == 4:
        ink &= img[..., 3] > 0.03
    r, c = np.where(ink)
    if not r.size:
        return img.shape[1] / img.shape[0], img
    return (c.max() - c.min() + 1) / (r.max() - r.min() + 1), img


def build(stem="Figure_1", schematic=True, formats=("svg", "pdf", "png")):
    avail_w = FIG_W - MARGIN["left"] - MARGIN["right"]
    panel_w = (avail_w - GAP_W) / 2.0
    row_h = panel_w / ROW_ASPECT[0]

    strip_h, img, asp = 0.0, None, None
    if schematic:
        if not os.path.exists(SCHEMATIC):
            raise SystemExit("schematic not found: %s" % SCHEMATIC)
        asp, img = schematic_ink()
        strip_h = avail_w / asp

    fig_h = (TRIM + (strip_h + SCHEM_GAP if schematic else 0.0)
             + MARGIN["top"] + row_h + GAP_H + row_h + MARGIN["bottom"])

    style = replace(_NAT_STYLE, figsize=(FIG_W, fig_h))
    fx = lambda v: v / FIG_W
    fy = lambda v: v / fig_h
    info = {}

    with gs.apply(style, **RC_EXTRA):
        fig = plt.figure(figsize=(FIG_W, fig_h), dpi=style.dpi)

        if schematic:
            ax = fig.add_axes([fx(MARGIN["left"]), 1.0 - fy(TRIM + strip_h),
                               fx(avail_w), fy(strip_h)])
            ax.imshow(img, interpolation="antialiased", aspect="auto")
            ax.set_axis_off()
            fig.text(fx(TRIM), 1.0 - fy(TRIM), "a",
                     fontsize=style.panel_letter_size,
                     fontweight=style.panel_letter_weight, ha="left", va="top")

        y0 = TRIM + (strip_h + SCHEM_GAP if schematic else 0.0) + MARGIN["top"]
        for letter, (kind, pid), r, c in DATA_PANELS:
            x = MARGIN["left"] + c * (panel_w + GAP_W)
            y = y0 + r * (row_h + GAP_H)
            ax = fig.add_axes([fx(x), 1.0 - fy(y + row_h), fx(panel_w), fy(row_h)])
            if kind == "registry":
                info[letter] = p12.draw_panel(ax, pid, style, letter=False)
            elif kind == "sec24":
                info[letter] = draw_sec_24h(ax, style)
            elif kind == "tht":
                info[letter] = draw_tht_bar(ax, style)
            gs.panel_letter(ax, letter, style)

        paths = []
        for ext in formats:
            q = os.path.join(OUT, f"{stem}.{ext}")
            fig.savefig(q, dpi=600 if ext == "png" else None, facecolor="white")
            paths.append(q)
    plt.close(fig)

    print(f"=== Figure 1   {FIG_W/MM:.4f} x {fig_h/MM:.4f} mm")
    if schematic:
        print(f"   a  schematic (slide 2)  {avail_w/MM:6.2f} x {strip_h/MM:5.2f} mm  "
              f"aspect {asp:.3f}:1   {img.shape[1]/avail_w:.0f} dpi")
    for o in info.get("b", []):
        print(f"   b  {o['label']:6} loop {o['loop']:3d} µL ×{o['k']:.0f}  "
              f"peak {o['peak']:6.2f} mAU  monomer AUC {o['mono']:7.3f}")
    print(f"   c  2A  CD, 24 h at 6/25/37 °C")
    print(f"   d  1D  CD, SEC-isolated oligomers, 37 °C")
    e = info.get("e", {})
    if e:
        lo, hi = e["ci_as_measured"]
        print(f"   e  ThT  {e['fold_as_measured']:.1f}×  p = {e['p_as_measured']:.4f}  "
              f"95% CI {lo:.1f}–{hi:.1f}×")
    tier = 185.0 if fig_h/MM <= 185.0 else (210.0 if fig_h/MM <= 210.0 else 225.0)
    print(f"   page {fig_h/MM:.2f} mm -> {tier:.0f} mm tier")
    for q in paths:
        print("   wrote", q)
    return paths


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-schematic", action="store_true")
    ap.add_argument("--stem", default="Figure_1")
    a = ap.parse_args(argv)
    build(stem=a.stem, schematic=not a.no_schematic)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
