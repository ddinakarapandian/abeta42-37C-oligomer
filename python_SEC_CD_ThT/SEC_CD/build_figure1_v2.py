"""build_figure1_v2.py — Figure 1, refined to the Nature type contract.

Same data as build_figure1.py; this version fixes three defects that a font census
of that file's PDF made explicit.  It writes ``Figure_1_v2.*`` and never touches
``Figure_1.*``.

    a   workflow schematic, full width, LABELS RESET AS LIVE ARIAL TEXT
    b   SEC, 100 µM, 24 h reinjection at 6 / 25 / 37 °C, raw mAU
    c   CD after 24 h at 6 / 25 / 37 °C            (manuscript Fig 2, registry 2A)
    d   CD, SEC-isolated oligomers, 37 °C          (registry 1D)
    e   ThT, fibril against oligomer

WHAT WAS WRONG, MEASURED
------------------------
1. PANEL a's TYPE WAS NEITHER THE RIGHT FONT NOR THE RIGHT SIZE.  The BioRender
   slide sets its labels in **Roboto** at 12.9-13.13 pt bold (headers), 10.75-11.06
   pt (device labels) and 7 pt ("Aggregates").  Placed at 160.39 mm on a 254 mm
   slide the scale factor is 0.6315, so on the page those became:
       headers      8.15 pt   -- OVER the 7 pt ceiling
       labels       6.79 pt   -- legal
       Aggregates   4.42 pt   -- UNDER the 5 pt floor
   Both ends of the contract broken, in a typeface the rest of the page does not
   use.  Rescaling the strip cannot fix this: shrinking it to bring the headers
   under 7 pt pushes "Aggregates" further below 5.
   FIX: export the slide with all 12 text shapes HIDDEN, then re-set every label
   in matplotlib as live Arial at this page's sizes.  The crop is taken at the
   union bounding box of the slide's shapes IN SLIDE COORDINATES, not at the ink,
   so mapping each label's recorded (x, y, w, h) back into the placed image is
   exact rather than inferred.

2. PANEL LETTERS b AND d FELL OFF THE PAGE.  ``gs.panel_letter`` offsets by
   ``panel_letter_xy = (-0.20, 1.13)`` in AXES fractions; -0.20 of a 71.74 mm panel
   is -14.35 mm, and the left column starts at 14.22 mm, so b and d landed at
   x = -0.1 mm -- outside the trim.  Letters also sat at three different x offsets
   (3.3, -0.1, 88.5 mm) because a was placed in figure coordinates and the rest in
   axes fractions.
   FIX: every letter, a included, is placed in FIGURE coordinates at one constant
   inset from its own panel box.

3. TITLE SIZES.  These were already uniform at 7.0 pt -- the census showed 7.0 pt
   x69 and 8.0 pt bold x5 with no third size -- so the apparent mismatch was panel
   a's Roboto, item 1.  Nothing to change, recorded so it is not "fixed" twice.

    python build_figure1_v2.py
    python build_figure1_v2.py --verify      # font census of the delivered PDF
"""
from __future__ import annotations

import argparse
import json
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
from build_figure1 import draw_sec_24h, draw_tht_bar, DATA_PANELS

OUT = os.path.join(HERE, "figures")
ASSETS = os.path.join(HERE, "assets")
ART = os.path.join(ASSETS, "fig1a_artonly_crop.png")
LABELS = os.path.join(ASSETS, "slide2_labels.json")

SCHEM_GAP = 0.150

#: Panel titles set here rather than inherited, because the registry's wording for
#: 2A ("CD of Abeta42 solutions after 24 h at three temperatures") does not say the
#: thing that separates it from d: c is the UNFRACTIONATED mixture, d is the SEC
#: ISOLATE. Read together without that, the two CD panels look like the same
#: measurement twice.
TITLES = {
    "c": u"Unfractionated mixture, 24 h",
    "d": u"SEC-isolated oligomers, 37 °C",
}
#: Panel a's own type scale, chosen inside the contract rather than inherited.
#: The slide's own hierarchy is preserved (headers bold and largest, device labels
#: smaller, "Aggregates" smallest) but every size is now legal Arial.
A_HEADER_PT = 7.0
A_LABEL_PT = 6.0
A_SMALL_PT = 5.5
#: one constant inset for every panel letter, in inches from its panel's top-left
LETTER_DX = 0.285      # horizontal inset, clear of the y-axis block
LETTER_DY = 0.075      # air between the letter's baseline and its panel


def load_labels(path=LABELS):
    d = json.load(open(path, encoding="utf-8"))
    return d["slide"], d["bbox"], d["labels"]


#: The slide writes "37 ° C" with a space before the degree sign; now that the
#: labels are live text that is ours to correct rather than bake into a raster.
TEXT_FIXES = ((u"° C", u"°C"),)


def clean_label(t: str) -> str:
    for bad, good in TEXT_FIXES:
        t = t.replace(bad, good)
    return t


def label_size(pt_native):
    """Slide point size -> this page's size, preserving the slide's hierarchy."""
    if pt_native >= 12.0:
        return A_HEADER_PT, True
    if pt_native >= 9.0:
        return A_LABEL_PT, False
    return A_SMALL_PT, False


def draw_schematic(fig, x0, y_top, w, h, fx, fy):
    """Artwork as a raster, every label re-set as live Arial text on top of it."""
    (SW, SH), (L, T, R, B), labels = load_labels()
    img = mpimg.imread(ART)
    ax = fig.add_axes([fx(x0), 1.0 - fy(y_top + h), fx(w), fy(h)])
    # extent puts the raster in SLIDE INCHES, the same coordinates the recorded
    # label positions are in, so a label drops straight onto the artwork it names.
    # Without it imshow lays the image out in pixels while the limits below are in
    # inches, and the artwork lands entirely outside the view -- which is exactly
    # how the first v2 build came out: labels present, schematic blank.
    ax.imshow(img, extent=(L, R, B, T), interpolation="antialiased", aspect="auto")
    ax.set_xlim(L, R)
    ax.set_ylim(B, T)                 # y grows downward, as in slide coordinates
    ax.set_axis_off()

    used = []
    for lb in labels:
        size, bold = label_size(lb["pt"])
        cx = lb["x"] + lb["w"] / 2.0                      # shapes are centre-aligned
        cy = lb["y"] + lb["h"] / 2.0
        ax.text(cx, cy, clean_label(lb["text"]), ha="center", va="center",
                fontsize=size, fontweight="bold" if bold else "normal",
                family="Arial", color="#000000", linespacing=1.15, zorder=5)
        used.append((lb["text"].replace("\n", " "), lb["pt"], size))
    return ax, used


def build(stem="Figure_1_v2", formats=("svg", "pdf", "png")):
    for p in (ART, LABELS):
        if not os.path.exists(p):
            raise SystemExit("missing asset: %s" % p)

    avail_w = FIG_W - MARGIN["left"] - MARGIN["right"]
    panel_w = (avail_w - GAP_W) / 2.0
    row_h = panel_w / ROW_ASPECT[0]

    (SW, SH), (L, T, R, B), _ = load_labels()
    strip_h = avail_w / ((R - L) / (B - T))

    # HEADROOM FOR PANEL a's LETTER.  b-e sit under MARGIN["top"], which already
    # holds their titles, so their letters have somewhere to go.  The strip has no
    # title, so without this its letter is drawn above the page: the first build
    # put "a" at y = -0.32 mm and left a 0.741 mm top margin against a 3.3 mm
    # budget.  Reserve exactly the letter's cap height plus the same clear air the
    # other letters get.
    cap_in = _NAT_STYLE.panel_letter_size * 0.716 / 72.0
    letter_h = cap_in + LETTER_DY
    strip_top = TRIM + letter_h

    fig_h = (strip_top + strip_h + SCHEM_GAP + MARGIN["top"]
             + row_h + GAP_H + row_h + MARGIN["bottom"])
    style = replace(_NAT_STYLE, figsize=(FIG_W, fig_h))
    fx = lambda v: v / FIG_W
    fy = lambda v: v / fig_h

    def letter(ax_box_x, ax_box_y, ch):
        """One rule for every panel letter, in FIGURE coordinates: the glyph's TOP
        sits ``LETTER_DY`` above its panel, so no letter can overhang the trim."""
        fig.text(fx(max(ax_box_x - LETTER_DX, TRIM)),
                 1.0 - fy(ax_box_y - LETTER_DY), ch,
                 fontsize=style.panel_letter_size,
                 fontweight=style.panel_letter_weight, ha="left", va="baseline")

    info = {}
    with gs.apply(style, **RC_EXTRA):
        fig = plt.figure(figsize=(FIG_W, fig_h), dpi=style.dpi)

        _, used = draw_schematic(fig, MARGIN["left"], strip_top, avail_w,
                                 strip_h, fx, fy)
        letter(MARGIN["left"], strip_top, "a")

        y0 = strip_top + strip_h + SCHEM_GAP + MARGIN["top"]
        for ch, (kind, pid), r, c in DATA_PANELS:
            x = MARGIN["left"] + c * (panel_w + GAP_W)
            y = y0 + r * (row_h + GAP_H)
            ax = fig.add_axes([fx(x), 1.0 - fy(y + row_h), fx(panel_w), fy(row_h)])
            if kind == "registry":
                # d holds one trace, so its legend only restates its title --
                # figure_fig12 suppresses it for the same reason.
                info[ch] = p12.draw_panel(ax, pid, style, letter=False,
                                          title=ch not in TITLES,
                                          legend=(pid != "1D"))
                if ch in TITLES:
                    gs.set_title(ax, TITLES[ch], style)
            elif kind == "sec24":
                info[ch] = draw_sec_24h(ax, style)
            elif kind == "tht":
                info[ch] = draw_tht_bar(ax, style)
            letter(x, y, ch)

        paths = []
        for ext in formats:
            q = os.path.join(OUT, f"{stem}.{ext}")
            fig.savefig(q, dpi=600 if ext == "png" else None, facecolor="white")
            paths.append(q)
    plt.close(fig)

    print(f"=== Figure 1 v2   {FIG_W/MM:.4f} x {fig_h/MM:.4f} mm")
    print(f"   a  schematic {avail_w/MM:.2f} x {strip_h/MM:.2f} mm, "
          f"{mpimg.imread(ART).shape[1]/avail_w:.0f} dpi, labels re-set as live Arial:")
    for t, native, now in sorted(used, key=lambda u: -u[1]):
        print(f"        {native:6.2f} pt Roboto -> {now:.1f} pt Arial   {t[:38]}")
    # Panel e now returns one entry per bracket, not the single fibril/oligomer
    # statistic it used to: the bar panel gained a third state (monomer at t = 0).
    e = info.get("e", {})
    for name, r in e.items():
        a, b = name.split("_over_")
        print(f"   e  ThT {a:8s}/ {b:8s} = {r['fold']:5.1f}×  p = {r['p']:.4f}   "
              f"({r['a_mean']:.1f} ± {r['a_sd']:.1f} vs {r['b_mean']:.1f} ± {r['b_sd']:.1f} mAU)")
    tier = 185.0 if fig_h/MM <= 185.0 else (210.0 if fig_h/MM <= 210.0 else 225.0)
    print(f"   page {fig_h/MM:.2f} mm -> {tier:.0f} mm tier")
    for q in paths:
        print("   wrote", q)
    return paths


def verify(stem="Figure_1_v2"):
    import collections, fitz
    p = fitz.open(os.path.join(OUT, f"{stem}.pdf"))[0]
    sizes, letters = collections.Counter(), []
    for b in p.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                sizes[round(s["size"], 2)] += 1
                if s["text"].strip() in list("abcde"):
                    letters.append((s["text"].strip(), round(s["size"], 2),
                                    round(s["bbox"][0] / 72 * 25.4, 2),
                                    round(s["bbox"][1] / 72 * 25.4, 2), s["font"]))
    pix = p.get_pixmap(dpi=1200)
    a = np.frombuffer(pix.samples, np.uint8).reshape(pix.h, pix.w, pix.n)
    ink = a[:, :, :3].min(2) < 245
    r, c = np.where(ink); k = 25.4 / 1200
    print(f"page {p.rect.width/72*25.4:.4f} x {p.rect.height/72*25.4:.4f} mm")
    print("font sizes:", dict(sorted(sizes.items())))
    bad = [s for s in sizes if (s > 7.01 and s != 8.0) or s < 5.0]
    print("sizes outside 5–7 pt (8 pt letters excepted):", bad or "none")
    print("panel letters:")
    for ch, sz, x, y, fn in sorted(letters, key=lambda t: t[0]):
        print(f"   {ch}  {sz:.1f} pt {fn:16} x={x:7.2f} y={y:7.2f} mm")
    xs = sorted({round(t[2], 2) for t in letters})
    print(f"   distinct letter x positions: {xs}  (expect two, one per column)")
    print(f"margins L {c.min()*k:.3f} R {(pix.w-1-c.max())*k:.3f} "
          f"T {r.min()*k:.3f} B {(pix.h-1-r.max())*k:.3f} mm")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--stem", default="Figure_1_v2")
    a = ap.parse_args(argv)
    build(stem=a.stem)
    if a.verify:
        print()
        verify(a.stem)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
