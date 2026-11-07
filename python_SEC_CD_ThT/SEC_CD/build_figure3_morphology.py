"""build_figure3_morphology.py — Figure 3: silver-stained SDS-PAGE + cryo-EM (+ negative
stain), built to the Nature Communications contract.

    python build_figure3_morphology.py --equal  # *** DELIVERED ***  a gel | b neg | c cryo
                                                #     three side by side, EQUAL HEIGHT, 180 x 59.8 mm
    python build_figure3_morphology.py          #     a gel | b cryo-EM            180 x 81.9 mm
    python build_figure3_morphology.py --row    #     a gel | b neg | c cryo       180 x 84.3 mm
    python build_figure3_morphology.py --three  #     a gel | b neg / c cryo       180 x 107.6 mm

THE GEL'S ANNOTATION IS REDRAWN AS LIVE VECTOR TEXT, and that is what makes --equal
possible.  In the source the ladder and lane labels are BURNT INTO THE RASTER, so their
printed size is set by panel width -- 90 px of glyph in a 4498 px image.  An equal-height
row of three gives the gel only ~62 mm, which renders them at 3.5 pt, under Nature's 5 pt
floor.  (Burnt-in, the gel needs >= 88.2 mm for 5 pt and >= 123.4 mm for 7 pt.)

So the slab is cropped out and the ladder, lane labels and load-group rules are redrawn at
7 pt in SOURCE PIXEL COORDINATES -- see GEL_SLAB / GEL_LADDER / GEL_LANES / GEL_GROUPS and
place_gel_vector().  Panel width then stops constraining legibility.  Measured positions:

    slab      cols 445-4424, rows 133-3182  (3980 x 3050, aspect 1.3049)
    ladder    y = 312 476 619 969 1172 1496 1751 2061 2507  ->  200 150 100 52 37 25 20 15 10
    lanes     x = 658 1526 1978 2430 3275 3726 4178          ->  M SDS S1 S2 SDS S1 S2
    groups    x 1307-2625 and 3055-4374                      ->  "1 uL loaded" / "2 uL loaded"

Finding the slab needed care: a mid-grey test locks onto the darkest lane, so it is found
with a not-white test (a < 248) at 0.60 column/row occupancy.

The other three layouts are kept for comparison; they leave the gel annotation as raster
and so trade legibility against cryo resolution:

    --row     gel pinned at GEL_ROW_MM (92 mm, 5.2 pt); cryo small, so 335 dpi
    --three   columns solved to end level; gel 121 mm (6.9 pt), cryo 258 dpi
    default   gel 89.12 mm (5.1 pt), cryo 82.37 mm at 158 dpi

SOURCES
    a  silver-stained SDS-PAGE   Dropbox/.../Data/SDS/SDS_gel.png            4498 x 3654
    b  negative-stain TEM        Writing/Asset 1.png (Shirin's prepared asset) 1821 x 1820
    c  cryo-EM                   extracted from Shirin's Illustrator file      512 x  512
                                 FromShirin/37Cpaperimage.ai, embedded xref 44

WHY THE CRYO IMAGE COMES OUT OF THE .ai.  Every loose cryo PNG in the project is smaller
or dirtier than the copy embedded in her Illustrator file:

    FromShirin/Finalized cryo-EM Image.png   382 x 402  + burnt-in caption sentence
    FromShirin/image2.png                    346 x 324
    FromShirin/image110nm.png                339 x 339
    Inkscape/CryoEM.png, CryoEMv2.png       1188 x 488  two-panel montage WITH the gridded
                                                        zoom Shirin asked to drop
    .ai embedded xref 44                     512 x 512  cleanest single frame  <-- used

*** RESOLUTION WARNING — READ BEFORE SUBMITTING ***
Nature Communications wants >= 300 dpi for photographic images at final size.  512 px is
not much, so the cryo panel's compliance depends entirely on how wide it is drawn:

    --row      38.79 mm  ->  335 dpi   OK
    --three    50.47 mm  ->  258 dpi   under
    --equal    57.24 mm  ->  227 dpi   under   <-- the delivered layout
    default    82.37 mm  ->  158 dpi   well under

Every build prints its achieved dpi, so the choice is explicit rather than silent.

THE REAL FIX is to re-export the frame from the original micrograph.  Its filename is
burnt into the top of the embedded image and is preserved here:

    25oct08a_00001sq_v01_00002h1_00003ex.mrc

a Leginon/Appion name, so the original lives on the microscope facility store, not in
Dropbox (searched: no .mrc anywhere under Dropbox).  A 2048 px export of that frame would
let the cryo panel run full width at 600 dpi.

WHAT IS CLEANED UP
  * the burnt-in .mrc filename strip is cropped off the top of the cryo frame (internal
    provenance, not for publication)
  * the burnt-in raster "100nm" bar is cropped off the bottom and REDRAWN as vector at
    7 pt Arial.  Scale measured from the original bar: 161 px = 100 nm, i.e. 1.61 px/nm.
  * in --equal the gel's burnt-in annotation is CROPPED OFF and redrawn as live 7 pt
    vector text; in the other layouts it is kept as raster and its printed size is
    reported per build
  * the negative stain keeps Shirin's own 100 nm bar

CONTRACT: exactly 180.0000 mm wide, Arial only, panel letters 8 pt bold in FIGURE
coordinates, no bbox_inches="tight" (it silently changes the width).
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
from matplotlib.patches import Rectangle
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from secfigs.figure_fig12 import FIG_W, MM, TRIM                      # noqa: E402

WRITING = os.path.dirname(HERE)
SRC = r"/path/to/Ab42oligomers37C/lab_main_folder"
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

GEL = f"{SRC}/Data/SDS/SDS_gel.png"
NEG = os.path.join(WRITING, "Asset 1.png")
AI = f"{SRC}/Writing/Inkscape/FromShirin/37Cpaperimage.ai"

# Measured on the 512 x 512 embedded frame:
#   rows 13-24   white band carrying the burnt-in .mrc filename
#   rows 482-489 the raster "100nm" label
#   rows 492-494 the bar itself (row 494 is 161 px of black)
# so the crops must clear row 24 at the top and row 482 at the bottom.  A previous
# 30-row bottom crop stopped at row 481 and left the top pixels of the label showing
# along the panel edge.
CRYO_TOP_CROP = 26
CRYO_BOT_CROP = 36
CRYO_PX_PER_NM = 1.61   # measured: the 100 nm bar spans 161 px (row 494)

LETTER_PT = 8.0
LABEL_PT = 7.0
GAP = 0.075             # inches between image panels
LETTER_PAD = 0.035      # inches between letter baseline and image top
GEL_ROW_MM = 92.0       # gel width in --row mode; 88.2 mm is the 5 pt floor


#: --------------------------------------------------------------------------------
#: Gel geometry, measured off SDS_gel.png (4498 x 3654) so its burnt-in annotation can
#: be cropped away and redrawn as live vector text.  Auto-detection needed care: a
#: mid-grey test locks onto the darkest lane rather than the slab, so the slab is found
#: with a not-white test (a < 248) at a 0.60 column/row occupancy threshold.
GEL_SLAB = (445, 4424, 133, 3182)              # x0, x1, y0, y1 -> 3980 x 3050, aspect 1.3049
GEL_LADDER = [(312, "200"), (476, "150"), (619, "100"), (969, "52"), (1172, "37"),
              (1496, "25"), (2061, "15"), (1751, "20"), (2507, "10")]
GEL_LANES = [(658, "M"), (1526, "SDS"), (1978, "S1"), (2430, "S2"),
             (3275, "SDS"), (3726, "S1"), (4178, "S2")]
GEL_GROUPS = [((1307, 2625), u"1 µL loaded"), ((3055, 4374), u"2 µL loaded")]


def load_rgb(path):
    return np.asarray(Image.open(path).convert("RGB"))


def load_gel_slab():
    """The gel image with its burnt-in ladder/lane annotation cropped off."""
    x0, x1, y0, y1 = GEL_SLAB
    return load_rgb(GEL)[y0:y1 + 1, x0:x1 + 1, :]


def place_gel_vector(fig, x_in, y_top_in, w_in, fh_in, letter="a",
                     gutter_l_in=0.30, gutter_b_in=0.34):
    """Gel panel with the ladder, lane labels and load-group rules redrawn as LIVE TEXT.

    Everything is drawn in SOURCE PIXEL coordinates -- the axes covers the slab plus its
    two gutters, so the measured positions above can be used directly with no rescaling
    arithmetic, and the annotation stays locked to the bands whatever size the panel is.

    This is the point of the exercise: with the labels burnt into the raster their printed
    size is fixed by panel width (90 px of glyph in a 4498 px image), so an equal-height
    row of three would render them at 3.5 pt.  Redrawn as vector they are 7 pt at any width.
    """
    x0, x1, y0, y1 = GEL_SLAB
    slab_w, slab_h = x1 - x0 + 1, y1 - y0 + 1
    img_w_in = w_in - gutter_l_in
    ppi = slab_w / img_w_in                      # source px per inch on the page
    gl, gb = gutter_l_in * ppi, gutter_b_in * ppi
    h_in = slab_h / ppi + gutter_b_in

    fx, fy = (lambda v: v / FIG_W), (lambda v: v / fh_in)
    ax = fig.add_axes([fx(x_in), 1.0 - fy(y_top_in + h_in), fx(w_in), fy(h_in)])
    ax.set_xlim(x0 - gl, x1)
    ax.set_ylim(y1 + gb, y0)                     # inverted: image row order
    ax.imshow(load_gel_slab(), extent=(x0, x1, y1, y0), aspect="auto",
              interpolation="lanczos", zorder=1)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    tick = 0.055 * ppi                           # tick length, inches -> source px
    pad = 0.022 * ppi
    for ty, lab in GEL_LADDER:
        ax.plot([x0 - tick, x0], [ty, ty], color="black", lw=0.5,
                solid_capstyle="butt", zorder=3, clip_on=False)
        ax.text(x0 - tick - pad, ty, lab, fontsize=LABEL_PT, family="Arial",
                ha="right", va="center", zorder=3)
    # Vertical budget below the slab, in inches: a 7 pt line is ~0.097 in tall, so the
    # lane labels occupy 0.045-0.142 and the rule has to clear that.  An earlier version
    # put the labels at 0.10 and the rule at 0.175, which drew the rule THROUGH the
    # descender line of "SDS/S1/S2" and made each lane look individually underlined.
    for lx, lab in GEL_LANES:
        ax.text(lx, y1 + 0.045 * ppi, lab, fontsize=LABEL_PT, family="Arial",
                ha="center", va="top", zorder=3)
    for (gx0, gx1), lab in GEL_GROUPS:
        ry = y1 + 0.170 * ppi
        ax.plot([gx0, gx1], [ry, ry], color="black", lw=0.6,
                solid_capstyle="butt", zorder=3, clip_on=False)
        ax.text((gx0 + gx1) / 2.0, ry + 0.022 * ppi, lab, fontsize=LABEL_PT,
                family="Arial", ha="center", va="top", zorder=3)

    fig.text(fx(x_in), 1.0 - fy(y_top_in - LETTER_PAD), letter,
             fontsize=LETTER_PT, fontweight="bold", family="Arial",
             ha="left", va="baseline")
    return h_in


def load_cryo():
    """Cleanest cryo frame = the raster embedded in Shirin's .ai, minus its burnt-in strips."""
    import fitz
    d = fitz.open(AI)
    best = None
    for pno in range(d.page_count):
        for x in d[pno].get_images(full=True):
            info = d.extract_image(x[0])
            if best is None or info["width"] * info["height"] > best[0] * best[1]:
                best = (info["width"], info["height"], info["image"])
    d.close()
    if best is None:
        raise SystemExit(f"no embedded raster found in {AI}")
    import io
    a = np.asarray(Image.open(io.BytesIO(best[2])).convert("RGB"))
    return a[CRYO_TOP_CROP:a.shape[0] - CRYO_BOT_CROP, :, :]


def place(fig, x_in, y_top_in, w_in, h_in, arr, letter, fh_in, scale_nm=None):
    """Add an image panel at an absolute position, letter in FIGURE coordinates."""
    fx, fy = (lambda v: v / FIG_W), (lambda v: v / fh_in)
    ax = fig.add_axes([fx(x_in), 1.0 - fy(y_top_in + h_in), fx(w_in), fy(h_in)])
    ax.imshow(arr, interpolation="lanczos", aspect="auto")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    fig.text(fx(x_in), 1.0 - fy(y_top_in - LETTER_PAD), letter,
             fontsize=LETTER_PT, fontweight="bold", family="Arial",
             ha="left", va="baseline")
    if scale_nm:
        # Vector scale bar, redrawn from the measured px/nm of the source frame.
        # Black text on an opaque white box, NOT white text with a stroke outline:
        # matplotlib rasterises stroked text into paths, which silently drops the
        # label out of the PDF's text layer (checked with get_text()) and breaks the
        # live-vector-text requirement.  A bbox keeps the glyphs as real text.
        frac = (scale_nm * CRYO_PX_PER_NM) / arr.shape[1]
        x0, y0 = 0.05, 0.06
        ax.add_patch(Rectangle((x0, y0), frac, 0.022, transform=ax.transAxes,
                               facecolor="white", edgecolor="none", zorder=5))
        ax.text(x0 + frac / 2.0, y0 + 0.040, f"{scale_nm} nm", transform=ax.transAxes,
                fontsize=LABEL_PT, family="Arial", color="black", ha="center",
                va="bottom", zorder=6,
                bbox=dict(facecolor="white", edgecolor="none", pad=0.8))
    return ax


def build(three_panel=False):
    gel, cryo = load_rgb(GEL), load_cryo()
    ar = lambda a: a.shape[1] / a.shape[0]           # noqa: E731  width/height

    avail = FIG_W - 2 * TRIM
    lp = LETTER_PT * 0.716 / 72.0 + LETTER_PAD
    if not three_panel:
        # DEFAULT: one row, matched heights — gel and cryo-EM, as requested.
        h = (avail - GAP) / (ar(gel) + ar(cryo))
        wg, wc = ar(gel) * h, ar(cryo) * h
        top = TRIM + lp
        fh = top + h + TRIM
        fig = plt.figure(figsize=(FIG_W, fh))
        place(fig, TRIM, top, wg, h, gel, "a", fh)
        place(fig, TRIM + wg + GAP, top, wc, h, cryo, "b", fh, scale_nm=100)
        panels = [("a gel", gel, wg), ("b cryo-EM", cryo, wc)]
        stem = "Figure3_morphology"
    elif three_panel == "equal":
        # THREE SIDE BY SIDE, EQUAL HEIGHT, tops AND bottoms flush -- no dead space.
        # Only possible because the gel's annotation is now vector (place_gel_vector):
        # with it burnt into the raster this layout put the labels at 3.5 pt.
        #
        # Solve for the common block height H.  The gel panel is the slab plus a left
        # gutter (ladder) and a bottom gutter (lane labels + rules), so its width is
        #     W_gel = gutter_l + ar_slab * (H - gutter_b)
        # and the micrographs are simply W = ar * H.  Widths therefore differ -- the
        # sources have different aspects -- but every panel is the same height, which is
        # what reads as "equal" across a row.
        neg = load_rgb(NEG)
        slab_ar = (GEL_SLAB[1] - GEL_SLAB[0] + 1) / (GEL_SLAB[3] - GEL_SLAB[2] + 1)
        GL, GB = 0.30, 0.34
        row = avail - 2 * GAP
        H = (row - GL + slab_ar * GB) / (slab_ar + ar(neg) + ar(cryo))
        wg = GL + slab_ar * (H - GB)
        wn, wc = ar(neg) * H, ar(cryo) * H
        top = TRIM + lp
        fh = top + H + TRIM
        fig = plt.figure(figsize=(FIG_W, fh))
        place_gel_vector(fig, TRIM, top, wg, fh, "a", GL, GB)
        xn = TRIM + wg + GAP
        place(fig, xn, top, wn, H, neg, "b", fh)
        place(fig, xn + wn + GAP, top, wc, H, cryo, "c", fh, scale_nm=100)
        panels = [("a gel (vector labels)", load_gel_slab(), wg - GL),
                  ("b negative stain", neg, wn), ("c cryo-EM", cryo, wc)]
        stem = "Figure3_morphology_equal"
        print("!! negative-stain panel: its '100 nm' annotation is a LABEL BOX WITH NO BAR")
        print("   in the source asset, so no scale can be verified or redrawn.")
    elif three_panel == "row":
        # ALL THREE IN ONE HORIZONTAL ROW.
        #
        # Matched heights do NOT work here: the row would give the gel only 61.96 mm,
        # and its annotation is burnt-in raster, so the ladder and lane labels would
        # render at 3.5 pt -- below Nature's 5 pt floor.  (The gel needs >= 88.2 mm for
        # 5 pt, >= 123.4 mm for 7 pt.)  So the gel is pinned at GEL_ROW_MM and the two
        # micrographs share what is left, vertically centred against it.  Ragged bottoms,
        # but every panel is legible and the cryo frame actually gains resolution over
        # the two-panel form (335 dpi against 158) because its panel is smaller.
        neg = load_rgb(NEG)
        wg = GEL_ROW_MM * MM
        wm = (avail - 2 * GAP - wg) / 2.0
        hg, hn, hc = wg / ar(gel), wm / ar(neg), wm / ar(cryo)
        h = max(hg, hn, hc)
        top = TRIM + lp
        fh = top + h + TRIM
        fig = plt.figure(figsize=(FIG_W, fh))
        # TOP-aligned, not vertically centred: centring drops the b and c letters below
        # a's baseline, which reads as a misalignment rather than a design choice.  With
        # tops flush all three letters sit on one line and the slack collects as a single
        # block of whitespace under the two micrographs.
        place(fig, TRIM, top, wg, hg, gel, "a", fh)
        xn = TRIM + wg + GAP
        place(fig, xn, top, wm, hn, neg, "b", fh)
        place(fig, xn + wm + GAP, top, wm, hc, cryo, "c", fh, scale_nm=100)
        panels = [("a gel", gel, wg), ("b negative stain", neg, wm), ("c cryo-EM", cryo, wm)]
        stem = "Figure3_morphology_row"
        print("!! negative-stain panel: its '100 nm' annotation is a LABEL BOX WITH NO BAR")
        print("   in the source asset, so no scale can be verified or redrawn.")
    else:
        # gel left, negative stain + cryo stacked right.  The right-column width is
        # SOLVED so the two columns finish at the same height -- an earlier version
        # fixed it at 35.5 % of the row and left ~28 mm of dead white under the gel.
        #   (avail - GAP - wr)/ar_gel = wr*(1/ar_neg + 1/ar_cryo) + lp
        neg = load_rgb(NEG)
        K = 1.0 / ar(neg) + 1.0 / ar(cryo)
        wr = (avail - GAP - ar(gel) * lp) / (1.0 + ar(gel) * K)
        wg = avail - GAP - wr
        hg = wg / ar(gel)
        hn, hc = wr / ar(neg), wr / ar(cryo)
        top = TRIM + lp
        fh = top + hg + TRIM
        fig = plt.figure(figsize=(FIG_W, fh))
        place(fig, TRIM, top, wg, hg, gel, "a", fh)
        xr = TRIM + wg + GAP
        place(fig, xr, top, wr, hn, neg, "b", fh)
        place(fig, xr, top + hn + lp, wr, hc, cryo, "c", fh, scale_nm=100)
        panels = [("a gel", gel, wg), ("b negative stain", neg, wr), ("c cryo-EM", cryo, wr)]
        stem = "Figure3_morphology_3panel"
        print("!! negative-stain panel: its '100 nm' annotation is a LABEL BOX WITH NO BAR")
        print("   in the source asset, so no scale can be verified or redrawn. Shirin needs")
        print("   to supply the bar length before this panel is publishable.")

    paths = []
    for ext in ("pdf", "svg", "png"):
        p = os.path.join(OUT, f"{stem}.{ext}")
        fig.savefig(p, dpi=600 if ext == "png" else None, facecolor="white")
        paths.append(p)
    plt.close(fig)

    print(f"=== {stem}   {FIG_W/MM:.4f} x {fh/MM:.4f} mm")
    for name, arr, w in panels:
        dpi = arr.shape[1] / w
        flag = "OK" if dpi >= 300 else "*** BELOW NATURE'S 300 dpi MINIMUM ***"
        print(f"   {name:20s} {arr.shape[1]:5d} x {arr.shape[0]:5d} px  "
              f"panel {w/MM:6.2f} mm  ->  {dpi:6.0f} dpi   {flag}")
        if name.startswith("a gel (vector"):
            print(f"   {'':20s} ladder/lane/group labels are LIVE VECTOR TEXT at "
                  f"{LABEL_PT:.0f} pt — panel width no longer constrains legibility")
        elif name.startswith("a gel"):
            # burnt-in annotation: its printed size is set by the panel width
            pt = 90.0 / arr.shape[1] * w * 72.0
            print(f"   {'':20s} burnt-in gel labels render at ~{pt:.1f} pt "
                  f"({'ok, 5-7 pt window' if 5 <= pt <= 7.6 else 'OUTSIDE the 5-7 pt window'})")
    cw = [w for n, _a, w in panels if "cryo" in n][0]
    print(f"   cryo frame: filename strip and raster bar cropped; 100 nm bar redrawn as vector")
    print(f"   to reach 300 dpi at this panel width the cryo source needs "
          f"{int(round(300 * cw))} px (it has {cryo.shape[1]})")
    print("   original micrograph: 25oct08a_00001sq_v01_00002h1_00003ex.mrc — a Leginon name,")
    print("   so the full-resolution frame is on the facility store, not in Dropbox.")
    for p in paths:
        print("   wrote", p)
    return paths


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--three", action="store_true",
                    help="add the negative-stain panel, gel left / micrographs stacked right")
    ap.add_argument("--equal", action="store_true",
                    help="THREE SIDE BY SIDE, equal height, vector gel labels")
    ap.add_argument("--row", action="store_true",
                    help="all three in ONE HORIZONTAL ROW (gel pinned wide enough to stay legible)")
    a = ap.parse_args(argv)
    build(three_panel=("equal" if a.equal else "row" if a.row else a.three))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
