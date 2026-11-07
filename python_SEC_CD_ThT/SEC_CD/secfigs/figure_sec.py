# -*- coding: utf-8 -*-
r"""
secfigs.figure_sec -- the combined SEC display item for the 37 C Abeta42
oligomer manuscript.

Merges the former Figure 3 (A/B/C) and Figure 4 (A/B) into ONE Nature
Communications figure, drawn as a single native matplotlib canvas via a
GridSpec.  No PNGs are pasted; every panel is drawn by the same
``secfigs.style`` primitives (``gs.axes`` / ``gs.set_labels`` / ``gs.set_title``
/ ``gs.legend`` / ``gs.annotate_peak``) that ``build_fig3_panels.py`` and
``build_fig4_sec_panels.py`` call for the standalone panels, and the CSVs are
read through the identical 4-step contract.

Panel order and rationale  (2 x 2, revised 2026-08-07)
------------------------------------------------------
    a   250 uM, 37 C   (was Fig 3A)   3 traces   2, 4, 6 h
    b    25 uM, 37 C   (was Fig 3C)   7 traces   0, 3, 6, 12, 24, 48, 72 h
    c   100 uM, 37 C   (was Fig 4B)   7 traces   0, 2, 4, 8, 12, 24, 48 h
    d   100 uM, 25 C   (was Fig 4A)   9 traces   2, 4, 6, 8, 10, 12, 24, 48, 72 h

The previous five-panel version carried 100 uM / 37 C TWICE -- old panel b
(raw mAU, from Fig 3B) and old panel e (normalised, from Fig 4B) are the SAME
SEVEN CSV FILES, verified 7 of 7.  The duplicate is dropped and the single
surviving 100 uM / 37 C panel is used as a PIVOT that belongs to both stories:

    a + b + c   concentration series at one temperature (250 / 25 / 100 uM)
    c + d       temperature pair at one concentration (37 C vs 25 C),
                side by side in the bottom row, identically scaled

Four equal tiles, not the uneven 3 + 2 of the old six-column grid: with only
four panels there is no block structure left to signal, and the bottom-row
temperature pair must be the same size for a reader to compare peak heights
directly.

EVERY PANEL IS MIN-MAX NORMALISED
---------------------------------
Forced by the layout: c and d must share a y axis for the temperature
comparison to be like-for-like, and c is also a member of the concentration
series, so a and b must match it.  It is also the better readout -- raw
absorbance across 250 / 100 / 25 uM mostly encodes how much material was
loaded, whereas the claim is about the oligomer FRACTION.

``normalized=True`` is a property of THIS FIGURE and is passed at draw time.
``secfigs.panel_data.PANELS`` is left untouched: it remains the record of what
the user's own originals were (3A/3C raw, 4A/4B normalised), and ``load()``
carries that original flag through to ``orig_normalized`` for the report.

Shared axes
-----------
x  (6, 20) ml on all four.  Every trace's maximum lies in 8.49 - 14.57 ml, so
   the window clips no peak in any panel; outside it every trace is <= 0.106 of
   full scale (the sub-6 ml void region, excluded in the originals too) and
   <= 0.06 beyond 20 ml.  20 ml rather than the old 18 ml of the Fig-4 pair
   because panel b grows a late-eluting species at ~17 ml that reaches 0.295 at
   72 h and only decays to baseline by ~19 ml; an 18 ml right edge would cut its
   shoulder.  Major ticks land on 6, 8, ..., 20.
y  (-0.05, 1.05) on all four, majors 0.0 - 1.0 by 0.2, minors 0.1 -- the Griffin
   normalised convention, identical everywhere, so the four panels are
   comparable at a glance without reading the tick labels.

Because the y axis is now literally the same axis four times, the y AXIS LABEL
is drawn once per row, on the left column only (panels a and c).  The y tick
labels are kept on all four so no panel has to be read across the gutter.
The Oligomer / Monomer pointers are likewise drawn ONCE, on panel a.

Nature Communications geometry
------------------------------
    page width  180.00 mm  exactly  (7.08661 in)
    page height 159.51 mm           (<= 230 mm)
    raster dpi  600                 (>= 300)
    fonts       Arial, 5.5 - 8.0 pt
    SVG text    live (svg.fonttype = 'none'); PDF text live (pdf.fonttype = 42)
    NO bbox_inches='tight' anywhere -- the margins below are absolute inches and
    the page stays exactly 180 mm wide.

Legend decision
---------------
PER-PANEL legends, not a shared one.  The four trace sets do not coincide:
a is {2,4,6 h}; b is {0,3,6,12,24,48,72 h}; c is {0,2,4,8,12,24,48 h};
d is {2,4,6,8,10,12,24,48,72 h}.  No two panels share a set, so a figure-level
legend would have to list the union (11 entries, most absent from any given
panel).  Four small in-axes legends cost ~19% of each panel's width in a region
(17.5-20 ml) where every trace is on baseline.

Colour
------
One ordered ramp for the whole figure, light -> dark with increasing incubation
time, built from Griffin Figure 1A's three hue families (coral / teal / gold),
defined once in ``secfigs.style`` (RAMP9 / RAMP3 / ramp).  Colour tracks
POSITION IN THE TIME SERIES, not absolute time: b and c both take the first
seven of RAMP9 even though their time points differ, so the same colour means
"nth time point", and every panel reads light -> dark as incubation proceeds.

Usage
-----
    from secfigs.figure_sec import CombinedSECFigure
    CombinedSECFigure().build()          # draws + saves svg/png/pdf
    python -m secfigs.figure_sec         # same, plus the geometry report
"""

from __future__ import annotations

import os
import sys
from dataclasses import replace

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from secfigs import style as gs                       # noqa: E402
from secfigs.panel_data import PANELS                 # noqa: E402

# The loader is imported from the Figure-3 builder so that the combined figure
# provably consumes the CSVs through the very same code path as the standalone
# panels (pd.read_csv utf-16/tab/skiprows=3 -> first two numeric columns ->
# truncate 0..25 ml -> subtract min only if min < 0 -> optional min-max).
from build_fig3_panels import load_csv                # noqa: E402


# ==========================================================================
# 1.  page geometry -- every number in INCHES, converted to figure fractions
# ==========================================================================
MM = 1.0 / 25.4

FIG_W_MM = 180.0                       # Nature Communications double column
FIG_W = FIG_W_MM * MM                  # 7.086614 in

# right=0.17 in is NOT slack: the last x tick label ("20") is centred on the
# right spine and overhangs it by half its width, 0.046 in / 28 px at 600 dpi.
# With right=0.05 the page had only 3 px (0.127 mm) of clear trim margin;
# 0.17 in restores ~3.1 mm.  left=0.52 in is set by the panel letter, which sits
# 0.42 in left of the axes box (the y label needs only 0.30 in).
MARGIN = dict(left=0.52, right=0.17, top=0.26, bottom=0.40)   # in

# GAP_W must clear, right to left: the right panel's letter (0.42 in out from
# its axes edge, ~0.09 in wide at 8 pt bold) and then trim.  0.62 in leaves
# 0.11 in between the left panel's right spine and the right panel's letter.
GAP_W = 0.62                           # in, between the two columns
# GAP_H must clear the top row's x tick labels + x axis label (~0.36 in) and the
# bottom row's title (~0.16 in).  0.62 in leaves ~0.10 in clear.
GAP_H = 0.62                           # in, between the two rows
ROW_H = 2.50                           # in, axes height -- SAME for both rows

FIG_H = MARGIN["top"] + 2 * ROW_H + GAP_H + MARGIN["bottom"]
FIG_H_MM = FIG_H / MM                  # 159.512 mm

NCOL = 2
NROW = 2
_AVAIL_W = FIG_W - MARGIN["left"] - MARGIN["right"]
PANEL_W = (_AVAIL_W - (NCOL - 1) * GAP_W) / NCOL       # 2.8883 in
WSPACE = GAP_W / PANEL_W
HSPACE = GAP_H / ROW_H


# ==========================================================================
# 2.  the small-format Griffin style
# ==========================================================================
# gs.STYLE is calibrated for a 12 x 9.53 in single panel.  scaled() takes the
# whole dataclass down uniformly; the overrides below then lift the strokes and
# the tick pad back off the floor, because a uniform 0.279x would put the spines
# at 0.17 pt (below the 0.25 pt print minimum) and the tick pad at 1.0 pt.
_SCALE = 6.0 / gs.STYLE.tick_labelsize          # 6 pt tick labels -> 0.2791

SMALL = replace(
    gs.STYLE.scaled(_SCALE),
    figsize=(FIG_W, FIG_H),
    dpi=600,
    # -- typography, all inside Nature's 5-8 pt body-text band ---------------
    tick_labelsize=6.0,
    tick_pad=1.8,
    axis_label_size=7.0,
    axis_label_pad=3.0,
    title_size=7.0,
    title_pad=3.5,
    legend_fontsize=5.5,
    annotation_size=6.0,
    panel_letter_size=8.0,
    # -- strokes ------------------------------------------------------------
    spine_linewidth=0.50,
    trace_linewidth=0.85,
    major_tick_length=2.40,
    minor_tick_length=1.40,
    major_tick_width=0.50,
    minor_tick_width=0.40,
    # -- legend -------------------------------------------------------------
    legend_handlelength=1.50,
    legend_handletextpad=0.45,
    legend_labelspacing=0.25,
    legend_borderpad=0.30,
    # -- peak pointers ------------------------------------------------------
    # offsets scale with the annotation size: (2.5, 40.8) pt at 24 pt text.
    annot_marker_size=2.2,
    annot_offset_right=(0.63, 10.20),
    annot_offset_left=(-5.08, 8.58),
)

RC_EXTRA = {
    "svg.fonttype": "none",            # live text in the SVG
    "pdf.fonttype": 42,                # live TrueType text in the PDF
    "ps.fonttype": 42,
    "svg.hashsalt": "figure_sec",      # deterministic element ids
    "legend.columnspacing": 0.8,
    "path.simplify": False,            # keep every raw vertex
}


# ==========================================================================
# 3.  colour ramp
# ==========================================================================
#: The ramp lives in secfigs.style so the standalone panel builders and this
#: combined figure cannot drift apart.  They previously each carried their own
#: copy, and the same time point ended up a different colour in
#: figures/panels/3B.svg and in the combined figure.
RAMP3 = gs.RAMP3
RAMP9 = gs.RAMP9
ramp = gs.ramp


# ==========================================================================
# 4.  peak-pointer geometry (fractions, so it transfers to any panel)
# ==========================================================================
# Griffin Figure 1A, xlim (6,18) span 12:
#   oligomer apex (8.60, 0.94) -> triangle (9.33, 0.45)   dx=+0.0608*span y=0.478*h
#   monomer  apex (14.30, 1.00) -> triangle (13.62, 0.77) dx=-0.0567*span y=0.770*h
#
# With every panel normalised, the oligomer apex and the monomer apex are BOTH
# exactly 1.000 in panel a (the 4 h / 6 h traces peak on the oligomer, the 2 h
# trace peaks on the monomer), so 0.478 and 0.770 of h separate the two labels
# by 0.29 of full scale on their own.  The y_frac override that the raw-mAU
# version of panel a needed (its oligomer apex was 1.4x its monomer apex, which
# collapsed 0.478*h and 0.770*h to within 1.5 mAU) is therefore GONE.
# Griffin Figure 1A places BOTH triangles as ">", and its oligomer triangle sits
# 0.72 mL to the RIGHT of the oligomer apex (measured: pointer 9.32 mL, apex
# 8.60 mL).  Copied literally that makes both arrows point right and both land on
# the monomer -- an outright mislabel, not a style quirk.  The oligomer pointer is
# therefore mirrored: placed to the LEFT of its apex... no -- placed to the RIGHT
# and flipped to "<", so it still occupies Griffin's position in the inter-peak
# gap but points back at the peak it names.
POINTER = {
    "Oligomer": dict(dx_frac=+0.0608, y_frac=0.478, side="right", marker="<"),
    "Monomer":  dict(dx_frac=-0.0567, y_frac=0.770, side="left",  marker=">"),
}
OLIG_WIN = (7.0, 11.5)
MONO_WIN = (12.0, 18.0)

POINTER_OVERRIDE = {}


def _apex(traces, win):
    """(volume, height) of the tallest peak in `win` over all traces."""
    best = (float("nan"), -np.inf)
    for V, A in traces:
        m = (V >= win[0]) & (V <= win[1])
        if not m.any():
            continue
        i = int(np.argmax(A[m]))
        if A[m][i] > best[1]:
            best = (float(V[m][i]), float(A[m][i]))
    return best


def _short(label: str) -> str:
    """'t = 24 h' -> '24 h'; 't = 0' -> '0 h'.  Legends must be narrow here."""
    s = label.replace("t =", "").strip()
    return "0 h" if s == "0" else s


# ==========================================================================
# 5.  panel specifications
# ==========================================================================
MU = u"µ"
DEG = u"°"

#: Shared axes -- see the module docstring for the evidence behind both.
XLIM = (6.0, 20.0)
YLIM = (-0.05, 1.05)
X_MAJOR = 2.0
Y_MAJOR = 0.2
Y_MINOR = 0.1

YLABEL = u"Normalized absorbance"
XLABEL = u"Elution volume (mL)"

PANEL_SPEC = [
    dict(letter="a", pid="3A", grid=(0, 0),
         title=u"250 %sM, 37 %sC" % (MU, DEG),
         ylabel=YLABEL, annotate=True),
    dict(letter="b", pid="3C", grid=(0, 1),
         title=u"25 %sM, 37 %sC" % (MU, DEG),
         ylabel=None, annotate=False),
    dict(letter="c", pid="4B", grid=(1, 0),
         title=u"100 %sM, 37 %sC" % (MU, DEG),
         ylabel=YLABEL, annotate=False),
    dict(letter="d", pid="4A", grid=(1, 1),
         title=u"100 %sM, 25 %sC" % (MU, DEG),
         ylabel=None, annotate=False),
]

# Panel-letter placement, in INCHES relative to the axes box corner.  Anchored
# in FIGURE coordinates at draw time (see _place_letters) so a narrow panel does
# not push its letter off the canvas the way an axes-fraction offset would.
LETTER_DX_IN = 0.42        # left of the axes left edge
LETTER_DY_IN = 0.045       # above the axes top edge (baseline)


# ==========================================================================
# 6.  the figure
# ==========================================================================
class CombinedSECFigure(object):
    """One Nature-format display item carrying all four SEC panels."""

    def __init__(self, style=SMALL, outdir=None, stem="Figure_SEC_Combined"):
        self.style = style
        self.outdir = outdir or os.path.join(_REPO, "figures")
        self.stem = stem
        self.fig = None
        self.axes = {}
        self.data = {}
        self.pointers = {}

    # -- data ------------------------------------------------------------
    def load(self):
        """Read every CSV through the documented 4-step contract.

        ``normalize=True`` unconditionally: this figure normalises all four
        panels (see module docstring).  ``PANELS[pid]['normalized']`` is NOT
        edited -- it is preserved here as ``orig_normalized`` so the report can
        still say which of the user's originals were raw.
        """
        for spec in PANEL_SPEC:
            P = PANELS[spec["pid"]]
            traces = [load_csv(f, normalize=True) for f in P["csvs"]]
            labels = [_short(l) for l in P["labels"]]
            colors = ramp(len(traces), self.style)
            assert len(set(colors)) == len(traces), "colours must be distinct"
            self.data[spec["letter"]] = dict(
                traces=traces, labels=labels, colors=colors,
                csvs=list(P["csvs"]), raw_labels=list(P["labels"]),
                orig_normalized=P["normalized"])
        return self

    # -- one panel -------------------------------------------------------
    def _draw_panel(self, ax, spec):
        d = self.data[spec["letter"]]
        st = self.style

        for (V, A), c, lab in zip(d["traces"], d["colors"], d["labels"]):
            ax.plot(V, A, color=c, lw=st.trace_linewidth, label=lab,
                    solid_capstyle="round", solid_joinstyle="round")

        gs.axes(ax, style=st, normalized=True, x_major=X_MAJOR,
                y_major=Y_MAJOR, y_minor=Y_MINOR,
                x_range=XLIM, y_range=YLIM)
        gs.set_labels(ax, XLABEL, spec["ylabel"], style=st)
        gs.set_title(ax, spec["title"], style=st)
        gs.legend(ax, style=st)

        if spec["annotate"]:
            span = XLIM[1] - XLIM[0]
            for text, win in (("Oligomer", OLIG_WIN), ("Monomer", MONO_WIN)):
                q = dict(POINTER[text])
                q.update(POINTER_OVERRIDE.get((spec["letter"], text), {}))
                pv, ph = _apex(d["traces"], win)
                x = pv + q["dx_frac"] * span
                y = q["y_frac"] * ph
                # marker comes from POINTER, not hardcoded: the two pointers sit on
                # opposite sides of their peaks, so they cannot share a glyph without
                # one of them aiming at the wrong species.
                gs.annotate_peak(ax, x, y, text, side=q["side"], style=st,
                                 marker=q["marker"])
                self.pointers[(spec["letter"], text)] = (x, y, pv, ph)

    # -- panel letters, anchored in FIGURE coordinates --------------------
    def _place_letters(self):
        fw, fh = self.style.figsize
        for spec in PANEL_SPEC:
            pos = self.axes[spec["letter"]].get_position()
            self.fig.text(pos.x0 - LETTER_DX_IN / fw,
                          pos.y1 + LETTER_DY_IN / fh,
                          spec["letter"],
                          transform=self.fig.transFigure,
                          fontsize=self.style.panel_letter_size,
                          fontweight=self.style.panel_letter_weight,
                          ha="left", va="baseline",
                          color=self.style.foreground)

    # -- draw ------------------------------------------------------------
    def draw(self):
        if not self.data:
            self.load()
        st = self.style
        fw, fh = st.figsize

        with gs.apply(style=st, **RC_EXTRA):
            fig = plt.figure(figsize=(fw, fh), dpi=st.dpi)
            grid = GridSpec(
                NROW, NCOL, figure=fig,
                left=MARGIN["left"] / fw,
                right=1.0 - MARGIN["right"] / fw,
                top=1.0 - MARGIN["top"] / fh,
                bottom=MARGIN["bottom"] / fh,
                wspace=WSPACE, hspace=HSPACE,
            )
            for spec in PANEL_SPEC:
                r, c = spec["grid"]
                ax = fig.add_subplot(grid[r, c])
                self.axes[spec["letter"]] = ax
                self._draw_panel(ax, spec)

            self.fig = fig
            self._place_letters()
        return fig

    # -- save ------------------------------------------------------------
    def save(self, formats=("svg", "png", "pdf")):
        if self.fig is None:
            self.draw()
        os.makedirs(self.outdir, exist_ok=True)
        out = {}
        with gs.apply(style=self.style, **RC_EXTRA):
            for ext in formats:
                p = os.path.join(self.outdir, "%s.%s" % (self.stem, ext))
                # NO bbox_inches -- the page must stay exactly 180 mm wide.
                self.fig.savefig(p, dpi=self.style.dpi, format=ext,
                                 facecolor="white", edgecolor="none")
                out[ext] = p
        return out

    def build(self, formats=("svg", "png", "pdf")):
        self.load()
        self.draw()
        files = self.save(formats)
        plt.close(self.fig)
        return files

    # -- geometry report --------------------------------------------------
    def geometry(self):
        fw, fh = self.style.figsize
        g = dict(page_mm=(fw / MM, fh / MM), dpi=self.style.dpi,
                 px=(int(round(fw * self.style.dpi)),
                     int(round(fh * self.style.dpi))),
                 margin_in=dict(MARGIN), gap_in=(GAP_W, GAP_H),
                 xlim=XLIM, ylim=YLIM, panels={})
        for spec in PANEL_SPEC:
            p = self.axes[spec["letter"]].get_position()
            g["panels"][spec["letter"]] = dict(
                x0_in=p.x0 * fw, y0_in=p.y0 * fh,
                w_in=p.width * fw, h_in=p.height * fh,
                aspect=(p.width * fw) / (p.height * fh),
                n_traces=len(self.data[spec["letter"]]["traces"]))
        return g


# ==========================================================================
# 7.  cli
# ==========================================================================
def main():
    fig = CombinedSECFigure()
    fig.load()
    fig.draw()
    files = fig.save()
    g = fig.geometry()
    print("page  %.3f x %.3f mm  @ %d dpi  = %d x %d px"
          % (g["page_mm"][0], g["page_mm"][1], g["dpi"], g["px"][0], g["px"][1]))
    print("x %s   y %s" % (g["xlim"], g["ylim"]))
    for k, v in g["panels"].items():
        print("  %s  %.4f x %.4f in  aspect %.3f  at (%.4f, %.4f) in  n=%d"
              % (k, v["w_in"], v["h_in"], v["aspect"], v["x0_in"], v["y0_in"],
                 v["n_traces"]))
    for k, v in files.items():
        print("  %-4s %s" % (k, v))
    plt.close(fig.fig)
    return files


if __name__ == "__main__":
    main()
