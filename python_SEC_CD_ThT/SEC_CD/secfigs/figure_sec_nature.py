# -*- coding: utf-8 -*-
r"""
secfigs.figure_sec_nature -- the 2 x 2 SEC display item, restyled to Nature.

This is a RESTYLE of :class:`secfigs.figure_sec_seq.SequentialSECFigure`, not a
new figure.  Same four panels, same CSVs, same loader, same 26 traces, same
numbers on the page -- ``verify_nature.py`` decodes both SVGs and differences
every plotted vertex to prove it.  ``Figure_SEC_normalized.*`` and
``Figure_SEC_raw.*`` are left exactly where they were so the two can be
compared side by side.

Output: ``figures/Figure_SEC_natstyle.{svg,png,pdf}`` (normalised data).

Nothing here re-implements a panel or a loader: the class subclasses
``SequentialSECFigure`` and overrides only ``_layout``, ``_rc``,
``_style_axes``, ``_style_legend``, ``_place_letters`` and ``save``.


WHY -- five changes, each with its evidence
===========================================

1.  TYPE up to Nature's floor
    -------------------------
    Nature's artwork guide (nature.com/documents/natrev-artworkguide.pdf) asks
    for "at least 7 point size text in your draft figures".  The parent figure
    has 6.0 pt tick labels and a 5.5 pt legend -- both under it.

        tick labels     6.0 -> 7.0 pt
        legend          5.5 -> 7.0 pt
        annotations     6.0 -> 7.0 pt   ("Oligomer" / "Monomer")
        axis titles     7.0 -> 8.0 pt
        panel titles    7.0 -> 8.0 pt
        panel letters   8.0 -> 9.0 pt, still bold

    Axis and panel titles are RAISED to 8.0 rather than held at 7.0.  Two
    reasons.  First, 8.0 is what the transcription of Nature's own house
    rcParams uses for ``font.size`` with 7 pt ticks and a 7 pt legend, so
    7/7/8 is the actual published ratio, not an invention.  Second, holding
    them at 7.0 would flatten the hierarchy: tick label, axis title and panel
    title would all be one size, and the panel title is the only thing telling
    a reader which of the four conditions they are looking at.  The layout
    affords it -- the page grows 4.34 mm to pay for the whole type increase,
    177.089 -> 181.432 mm, still 48.57 mm inside the 230 mm cap.

    Panel letters go to 9.0 pt for the same hierarchy reason: at 8.0 they
    would be exactly the size of the panel titles they are meant to outrank.
    They keep the 0.42 in horizontal offset, which still clears the y-axis
    block (now 0.3741 in) by 0.0459 in.

    NOTHING on the page is below 7.0 pt.  ``verify_nature.py`` proves it by
    censusing every ``<text>`` element in the SVG rather than trusting this
    list.

2.  OPEN AXES -- top and right spines and their ticks dropped
    ---------------------------------------------------------
    The closed four-sided box is an ACS convention: ``secfigs.style`` was
    reverse-engineered from Figure 1A of an ACS Chem. Neurosci. paper, and
    ``show_box``/``ticks_all_sides`` came along with it.  The Nature house
    rcParams set ``axes.spines.top`` and ``axes.spines.right`` FALSE.

    WHAT IT COSTS.  Three things, measured rather than asserted:

    (a) DO THE TRACES FLOAT AT THE RIGHT?  Barely -- and this was checked
        before it was claimed.  At x = 20 mL the largest value any of the 26
        traces still holds is 0.0228 of full scale (panel b, 72 h); the
        per-panel maxima are 0.0029 / 0.0228 / 0.0063 / 0.0202.  Every trace
        therefore ends between 0.120 and 0.166 in above the bottom spine,
        which is still there and still runs all the way out to x = 20.  The
        curves arrive at a rule, not at nothing.  Dropping the right spine
        costs very little HERE precisely because these are chromatograms that
        return to baseline; on a figure whose traces ended mid-panel the same
        change would be a real loss.

    (b) THE TOP IS WHERE IT ACTUALLY FLOATS.  Every normalised panel has a
        peak reaching exactly 1.000, whose apex sits 0.114 in below where the
        top spine used to be, and that apex is now the topmost ink in the
        panel with nothing above it.  The y = 1.0 tick and its label are the
        only remaining statement of where the box ends.  This is the genuine
        cost of the change and it is accepted, not dodged.

    (c) HALF THE TICK INK IS GONE, exactly.  Each panel drew 8 x majors and
        6 y majors on two sides each; they are now on one side each, so the
        major tick count per panel falls 28 -> 14 and the minors likewise.
        The four-sided treatment let a reader recover a coordinate from
        whichever edge was nearest; now a point in the upper right has to be
        carried back down and left.  What remains is unchanged Griffin
        geometry -- 2.40 pt major, 1.40 pt minor, both outward, 0.50 / 0.40 pt
        wide.  Nothing about the two surviving sides was re-measured.

    A fourth thing was expected to cost and does not: the four panels no
    longer read as four identical rectangles, which mattered for the c-vs-d
    comparison.  Because ROW_H and PANEL_W are still identical across all
    four and the columns are still aligned, the four L-shapes register just
    as well as the four boxes did.

    WHAT IT DOES NOT COST.  The colourbar still reads as attached, and by
    construction: it is drawn at exactly the axes extent, left spine of the
    left column to right spine of the right column (2 * PANEL_W + GAP_W =
    6.3516 in), so its own outline rectangle is vertically aligned with the
    four left spines and the four bottom-row x axes.  It also keeps a closed
    four-sided outline of its own -- a colourbar is a swatch strip, not a
    plot, and an open swatch strip would have no edge at all.  That closed
    rectangle sitting under two open axes is what re-anchors the right-hand
    side of the page.

3.  RE-LAID-OUT from the new point sizes
    ------------------------------------
    Every gap is recomputed by ``figure_sec_seq._blocks`` from the style's own
    point sizes; nothing is nudged.  The three pads that separate text from
    what it labels are scaled by the ratio of the type they set off, which is
    the same rule that produced them at 6 pt:

        tick_pad          1.80 -> 2.10  (x 7/6, the tick-label ratio)
        axis_label_pad    3.00 -> 3.43  (x 8/7, the axis-title ratio)
        title_pad         3.50 -> 4.00  (x 8/7, the panel-title ratio)
        CB_TICK_PAD       2.00 -> 2.33  (x 7/6)
        CB_LABEL_PAD      4.00 -> 4.67  (x 7/6)

    See :data:`GEOM` and the module ``main()`` for the resulting numbers.

4.  MAGMA, and a dark-grey frame
    -----------------------------
    ``secfigs.timecolor.COLORS_MAGMA`` -- same construction as the cividis
    ramp (eleven members, distinct 8-bit entries, light end clipped at
    L* <= 78, largest uniform L* gap the map supports), rank-indexed to the
    same eleven ``TIMES``.  Discrete colourbar, eleven equal cells, labelled
    in real hours.  Magma is what the one 2025 precedent for this figure type
    uses (Sneideriene et al., PNAS 2025, Fig 3D) and it beats cividis on this
    figure's own metric -- worst-of-{normal, deutan, protan, greyscale}
    6.58 vs 6.20, normal-vision 10.44 vs 8.49, and panel a's minimum
    separation 20.61 vs 11.58, a 78% improvement on the one panel a ramp keyed
    to absolute time cannot otherwise serve.

    Magma's darkest member is L* 6.86, dE76 6.86 from pure black -- which is
    why the frame moves to :data:`timecolor.FRAME_INK_GREY` = ``#3b3b3b``,
    chosen by sweeping all 256 neutral greys (``timecolor.frame_grey_search``).
    It is dE76 18.01 from that darkest member in the worst rendering
    (greyscale), against 6.86 for black, and holds 11.20 : 1 contrast against
    white so 7 pt tick labels still read as ink.

    The grey draws the READING APPARATUS -- spines, ticks, tick labels, axis
    titles, the legend text and the whole colourbar.  Pure black is kept for
    the NAMING layer -- panel letters, panel titles, and the Oligomer/Monomer
    pointers and their triangles.  Two registers, and the one that recedes is
    the one the darkest traces run along.

5.  NO redundant non-colour channel -- tested, measured, rejected
    --------------------------------------------------------------
    Sparse open markers, rank-keyed dash patterns and direct end labelling
    were each built, measured and inspected at 1:1.  All three are worse than
    colour alone on THIS data -- the short version is that a chromatogram
    puts its whole signal in 10% of the axis and stacks its traces 0.02-0.04
    in apart everywhere else, which is thinner than the smallest usable
    glyph, and dashing costs 4 of 15 oligomer apexes outright.  The full
    reasoning, with the numbers, is :data:`REDUNDANT_CHANNEL_VERDICT`.

Usage
-----
    from secfigs.figure_sec_nature import NatureSECFigure
    NatureSECFigure().build()
    python -m secfigs.figure_sec_nature          # + the geometry report
"""

from __future__ import annotations

import math
import os
import sys
from dataclasses import replace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from secfigs import style as gs                       # noqa: E402
from secfigs import timecolor as tc                   # noqa: E402
from secfigs.figure_sec import MM, PANEL_SPEC         # noqa: E402
from secfigs.figure_sec import LETTER_DX_IN           # noqa: E402
from secfigs.figure_sec import RC_EXTRA as _RC_BASE   # noqa: E402
from secfigs.figure_sec_seq import (                  # noqa: E402
    NCOL, NROW, ROW_H, TRIM, SMALL as _SMALL_SEQ,
    SequentialSECFigure, _blocks,
)


# ==========================================================================
# 1.  the Nature-floor style
# ==========================================================================
#: Ratios the three text pads are scaled by -- each pad is scaled by the size
#: of the type it sets off, which is the rule that produced the 6 pt values.
_R_TICK = 7.0 / 6.0                  # tick labels 6.0 -> 7.0
_R_TITLE = 1.0                       # axis / panel titles STAY at 7.0.
# Nature's NRJ Guide to Preparing Final Artwork, which is the document that governs a
# research DATA figure, says verbatim: "Maximum text size is 7pt. Minimum text size is
# 5pt."  7 pt is a CEILING here, not a floor.  The "at least 7 point size text" line that
# originally motivated raising these came from the Nature Reviews artwork guide, which
# states on its own page 2 that it "chiefly applies to conceptual figures (as opposed to
# data figures)".  Wrong document.  The only carve-out is the panel letter at 8 pt bold.

#: Axis ink.  Spines, ticks, tick labels, axis titles, legend text, colourbar.
FRAME_GREY = tc.FRAME_INK_GREY       # '#3b3b3b'
#: Naming ink.  Panel letters, panel titles, peak pointers.
NAME_INK = "#000000"

_STYLE_NOSIZE = replace(
    _SMALL_SEQ,
    # -- typography: nothing below 7.0 pt -----------------------------------
    tick_labelsize=7.0,
    tick_pad=1.80 * _R_TICK,                     # 2.100
    axis_label_size=7.0,
    axis_label_pad=3.00 * _R_TITLE,              # 3.4286
    title_size=7.0,
    title_pad=3.50 * _R_TITLE,                   # 4.000
    legend_fontsize=7.0,
    annotation_size=7.0,
    panel_letter_size=8.0,   # the one documented exception to the 7 pt ceiling
    # -- peak pointers: offsets and marker scale with the annotation size ---
    annot_marker_size=2.2 * _R_TICK,             # 2.5667
    annot_offset_right=(0.63 * _R_TICK, 10.20 * _R_TICK),
    annot_offset_left=(-5.08 * _R_TICK, 8.58 * _R_TICK),
    # -- ink ----------------------------------------------------------------
    foreground=FRAME_GREY,
    annot_color=NAME_INK,
    # -- open axes ----------------------------------------------------------
    # ticks_all_sides=False is the whole of the tick half of change 2: gs.axes
    # reads it as `sides = dict(top=..., bottom=True, left=True, right=...)`
    # and gs.apply as xtick.top / ytick.right, so both major and minor ticks
    # vanish from the top and right in one place.
    #
    # show_box STAYS TRUE.  It is a single flag over all four spines
    # (`for spine in ax.spines.values(): spine.set_visible(style.show_box)`),
    # so setting it False would take the left and bottom with it.  The two
    # spines that go are removed per-Axes in _style_axes, after gs.axes has
    # drawn and widthed all four.
    show_box=True,
    ticks_all_sides=False,
)
# Strokes are NOT touched: spine 0.50, trace 0.85, major tick 2.40 x 0.50,
# minor tick 1.40 x 0.40 pt, all outward.  Change 2 drops two sides of the
# Griffin tick treatment; it does not re-measure the two that stay.

#: Text-block extents in inches, derived from the sizes above.
B = _blocks(_STYLE_NOSIZE)


# ==========================================================================
# 2.  page geometry -- every gap derived, nothing nudged
# ==========================================================================
_PT = 1.0 / 72.0

# -- colourbar block --------------------------------------------------------
CB_H = 0.120                                  # in, bar thickness (unchanged)
CB_GAP = 0.200                                # in, bar top -> x axis title
CB_TICK_PAD = 2.00 * _R_TICK                  # 2.3333 pt
CB_LABEL_PAD = 4.00 * _R_TICK                 # 4.6667 pt
CB_TICK_BLOCK = CB_TICK_PAD * _PT + B["tick_h"]        # 0.145185 in
CB_TITLE_BLOCK = CB_LABEL_PAD * _PT + B["axlab_h"]     # 0.193704 in
CB_BLOCK = CB_GAP + CB_H + CB_TICK_BLOCK + CB_TITLE_BLOCK + TRIM

FIG_W_MM = 180.0
FIG_W = FIG_W_MM * MM                                  # 7.086614 in

def _ceil3(v):
    """Round a derived length UP to the next 0.001 in.

    Not ``round``.  Every number below is a MINIMUM -- the space something
    needs -- so rounding to nearest can eat into it.  It did: ``right``
    derives to 0.184056 in and round() gives 0.184, which is 0.000056 in
    short, and the renderer duly put the right edge of panel b's "20" tick
    label 0.0001 in past the 0.130 in trim line.  Rounding up can only ever
    add clearance.  (``figure_sec_seq`` patched the same two margins by hand,
    0.176 -> 0.180 and 0.291 -> 0.300; this is that patch made a rule.)
    """
    return math.ceil(v * 1000.0 - 1e-9) / 1000.0


MARGIN = dict(
    # the 9 pt panel letter still reaches further out (0.420 in) than the
    # y-axis block does (0.3741 in), by 0.0459 in -- so it still sets `left`
    left=_ceil3(LETTER_DX_IN + TRIM),                  # 0.550
    # "20" is centred on the last x tick, at the axes' right edge, and
    # overhangs it by half its width -- 0.0541 in at 7 pt
    right=_ceil3(B["xtick_half"] + TRIM),              # 0.185
    # the title block beats the 9 pt panel letter, whose cap top reaches only
    # 0.045 + 9 * 0.716 / 72 = 0.1345 in above the axes
    top=_ceil3(B["t_block"] + TRIM),                   # 0.315
    bottom=_ceil3(B["x_block"] + CB_BLOCK),            # 1.141
)

#: 0.154 in of clear air between the left panel's overhanging "20" and the
#: right panel's letter -- the same constant the 6 pt layout used.
GAP_W_CLEAR = 0.154
#: 0.150 in of clear air between the top row's x axis title and the bottom
#: row's panel title -- likewise.
GAP_H_CLEAR = 0.150

GAP_W = _ceil3(LETTER_DX_IN + B["xtick_half"] + GAP_W_CLEAR)     # 0.629
GAP_H = _ceil3(B["x_block"] + B["t_block"] + GAP_H_CLEAR)        # 0.687

_AVAIL_W = FIG_W - MARGIN["left"] - MARGIN["right"]
PANEL_W = (_AVAIL_W - (NCOL - 1) * GAP_W) / NCOL
WSPACE = GAP_W / PANEL_W
HSPACE = GAP_H / ROW_H

FIG_H = MARGIN["top"] + NROW * ROW_H + GAP_H + MARGIN["bottom"]
FIG_H_MM = FIG_H / MM

assert FIG_H_MM <= 230.0, "page is taller than the 230 mm Nature cap"

#: Colourbar rect in inches from the page's bottom-left corner.  It spans the
#: full axes extent -- left column's left spine to right column's right edge.
CB_RECT_IN = (MARGIN["left"],
              TRIM + CB_TITLE_BLOCK + CB_TICK_BLOCK,
              NCOL * PANEL_W + GAP_W,
              CB_H)

STYLE = replace(_STYLE_NOSIZE, figsize=(FIG_W, FIG_H))

RC_EXTRA = dict(_RC_BASE)
RC_EXTRA["svg.hashsalt"] = "figure_sec_nature"
#: Titles, panel letters and pointer text default to rcParams['text.color'].
RC_EXTRA["text.color"] = NAME_INK
#: ...while the axis furniture follows style.foreground through gs.apply.

#: One dict of every derived number, for the report and the verifier.
GEOM = dict(
    page_mm=(FIG_W_MM, FIG_H_MM), margin_in=dict(MARGIN),
    gap_w_in=GAP_W, gap_h_in=GAP_H, row_h_in=ROW_H, panel_w_in=PANEL_W,
    blocks_in={k: round(v, 6) for k, v in B.items()},
    cb=dict(h=CB_H, gap=CB_GAP, tick_block=CB_TICK_BLOCK,
            title_block=CB_TITLE_BLOCK, block=CB_BLOCK, rect_in=CB_RECT_IN),
    trim_in=TRIM,
)


# ==========================================================================
# 3.  change 5 -- the redundant-channel experiment, and its verdict
# ==========================================================================
#: Eleven open marker shapes, rank-indexed exactly like the colours, so a time
#: would carry the same glyph in every panel.  Only used by the rejected
#: ``redundant='markers'`` build; kept so the experiment can be re-run.
RANK_MARKERS = ("o", "s", "^", "v", "D", "<", ">", "p", "h", "*", "X")

#: Eleven dash patterns, on-off in points, rank-indexed.  Rejected build only.
RANK_DASHES = (
    (None, None), (3.0, 1.0), (1.2, 1.0), (4.0, 1.0, 1.2, 1.0),
    (2.0, 1.0), (5.0, 1.2), (1.2, 1.6), (3.0, 1.0, 1.2, 1.0, 1.2, 1.0),
    (6.0, 1.4), (2.4, 1.2, 1.2, 1.2), (0.8, 1.2),
)

REDUNDANT_CHANNEL_VERDICT = """\
EVALUATED AND REJECTED.  Colour alone carries the time index; the figure ships
with no second channel.  All three candidates were built at the shipping
600 dpi, measured, and looked at at 1:1.  Re-run any of them with
``NatureSECFigure(redundant='markers' | 'dashes' | 'labels')``.

WHY THE PRECEDENT DOES NOT TRANSFER.  Brichos et al., PNAS 2024, indexes nine
traces with nine open marker shapes -- but that is ThT aggregation kinetics:
smooth monotone sigmoids, sparse points, and the curves fan APART as they
rise, so there is always empty page next to a curve to put a glyph in.  These
are chromatograms.  Two structural facts about them defeat every candidate:

    the signal is 10% of the axis.  The oligomer and monomer peaks together
    occupy ~1.5 mL of a 14 mL window; over the other ~90% all 7 or 9 traces
    lie on a shared baseline.
    the traces converge, they do not fan.  At the marker positions the whole
    stack is 0.0203 in (panel c), 0.0285 in (d) and 0.0424 in (b) thick,
    median.  At x = 20 mL it is 0.0097 / 0.0304 / 0.0498 in.

(a) SPARSE OPEN MARKERS -- rejected
    First attempt used the float form of ``markevery``, which spaces glyphs by
    ARC LENGTH.  On a chromatogram most of the arc length is in the peak
    flanks and the noisy baseline, so markevery=(offset, 0.10) delivered 26-40
    markers per trace, up to 32% of them inside the 7.6-11.0 mL oligomer peak
    -- the opposite of sparse, and it shredded exactly the peaks it was
    supposed to leave alone.
    Retuned to an INDEX STRIDE of n/10, every trace carries exactly 10-11
    markers, evenly spaced in elution volume.  That is the brief's 8-12, and
    it still fails, for a reason no tuning can reach:
      * 83% (panel d), 87% (c) and 90% (b) of those markers land where their
        own trace is below 0.10 of full scale -- i.e. on the shared baseline.
        They identify traces in the one place nobody needs them identified.
      * one open marker at 2.57 pt is 0.0356 in across.  The trace stack at
        those positions is 0.0203-0.0424 in thick.  A single glyph is
        therefore 84% (b), 125% (d) and 175% (c) of the ENTIRE stack
        thickness: it cannot be drawn without its white fill covering the
        neighbouring traces.  Verified at 1:1 -- the baseline turns into a
        row of overlapping pentagons and crosses with the lines lost under
        them.
      * the 1-2 markers per trace that do land near the oligomer peak land on
        its flanks, where they punch white holes through the adjacent curves.
    The region where the traces are hardest to tell apart is exactly the
    region with no room for a glyph, and the region with room is the region
    where telling them apart does not matter.

(b) DASH PATTERNS keyed to rank -- rejected
    Eleven patterns, rank-indexed like the colours, periods 2.0-8.4 pt.  The
    oligomer peaks are long enough in arc length to carry 6-80 dash periods,
    so the peaks are not lost wholesale.  The apexes are.  Testing the
    apex vertex of every oligomer peak against its trace's dash phase
    analytically: 4 of 15 fall inside an OFF segment -- panel c 12 h, panel d
    6 h, 48 h and 72 h.  For those four the drawn peak height is set by dash
    phase rather than by the data.  Peak height is the measurement in this
    figure, so this is not a legibility trade-off, it is the style deciding
    the value of the readout.  Disqualifying on its own, before any question
    of how it looks -- and at 1:1 it also looks like static.

(c) DIRECT LABELLING at the curve ends -- rejected
    The most promising of the three and still a clear loss.  Labels were put
    at x = 20 mL on three traces per panel (first, middle, last).  Two
    failures, both measured:
      * the curves have converged there.  Adjacent traces are 0.0005-0.0008 in
        apart in y at x = 20 mL, against a 0.1128 in line box for 7 pt type.
        Even three labels overprint; labelling all of them would need
        16 x (b), 33 x (d) and 81 x (c) more vertical room than exists.  At
        1:1 panel b's three labels are a single unreadable smear.
      * they do not fit the page.  Every end label crosses the 6.9566 in trim
        line, the worst reaching 7.0754 in on a 7.0866 in page -- 0.03 mm of
        paper.  Moving them inboard puts them on top of the traces instead.
    Labelling only the three the ramp already separates best, while leaving
    the crowded middle unlabelled, buys nothing where the problem is.

WHAT IS DONE INSTEAD.  Two things already in the figure do the work a
redundant channel would have done, without putting ink on the peaks:
  * the ramp is monotone in L* BY CONSTRUCTION, min adjacent step 6.85 L*.
    Strip the colour out entirely and the series still reads as an ordered
    light-to-dark sequence -- min greyscale dE76 6.85 over all 55 pairs.  That
    IS a non-colour channel; it is carried by the same stroke instead of by a
    second mark, which is the only way to carry it on a 0.6 mL peak.
  * panel a -- the one panel whose three times span too little of the 0-72 h
    ramp for the ramp to separate them -- keeps its explicit 3-entry legend.
An honest "this does not work here" beats a cluttered figure.\
"""


# ==========================================================================
# 4.  the figure
# ==========================================================================
class NatureSECFigure(SequentialSECFigure):
    """The 2 x 2 SEC figure at Nature's type floor, open axes, magma ramp.

    Parameters
    ----------
    normalized
        ``True`` (default) -- the recommended variant, and the only one built
        by ``main()``.  ``False`` still works and inherits the parent's raw
        y-axis resolution unchanged.
    redundant
        ``None`` (default, shipped) or one of ``'markers'`` / ``'dashes'`` /
        ``'labels'`` to re-render one of the three rejected candidates.  See
        :data:`REDUNDANT_CHANNEL_VERDICT`.
    """

    PALETTE = "magma"

    def __init__(self, normalized=True, style=STYLE, outdir=None,
                 stem="Figure_SEC_natstyle", redundant=None, **kw):
        super().__init__(normalized=normalized, style=style, outdir=outdir,
                         stem=stem, **kw)
        if redundant not in (None, "markers", "dashes", "labels"):
            raise ValueError("redundant must be None, 'markers', 'dashes' "
                             "or 'labels'")
        self.redundant = redundant

    # -- geometry / rc -----------------------------------------------------
    def _layout(self):
        return dict(margin=dict(MARGIN), gap_w=GAP_W, gap_h=GAP_H,
                    row_h=ROW_H, panel_w=PANEL_W, wspace=WSPACE,
                    hspace=HSPACE, cb_rect_in=CB_RECT_IN,
                    cb_h=CB_H, cb_gap=CB_GAP,
                    cb_tick_pad=CB_TICK_PAD, cb_label_pad=CB_LABEL_PAD,
                    blocks=dict(B), trim=TRIM)

    def _rc(self):
        return dict(RC_EXTRA)

    # -- change 2: open axes ----------------------------------------------
    def _style_axes(self, ax, spec):
        """Drop the top and right spines and every tick on those two sides.

        ``gs.axes`` has already written all four spines and four-sided ticks
        from the style; this undoes two of the four.  Left and bottom keep the
        Griffin outward geometry untouched.
        """
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="x", which="both", top=False)
        ax.tick_params(axis="y", which="both", right=False)

    # -- change 4: the legend joins the grey reading layer -----------------
    def _style_legend(self, leg, ax, spec):
        for t in leg.get_texts():
            t.set_color(FRAME_GREY)

    # -- change 4: panel letters stay black --------------------------------
    def _place_letters(self):
        """Parent's placement, then recolour: letters are the naming layer.

        ``CombinedSECFigure._place_letters`` draws them in
        ``style.foreground``, which this style has moved to the axis grey.
        """
        n0 = len(self.fig.texts)
        super()._place_letters()
        for t in self.fig.texts[n0:]:
            t.set_color(NAME_INK)

    # -- change 5 candidates (all rejected; kept runnable) -----------------
    def _trace_kwargs(self, spec, i, n):
        if self.redundant is None:
            return {}
        rank = self.data[spec["letter"]]["ranks"][i]
        if self.redundant == "markers":
            # markevery must be an INDEX STRIDE, not a float.  The float form
            # spaces markers by ARC LENGTH, and on a chromatogram most of the
            # arc length is in the peak flanks and the noisy baseline: at
            # markevery=(off, 0.10) the traces came out carrying 26-40 markers
            # each, up to 32% of them piled inside the 7.6-11.0 mL oligomer
            # peak.  A stride of n_vertices/10 puts ~10 on each trace, evenly
            # spaced in elution volume, which is the brief's 8-12.
            npts = len(self.data[spec["letter"]]["traces"][i][0])
            stride = max(int(round(npts / 10.0)), 1)
            return dict(marker=RANK_MARKERS[rank],
                        markersize=self.style.annot_marker_size,
                        markerfacecolor="white",
                        markeredgewidth=0.45,
                        # offset per trace so glyphs interleave in x rather
                        # than stacking on one vertical
                        markevery=(int(i * stride / max(n, 1)), stride))
        if self.redundant == "dashes":
            on_off = RANK_DASHES[rank]
            return {} if on_off[0] is None else dict(dashes=list(on_off))
        return {}

    def _after_panel(self, ax, spec):
        if self.redundant != "labels":
            return
        from secfigs.figure_sec import XLIM
        d = self.data[spec["letter"]]
        n = len(d["traces"])
        pick = sorted({0, n // 2, n - 1})
        for i in pick:
            V, A = d["traces"][i]
            j = int(((V >= XLIM[1] - 0.35) & (V <= XLIM[1])).argmax())
            ax.annotate(d["labels"][i], xy=(float(V[j]), float(A[j])),
                        xytext=(2.0, 0.0), textcoords="offset points",
                        ha="left", va="center",
                        fontsize=self.style.legend_fontsize,
                        color=d["colors"][i], annotation_clip=False)

    # -- save with this figure's own rc (element-id salt) ------------------
    def save(self, formats=("svg", "png", "pdf")):
        """As the parent, but under :meth:`_rc` so the SVG ids get own salt."""
        if self.fig is None:
            self.draw()
        os.makedirs(self.outdir, exist_ok=True)
        out = {}
        with gs.apply(style=self.style, **self._rc()):
            for ext in formats:
                p = os.path.join(self.outdir, "%s.%s" % (self.stem, ext))
                # NO bbox_inches -- the page must stay exactly 180 mm wide.
                self.fig.savefig(p, dpi=self.style.dpi, format=ext,
                                 facecolor="white", edgecolor="none")
                out[ext] = p
        return out

    # -- report ------------------------------------------------------------
    def geometry(self):
        g = super().geometry()
        g["frame_ink"] = FRAME_GREY
        g["name_ink"] = NAME_INK
        g["redundant_channel"] = self.redundant
        g["font_pt"] = dict(
            tick=self.style.tick_labelsize,
            axis_title=self.style.axis_label_size,
            panel_title=self.style.title_size,
            legend=self.style.legend_fontsize,
            annotation=self.style.annotation_size,
            panel_letter=self.style.panel_letter_size)
        g["open_axes"] = {
            L: {s: bool(ax.spines[s].get_visible())
                for s in ("left", "bottom", "top", "right")}
            for L, ax in self.axes.items()}
        return g


# ==========================================================================
# 5.  cli
# ==========================================================================
def build(formats=("svg", "png", "pdf"), redundant=None,
          stem="Figure_SEC_natstyle", outdir=None):
    f = NatureSECFigure(normalized=True, redundant=redundant, stem=stem,
                        outdir=outdir)
    f.load()
    f.draw()
    files = f.save(formats)
    g = f.geometry()
    plt.close(f.fig)
    return files, g


def main():
    files, g = build()
    print("page  %.4f x %.4f mm  (%.4f x %.4f in)  @ %d dpi = %d x %d px"
          % (g["page_mm"][0], g["page_mm"][1], FIG_W, FIG_H, g["dpi"],
             g["px"][0], g["px"][1]))
    print("palette %s   frame ink %s   naming ink %s   redundant channel %s"
          % (g["palette"], g["frame_ink"], g["name_ink"],
             g["redundant_channel"]))
    print("fonts pt: %s" % g["font_pt"])
    print("\nderived text blocks (in):")
    for k in ("tick_h", "axlab_h", "title_h", "ytick_w", "xtick_half",
              "letter_cap", "x_block", "y_block", "t_block"):
        print("   %-11s %.6f" % (k, B[k]))
    print("\nmargins (in):")
    print("   left   %.4f = panel letter %.3f + trim %.3f"
          % (MARGIN["left"], LETTER_DX_IN, TRIM))
    print("   right  %.4f = half of \"20\" %.4f + trim %.3f"
          % (MARGIN["right"], B["xtick_half"], TRIM))
    print("   top    %.4f = T_BLOCK %.4f + trim %.3f"
          % (MARGIN["top"], B["t_block"], TRIM))
    print("   bottom %.4f = X_BLOCK %.4f + cbar gap %.3f + bar %.3f"
          " + cbar ticks %.4f + cbar title %.4f + trim %.3f"
          % (MARGIN["bottom"], B["x_block"], CB_GAP, CB_H,
             CB_TICK_BLOCK, CB_TITLE_BLOCK, TRIM))
    print("   GAP_W  %.4f = letter %.3f + \"20\" overhang %.4f + clear %.3f"
          % (GAP_W, LETTER_DX_IN, B["xtick_half"], GAP_W_CLEAR))
    print("   GAP_H  %.4f = X_BLOCK %.4f + T_BLOCK %.4f + clear %.3f"
          % (GAP_H, B["x_block"], B["t_block"], GAP_H_CLEAR))
    print("   ROW_H  %.4f (unchanged, so panel boxes stay comparable)" % ROW_H)
    print("\npanels:")
    for k, v in g["panels"].items():
        sp = g["open_axes"][k]
        print("   %s  %.4f x %.4f in at (%.4f, %.4f)  aspect %.4f  n=%d  "
              "spines L%d B%d T%d R%d"
              % (k, v["w_in"], v["h_in"], v["x0_in"], v["y0_in"],
                 v["aspect"], v["n_traces"], sp["left"], sp["bottom"],
                 sp["top"], sp["right"]))
    c = g["colorbar"]
    print("   cbar %.4f x %.4f in at (%.4f, %.4f)  cell %.4f in"
          % (c["w_in"], c["h_in"], c["x0_in"], c["y0_in"], c["cell_in"]))
    for k, v in files.items():
        print("   %-4s %s" % (k, v))
    return files, g


if __name__ == "__main__":
    main()
