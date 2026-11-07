"""
secfigs.style -- the "Griffin Figure 1A" look, as a reusable matplotlib style.

Every number in :class:`GriffinStyle` was MEASURED off the rasterised panel A of
Figure 1 of

    Griffin, Dinakarapandian, Paravastu & Wood,
    ACS Chem. Neurosci. 2026 (submission 2, 05/19/26)

i.e. the file ``griffin_fig1.png`` extracted from the submitted .docx
(3047 x 2404 px RGBA).  Nothing here is eyeballed; see the module-level
``MEASUREMENTS`` dict for the raw pixel evidence behind each value.

Reference geometry
------------------
Panel A's axes box occupies columns 271..1466 and rows 145..1085 of the PNG,
i.e. 1195 x 940 px, with ``xlim=(6, 18)`` and ``ylim=(-0.05, 1.05)``.
Tick-label digits are 28 px tall; Arial digit height is 0.716 em, so
1 em of tick text = 39.1 px.  Reading the tick labels as 22 pt fixes the
raster scale at

    k = 39.1 px / 22 pt = 1.777 px per point

and every other measured pixel length below is divided by k to give points.
At that scale the axes box is 672 x 529 pt = 9.33 x 7.34 in, which is what a
``figsize=(12, 9.53)`` figure gives you after ``tight_layout`` -- so the point
values in this module are calibrated for ``figsize=(12, 9.53)``,
which reproduces the reference axes-box aspect ratio of 1.270 to 0.3%.  Scale the whole
dataclass with :meth:`GriffinStyle.scaled` if you need a different canvas.

Provenance of the drawing code
------------------------------
The four-sided outward tick treatment comes from ``_configure_axes``, which is
BYTE-IDENTICAL in the lab's ``sec1.py`` and ``sec3.py``.  The only difference
between those two files that touches panel A is that ``sec1.py`` force-bolds the
tick labels and the legend text; the published Figure 1A has REGULAR-weight tick
labels and legend text (measured stem/em 0.087-0.096 for ticks and 0.072 for the
legend, vs 0.145 for the genuinely bold title).  So panel A follows the
``sec3.py`` variant, plus a bold title, plus hand-placed peak annotations that
exist in neither library.  This module reproduces the published panel, not
``sec1.py``'s defaults.

Usage
-----
    from secfigs import style as gs

    with gs.apply():
        fig, ax = plt.subplots(figsize=gs.STYLE.figsize, dpi=gs.STYLE.dpi)
        for (x, y), c, lab in zip(traces, gs.STYLE.palette, labels):
            ax.plot(x, y, color=c, lw=gs.STYLE.trace_linewidth, label=lab)
        gs.axes(ax, normalized=True, x_major=2.0)
        gs.set_labels(ax, "Elution Volume (ml)", "Normalized Absorbance")
        gs.set_title(ax, "Normalized SEC Profile ...")
        gs.legend(ax)
        gs.annotate_peak(ax, 9.33, 0.45, "Oligomer", side="right")
        gs.annotate_peak(ax, 13.62, 0.77, "Monomer", side="left")
        gs.panel_letter(ax, "A")
        fig.savefig(out, dpi=gs.STYLE.dpi, bbox_inches="tight")
"""

from __future__ import annotations

import warnings
from contextlib import contextmanager
from dataclasses import dataclass, replace
from typing import Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MultipleLocator


# --------------------------------------------------------------------------
# Raw evidence.  Pixel measurements on griffin_fig1.png, panel A.
# k = 1.777 px/pt (see module docstring).
# --------------------------------------------------------------------------
PX_PER_PT = 1.777

MEASUREMENTS = {
    "axes_box_px": {"left": 271, "right": 1466, "top": 145, "bottom": 1085,
                    "width": 1195, "height": 940},
    "axes_limits": {"xlim": (6.0, 18.0), "ylim": (-0.05, 1.05)},
    # integrated ink (alpha-weighted) across a stroke, in px
    "spine_ink_px": {"top": 1.156, "bottom": 1.176, "left": 1.143, "right": 1.147},
    "tick_width_ink_px": {"major": 1.070, "minor": 0.797},
    "tick_length_px": {"major": 13.3, "minor": 6.7},   # both point OUTWARD
    "tick_sides": "all four; direction='out' verified above top spine, "
                  "right of right spine, left of left spine, below bottom spine",
    "trace_ink_px": {"WT": 4.74, "N27D": 4.80, "Q15E": 4.83},
    "cap_height_px": {"tick_labels": 28, "xlabel": 31, "ylabel": 31,
                      "title": 31, "legend": 30, "annotation": 31,
                      "panel_letter": 69},
    # stem width / em -- Arial regular ~0.09, Arial bold ~0.14
    "stem_over_em": {"tick_labels": 0.090, "xlabel": 0.083, "legend": 0.072,
                     "annotation": 0.082, "title": 0.145, "panel_letter": 0.146},
    "tick_label_pad_px": 6,        # end of major tick -> top of tick label
    "axis_label_gap_px": 32,       # tick-label bbox -> xlabel bbox
    "title_gap_px": 50,            # title bbox bottom -> top spine
    "legend_swatch_px": {"length": 71, "gap_to_text": 24, "row_pitch": 45.5},
    "legend_anchor_px": {"content_top": 165, "text_right": 1469},  # vs axes 145 / 1466
    "annotation_triangle_px": {"oligomer": {"rows": (653, 666), "cols": (597, 608)},
                               "monomer":  {"rows": (382, 394), "cols": (1025, 1035)}},
    "annotation_offsets_px": {"oligomer": {"dx": -4.5, "dy_above": 72.5},
                              "monomer":  {"dx": +36.0, "dy_above": 61.0}},
    "swatch_rgb": {"WT": (251, 106, 74), "N27D": (102, 194, 164),
                   "Q15E": (254, 227, 145)},
    "x_ticks": {"major": 2.0, "minor": 0.5},
    "y_ticks": {"major": 0.2, "minor": 0.1},
    "legend_row_pitch_per1000_boxW": 38.11,
    "legend_content_top_per1000_boxH": 39.36,
}

# The point values in GriffinStyle were seeded from PX_PER_PT and then closed
# by rendering figures/style_probe.png and re-measuring it with the SAME code
# used on the reference.  Final agreement, expressed per 1000 of the axes-box
# width so it is scale invariant (probe / reference):
CALIBRATION = {
    "axes_box_aspect":      (1.270, 1.266, 0.997),
    "spine_ink":            (0.990, 0.998, 1.008),
    "major_tick_length":    (10.888, 10.753, 0.988),
    "minor_tick_length":    (5.025, 5.018, 0.999),
    "major_tick_width":     (0.896, 0.898, 1.002),
    "minor_tick_width":     (0.667, 0.666, 0.998),
    "trace_linewidth":      (3.969, 4.012, 1.011),
    "legend_swatch_length": (57.789, 57.706, 0.999),
    "legend_row_pitch":     (38.11, 37.63, 0.987),
    "legend_content_top":   (39.36, 39.95, 1.015),
    "cap_tick_label":       (23.451, 23.656, 1.009),
    "cap_legend":           (24.288, 23.656, 0.974),
    "cap_title":            (25.96, 25.81, 0.994),
    "cap_panel_letter":     (57.789, 57.706, 0.999),
}


# ---------------------------------------------------------------------------
# The single trace palette for every SEC panel, standalone or combined.
#
# Griffin's Figure 1A uses three ColorBrewer families (Reds / BuGn / YlOrBr) and
# needs only three colours.  These panels need up to nine, and the obvious
# extension -- three shades per family at adjacent ColorBrewer stops -- put
# #8FD4BC next to #66C2A4 at dE 11.0, which is not separable as a 0.85 pt line.
# Stops were re-chosen by maximising the MINIMUM pairwise CIE76 dE over all 36
# pairs, subject to every member being darker than L* = 78 so it reads on white:
# min dE 11.0 -> 21.7.  Order within each family is light -> dark, so a time
# series still reads as a progression.
#
# These are defined ABOVE GriffinStyle because they are the DEFAULTS of its
# ``palette`` / ``palette_extended`` fields.  They used to sit at the bottom of
# the module as independent constants, and ``palette`` carried a third, different
# set of hexes -- so the dataclass and the ramp disagreed about what the house
# colours were, and ``ramp()`` read the constants while ``apply()`` read the
# dataclass.  One source now.
RAMP9 = ("#FC9272", "#FB6A4A", "#CB181D",     # coral  light -> dark
         "#66C2A4", "#238B45", "#00441B",     # teal   light -> dark
         "#FE9929", "#993404", "#662506")     # gold   light -> dark
RAMP3 = ("#FB6A4A", "#66C2A4", "#FE9929")     # Griffin's own three; his pale
                                              # #FEE391 replaced -- at L* 90 it
                                              # is invisible as a thin stroke.

#: The three swatches as MEASURED off the legend rows of griffin_fig1.png
#: (see ``MEASUREMENTS['swatch_rgb']``).  Kept for provenance only: the house
#: palette is :data:`RAMP3`, which replaces the pale #FEE391.  This tuple is
#: what ``GriffinStyle.palette`` used to default to, and the mismatch with
#: RAMP3 is the disagreement that has now been resolved in RAMP3's favour.
GRIFFIN_SWATCHES = ("#FB6A4A", "#66C2A4", "#FEE391")


# --------------------------------------------------------------------------
# The style
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class GriffinStyle:
    """Every measured parameter of Griffin Figure 1A, in points."""

    # -- canvas ------------------------------------------------------------
    # 9.33 x 7.34 in axes box measured; (12, 9.53) + bbox_inches=tight gives
    # 9.30 x 7.34 in -> aspect 1.266 vs the reference 1.270.
    figsize: tuple = (12.0, 9.53)
    dpi: int = 300

    # -- typography --------------------------------------------------------
    font_family: str = "Arial"
    fallback_family: str = "DejaVu Sans"

    tick_labelsize: float = 21.5       # 28 px cap / 0.716 / 1.777 (23.5 per 1000 boxW)
    tick_label_weight: str = "normal"  # stem/em 0.090 -> Arial Regular
    tick_pad: float = 3.5              # 6 px / 1.777 = 3.4 (mpl default)

    axis_label_size: float = 24.0      # 31 px cap / 0.716 / 1.777 = 24.4
    axis_label_weight: str = "normal"  # stem/em 0.083
    axis_label_pad: float = 18.0       # 32 px / 1.777

    title_size: float = 24.0           # 31 px cap, same as axis labels
    title_weight: str = "bold"         # stem/em 0.145 -> Arial Bold
    title_pad: float = 28.0            # 50 px / 1.777
    title_loc: str = "center"

    legend_fontsize: float = 22.4      # 30 px cap = 24.3 per 1000 of box width
    legend_weight: str = "normal"      # stem/em 0.072

    annotation_size: float = 24.0      # 31 px cap
    annotation_weight: str = "normal"  # stem/em 0.082

    panel_letter_size: float = 54.0    # 69 px cap / 0.716 / 1.777 = 54.2
    panel_letter_weight: str = "bold"  # stem/em 0.146
    # axes-fraction position of the letter's top-left corner, from the
    # composite: (28-271)/1195 = -0.203 ; (1085-21)/940 = 1.132
    panel_letter_xy: tuple = (-0.20, 1.13)

    # -- strokes -----------------------------------------------------------
    spine_linewidth: float = 0.60      # 1.155 px ink = 0.97 per 1000 of box width
    trace_linewidth: float = 2.65      # 4.79 px ink = 3.97 per 1000 of box width

    major_tick_length: float = 7.46    # 13.3 px ink = 11.14 per 1000 of box width
    minor_tick_length: float = 3.76    #  6.7 px ink =  5.61 per 1000 of box width
    major_tick_width: float = 0.60     # 1.070 px = 0.90 per 1000 of box width
    minor_tick_width: float = 0.445    # 0.797 px = 0.67 per 1000 of box width
    tick_direction: str = "out"
    ticks_all_sides: bool = True
    show_box: bool = True              # all four spines visible

    # -- locators ----------------------------------------------------------
    x_major: float | None = 2.0
    x_minor: float | None = None       # None -> AutoMinorLocator() (4 subdiv)
    y_major_normalized: float = 0.2
    y_minor_normalized: float = 0.1
    y_pad_normalized: float = 0.05     # ylim = (-0.05, 1.05)

    # -- colour ------------------------------------------------------------
    # THE trace palette.  ``ramp(n, style)`` reads these two fields and nothing
    # else, so replacing them on the style object is the supported way to
    # recolour every panel at once.
    #
    # Provenance: the swatch rows of the PNG sample to
    #   row 183 cols 1175-1240 -> (251,106,74)   #FB6A4A
    #   row 229                -> (102,194,164)  #66C2A4
    #   row 274                -> (254,227,145)  #FEE391
    # and that measured triple is kept as GRIFFIN_SWATCHES.  The house palette
    # is RAMP3, which replaces the third: at L* 90, #FEE391 is invisible as a
    # thin stroke on white.  These fields used to default to the measured
    # triple (and to a 9-colour extension of it) while ``ramp()`` handed out
    # RAMP3/RAMP9 -- two disagreeing definitions of the same thing.
    palette: tuple = RAMP3
    # the same three hue families at three shades each, for 4-9 trace panels
    palette_extended: tuple = RAMP9
    foreground: str = "#000000"

    # -- legend ------------------------------------------------------------
    legend_loc: str = "upper right"
    legend_frameon: bool = False       # frameless -- no box in the reference
    legend_handlelength: float = 1.72  # 69 px swatch = 57.8 per 1000 of box width
    legend_handletextpad: float = 0.60  # 24 px / (24 pt * 1.777)
    legend_labelspacing: float = 0.20  # 45.5 px row pitch = 38.1 per 1000 of box width
    legend_borderaxespad: float = 0.0  # text right edge overshoots the spine 3 px
    legend_borderpad: float = 0.61     # legend content top 39 per 1000 of box height

    # -- peak annotations ("Oligomer" / "Monomer") -------------------------
    # filled triangle 12 x 11 px -> markersize 6.8 pt; matplotlib '>' / '<'
    annot_marker_size: float = 6.8
    annot_marker_right: str = ">"
    annot_marker_left: str = "<"
    annot_color: str = "#000000"
    # text offset from the marker, in points. Measured:
    #   Oligomer: text left edge 4.5 px right of marker centre, baseline
    #             72.5 px above it  -> (+2.5, +40.8) pt, ha='left'
    #   Monomer:  text right edge 36 px left of marker centre, baseline
    #             61 px above it    -> (-20.3, +34.3) pt, ha='right'
    annot_offset_right: tuple = (2.5, 40.8)
    annot_offset_left: tuple = (-20.3, 34.3)

    def scaled(self, factor: float) -> "GriffinStyle":
        """Uniformly rescale every point-valued field (fonts and strokes).

        ``STYLE.scaled(0.5)`` with ``figsize=(6, 5)`` gives the identical
        printed appearance at half the canvas size.
        """
        pt_fields = (
            "tick_labelsize", "tick_pad", "axis_label_size", "axis_label_pad",
            "title_size", "title_pad", "legend_fontsize", "annotation_size",
            "panel_letter_size", "spine_linewidth", "trace_linewidth",
            "major_tick_length", "minor_tick_length", "major_tick_width",
            "minor_tick_width", "annot_marker_size",
        )
        kw = {f: getattr(self, f) * factor for f in pt_fields}
        kw["figsize"] = (self.figsize[0] * factor, self.figsize[1] * factor)
        kw["annot_offset_right"] = tuple(v * factor for v in self.annot_offset_right)
        kw["annot_offset_left"] = tuple(v * factor for v in self.annot_offset_left)
        return replace(self, **kw)


STYLE = GriffinStyle()


#: Ratio between Arial Bold and Arial Regular stem weight, taken from the two
#: measurements already recorded on GriffinStyle rather than chosen: the reference's
#: bold roles run stem/em 0.145 (title) and 0.146 (panel letter), its regular roles
#: 0.090 (tick labels) and 0.083 (axis labels).  0.1455 / 0.0865 = 1.68.
BOLD_STEM_RATIO = 1.68


def bold_variant(base: GriffinStyle = None, frame: bool = True,
                 ratio: float = BOLD_STEM_RATIO) -> GriffinStyle:
    """The reference style with everything that carries weight set bold.

    THIS IS A DEPARTURE FROM THE REFERENCE, NOT A CORRECTION TO IT.  Griffin Fig. 1
    draws its tick labels and axis labels in Arial REGULAR -- that is what the
    stem/em measurements on GriffinStyle record -- and only its title and panel
    letter in Bold.  Use this when the house preference is bold throughout.

    ``frame=True`` also thickens the box and both tick trains by the same 1.68
    factor, because a bold label inside a hairline frame reads as a mismatch
    rather than as emphasis.  The factor is the text ratio, so the frame gains
    exactly the weight the type does:

        spine        0.600 -> 1.008 pt
        major tick   0.600 -> 1.008 pt
        minor tick   0.445 -> 0.748 pt

    Everything else -- sizes, tick lengths, spacings, colours, figsize -- is left
    exactly as calibrated.
    """
    from dataclasses import replace
    base = STYLE if base is None else base
    kw = dict(tick_label_weight="bold", axis_label_weight="bold",
              legend_weight="bold", annotation_weight="bold")
    if frame:
        kw.update(spine_linewidth=base.spine_linewidth * ratio,
                  major_tick_width=base.major_tick_width * ratio,
                  minor_tick_width=base.minor_tick_width * ratio)
    return replace(base, **kw)


# --------------------------------------------------------------------------
# font resolution
# --------------------------------------------------------------------------
def _resolve_family(style: GriffinStyle) -> str:
    from matplotlib import font_manager
    names = {f.name for f in font_manager.fontManager.ttflist}
    if style.font_family in names:
        return style.font_family
    warnings.warn(
        f"secfigs.style: {style.font_family!r} not installed; "
        f"falling back to {style.fallback_family!r}. Metrics will differ slightly."
    )
    return style.fallback_family


# --------------------------------------------------------------------------
# apply()
# --------------------------------------------------------------------------
@contextmanager
def apply(style: GriffinStyle = STYLE, **rc_overrides):
    """Context manager installing the Griffin rcParams.

    Only rcParams that are genuinely global live here; anything that has to be
    set per-Axes (four-sided ticks, locators, spine widths on an existing axes)
    is applied by :func:`axes`.
    """
    family = _resolve_family(style)
    rc = {
        "font.family": "sans-serif",
        "font.sans-serif": [family, "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": family,
        "mathtext.it": f"{family}:italic",
        "mathtext.bf": f"{family}:bold",

        "figure.figsize": style.figsize,
        "figure.dpi": style.dpi,
        "savefig.dpi": style.dpi,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.facecolor": "white",

        "axes.linewidth": style.spine_linewidth,
        "axes.edgecolor": style.foreground,
        "axes.labelsize": style.axis_label_size,
        "axes.labelweight": style.axis_label_weight,
        "axes.labelpad": style.axis_label_pad,
        "axes.labelcolor": style.foreground,
        "axes.titlesize": style.title_size,
        "axes.titleweight": style.title_weight,
        "axes.titlepad": style.title_pad,
        "axes.titlelocation": style.title_loc,
        "axes.grid": False,
        "axes.spines.top": style.show_box,
        "axes.spines.right": style.show_box,
        "axes.spines.left": style.show_box,
        "axes.spines.bottom": style.show_box,
        "axes.prop_cycle": mpl.cycler(color=list(style.palette_extended)),

        "lines.linewidth": style.trace_linewidth,
        "lines.solid_capstyle": "round",
        "lines.solid_joinstyle": "round",

        "xtick.labelsize": style.tick_labelsize,
        "ytick.labelsize": style.tick_labelsize,
        "xtick.color": style.foreground,
        "ytick.color": style.foreground,
        "xtick.direction": style.tick_direction,
        "ytick.direction": style.tick_direction,
        "xtick.top": style.ticks_all_sides,
        "xtick.bottom": True,
        "ytick.left": True,
        "ytick.right": style.ticks_all_sides,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.major.size": style.major_tick_length,
        "ytick.major.size": style.major_tick_length,
        "xtick.minor.size": style.minor_tick_length,
        "ytick.minor.size": style.minor_tick_length,
        "xtick.major.width": style.major_tick_width,
        "ytick.major.width": style.major_tick_width,
        "xtick.minor.width": style.minor_tick_width,
        "ytick.minor.width": style.minor_tick_width,
        "xtick.major.pad": style.tick_pad,
        "ytick.major.pad": style.tick_pad,

        "legend.fontsize": style.legend_fontsize,
        "legend.frameon": style.legend_frameon,
        "legend.handlelength": style.legend_handlelength,
        "legend.handletextpad": style.legend_handletextpad,
        "legend.labelspacing": style.legend_labelspacing,
        "legend.borderaxespad": style.legend_borderaxespad,
        "legend.borderpad": style.legend_borderpad,
        "legend.loc": style.legend_loc,
    }
    rc.update(rc_overrides)
    with mpl.rc_context(rc):
        yield style


def figure(style: GriffinStyle = STYLE, **kw):
    """``plt.subplots`` with the style's figsize/dpi. Call inside ``apply()``."""
    kw.setdefault("figsize", style.figsize)
    kw.setdefault("dpi", style.dpi)
    return plt.subplots(**kw)


# --------------------------------------------------------------------------
# axes(ax) -- four-sided outward ticks + full box frame
# --------------------------------------------------------------------------
def axes(ax,
         style: GriffinStyle = STYLE,
         normalized: bool = True,
         x_major: float | None = "inherit",
         x_minor: float | None = "inherit",
         y_major: float | None = None,
         y_minor: float | None = None,
         x_range: Sequence[float] | None = None,
         y_range: Sequence[float] | None = None):
    """Install the Griffin box frame and four-sided outward ticks on ``ax``.

    Works whether or not you are inside :func:`apply` -- every property is set
    explicitly here, so this is the one call that guarantees the look.

    Parameters
    ----------
    normalized
        True for min-max normalised traces: y major 0.2 / minor 0.1, labels
        formatted ``%.1f``, ``ylim = (-0.05, 1.05)`` (exactly the reference).
        False for raw mAU: y locators are left to ``y_major``/``y_minor``
        or AutoMinorLocator.
    x_major, x_minor, y_major, y_minor
        Tick spacings. ``x_major``/``x_minor`` default to the style's values
        (2.0 and AutoMinorLocator -> 0.5 at major 2).
    """
    if x_major == "inherit":
        x_major = style.x_major
    if x_minor == "inherit":
        x_minor = style.x_minor

    # --- locators ---------------------------------------------------------
    if x_major:
        ax.xaxis.set_major_locator(MultipleLocator(x_major))
    if x_minor:
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor))
    else:
        ax.xaxis.set_minor_locator(AutoMinorLocator())

    if normalized:
        step = y_major or style.y_major_normalized
        ticks = np.arange(0.0, 1.0 + step / 2.0, step)
        ax.yaxis.set_major_locator(MultipleLocator(step))
        ax.set_yticks(ticks)
        ax.set_yticklabels([f"{t:0.1f}" for t in ticks])
        ax.yaxis.set_minor_locator(
            MultipleLocator(y_minor or style.y_minor_normalized))
        pad = style.y_pad_normalized
        ax.set_ylim(-pad, 1.0 + pad)
    else:
        if y_major:
            ax.yaxis.set_major_locator(MultipleLocator(y_major))
        if y_minor:
            ax.yaxis.set_minor_locator(MultipleLocator(y_minor))
        else:
            ax.yaxis.set_minor_locator(AutoMinorLocator())

    # --- ticks: all four sides, pointing out ------------------------------
    sides = dict(top=style.ticks_all_sides, bottom=True,
                 left=True, right=style.ticks_all_sides)
    ax.tick_params(axis="both", which="major",
                   labelsize=style.tick_labelsize,
                   length=style.major_tick_length,
                   width=style.major_tick_width,
                   direction=style.tick_direction,
                   color=style.foreground,
                   pad=style.tick_pad,
                   labeltop=False, labelright=False, **sides)
    ax.tick_params(axis="both", which="minor",
                   length=style.minor_tick_length,
                   width=style.minor_tick_width,
                   direction=style.tick_direction,
                   color=style.foreground, **sides)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight(style.tick_label_weight)

    # --- box frame --------------------------------------------------------
    for spine in ax.spines.values():
        spine.set_visible(style.show_box)
        spine.set_linewidth(style.spine_linewidth)
        spine.set_color(style.foreground)

    if x_range is not None:
        ax.set_xlim(x_range)
    if y_range is not None:
        ax.set_ylim(y_range)
    return ax


# --------------------------------------------------------------------------
# text helpers
# --------------------------------------------------------------------------
def set_labels(ax, xlabel: str | None = None, ylabel: str | None = None,
               style: GriffinStyle = STYLE):
    """Axis titles at the reference size/weight/pad."""
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=style.axis_label_size,
                      fontweight=style.axis_label_weight,
                      labelpad=style.axis_label_pad)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=style.axis_label_size,
                      fontweight=style.axis_label_weight,
                      labelpad=style.axis_label_pad)
    return ax


def set_title(ax, text: str, style: GriffinStyle = STYLE):
    """Bold sans title above the axes box (measured gap 50 px = 28 pt)."""
    ax.set_title(text, fontsize=style.title_size,
                 fontweight=style.title_weight,
                 pad=style.title_pad, loc=style.title_loc)
    return ax


def legend(ax, style: GriffinStyle = STYLE, handles=None, labels=None, **kw):
    """Frameless legend inside the axes, upper-right, regular weight."""
    opts = dict(loc=style.legend_loc,
                frameon=style.legend_frameon,
                fontsize=style.legend_fontsize,
                handlelength=style.legend_handlelength,
                handletextpad=style.legend_handletextpad,
                labelspacing=style.legend_labelspacing,
                borderaxespad=style.legend_borderaxespad,
                borderpad=style.legend_borderpad)
    opts.update(kw)
    leg = ax.legend(handles, labels, **opts) if handles is not None \
        else ax.legend(**opts)
    for t in leg.get_texts():
        t.set_fontweight(style.legend_weight)
        t.set_fontsize(opts["fontsize"])
    return leg


def annotate_peak(ax, x: float, y: float, text: str, side: str = "right",
                  style: GriffinStyle = STYLE,
                  dx_pt: float | None = None, dy_pt: float | None = None,
                  marker: str | None = None, **kw):
    """The "Oligomer" / "Monomer" pointer: small filled triangle plus text.

    ``(x, y)`` are DATA coordinates and locate the TRIANGLE, not the text.
    The text sits above the triangle, offset in points.

    side='right'
        triangle '>' , text left-aligned and up-and-right of it.
        Measured offset (+2.5, +40.8) pt -- this is the "Oligomer" placement.
    side='left'
        triangle '<' , text right-aligned and up-and-left of it.
        Measured offset (-20.3, +34.3) pt -- the "Monomer" placement
        (its triangle in the reference is actually '>', pass marker='>'
        to reproduce that exactly).

    Returns ``(line, annotation)``.
    """
    if side not in ("right", "left"):
        raise ValueError("side must be 'right' or 'left'")

    if marker is None:
        marker = style.annot_marker_right if side == "right" \
            else style.annot_marker_left
    if side == "right":
        odx, ody = style.annot_offset_right
        ha = "left"
    else:
        odx, ody = style.annot_offset_left
        ha = "right"
    if dx_pt is not None:
        odx = dx_pt
    if dy_pt is not None:
        ody = dy_pt

    (line,) = ax.plot([x], [y], marker=marker, linestyle="none",
                      markersize=style.annot_marker_size,
                      color=style.annot_color,
                      markeredgewidth=0.0, clip_on=False, zorder=5)
    ann = ax.annotate(
        text, xy=(x, y), xycoords="data",
        xytext=(odx, ody), textcoords="offset points",
        ha=ha, va="baseline",
        fontsize=kw.pop("fontsize", style.annotation_size),
        fontweight=kw.pop("fontweight", style.annotation_weight),
        color=kw.pop("color", style.annot_color),
        annotation_clip=False, zorder=5, **kw)
    return line, ann


def panel_letter(ax, letter: str, style: GriffinStyle = STYLE,
                 xy: tuple | None = None, **kw):
    """Large bold capital outside the top-left corner of the axes box."""
    x, y = xy if xy is not None else style.panel_letter_xy
    return ax.text(x, y, letter, transform=ax.transAxes,
                   fontsize=kw.pop("fontsize", style.panel_letter_size),
                   fontweight=kw.pop("fontweight", style.panel_letter_weight),
                   ha="left", va="top",
                   color=kw.pop("color", style.foreground),
                   clip_on=False, **kw)


def ramp(n: int, style: GriffinStyle | None = None) -> list:
    """``n`` trace colours, taken FROM THE STYLE OBJECT.

    ``style.palette`` supplies a 3-trace panel (this is the pair 1A and 1B draw
    with, via ``ramp(3, style)[:2]``); ``style.palette_extended`` supplies every
    other length.

    ``style=None`` falls back to the module constants :data:`RAMP3` /
    :data:`RAMP9` -- which are exactly what those two fields default to, so
    ``ramp(n)`` returns what it always did and no existing caller breaks.  What
    changes is that ``ramp(n, dataclasses.replace(STYLE, palette=...))`` now
    honours the replacement.  It could not before: this function read the module
    constants unconditionally and took no style argument at all, so colour was
    the one house-style property that did not propagate through the dataclass.
    """
    pal3 = RAMP3 if style is None else tuple(style.palette)
    pal9 = RAMP9 if style is None else tuple(style.palette_extended)
    if n == 3:
        if len(pal3) < 3:
            raise ValueError(
                f"style.palette has {len(pal3)} colours; a 3-trace panel needs 3")
        return list(pal3[:3])
    if n > len(pal9):
        raise ValueError(
            f"no ramp for {n} traces (style.palette_extended has {len(pal9)})")
    return list(pal9[:n])


__all__ = ["GriffinStyle", "STYLE", "MEASUREMENTS", "PX_PER_PT",
           "RAMP3", "RAMP9", "GRIFFIN_SWATCHES", "ramp",
           "apply", "figure", "axes", "set_labels", "set_title",
           "legend", "annotate_peak", "panel_letter"]
