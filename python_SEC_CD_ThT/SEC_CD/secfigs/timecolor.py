# -*- coding: utf-8 -*-
r"""
secfigs.timecolor -- ONE colour per incubation time, shared by every SEC panel.

The problem this replaces
-------------------------
``secfigs.style.ramp(n)`` hands out colours by POSITION in each panel's own
series, so the same hex means a different time in a different panel:
``#fe9929`` is 72 h in panel b, 48 h in panel c and 24 h in panel d.  Since the
whole point of the 2 x 2 is to let a reader compare panel c (100 uM, 37 C) with
panel d (100 uM, 25 C), that is a readability defect, not a nicety.

Here colour is a pure function of t.  ``color_for(24)`` returns ``#263d6e``
whatever panel it is asked from, and the four panels become directly
comparable.

The eleven times
----------------
The union over the four panels is

    a  250 uM 37 C        2  4  6
    b   25 uM 37 C   0  3    6      12  24  48  72
    c  100 uM 37 C   0  2  4     8  12  24  48
    d  100 uM 25 C      2  4  6  8  10  12  24  48  72
    ------------------------------------------------
    union            0  2  3  4  6  8  10  12  24  48  72     (11 values)


WHY RANK AND NOT log(t + 1)
===========================
t runs 0 -> 72 h but eight of the eleven samples are at or below 12 h, so any
map monotone in t itself piles those eight into one end of the ramp.  Measured
on the ramp below -- minimum pairwise CIE76 dE over all 55 pairs, worst of
{normal vision, deuteranopia, protanopia, greyscale}:

    u = t / 72                     WORST dE   0.85    8 pairs below dE 5
    u = sqrt(t) / sqrt(72)         WORST dE   2.09    2 pairs below dE 5
    u = log10(t+1) / log10(73)     WORST dE   2.39    0 pairs below dE 5 in
                                                      colour, but greyscale
                                                      collapses: the 8/10/12 h
                                                      run lands 2.4-2.8 L*
                                                      units apart
    u = i / 10   (rank)            WORST dE   5.86    0 pairs below dE 5 in
                                                      any rendering

log(t + 1) is the best of the continuous transforms and still fails
requirement 2: it survives in colour and dies the moment the figure is printed
in greyscale.  Rank spacing is the only mapping that keeps every pair separated
in EVERY rendering, and it is the honest one here -- these eleven values are
the experiment's sampling schedule, an ordered set of design points, not eleven
draws from a continuous field.

The cost of rank spacing is that the colourbar is not a linear time axis.
:func:`add_colorbar` therefore draws it as eleven EQUAL CELLS, every one TICKED
AND LABELLED WITH THE ACTUAL MEASURED HOURS, so cell width visibly carries no
information and elapsed time has to be read off the labels.


WHY CIVIDIS
===========
Selection experiment
--------------------
Search space for each candidate map: its distinct 8-bit entries (the palette
has to survive being written as hex -- 255 for cividis).

Constraint: every member L* <= 78, so the palest still reads as a 0.85 pt
stroke on white.  This is the same threshold ``style.RAMP9`` was closed
against.

Objective: maximise the minimum pairwise CIE76 dE over the TWELVE inks that
actually land on the page -- the 11 trace colours PLUS ``#000000``, which is
what ``style.foreground`` draws the spines, ticks, tick labels, axis labels and
the Oligomer/Monomer pointers in -- taken as the worst of {normal,
deuteranopia, protanopia, greyscale}.

Members are placed at uniform L* spacing rather than uniform colormap
coordinate, because greyscale dE between two colours is EXACTLY |dL*|, so the
greyscale floor is min adjacent |dL*| and uniform L* is what maximises it.
Including black in the ink set is what stops the optimiser cheating: without
it, magma and inferno win by driving the 72 h member to L* 5, which is a trace
no reader can tell from the frame.

    cmap      L*hi  L*lo   normal  deutan  protan  grey  tritan   WORST
    magma     77.7   6.9    10.44    6.58    9.27  6.85    5.35    6.58
    YlGnBu    78.0  13.4    10.74    7.49    8.67  6.22    5.18    6.22
    cividis   77.9  13.8     8.49    8.90    8.38  6.20    9.37    6.20
    viridis   78.0  15.2     9.66    6.19    8.23  6.16    4.53    6.16
    inferno   77.6   7.4     9.85    8.25    5.79  6.82    5.50    5.79
    Blues     77.9  20.9     6.48    5.63    5.98  5.54    5.87    5.54
    plasma    77.8  15.5     9.93    7.51    4.52  6.10    5.44    4.52
    YlOrRd    77.7  26.0     6.50    4.42    6.70  5.03    3.45    4.42

The single-hue ColorBrewer maps lose outright, and for a structural reason:
clipped at L* <= 78 a single-hue ramp has only its own dark end left to spend,
and neither reaches far enough -- Blues bottoms out at L* 20.9, YlOrRd at
L* 26.0, so their L* budget is 57 and 52 units against 64 for cividis.  Chroma
is the only axis they vary on, and chroma is exactly what a colour-blind or
greyscale reader loses.  Blues 5.54 and YlOrRd 4.42 are the consequence.

Tiebreak
--------
The top four are a tie: 6.58 / 6.22 / 6.20 / 6.16 is a 6% spread.  What
separates them is WHERE the worst case lives.

cividis is the only candidate whose limit is the GREYSCALE channel -- the L*
budget, which is physics and which no colormap can beat.  Its three
colour-vision channels all sit ABOVE that floor (deutan 8.90, protan 8.38,
tritan 9.37), so a dichromat sees the series at least as well separated as a
greyscale print does, and very nearly as well as a trichromat.  That is the
property cividis was built for, and for a figure whose entire job is to convey
an ORDER it is the property that matters: every reader gets the same ordering
cue, not merely a legible one.

Every other candidate is limited by a colour-vision deficiency instead, i.e. a
CVD reader sees a strictly worse figure than a normal reader:

    magma    headline 6.58 is deuteranopia; and under tritanopia its 0 h / 2 h
             pair falls to dE76 5.35 / dE2000 2.83, which is at the JND, while
             its 72 h member sits at L* 6.9 (contrast 18 : 1) -- effectively
             ink.  Its apparent win is bought at both ends.
    viridis  6.16 deuteranopia, and 4.53 under tritanopia
    inferno  5.79 protanopia
    plasma   4.52 protanopia
    YlGnBu   6.22 is greyscale, but tritanopia drops to 5.18

cividis's known cost is low chroma at mid-ramp (C* 2.4 at t = 6 h): the middle
of the series is near-neutral.  That is bought deliberately -- it is the price
of routing all of the discriminability through luminance, the one channel every
reader and every printer has.


MEASURED PROPERTIES OF THE DELIVERED ELEVEN
===========================================
    lightest member    #d4c15f   L* 77.85   contrast  1.81 : 1 vs white
    darkest member     #00224e   L* 13.79   contrast 15.69 : 1 vs white,
                                            and dE76 34.4 from the black frame
    L* strictly decreasing with t, adjacent step 6.20 - 6.56
    relative luminance strictly decreasing, min adjacent |dY| 0.0130
    L* budget 64.06 units over 10 steps -> greyscale floor 6.20 is 97% of the
    6.41 theoretical maximum for this colormap under the L* <= 78 light clip

    rendering        min dE76              min dE2000    median dE76
    normal            8.49  (6 h / 8 h)      5.06           41.58
    deuteranopia      8.90  (6 h / 8 h)      5.38           43.92
    protanopia        8.38  (6 h / 8 h)      5.31           41.61
    tritanopia        9.37  (0 h / 2 h)      5.34           64.33
    greyscale         6.20  (12 h / 24 h)    4.18           25.34

Dichromacy is simulated with the Vienot, Brettel & Mollon (1999) single-plane
projection in LMS applied to LINEAR sRGB -- see :func:`simulate_cvd`, which is
implemented here rather than imported so every number above can be re-derived
by running ``python -m secfigs.timecolor``.  The visual proof is
``figures/timecolor_probe.png``.

Usage
-----
    from secfigs import timecolor as tc

    for t, (x, y) in zip(panel_times, traces):
        ax.plot(x, y, color=tc.color_for(t), lw=0.85, label=f"{t} h")

    tc.add_colorbar(fig, [ax_a, ax_b, ax_c, ax_d])
"""

from __future__ import annotations

import itertools

import numpy as np

__all__ = [
    "TIMES", "COLORS", "CMAP", "NORM", "BOUNDARIES", "FRAME_INK",
    "COLORS_MAGMA", "CMAP_MAGMA", "NORM_MAGMA", "FRAME_INK_GREY",
    "PALETTES", "DEFAULT_PALETTE", "palette_colors", "cmap_norm_for",
    "color_for", "colors_for", "rank_for", "add_colorbar",
    "srgb_to_linear", "linear_to_srgb", "rgb_to_lab", "relative_luminance",
    "delta_e76", "delta_e2000", "simulate_cvd", "to_greyscale",
    "derive_from_cmap", "derive_from_cividis", "derive_from_magma",
    "frame_grey_search", "validate", "report",
]


# ==========================================================================
# 1. The palette
# ==========================================================================
#: The union of incubation times over panels a-d, in hours, ascending.
TIMES: tuple[int, ...] = (0, 2, 3, 4, 6, 8, 10, 12, 24, 48, 72)

#: What style.foreground draws the frame, ticks and text in.  Counted as a
#: twelfth ink during selection; see the module docstring.
FRAME_INK: str = "#000000"

#: cividis coordinates the eleven members were taken from.  Reproduced exactly
#: by derive_from_cividis().
_CIVIDIS_U: tuple[float, ...] = (
    0.83986, 0.75782, 0.67579, 0.58985, 0.50391, 0.41407,
    0.32813, 0.23828, 0.15235, 0.06250, 0.00000,
)

#: One colour per entry of :data:`TIMES`, light (t = 0) -> dark (t = 72).
COLORS: tuple[str, ...] = (
    "#d4c15f",   #  0 h   L* 77.85   C* 51.3
    "#beaf6b",   #  2 h   L* 71.33   C* 37.0
    "#a89e73",   #  3 h   L* 64.98   C* 24.0
    "#928d78",   #  4 h   L* 58.52   C* 12.0
    "#7e7c78",   #  6 h   L* 52.07   C*  2.4
    "#6a6c71",   #  8 h   L* 45.61   C*  3.1
    "#565c6d",   # 10 h   L* 39.13   C* 10.6
    "#404c6c",   # 12 h   L* 32.57   C* 20.5
    "#263d6e",   # 24 h   L* 26.38   C* 32.2
    "#002e6a",   # 48 h   L* 20.02   C* 40.8
    "#00224e",   # 72 h   L* 13.79   C* 31.5
)

assert len(TIMES) == len(COLORS) == len(_CIVIDIS_U) == 11


# --------------------------------------------------------------------------
# 1b.  The MAGMA alternative -- same construction, different tiebreak
# --------------------------------------------------------------------------
# Added 2026-08-08 for the Nature-style revision.  NOTHING above is touched:
# COLORS / CMAP / NORM / FRAME_INK still mean cividis-on-black, and every
# existing caller keeps the palette it was validated against.
#
# WHY A SECOND PALETTE EXISTS
# ---------------------------
# The cividis tiebreak (see above) assumed a PURE BLACK frame and therefore
# counted #000000 as a twelfth ink.  Under that constraint magma's win in the
# headline metric was declared "bought at both ends", because its darkest
# member lands at L* 6.86 -- dE76 6.86 from black in greyscale, i.e. a trace a
# reader cannot tell from the spine it runs along.
#
# That objection is an objection to the FRAME, not to the colormap.  Draw the
# spines, ticks and tick labels in the dark grey :data:`FRAME_INK_GREY` instead
# and it evaporates: the darkest member is then dE76 18.01 from the frame ink
# (worst of {normal, deutan, protan, greyscale}), against 6.86 for black.  With
# that repaired, magma is strictly the better ramp on this figure's own metric:
#
#     worst-of-{normal, deutan, protan, greyscale} min dE76
#         magma    6.58        cividis  6.20
#     normal-vision min dE76
#         magma   10.44        cividis  8.49
#     panel a (2 / 4 / 6 h only -- one sixth of the 0-72 h ramp)
#         magma   20.61        cividis 11.58        (+78%)
#
# and it is what the one 2025 precedent for this exact figure type uses
# (Sneideriene et al., PNAS 2025, Fig 3D: traces indexed by time in hours under
# an unmodified matplotlib magma).
#
# Magma's remaining weakness is unchanged and is recorded honestly: under
# TRITANOPIA the 0 h / 2 h pair falls to dE76 5.35 / dE2000 2.83, at the JND.
# Cividis does not have that failure (9.37).  Tritanopia is ~0.01% of the
# population against ~6% of males for deutan/protan, and panel a -- the only
# panel where 0-2 h separation is load bearing -- carries an explicit legend,
# so the trade is taken deliberately, not by accident.

#: magma coordinates the eleven members were taken from.  Reproduced exactly by
#: ``derive_from_magma()``.
_MAGMA_U: tuple[float, ...] = (
    0.82814, 0.76564, 0.69923, 0.62892, 0.55470, 0.47657,
    0.39844, 0.31641, 0.23828, 0.16407, 0.10156,
)

#: One colour per entry of :data:`TIMES`, light (t = 0) -> dark (t = 72).
#: Same construction as :data:`COLORS`: distinct 8-bit magma entries, light end
#: clipped at L* <= 78 so the palest still reads as a 0.85 pt stroke on white,
#: eleven members at the largest uniform L* gap the map can support.
COLORS_MAGMA: tuple[str, ...] = (
    "#feae77",   #  0 h   L* 77.69   C* 46.3
    "#fc9065",   #  2 h   L* 70.59   C* 54.7
    "#f7705c",   #  3 h   L* 63.46   C* 62.0
    "#e85362",   #  4 h   L* 56.32   C* 63.0
    "#cd4071",   #  6 h   L* 49.30   C* 59.1
    "#ad347c",   #  8 h   L* 42.46   C* 57.1
    "#8c2981",   # 10 h   L* 35.49   C* 58.7
    "#6b1d81",   # 12 h   L* 28.35   C* 62.6
    "#4c117a",   # 24 h   L* 21.37   C* 65.2
    "#2c115f",   # 48 h   L* 14.19   C* 52.7
    "#150e38",   # 72 h   L*  6.86   C* 31.3
)

#: Axis ink for the magma build: spines, ticks, tick labels, axis labels.
#:
#: CHOSEN BY MEASUREMENT, not taste -- see :func:`frame_grey_search`, which
#: sweeps all 8-bit neutral greys and is what produced this value.  Two
#: constraints, then one objective:
#:
#:   1. worst-of-{normal, deutan, protan, greyscale} dE76 from the DARKEST
#:      member (#150e38) >= 10.  Satisfied from #2a2a2a (L* 17.0) upward;
#:      pure black gives only 6.86 and is what this replaces.
#:   2. contrast >= 7 : 1 against white, so 7 pt tick labels still read as ink
#:      (WCAG AAA for normal-size text).  Satisfied up to #595959 (L* 37.9).
#:
#:   objective: maximise the minimum |dL*| to ANY of the eleven members, since
#:   greyscale dE is exactly |dL*| and the frame must not share a lightness
#:   with a trace in a greyscale print.
#:
#: 48 greys clear both constraints, #2a2a2a (L* 17.06) through #595959
#: (L* 37.82).  The objective tops out at 3.54 and TWO greys are inside 0.10
#: L* of that, in different inter-member gaps -- #4b4b4b (3.54) and #3b3b3b
#: (3.48).  The objective cannot separate them, so contrast does: #4b4b4b
#: holds only 8.72 : 1 and reads washed out at 7 pt, #3b3b3b holds 11.20 : 1.
#: #3b3b3b sits at L* 24.87, 3.48 off the nearest member, in a gap whose own
#: midpoint would give 3.49 -- 99.7% of the best that gap allows.  It lands
#: almost exactly between the 12 h member (L* 28.35) and the 24 h (L* 21.37).
#:
#:     #3b3b3b   L* 24.87   contrast 11.20 : 1 vs white
#:               dE76 to the darkest member #150e38 :
#:                   normal 36.10 / deutan 34.61 / protan 35.69 /
#:                   greyscale 18.01        (pure black scores 6.86)
FRAME_INK_GREY: str = "#3b3b3b"

#: ``name -> eleven hexes``.  ``palette_colors`` and friends resolve through it.
PALETTES: dict[str, tuple[str, ...]] = {
    "cividis": COLORS,
    "magma": COLORS_MAGMA,
}

#: What a bare ``color_for(t)`` means.  Unchanged, so no existing caller moves.
DEFAULT_PALETTE: str = "cividis"

#: The frame ink each palette was closed against.
FRAME_INKS: dict[str, str] = {
    "cividis": FRAME_INK,
    "magma": FRAME_INK_GREY,
}

assert len(COLORS_MAGMA) == len(_MAGMA_U) == len(TIMES) == 11

_INDEX = {float(t): i for i, t in enumerate(TIMES)}


def palette_colors(palette: str = DEFAULT_PALETTE) -> tuple[str, ...]:
    """The eleven hexes of a named palette. ``'cividis'`` or ``'magma'``."""
    try:
        return PALETTES[palette]
    except KeyError:
        raise KeyError(f"unknown palette {palette!r}; have {sorted(PALETTES)}")


def rank_for(t) -> int:
    """Ordinal position of ``t`` in :data:`TIMES` (0-10).

    This is the quantity :data:`NORM` is defined over.
    """
    key = float(t)
    if key in _INDEX:
        return _INDEX[key]
    diffs = [abs(key - float(v)) for v in TIMES]      # tolerate 24.0000001
    k = int(np.argmin(diffs))
    if diffs[k] <= 1e-6:
        return k
    raise KeyError(
        f"t = {t!r} h is not one of the measured incubation times {TIMES}. "
        "timecolor is a lookup, not an interpolation: if a new time point was "
        "collected, add it to TIMES and re-run derive_from_cividis()."
    )


def color_for(t, palette: str = DEFAULT_PALETTE) -> str:
    """Hex colour for incubation time ``t`` in hours. Identical in every panel."""
    return palette_colors(palette)[rank_for(t)]


def colors_for(ts, palette: str = DEFAULT_PALETTE) -> list[str]:
    """``[color_for(t) for t in ts]``."""
    return [color_for(t, palette) for t in ts]


# ==========================================================================
# 2. matplotlib objects for the colourbar
# ==========================================================================
def _build_cmap_norm(colors=None, name: str = "sec_incubation_time"):
    from matplotlib.colors import BoundaryNorm, ListedColormap
    cmap = ListedColormap(list(COLORS if colors is None else colors), name=name)
    cmap.set_bad("#ffffff")
    bounds = np.arange(len(TIMES) + 1, dtype=float) - 0.5   # -0.5 .. 10.5
    return cmap, BoundaryNorm(bounds, cmap.N), bounds


#: ``NORM`` is defined over RANK (0-10), not over hours -- the mapping is
#: ordinal, see the module docstring.  Feed it ``rank_for(t)``;
#: :func:`add_colorbar` relabels those ranks with the measured hours so the
#: reader never sees a rank.
CMAP, NORM, BOUNDARIES = _build_cmap_norm()

#: The same three objects for the magma build.  BOUNDARIES is shared -- the
#: eleven equal cells are a property of the rank mapping, not of the colours.
CMAP_MAGMA, NORM_MAGMA, _ = _build_cmap_norm(
    COLORS_MAGMA, name="sec_incubation_time_magma")

_CMAP_NORM = {"cividis": (CMAP, NORM), "magma": (CMAP_MAGMA, NORM_MAGMA)}


def cmap_norm_for(palette: str = DEFAULT_PALETTE):
    """``(ListedColormap, BoundaryNorm)`` for a named palette."""
    palette_colors(palette)                     # validates the name
    return _CMAP_NORM[palette]


def add_colorbar(fig, ax_or_axes=None, *, cax=None,
                 orientation: str = "vertical",
                 label: str = "Incubation time (h)",
                 fraction: float = 0.030, pad: float = 0.018,
                 aspect: float = 34.0, labelsize: float | None = None,
                 titlesize: float | None = None, linewidth: float = 0.6,
                 labelpad: float | None = None, tickpad: float = 2.0,
                 palette: str = DEFAULT_PALETTE,
                 framecolor: str | None = None,
                 **kw):
    """Discrete colourbar: eleven equal cells labelled with the measured hours.

    Parameters
    ----------
    fig
        Figure to attach to.
    ax_or_axes
        A single Axes or a sequence of Axes to steal space from -- pass all
        four panels to get a bar spanning the whole 2 x 2.  Ignored when
        ``cax`` is given.
    cax
        An Axes to draw the bar INTO, at a rect the caller has already derived.
        Mutually exclusive with ``ax_or_axes``/``fraction``/``pad``/``aspect``:
        when it is supplied the bar occupies exactly ``cax`` and steals no
        space from anything.  This is the path a figure with a hand-derived
        padding budget wants.
    orientation
        ``'vertical'`` (default) or ``'horizontal'``.
    label
        Colourbar title.  The unit lives in the label, so the tick labels are
        bare numbers.
    palette
        ``'cividis'`` (default, unchanged) or ``'magma'``.
    framecolor
        Ink for the bar's outline and its tick labels.  Defaults to whatever
        frame ink the chosen palette was closed against -- ``#000000`` for
        cividis, :data:`FRAME_INK_GREY` for magma.

    Notes
    -----
    The cells are EQUAL WIDTH because the mapping is by rank; the tick labels
    are the real times (0, 2, 3, 4, 6, 8, 10, 12, 24, 48, 72 h).  Cell width
    therefore carries no information and elapsed time must be read off the
    labels -- which is why they are all drawn.
    """
    from matplotlib.cm import ScalarMappable

    cmap, norm = cmap_norm_for(palette)
    ink = FRAME_INKS[palette] if framecolor is None else framecolor
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    common = dict(orientation=orientation, boundaries=BOUNDARIES,
                  values=np.arange(len(TIMES)), spacing="uniform",
                  ticks=np.arange(len(TIMES)))
    if cax is not None:
        # fraction / pad / aspect are make_axes parameters; with an explicit
        # cax there is no make_axes call and passing them is a TypeError.
        cbar = fig.colorbar(sm, cax=cax, **common, **kw)
    else:
        cbar = fig.colorbar(sm, ax=ax_or_axes, fraction=fraction, pad=pad,
                            aspect=aspect, **common, **kw)

    cbar.set_ticklabels([str(t) for t in TIMES])
    if titlesize is not None:
        cbar.set_label(label, size=titlesize)
    else:
        cbar.set_label(label)
    if labelpad is not None:
        cbar._long_axis().labelpad = labelpad
    if labelsize is not None:
        cbar.ax.tick_params(labelsize=labelsize)
    cbar.ax.tick_params(length=0.0, pad=tickpad, colors=ink)
    cbar.outline.set_linewidth(linewidth)
    cbar.outline.set_edgecolor(ink)
    cbar.dividers.set_linewidth(0.0)
    return cbar


# ==========================================================================
# 3. Colour science -- implemented, not assumed
# ==========================================================================
_M_RGB2XYZ = np.array([[0.4124564, 0.3575761, 0.1804375],
                       [0.2126729, 0.7151522, 0.0721750],
                       [0.0193339, 0.1191920, 0.9503041]])
_WHITE_D65 = np.array([0.95047, 1.00000, 1.08883])


def srgb_to_linear(c):
    """sRGB 0-1 -> linear-light 0-1 (IEC 61966-2-1)."""
    c = np.asarray(c, float)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def linear_to_srgb(c):
    """Linear-light 0-1 -> sRGB 0-1, clipped to gamut."""
    c = np.clip(np.asarray(c, float), 0.0, 1.0)
    return np.where(c <= 0.0031308, c * 12.92, 1.055 * c ** (1 / 2.4) - 0.055)


def relative_luminance(rgb):
    """CIE Y (0-1) of an sRGB triple -- the greyscale channel."""
    return (srgb_to_linear(rgb) @ _M_RGB2XYZ.T)[..., 1]


def _f_lab(t):
    d = 6.0 / 29.0
    return np.where(t > d ** 3, np.cbrt(t), t / (3 * d * d) + 4.0 / 29.0)


def rgb_to_lab(rgb):
    """sRGB 0-1 -> CIE L*a*b* under D65."""
    xyz = (srgb_to_linear(rgb) @ _M_RGB2XYZ.T) / _WHITE_D65
    fx, fy, fz = _f_lab(xyz[..., 0]), _f_lab(xyz[..., 1]), _f_lab(xyz[..., 2])
    return np.stack([116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz)], axis=-1)


def delta_e76(lab1, lab2) -> float:
    """CIE76 colour difference (Euclidean distance in L*a*b*)."""
    return float(np.linalg.norm(np.asarray(lab1) - np.asarray(lab2)))


def delta_e2000(lab1, lab2) -> float:
    """CIEDE2000 colour difference, kL = kC = kH = 1."""
    L1, a1, b1 = (float(v) for v in lab1)
    L2, a2, b2 = (float(v) for v in lab2)
    C1, C2 = np.hypot(a1, b1), np.hypot(a2, b2)
    Cbar = 0.5 * (C1 + C2)
    G = 0.5 * (1 - np.sqrt(Cbar ** 7 / (Cbar ** 7 + 25.0 ** 7))) if Cbar else 0.5
    a1p, a2p = (1 + G) * a1, (1 + G) * a2
    C1p, C2p = np.hypot(a1p, b1), np.hypot(a2p, b2)
    h1p = np.degrees(np.arctan2(b1, a1p)) % 360.0
    h2p = np.degrees(np.arctan2(b2, a2p)) % 360.0
    dLp, dCp = L2 - L1, C2p - C1p
    if C1p * C2p == 0:
        dhp = 0.0
    else:
        dhp = h2p - h1p
        dhp += -360.0 if dhp > 180 else (360.0 if dhp < -180 else 0.0)
    dHp = 2 * np.sqrt(C1p * C2p) * np.sin(np.radians(dhp) / 2)
    Lbp, Cbp = 0.5 * (L1 + L2), 0.5 * (C1p + C2p)
    if C1p * C2p == 0:
        hbp = h1p + h2p
    elif abs(h1p - h2p) <= 180:
        hbp = 0.5 * (h1p + h2p)
    elif h1p + h2p < 360:
        hbp = 0.5 * (h1p + h2p + 360)
    else:
        hbp = 0.5 * (h1p + h2p - 360)
    T = (1 - 0.17 * np.cos(np.radians(hbp - 30))
         + 0.24 * np.cos(np.radians(2 * hbp))
         + 0.32 * np.cos(np.radians(3 * hbp + 6))
         - 0.20 * np.cos(np.radians(4 * hbp - 63)))
    dtheta = 30 * np.exp(-(((hbp - 275) / 25.0) ** 2))
    Rc = 2 * np.sqrt(Cbp ** 7 / (Cbp ** 7 + 25.0 ** 7)) if Cbp else 0.0
    Sl = 1 + (0.015 * (Lbp - 50) ** 2) / np.sqrt(20 + (Lbp - 50) ** 2)
    Sc = 1 + 0.045 * Cbp
    Sh = 1 + 0.015 * Cbp * T
    Rt = -np.sin(np.radians(2 * dtheta)) * Rc
    return float(np.sqrt((dLp / Sl) ** 2 + (dCp / Sc) ** 2 + (dHp / Sh) ** 2
                         + Rt * (dCp / Sc) * (dHp / Sh)))


# -- Vienot, Brettel & Mollon (1999) dichromacy simulation ------------------
# LMS primaries as tabulated in that paper; the transform acts on LINEAR sRGB.
_M_RGB2LMS = np.array([[17.8824, 43.5161, 4.11935],
                       [3.45565, 27.1554, 3.86714],
                       [0.0299566, 0.184309, 1.46709]])
_M_LMS2RGB = np.linalg.inv(_M_RGB2LMS)

#: Projection onto the dichromat's reduced colour plane, in LMS.
_PLANE = {
    # protanope: the L cone is missing, L reconstructed from M and S
    "protan": np.array([[0.0, 2.02344, -2.52581],
                        [0.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0]]),
    # deuteranope: the M cone is missing, M reconstructed from L and S
    "deutan": np.array([[1.0, 0.0, 0.0],
                        [0.494207, 0.0, 1.24827],
                        [0.0, 0.0, 1.0]]),
    # tritanope: the S cone is missing (Brettel-form single-plane approximation)
    "tritan": np.array([[1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [-0.395913, 0.801109, 0.0]]),
}


def simulate_cvd(rgb, kind: str):
    """Simulate dichromatic vision. ``kind`` in {'protan', 'deutan', 'tritan'}.

    Vienot, Brettel & Mollon, *Color Res. Appl.* **24**, 243 (1999): convert
    linear sRGB to LMS, project onto the plane the dichromat's two surviving
    cone classes span, convert back.
    """
    if kind not in _PLANE:
        raise ValueError(f"kind must be one of {sorted(_PLANE)}, got {kind!r}")
    lin = srgb_to_linear(np.asarray(rgb, float))
    lms = lin @ _M_RGB2LMS.T
    return linear_to_srgb((lms @ _PLANE[kind].T) @ _M_LMS2RGB.T)


def to_greyscale(rgb):
    """Luminance-correct desaturation: keep CIE Y, drop chroma."""
    y = relative_luminance(rgb)
    g = linear_to_srgb(y)
    return np.stack([g, g, g], axis=-1)


# ==========================================================================
# 4. Derivation and self-validation
# ==========================================================================
def derive_from_cmap(cmap_name: str, l_max: float = 78.0, n: int = 11,
                     samples: int = 65536, frame_lstar: float = 0.0):
    """Re-derive an eleven-member ramp from scratch. Returns ``(hexes, u, L*)``.

    The recipe, exactly as run for the selection experiment:

    1. Enumerate the DISTINCT 8-BIT entries of ``cmap_name`` (255 for cividis,
       255 for magma) -- the palette has to survive being written as hex.
    2. Discard everything lighter than ``l_max``.
    3. Binary-search the largest ``g`` for which ``n`` entries can be chosen
       with L* descending, each at least ``g`` below the last AND the last
       itself at least ``g`` above ``frame_lstar`` -- the frame ink counts as
       an extra member of the set.  Greedily taking the highest L* that still
       clears the gap is exact for this problem.

    Step 3 maximises the greyscale floor, because greyscale dE is exactly
    |dL*|; every other rendering has more room and comes out above it.

    ``frame_lstar`` is 0 (pure black) for both delivered palettes.  The magma
    build keeps it at 0 even though it prints its frame in
    :data:`FRAME_INK_GREY`: the constraint's only job is to stop the optimiser
    running the dark end down to L* ~ 0, and the L* budget it produces is the
    same either way.  The grey is then chosen against the FINISHED ramp by
    :func:`frame_grey_search`, which is a separate, later decision.
    """
    import matplotlib as mpl
    from matplotlib.colors import to_hex
    u = np.linspace(0.0, 1.0, samples)
    raw = np.array([mpl.colormaps[cmap_name](v)[:3] for v in u])
    q8 = np.round(raw * 255).astype(int)
    uniq, first = np.unique(q8, axis=0, return_index=True)
    cols = uniq / 255.0
    us = u[first]
    L = rgb_to_lab(cols)[:, 0]
    o = np.argsort(-L)
    cols, us, L = cols[o], us[o], L[o]
    start = int(np.argmax(L <= l_max))

    def feasible(gap):
        picks, cur = [start], L[start]
        for _ in range(n - 1):
            c = np.where(L <= cur - gap)[0]
            if len(c) == 0:
                return None
            picks.append(int(c[0]))
            cur = L[picks[-1]]
        return picks if cur >= frame_lstar + gap else None

    lo, hi, best = 0.0, 12.0, None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        p = feasible(mid)
        if p is not None:
            best, lo = p, mid
        else:
            hi = mid
    return ([to_hex(cols[i]) for i in best],
            [round(float(us[i]), 5) for i in best],
            L[best])


def derive_from_cividis(l_max: float = 78.0, n: int = 11, samples: int = 65536):
    """Re-derive :data:`COLORS`. Thin wrapper on :func:`derive_from_cmap`."""
    return derive_from_cmap("cividis", l_max=l_max, n=n, samples=samples)


def derive_from_magma(l_max: float = 78.0, n: int = 11, samples: int = 65536):
    """Re-derive :data:`COLORS_MAGMA`. Wrapper on :func:`derive_from_cmap`."""
    return derive_from_cmap("magma", l_max=l_max, n=n, samples=samples)


def frame_grey_search(palette: str = "magma", min_de: float = 10.0,
                      min_contrast: float = 7.0, tie_tol: float = 0.10) -> dict:
    """Choose the axis-ink grey for ``palette`` by sweeping the 8-bit greys.

    This is the function that produced :data:`FRAME_INK_GREY`; running it is
    how that hex is audited.  See the constant's own docstring for the two
    constraints and the objective.  Returns the winner plus the whole
    feasible band, so the margin either side of the pick is visible.

    ``tie_tol`` closes the last 8-bit step honestly.  Two greys sit in
    different inter-member gaps and score within 0.06 L* of each other
    (#4b4b4b 3.54, #3b3b3b 3.48); at that separation the objective cannot tell
    them apart, so anything within ``tie_tol`` of the best is treated as tied
    and CONTRAST breaks it -- the darkest tied grey wins, because the whole
    point of the ink is to read as ink.
    """
    from matplotlib.colors import to_hex, to_rgb
    cols = np.array([to_rgb(c) for c in palette_colors(palette)])
    Ls = rgb_to_lab(cols)[:, 0]
    dark = cols[-1]

    rows = []
    for v in range(256):
        g = np.array([v / 255.0] * 3)
        Lg = float(rgb_to_lab(g)[0])
        Y = float(relative_luminance(g))
        des = {}
        for nm, fn in (("normal", lambda a: a),
                       ("deutan", lambda a: simulate_cvd(a, "deutan")),
                       ("protan", lambda a: simulate_cvd(a, "protan")),
                       ("greyscale", to_greyscale)):
            lab = rgb_to_lab(fn(np.stack([g, dark])))
            des[nm] = delta_e76(lab[0], lab[1])
        rows.append(dict(hex=to_hex(g), lstar=Lg,
                         contrast=1.05 / (Y + 0.05), de_vs_darkest=des,
                         worst_de=min(des.values()),
                         min_dlstar_any=float(np.min(np.abs(Ls - Lg)))))

    band = [r for r in rows
            if r["worst_de"] >= min_de and r["contrast"] >= min_contrast]
    if not band:
        raise RuntimeError("no grey satisfies both constraints")
    top = max(r["min_dlstar_any"] for r in band)
    tied = [r for r in band if r["min_dlstar_any"] >= top - tie_tol]
    best = max(tied, key=lambda r: r["contrast"])       # darkest of the tied
    # ceiling on "distance to the nearest member": half the widest L* gap the
    # ramp leaves anywhere.  A grey can only ever sit at the middle of a gap.
    ceiling = 0.5 * float(np.abs(np.diff(Ls)).max())
    # ...and the ceiling for the particular gap the winner landed in
    k = int(np.argmin(np.abs(Ls - best["lstar"])))
    nb = [Ls[i] for i in (k - 1, k + 1) if 0 <= i < len(Ls)]
    local = 0.5 * min(abs(best["lstar"] - x) + abs(best["lstar"] - Ls[k])
                      for x in nb)
    return dict(all=rows, band=band, best=best, tied=tied,
                band_hex=(band[0]["hex"], band[-1]["hex"]),
                band_lstar=(band[0]["lstar"], band[-1]["lstar"]),
                objective_top=top, objective_ceiling=ceiling,
                objective_local_ceiling=local,
                delivered=FRAME_INKS[palette])


def _rgbs(palette: str = DEFAULT_PALETTE):
    from matplotlib.colors import to_rgb
    return np.array([to_rgb(c) for c in palette_colors(palette)])


def report(palette: str = DEFAULT_PALETTE) -> dict:
    """Every measured number quoted in the module docstring, recomputed."""
    from matplotlib.colors import to_rgb
    rgbs = _rgbs(palette)
    lab = rgb_to_lab(rgbs)
    Y = relative_luminance(rgbs)
    pairs = list(itertools.combinations(range(len(TIMES)), 2))
    black_lab = rgb_to_lab(np.array(to_rgb(FRAME_INKS[palette])))
    out = {
        "lstar": lab[:, 0],
        "cstar": np.hypot(lab[:, 1], lab[:, 2]),
        "luminance": Y,
        "contrast_vs_white": 1.05 / (Y + 0.05),
        "lstar_monotone_decreasing": bool(np.all(np.diff(lab[:, 0]) < 0)),
        "luminance_monotone_decreasing": bool(np.all(np.diff(Y) < 0)),
        "min_adjacent_dLstar": float(np.abs(np.diff(lab[:, 0])).min()),
        "max_adjacent_dLstar": float(np.abs(np.diff(lab[:, 0])).max()),
        "min_adjacent_dY": float(np.abs(np.diff(Y)).min()),
        "darkest_to_frame_de76": delta_e76(lab[-1], black_lab),
        # greyscale dE is exactly |dL*|; for the black cividis frame this is
        # the darkest member's own L*, which is what the docstring quotes.
        "darkest_to_frame_grey": float(abs(lab[-1, 0] - black_lab[0])),
        "frame_ink": FRAME_INKS[palette],
        "palette": palette,
    }
    for name, arr in (("normal", rgbs),
                      ("deuteranopia", simulate_cvd(rgbs, "deutan")),
                      ("protanopia", simulate_cvd(rgbs, "protan")),
                      ("tritanopia", simulate_cvd(rgbs, "tritan")),
                      ("greyscale", to_greyscale(rgbs))):
        la = rgb_to_lab(arr)
        d76 = np.array([delta_e76(la[i], la[j]) for i, j in pairs])
        d00 = np.array([delta_e2000(la[i], la[j]) for i, j in pairs])
        i, j = pairs[int(d76.argmin())]
        out[name] = {
            "min_de76": float(d76.min()),
            "min_de76_pair_hours": (TIMES[i], TIMES[j]),
            "min_de2000": float(d00.min()),
            "median_de76": float(np.median(d76)),
        }
    return out


def validate(verbose: bool = True, palette: str = DEFAULT_PALETTE) -> dict:
    """Assert every claim in the docstring. Raises AssertionError on failure."""
    r = report(palette)
    assert r["lstar_monotone_decreasing"], "L* is not monotone in t"
    assert r["luminance_monotone_decreasing"], "luminance is not monotone in t"
    assert r["lstar"].max() <= 78.0, "lightest member is paler than L* 78"
    assert r["greyscale"]["min_de76"] > 6.0, "greyscale pairs collapse"
    assert r["darkest_to_frame_grey"] > r["greyscale"]["min_de76"], \
        "darkest trace is closer to the frame ink than traces are to each other"
    if palette == "cividis":
        # cividis was SELECTED on this property -- greyscale, i.e. the L*
        # budget, must be the binding channel.  magma is not: its binding
        # channel is deuteranopia (6.58) and it is chosen on the headline
        # figure instead.  Asserting it for magma would be asserting a claim
        # the magma docstring explicitly does not make.
        for k in ("normal", "deuteranopia", "protanopia", "tritanopia"):
            assert r[k]["min_de76"] > r["greyscale"]["min_de76"], \
                f"{k} is the binding channel, not greyscale -- tiebreak fails"
            assert r[k]["min_de2000"] > 5.0, f"{k} min dE2000 too small"
    for k in ("normal", "deuteranopia", "protanopia"):
        assert r[k]["min_de76"] > 6.0, f"{k} min dE76 too small"
    if verbose:
        print(f"PALETTE {palette}   frame ink {r['frame_ink']}")
        print(f"TIMES  {TIMES}")
        print(f"COLORS {list(palette_colors(palette))}")
        print(f"L*     {np.round(r['lstar'], 2).tolist()}")
        print(f"C*     {np.round(r['cstar'], 1).tolist()}")
        print(f"contrast vs white {np.round(r['contrast_vs_white'], 2).tolist()}")
        print(f"L* monotone decreasing : {r['lstar_monotone_decreasing']}  "
              f"(adjacent step {r['min_adjacent_dLstar']:.2f} - "
              f"{r['max_adjacent_dLstar']:.2f})")
        print(f"Y  monotone decreasing : {r['luminance_monotone_decreasing']}  "
              f"(min |dY| {r['min_adjacent_dY']:.4f})")
        print(f"darkest trace vs frame ink: dE76 "
              f"{r['darkest_to_frame_de76']:.2f}, greyscale "
              f"{r['darkest_to_frame_grey']:.2f}")
        for k in ("normal", "deuteranopia", "protanopia", "tritanopia",
                  "greyscale"):
            v = r[k]
            print(f"{k:13s} min dE76 {v['min_de76']:6.2f} "
                  f"({v['min_de76_pair_hours'][0]} h / "
                  f"{v['min_de76_pair_hours'][1]} h)   "
                  f"min dE2000 {v['min_de2000']:5.2f}   "
                  f"median dE76 {v['median_de76']:6.2f}")
    return r


if __name__ == "__main__":
    for _pal, _derive, _ref in (("cividis", derive_from_cividis, COLORS),
                                ("magma", derive_from_magma, COLORS_MAGMA)):
        print("=" * 74)
        validate(palette=_pal)
        hx, us, Ls = _derive()
        ok = [c.lower() for c in hx] == [c.lower() for c in _ref]
        print(f"\nderive_from_{_pal}() reproduces the delivered palette: {ok}")
        if not ok:
            print("  derived:", hx)
        print("  u =", us)
        print()

    print("=" * 74)
    fg = frame_grey_search("magma")
    b = fg["best"]
    print("FRAME GREY for magma -- 8-bit neutral sweep")
    print("  feasible band (worst dE76 >= 10 vs darkest AND contrast >= 7:1): "
          f"{len(fg['band'])} greys, {fg['band_hex'][0]} (L* "
          f"{fg['band_lstar'][0]:.2f}) .. {fg['band_hex'][1]} (L* "
          f"{fg['band_lstar'][1]:.2f})")
    print(f"  best objective in band {fg['objective_top']:.2f};  "
          f"{len(fg['tied'])} greys tie within 0.10 L*, contrast breaks it")
    print(f"  winner {b['hex']}  L* {b['lstar']:.2f}  "
          f"contrast {b['contrast']:.2f}:1  "
          f"min |dL*| to any member {b['min_dlstar_any']:.2f} "
          f"(its own gap allows {fg['objective_local_ceiling']:.2f}; "
          f"widest gap anywhere allows {fg['objective_ceiling']:.2f})")
    d = b["de_vs_darkest"]
    print("  dE76 vs darkest member: " + "  ".join(
        f"{k} {v:.2f}" for k, v in d.items()))
    print(f"  delivered FRAME_INK_GREY = {fg['delivered']}  "
          f"(matches winner: {fg['delivered'] == b['hex']})")
