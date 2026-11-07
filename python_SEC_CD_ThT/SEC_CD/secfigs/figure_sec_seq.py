# -*- coding: utf-8 -*-
r"""
secfigs.figure_sec_seq -- the 2 x 2 SEC display item with a TIME-KEYED ramp.

What this changes, and why
--------------------------
``secfigs.figure_sec.CombinedSECFigure`` colours a trace by its POSITION in its
own panel's series (``style.ramp(n)``).  The same hex therefore means different
incubation times in different panels -- ``#fe9929`` is 72 h in panel b, 48 h in
panel c and 24 h in panel d.  Panels c and d are the same 100 uM sample at 37 C
and 25 C, and that comparison is the entire reason the 2 x 2 exists, so a colour
that does not survive the trip from c to d is a readability defect.

Here colour is a pure function of t, taken from :mod:`secfigs.timecolor`
(cividis, eleven members at uniform L* spacing, rank-indexed).  ``color_for(24)``
is ``#263d6e`` in b, c and d alike, and a reader can follow one time point
across the figure.

Because colour now carries the same meaning everywhere, the four per-panel
legends (3 + 7 + 7 + 9 = 26 entries, all of them the string "N h") collapse into
ONE shared key: a discrete colourbar of eleven equal cells, every cell ticked
and labelled with the measured hours.

Two variants, one class
-----------------------
``SequentialSECFigure(normalized=True)``   -> figures/Figure_SEC_normalized.*
``SequentialSECFigure(normalized=False)``  -> figures/Figure_SEC_raw.*

Everything except the y axis is identical between them.  ``figure_sec.py`` is
imported, not copied: the panel definitions (:data:`figure_sec.PANEL_SPEC`), the
CSV loader (``build_fig3_panels.load_csv``, reached through the parent), the
pointer geometry (:data:`figure_sec.POINTER`), the apex finder
(``figure_sec._apex``) and the small-format Griffin style
(:data:`figure_sec.SMALL`) all come from there, and
``SequentialSECFigure`` subclasses ``CombinedSECFigure``.


COLOURBAR PLACEMENT -- horizontal, under the panels
===================================================
The page width is FIXED at 180.00 mm (Nature Communications double column); the
page height is not -- anything up to 230 mm is legal.  The two options cost
their space in different currencies:

    vertical, down the right edge
        bar 0.120 + tick pad 0.028 + tick labels ("72") 0.093 + label pad 0.056
        + rotated title 0.113 + clearance from the panels 0.200 + trim 0.130
        = 0.740 in of WIDTH.  The width is fixed, so that comes straight out of
        the panels: PANEL_W 2.8683 -> 2.5883 in, every panel 9.8% narrower and
        its aspect squeezed from 1.147 to 1.035.

    horizontal, under the panels
        clearance from the bottom row's x axis label 0.200 + bar 0.120 + tick
        pad 0.028 + tick labels 0.097 + label pad 0.056 + title 0.113 + trim
        0.130 = 0.744 in of HEIGHT, which the page has to spare: 159.51 mm ->
        177.09 mm, still 52.9 mm inside the 230 mm cap, and the panels keep
        their full 2.8683 in width.

Same ~0.74 in either way; only one of the two currencies is scarce.  The
horizontal bar also lands in the one band of the page that is already empty
across its full 7.09 in -- below the bottom row, where nothing lives but the x
axis label -- and spanning both columns is itself the statement that the key is
shared by all four panels.


PADDING -- every gap derived, nothing guessed
=============================================
Text blocks are computed from the style's own point sizes (:func:`_blocks`),
not measured off a render and hard-coded:

    tick label box   6.0 pt * 1.16 = 6.96 pt   (Arial ascent+descent)
    axis label box   7.0 pt * 1.16 = 8.12 pt
    title box        7.0 pt * 1.16 = 8.12 pt
    digit advance    0.556 em -> "20" is 6.67 pt at 6 pt

    X_BLOCK  space below an axes = major tick 2.40 + tick pad 1.80
             + tick label 6.96 + axis label pad 3.00 + axis label 8.12
             = 22.28 pt = 0.3094 in
    Y_BLOCK  space left of an axes = 2.40 + 1.80 + "0.0" 8.34
             + 3.00 + 8.12 = 23.66 pt = 0.3286 in
    T_BLOCK  space above an axes = title pad 3.50 + title 8.12
             = 11.62 pt = 0.1614 in

and then

    left    0.550 = panel letter offset 0.420 + TRIM 0.130.
            The letter (0.420 in out) is what sets this, not the y axis label
            (Y_BLOCK 0.329 in out) -- the letter is 0.091 in further out.
    right   0.180 = half of the last x tick label 0.046 + TRIM 0.130.
            "20" is centred on the right spine and overhangs it.
    top     0.300 = T_BLOCK 0.161 + TRIM 0.130 (rounded up 0.008).
            The title beats the panel letter, whose cap top reaches only
            0.045 + 8 pt * 0.716 = 0.125 in above the axes.
    GAP_W   0.620 = letter 0.420 + left panel's "20" overhang 0.046
            + 0.154 clear.  In the RAW variant the right column also carries a
            y axis label, but at Y_BLOCK 0.329 in it sits inside the 0.420 in
            the letter already reserves, so the gap does not have to grow.
    GAP_H   0.620 = X_BLOCK 0.309 + T_BLOCK 0.161 + 0.150 clear.
    bottom  1.052 = X_BLOCK 0.309 + CB_GAP 0.200 + CB_H 0.120
            + cbar ticks 0.124 + cbar title 0.169 + TRIM 0.130.

TRIM = 0.130 in = 3.30 mm: nothing on the page comes closer than that to the
trim edge, on any of the four sides.  (``figure_sec.py`` used left 0.52 / top
0.26, which leave 2.54 mm and 2.51 mm; both are lifted here.)


Y AXES
======
normalized
    Shared (-0.05, 1.05) on all four, majors 0.2, minors 0.1 -- the Griffin
    convention.  The axis is literally the same axis four times, so the y AXIS
    LABEL is drawn once per row, on the left column only.

raw
    The four panels do not share a scale and must not pretend to: over the
    6-20 mL window panel a reaches 20.577 mAU and panel b only 2.341 mAU, an
    8.8 : 1 ratio, so a shared axis would draw b as a 11% -tall smudge.  Each
    panel gets its own limits and its own y label.

    EXCEPT the bottom row.  c and d are the same 100 uM sample and their
    comparison is the point of that row, so they are checked for a shared
    range: c peaks at 9.739 mAU and d at 6.663, a ratio of only 1.46, and on
    the shared axis d still fills 68.4% of the data span.  That is not wasted
    space -- it is the measurement, and it is only visible if the two panels
    share a scale.  They do: both are (-0.487, 10.226) mAU, majors 2, minors 1,
    with identical tick labels as the visible signal that the scale is shared.

    Limits follow the normalised convention exactly: ``(-0.05 M, 1.05 M)``
    where M is the group maximum, so every raw panel fills 90.9% of its box,
    the same fraction a normalised panel does.  Major steps are chosen by
    :func:`_nice_step` to land ~5 majors on a 1 / 2 / 2.5 / 4 / 5 grid;
    minors are half of that.

POINTERS
    ``figure_sec`` puts the Oligomer/Monomer pointers at a fraction of their
    own peak's HEIGHT (0.478 and 0.770).  With both apexes at 1.000 that is
    fine, but in raw mAU panel a's oligomer apex is 20.577 and its monomer apex
    14.694, which drags the two labels to 9.84 and 11.31 mAU -- 6.8% of the
    axis apart, i.e. collapsed.  The fractions are therefore converted ONCE,
    from the normalised geometry, into AXES fractions
    (0.478 -> 0.4800, 0.770 -> 0.7455 of the box) and applied in every variant.
    The normalised figure is unchanged to the last digit; the raw figure puts
    its pointers at the same place on the page, 26.5% of the axis apart.

x axis
    Shared (6, 20) mL on all four, majors 2 -- unchanged.

Usage
-----
    from secfigs.figure_sec_seq import SequentialSECFigure
    SequentialSECFigure(normalized=True).build()
    SequentialSECFigure(normalized=False).build()
    python -m secfigs.figure_sec_seq          # both, plus the geometry report
"""

from __future__ import annotations

import math
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
from secfigs import timecolor as tc                   # noqa: E402
from secfigs.panel_data import PANELS                 # noqa: E402
from secfigs.figure_sec import (                      # noqa: E402
    MM, XLIM, X_MAJOR, Y_MAJOR, Y_MINOR, YLIM,
    LETTER_DX_IN, LETTER_DY_IN,
    MONO_WIN, OLIG_WIN, POINTER, POINTER_OVERRIDE,
    PANEL_SPEC, RC_EXTRA as _RC_BASE, SMALL as _SMALL_BASE,
    CombinedSECFigure, _apex, _short,
)
from build_fig3_panels import load_csv                # noqa: E402

# A panel keeps its own legend only if the shared colourbar cannot separate its
# traces -- in practice only panel a, whose three time points (2, 4, 6 h) span
# one sixth of a 0-72 h ramp.
LEGEND_MAX_ENTRIES = 3


# ==========================================================================
# 1.  text-block arithmetic -- the source of every gap below
# ==========================================================================
_PT = 1.0 / 72.0                 # inch per point
_LINE = 1.16                     # Arial ascent+descent as a multiple of em
_DIGIT = 0.556                   # Arial digit advance, em
_CAP = 0.716                     # Arial cap height, em

#: Trim-edge clearance.  Nothing on the page is closer than this to any edge.
TRIM = 0.130                     # in = 3.30 mm


def _blocks(st):
    """Derived text-block extents, in INCHES, for a given style."""
    tick_h = st.tick_labelsize * _LINE * _PT
    axlab_h = st.axis_label_size * _LINE * _PT
    title_h = st.title_size * _LINE * _PT
    b = dict(
        tick_h=tick_h, axlab_h=axlab_h, title_h=title_h,
        # widest y tick label is "0.0" (normalized) = 2 digits + a period
        ytick_w=(2 * _DIGIT + 0.278) * st.tick_labelsize * _PT,
        # last x tick label "20", centred on the right spine
        xtick_half=0.5 * (2 * _DIGIT) * st.tick_labelsize * _PT,
        letter_cap=st.panel_letter_size * _CAP * _PT,
    )
    b["x_block"] = ((st.major_tick_length + st.tick_pad) * _PT + tick_h
                    + st.axis_label_pad * _PT + axlab_h)
    b["y_block"] = ((st.major_tick_length + st.tick_pad) * _PT + b["ytick_w"]
                    + st.axis_label_pad * _PT + axlab_h)
    b["t_block"] = st.title_pad * _PT + title_h
    return b


B = _blocks(_SMALL_BASE)

# -- colourbar block --------------------------------------------------------
CB_H = 0.120                     # in, bar thickness
CB_GAP = 0.200                   # in, bar top -> bottom row's x axis label
CB_TICK_PAD = 2.0                # pt
# 4.0 pt, not 2.0: at 2.0 the rendered gap between the "8" tick label and the
# "Incubation time (h)" title measured 0.71 mm, the tightest text pair on the
# page.  4.0 pt puts it at 1.41 mm and costs 0.028 in of page height.
CB_LABEL_PAD = 4.0               # pt
CB_TICK_BLOCK = CB_TICK_PAD * _PT + B["tick_h"]          # 0.1244 in
CB_TITLE_BLOCK = CB_LABEL_PAD * _PT + B["axlab_h"]       # 0.1406 in
CB_BLOCK = CB_GAP + CB_H + CB_TICK_BLOCK + CB_TITLE_BLOCK + TRIM


# ==========================================================================
# 2.  page geometry
# ==========================================================================
FIG_W_MM = 180.0
FIG_W = FIG_W_MM * MM                                    # 7.086614 in

MARGIN = dict(
    left=round(LETTER_DX_IN + TRIM, 3),                  # 0.550
    right=round(B["xtick_half"] + TRIM, 3),              # 0.176 -> 0.180
    top=round(B["t_block"] + TRIM, 3),                   # 0.291 -> 0.300
    bottom=round(B["x_block"] + CB_BLOCK, 3),            # 1.052
)
MARGIN["right"] = 0.180
MARGIN["top"] = 0.300

GAP_W = round(LETTER_DX_IN + B["xtick_half"] + 0.154, 2)     # 0.62
GAP_H = round(B["x_block"] + B["t_block"] + 0.150, 2)        # 0.62
ROW_H = 2.50                     # in -- IDENTICAL to figure_sec, so the panel
                                 # boxes of the three figures are comparable

NCOL = NROW = 2
_AVAIL_W = FIG_W - MARGIN["left"] - MARGIN["right"]
PANEL_W = (_AVAIL_W - (NCOL - 1) * GAP_W) / NCOL             # 2.868307 in
WSPACE = GAP_W / PANEL_W
HSPACE = GAP_H / ROW_H

FIG_H = MARGIN["top"] + 2 * ROW_H + GAP_H + MARGIN["bottom"]
FIG_H_MM = FIG_H / MM                                        # 177.09 mm

#: Colourbar rect, in inches from the page's bottom-left corner.  It spans the
#: full axes extent -- left panel's left spine to right panel's right spine.
CB_RECT_IN = (MARGIN["left"],
              TRIM + CB_TITLE_BLOCK + CB_TICK_BLOCK,
              NCOL * PANEL_W + GAP_W,
              CB_H)

SMALL = replace(_SMALL_BASE, figsize=(FIG_W, FIG_H))

#: Same rc as figure_sec (live SVG/PDF text, no path simplification), with its
#: own hash salt so the two figures' SVG element ids cannot collide.
RC_EXTRA = dict(_RC_BASE)
RC_EXTRA["svg.hashsalt"] = "figure_sec_seq"


# ==========================================================================
# 3.  y axes
# ==========================================================================
YLABEL_NORM = u"Normalized absorbance"
YLABEL_RAW = u"Absorbance (mAU)"

#: Which panels share a raw y range.  a and b are alone (250 uM and 25 uM have
#: nothing to compare); c and d are the same 100 uM sample at two temperatures
#: and the bottom row exists to compare them, so they are one group.
RAW_GROUPS = (("a",), ("b",), ("c", "d"))

_NICE = (1.0, 2.0, 2.5, 4.0, 5.0, 10.0)


def _nice_step(m, target=5.0):
    """Major tick step giving ~`target` majors over [0, m], on a nice grid."""
    if m <= 0:
        raise ValueError("m must be positive")
    k = 10.0 ** math.floor(math.log10(m / target))
    cands = [c * k for c in _NICE]
    return min(cands, key=lambda s: abs(m / s - target))


def raw_ylim(m):
    """``(-0.05 m, 1.05 m)`` -- the normalised convention, scaled to the data."""
    return (-0.05 * m, 1.05 * m)


# -- pointer y, converted once from the normalised geometry -----------------
#: ``figure_sec.POINTER['Oligomer']['y_frac']`` is a fraction of the peak's own
#: HEIGHT.  Both apexes are 1.000 when normalised, and ylim is (-0.05, 1.05),
#: so the pointer's position as a fraction of the AXES BOX is fixed:
_span = YLIM[1] - YLIM[0]
POINTER_Y_AXFRAC = {k: (v["y_frac"] * 1.0 - YLIM[0]) / _span
                    for k, v in POINTER.items()}     # 0.480000, 0.745455


# ==========================================================================
# 4.  the figure
# ==========================================================================
class SequentialSECFigure(CombinedSECFigure):
    """The 2 x 2 SEC figure, coloured by incubation time, with one shared key.

    Parameters
    ----------
    normalized
        ``True``  -> min-max normalised traces, shared (-0.05, 1.05) y axis.
        ``False`` -> raw mAU, per-panel y limits (c and d shared).

    Subclassing
    -----------
    Everything a variant is likely to want to change is reachable without
    copying :meth:`load`, :meth:`_draw_panel` or :meth:`draw`:

    ``PALETTE``            class attribute, a key of :data:`timecolor.PALETTES`
    :meth:`_layout`        every page number in inches, as one dict
    :meth:`_rc`            the rcParams overlay
    :meth:`_trace_kwargs`  extra ``ax.plot`` kwargs, per trace
    :meth:`_style_axes`    runs right after ``gs.axes`` on each panel
    :meth:`_style_legend`  runs right after ``gs.legend``, when there is one
    :meth:`_after_panel`   runs last, once the panel is otherwise complete

    ``secfigs.figure_sec_nature.NatureSECFigure`` is the worked example.  The
    defaults below reproduce the shipped ``Figure_SEC_normalized`` /
    ``Figure_SEC_raw`` byte for byte.
    """

    #: Which :data:`timecolor.PALETTES` entry keys colour to time.
    PALETTE = "cividis"

    def __init__(self, normalized=True, style=SMALL, outdir=None, stem=None,
                 palette=None):
        if stem is None:
            stem = "Figure_SEC_normalized" if normalized else "Figure_SEC_raw"
        super().__init__(style=style, outdir=outdir, stem=stem)
        self.normalized = bool(normalized)
        self.palette = self.PALETTE if palette is None else palette
        self.cax = None
        self.cbar = None
        self.ylims = {}
        self.ysteps = {}

    # -- geometry / rc, as overridable single sources of truth -------------
    def _layout(self):
        """Every page number :meth:`draw` needs, in INCHES.

        The module-level constants are the default; a subclass that re-derives
        them for a different type size overrides this one method and inherits
        the drawing code unchanged.
        """
        return dict(margin=dict(MARGIN), gap_w=GAP_W, gap_h=GAP_H,
                    row_h=ROW_H, panel_w=PANEL_W, wspace=WSPACE,
                    hspace=HSPACE, cb_rect_in=CB_RECT_IN,
                    cb_h=CB_H, cb_gap=CB_GAP,
                    cb_tick_pad=CB_TICK_PAD, cb_label_pad=CB_LABEL_PAD,
                    blocks=dict(B), trim=TRIM)

    def _rc(self):
        """rcParams overlay handed to ``gs.apply``."""
        return dict(RC_EXTRA)

    # -- per-panel hooks, all no-ops by default ----------------------------
    def _trace_kwargs(self, spec, i, n):
        """Extra ``ax.plot`` kwargs for trace ``i`` of ``n`` in ``spec``."""
        return {}

    def _style_axes(self, ax, spec):
        """Runs immediately after ``gs.axes`` has configured ``ax``."""

    def _style_legend(self, leg, ax, spec):
        """Runs immediately after ``gs.legend``, only on panels that get one."""

    def _after_panel(self, ax, spec):
        """Runs once the panel is otherwise finished."""

    # -- data ------------------------------------------------------------
    def load(self):
        """Read every CSV through the documented 4-step contract.

        Colours come from :func:`secfigs.timecolor.color_for`, keyed on the
        incubation time parsed out of ``PANELS[pid]['labels']`` -- NOT on the
        trace's position in the panel.
        """
        for spec in PANEL_SPEC:
            P = PANELS[spec["pid"]]
            traces = [load_csv(f, normalize=self.normalized) for f in P["csvs"]]
            labels = [_short(l) for l in P["labels"]]
            hours = [float(l.split()[0]) for l in labels]
            colors = tc.colors_for(hours, self.palette)
            assert len(set(colors)) == len(traces), "colours must be distinct"
            self.data[spec["letter"]] = dict(
                traces=traces, labels=labels, colors=colors, hours=hours,
                ranks=[tc.rank_for(h) for h in hours],
                csvs=list(P["csvs"]), raw_labels=list(P["labels"]),
                orig_normalized=P["normalized"])
        self._resolve_y()
        return self

    def _resolve_y(self):
        """Fill ``self.ylims`` / ``self.ysteps``, one entry per panel letter."""
        if self.normalized:
            for spec in PANEL_SPEC:
                self.ylims[spec["letter"]] = YLIM
                self.ysteps[spec["letter"]] = (Y_MAJOR, Y_MINOR)
            return
        for group in RAW_GROUPS:
            m = max(float(A[(V >= XLIM[0]) & (V <= XLIM[1])].max())
                    for L in group for V, A in self.data[L]["traces"])
            lim = raw_ylim(m)
            step = _nice_step(m)
            for L in group:
                self.ylims[L] = lim
                self.ysteps[L] = (step, step / 2.0)

    def _ylabel(self, letter):
        if self.normalized:
            # one per ROW: the axis is the same axis four times, so the left
            # column speaks for the right.
            return YLABEL_NORM if letter in ("a", "c") else None
        return YLABEL_RAW               # every panel has its own scale

    # -- one panel -------------------------------------------------------
    def _draw_panel(self, ax, spec):
        d = self.data[spec["letter"]]
        st = self.style
        ylim = self.ylims[spec["letter"]]
        ymaj, ymin = self.ysteps[spec["letter"]]

        n = len(d["traces"])
        for i, ((V, A), c, lab) in enumerate(
                zip(d["traces"], d["colors"], d["labels"])):
            ax.plot(V, A, color=c, lw=st.trace_linewidth, label=lab,
                    solid_capstyle="round", solid_joinstyle="round",
                    **self._trace_kwargs(spec, i, n))

        gs.axes(ax, style=st, normalized=self.normalized, x_major=X_MAJOR,
                y_major=ymaj, y_minor=ymin, x_range=XLIM, y_range=ylim)
        self._style_axes(ax, spec)
        gs.set_labels(ax, u"Elution volume (mL)",
                      self._ylabel(spec["letter"]), style=st)
        gs.set_title(ax, spec["title"], style=st)

        # The colourbar replaces the legend on every panel EXCEPT a.
        #
        # a holds only t = 2, 4 and 6 h.  On a ramp spanning 0-72 h those three
        # land inside one sixth of the scale: #beaf6b / #928d78 / #7e7c78, min
        # pairwise dE 11.58 and 6.45 in greyscale, so 4 h and 6 h merge into one
        # line at 0.85 pt.  That is a real regression against the categorical
        # build, where the same three traces were dE 34.48 apart -- and it is
        # unavoidable for any ramp keyed to absolute time, because a short early
        # series simply does not span the scale.  Panel a is also the one panel
        # that never had a cross-panel colour collision to fix, so it pays the
        # whole cost of the change and gains nothing.  Giving it back its own
        # 3-entry key costs almost nothing (it has the emptiest upper right of
        # the four) and restores the discrimination the colourbar cannot.
        if len(d["labels"]) <= LEGEND_MAX_ENTRIES:
            self._style_legend(gs.legend(ax, style=st), ax, spec)

        if spec["annotate"]:
            span = XLIM[1] - XLIM[0]
            for text, win in ((u"Oligomer", OLIG_WIN), (u"Monomer", MONO_WIN)):
                q = dict(POINTER[text])
                q.update(POINTER_OVERRIDE.get((spec["letter"], text), {}))
                pv, ph = _apex(d["traces"], win)
                x = pv + q["dx_frac"] * span
                # axes fraction, not a fraction of this peak's height -- see
                # the module docstring.
                y = ylim[0] + POINTER_Y_AXFRAC[text] * (ylim[1] - ylim[0])
                gs.annotate_peak(ax, x, y, text, side=q["side"], style=st,
                                 marker=q["marker"])
                self.pointers[(spec["letter"], text)] = (x, y, pv, ph)

        self._after_panel(ax, spec)

    # -- the shared key ---------------------------------------------------
    def _add_colorbar(self):
        fw, fh = self.style.figsize
        L = self._layout()
        x0, y0, w, h = L["cb_rect_in"]
        self.cax = self.fig.add_axes([x0 / fw, y0 / fh, w / fw, h / fh])
        self.cbar = tc.add_colorbar(
            self.fig, cax=self.cax, orientation="horizontal",
            labelsize=self.style.tick_labelsize,
            titlesize=self.style.axis_label_size,
            linewidth=self.style.spine_linewidth,
            labelpad=L["cb_label_pad"], tickpad=L["cb_tick_pad"],
            palette=self.palette)
        return self.cbar

    # -- draw ------------------------------------------------------------
    def draw(self):
        if not self.data:
            self.load()
        st = self.style
        fw, fh = st.figsize
        L = self._layout()
        M = L["margin"]

        with gs.apply(style=st, **self._rc()):
            fig = plt.figure(figsize=(fw, fh), dpi=st.dpi)
            grid = GridSpec(
                NROW, NCOL, figure=fig,
                left=M["left"] / fw,
                right=1.0 - M["right"] / fw,
                top=1.0 - M["top"] / fh,
                bottom=M["bottom"] / fh,
                wspace=L["wspace"], hspace=L["hspace"],
            )
            for spec in PANEL_SPEC:
                r, c = spec["grid"]
                ax = fig.add_subplot(grid[r, c])
                self.axes[spec["letter"]] = ax
                self._draw_panel(ax, spec)

            self.fig = fig
            self._place_letters()
            self._add_colorbar()
        return fig

    # -- geometry report --------------------------------------------------
    def geometry(self):
        g = super().geometry()
        fw, fh = self.style.figsize
        L = self._layout()
        g["normalized"] = self.normalized
        g["palette"] = self.palette
        g["margin_in"] = dict(L["margin"])   # parent reports figure_sec's
        g["gap_in"] = dict(col=L["gap_w"], row=L["gap_h"], cbar=L["cb_gap"])
        g["trim_in"] = L["trim"]
        g["blocks_in"] = {k: round(v, 4) for k, v in L["blocks"].items()}
        p = self.cax.get_position()
        g["colorbar"] = dict(x0_in=p.x0 * fw, y0_in=p.y0 * fh,
                             w_in=p.width * fw, h_in=p.height * fh,
                             cell_in=(p.width * fw) / len(tc.TIMES))
        for spec in PANEL_SPEC:
            L = spec["letter"]
            g["panels"][L]["ylim"] = self.ylims[L]
            g["panels"][L]["y_major"] = self.ysteps[L][0]
            g["panels"][L]["hours"] = list(self.data[L]["hours"])
            g["panels"][L]["colors"] = list(self.data[L]["colors"])
            g["panels"][L]["ylabel"] = self._ylabel(L)
            g["panels"][L]["fill_frac"] = float(
                max(A[(V >= XLIM[0]) & (V <= XLIM[1])].max()
                    for V, A in self.data[L]["traces"]) / self.ylims[L][1])
        return g


# ==========================================================================
# 5.  cli
# ==========================================================================
def build_both(formats=("svg", "png", "pdf")):
    """Build and save both variants.  Returns ``{'normalized': ..., 'raw': ...}``."""
    out = {}
    for norm, key in ((True, "normalized"), (False, "raw")):
        f = SequentialSECFigure(normalized=norm)
        f.load()
        f.draw()
        files = f.save(formats)
        g = f.geometry()
        out[key] = dict(files=files, geometry=g)
        plt.close(f.fig)
    return out


def main():
    res = build_both()
    for key, r in res.items():
        g = r["geometry"]
        print("\n=== %s ===" % key)
        print("page  %.3f x %.3f mm  @ %d dpi  = %d x %d px"
              % (g["page_mm"][0], g["page_mm"][1], g["dpi"],
                 g["px"][0], g["px"][1]))
        print("margins in: %s" % g["margin_in"])
        print("gaps in:    %s   trim %.3f" % (g["gap_in"], g["trim_in"]))
        for k, v in g["panels"].items():
            print("  %s  %.4f x %.4f in  at (%.4f, %.4f)  n=%d  "
                  "ylim (%.4f, %.4f) maj %.4g  fill %.1f%%  ylab %s"
                  % (k, v["w_in"], v["h_in"], v["x0_in"], v["y0_in"],
                     v["n_traces"], v["ylim"][0], v["ylim"][1], v["y_major"],
                     100 * v["fill_frac"], v["ylabel"]))
        c = g["colorbar"]
        print("  cbar %.4f x %.4f in at (%.4f, %.4f)  cell %.4f in"
              % (c["w_in"], c["h_in"], c["x0_in"], c["y0_in"], c["cell_in"]))
        for k, v in r["files"].items():
            print("  %-4s %s" % (k, v))
    return res


if __name__ == "__main__":
    main()
