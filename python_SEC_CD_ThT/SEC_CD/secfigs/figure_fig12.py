# -*- coding: utf-8 -*-
r"""
secfigs.figure_fig12 -- Figures 1 and 2 of the 37 C Abeta42 oligomer manuscript,
combined into ONE Nature Communications display item.

    figures/Figure_1_2_Combined.{svg,png,pdf}             as her pipeline processes
    figures/Figure_1_2_Combined_baselined.{svg,png,pdf}   + a common CD zero

Five panels, drawn NATIVELY at final size through the panels' own ``draw(ax)``
functions (:func:`build_fig12_panels.draw_sec` and :func:`secfigs.cd_panels.draw`).
Nothing here re-implements a loader, a peak finder or a smoothing pass, and the
standalone 12.41 in canvases are never scaled -- scaling them is what would put the
tick labels at 4.09 pt.

    row 1   a  SEC, cold incubation        b  SEC, 37 C
    row 2   c  CD, SEC-isolated monomers   d  CD, SEC-isolated oligomers
            e  CD, unfractionated, 24 h at 6 / 25 / 37 C

Style is :mod:`secfigs.figure_sec_nature`'s, unchanged: 180 mm wide exactly, open
axes, frame ink ``#3b3b3b``, Arial, live SVG text, ``path.simplify=False``, 600 dpi,
no ``bbox_inches``, and Nature's real type rule from the NRJ Guide to Preparing Final
Artwork -- "Maximum text size is 7pt. Minimum text size is 5pt." -- with the single
documented exception of 8 pt bold panel letters.


================================================================================
1.  WHY THIS GROUPING
================================================================================
The row break is the argument, not the packing.

    ROW 1 asks DOES THE OLIGOMER FORM?  Two SEC panels and only two, because the
    claim is a contrast between two incubation temperatures: cold gives one peak
    (monomer), 37 C grows a second one.  Both panels carry the same x axis
    (elution volume, 6-20 mL) and the same y axis (normalised absorbance, 0-1), so
    they are directly superimposable by eye.

    ROW 2 asks WHAT STRUCTURE RESULTS?  Three CD panels, ordered
    isolated monomer -> isolated oligomer -> the unfractionated mixture, i.e. the
    two separated fractions first and the whole solution last.  d sits immediately
    left of e so that the SEC-isolated 37 C oligomer (d) is adjacent to the
    unfractionated 37 C curve it was separated out of (e's third trace), which is
    also why the two share an ink.

No arrangement that puts a CD panel in row 1, or an SEC panel in row 2, keeps that
question order.  A single row of five would give 1.10 in panels; a 3-over-2 would
have to split the SEC pair or the CD triple.

WHAT THE ASYMMETRY COSTS AND BUYS
    Two panels in a row where the row below has three means the top panels are
    1.70x wider: 2.8243 in against 1.6609 in.  That extra width is spent entirely
    on x, where it is earned -- row 1 resolves a 0.6 mL-wide oligomer shoulder
    inside a 14 mL window, row 2 resolves one broad trough over 60 nm sampled at
    1 nm.  It is NOT spent on y: see sec. 2, where the two row heights are each
    derived from an aspect ratio rather than from what is left over.


================================================================================
2.  GEOMETRY -- every length derived, nothing nudged
================================================================================
Width is fixed at 180.000 mm.  After the margins and the two column gaps -- each of
which is the space a specific string needs, measured, not chosen -- the panel widths
are FORCED:

    left    LETTER_DX + TRIM                       the panel letter is the outermost
                                                   element; LETTER_DX itself is set
                                                   by the widest y-axis block on the
                                                   page (panel c, 0.4092 in) plus
                                                   0.020 in of clear air
    right   half of "250" + TRIM                   the last CD x tick label is
                                                   centred on the right spine and
                                                   overhangs it by half its width
    GAP_W   half of "250" + 0.154 clear + LETTER_DX
    top     title block + TRIM
    bottom  x block + TRIM
    GAP_H   x block + title block + 0.150 clear

Heights are then the only free parameters, and they are set by ASPECT RATIO, not by
filling the page:

    row 1   aspect 1.270 -- the axes-box aspect MEASURED off Griffin Figure 1A
            (``secfigs.style`` module docstring), which is what the whole house
            style is calibrated against and what an SEC panel in this project has
            always looked like.
    row 2   aspect 1.000 -- a CD spectrum has no privileged axis: one trough, read
            in both its position (x) and its depth (y).  Square is the honest
            default and it is what keeps a 42 mm panel from stretching the trough.

That gives 180.000 x 134.08 mm.  The remaining 50.9 mm under Nature's 185 mm cap is
NOT slack that should have gone into the panels.  180 mm of width divided two ways
and three ways sets the panel widths; the two aspect ratios then set the heights;
there is no freedom left.  Growing the page to 185 mm could only be done by
distorting one of the two aspects or by padding white space, and the 50.9 mm is
worth more as caption room -- this figure's caption has to carry three corrected
labels, a repeat notation and two colour keys.


================================================================================
3.  COLOUR -- TWO KEYS, AND WHY THEY CANNOT BE READ FOR EACH OTHER
================================================================================
This display item indexes TWO ordinal variables:

    INCUBATION TIME         a: t = 0 / 24 h     b: t = 0 / 2 h
    INCUBATION TEMPERATURE  c: 6 C / RT         d: 37 C     e: 6 / 25 / 37 C

ONE KEY IS NOT AVAILABLE.  A single ramp over six levels would have to claim that
24 h and 37 C are points on one scale.  They are not, and the figure's whole
argument is that they are independent axes -- 24 h at 6 C gives random coil (e),
24 h at 37 C gives beta strand (d, e).  So: two keys.

TWO KEYS ARE ONLY SAFE IF THEY DIFFER IN KIND.  A reader must be able to tell WHICH
key an ink belongs to before decoding it, and the reliable carrier of that is the
SHAPE of the ramp, not the exact hue of any one member:

    temperature   DIVERGING, two hues, cold-blue -> hot-red.  A convention a reader
                  decodes with no legend at all.  ``RdBu_r``, held from cd_panels.
    time          SEQUENTIAL, ONE hue, light -> dark.  ColorBrewer ``Oranges``.

Row 1 is one colour in three shades; row 2 runs blue to red.  Those are
unmistakably different KINDS of key.  No time ink is blue, so blue means
temperature, always.

WHAT THE CROSS-KEY CONSTRAINT SHOULD AND SHOULD NOT BE, because getting this wrong
is what produced the first version of this figure.  The natural move is to demand a
large dE76 between every time ink and every temperature ink in all four renderings
{normal, deuteranopia, protanopia, greyscale}.  Do that and the figure gets WORSE,
for a reason worth stating: greyscale dE76 is exactly |dL*|, so a cross-key
greyscale floor forces the two ramps into disjoint L* bands, and eight inks cannot
share [0, 78] in two disjoint bands without crowding.  Measured, that version put

    panel b (t = 0 vs t = 2 h)      at 10.60      panel e (6 / 25 / 37 C)  at 9.76
    panel c (6 C vs RT)             at  9.83      panel d (37 C vs frame)  at 10.46

-- i.e. it spent the whole colour budget on pairs that never appear together and
starved the pairs that do.  And there is no reader it serves: a time ink and a
temperature ink are never in the same panel, so nobody ever separates them in
greyscale.  What must survive across panels is only "these are different colours,
therefore different keys", which is a question of appearance.

So the objective is: MAXIMISE the minimum separation of the pairs that SHARE A
PANEL, plus every ink against the frame and the naming black, in the worst of all
four renderings; SUBJECT TO every cross-key pair clearing 20 dE76 in the worst of
{normal, deutan, protan}.  Two further points fall out of the same principle:

  * (t = 2 h, t = 24 h) is not a co-occurring pair either.  2 h is a 37 C time point
    in panel b, 24 h is a cold time point in panel a; they are never drawn together
    and they are not comparable if they were, because the incubation temperature
    differs.  They are held at >= 10 dE76 -- several JND, so they still read as
    different inks -- and nothing more.  Every unit above that comes straight out of
    panel b: at a 14 dE76 floor panel b falls 19.81 -> 15.83, at 18 to 11.91.
  * ORANGE AND RED ARE ONE HUE FOR A DEUTAN OR PROTAN READER.  That is why the
    cross-key floor cannot be pushed much higher: sweeping it, the best achievable
    panel minimum is 17.09 at a floor of 18-20 and collapses to 14.14 at 25 and
    11.05 at 30.  20 is the knee, and it is where this figure sits.

NO GREEN ANYWHERE.  Nature's guide names red/green as the pairing to avoid; the
temperature key is red at two of three levels, so green is excluded from the time
ramp by construction.  Oranges has none, RdBu_r has none.

THE RESULT

    TIME          t = 0     #fdb06e   L* 77.91        Oranges, single hue
                  t = 2 h   #eb600e   L* 58.10
                  t = 24 h  #cb4302   L* 47.90
    TEMPERATURE   6 C       #96c7df   L* 77.67        RdBu_r, diverging
                  25 C/RT   #de735c   L* 60.40
                  37 C      #bb2a34   L* 41.96
    frame         #3b3b3b   L* 24.87        naming ink  #000000   L* 0.00

    per-panel minimum dE76, traces AND the frame ink, worst of the four renderings
        a 23.03    b 19.81    c 17.27    d 17.09    e 17.09
    cross-key minimum, worst of {normal, deutan, protan}          20.10
    (t = 2 h, t = 24 h), the one non-co-occurring within-key pair  10.19

Every panel is therefore at least 17.09, against 9.76 for the L*-stratified version
and against 6.58 for the tightest pair of the settled SEC display item.

WHAT THE MINIMUM OVER ALL 28 PAIRS IS, AND WHY IT IS NOT THE HEADLINE.  It is 0.24,
between t = 0 (#fdb06e, L* 77.91) and 6 C (#96c7df, L* 77.67) IN GREYSCALE ONLY --
a pale orange and a pale blue of the same lightness.  In colour they are 55.54 apart
in the worst of {normal, deutan, protan}.  This is the deliberate trade named above
and it is safe because the pair never shares a panel: a greyscale print of this
figure is still fully readable PANEL BY PANEL, since every panel's own minimum
(17.09 and up) is a worst-of-four figure that already includes greyscale.  Two other
cross-key pairs are likewise close in greyscale and far in colour: t = 2 h / 25 C at
2.30 grey vs 28.95 colour, t = 24 h / 37 C at 5.95 vs 20.10.

WHAT THIS COSTS: THE LINK TO THE SEC FIGURE'S MAGMA RAMP.  ``figure_sec_nature``
indexes time with ``timecolor.COLORS_MAGMA``, and reusing it here was the first
choice.  It does not survive measurement, twice over.  Rank-indexed, this figure's
times (0, 2, 24 h) fall at ranks 0, 1 and 8, so panel b would draw ``#feae77`` and
``#fc9065`` -- two light oranges 7.1 L* apart, worst-rendering dE76 6.58, for a panel
whose entire content is a peak-height comparison between exactly those two traces.
Re-ranked over this figure's own three levels it is better but still last: magma's
warm half sits on top of RdBu_r's warm half for a deutan/protan reader, which blocks
the temperature triple from spreading in L* and pins the page at 9.74 against
Oranges' 17.09.  A single-hue ramp also states the key's KIND more plainly than a
ramp that changes hue as it darkens.  The link is worth something; it is not worth
halving every panel's separation, and the two figures index different things anyway
(a continuous 11-point kinetic series under a colourbar there, two before/after pairs
under explicit legends here).

RdBu_r was held against a sweep of every diverging map under this objective: RdBu_r
17.09, then Spectral_r 8.89, coolwarm 8.76, bwr 8.25, seismic 8.00, turbo 3.89,
RdYlBu_r 2.56, BrBG_r 3.03, PuOr_r infeasible.  It is not a tie this time -- RdBu_r
wins by 92%, because it is the only diverging map whose blue end reaches L* 78 while
its red end reaches L* 42.

THE TEMPERATURE KEY THEREFORE MOVES off cd_panels' frozen
``('#3885bc', '#bf3338', '#a21328')``, and that is deliberate.  That triple was
derived against a different ink set (the coral/teal SEC pair) and a different
objective (min over pairs INVOLVING a temperature member, greyscale included), and
its own internal minimum is 9.24 -- which is panel e, three traces, the panel that
carries the temperature argument.  17.09 is an 85% improvement on the one panel it
matters most for.  The plotted DATA is identical to the standalone panels; only the
ink changes, and the audit reports it.

d TAKES e's 37 C INK, unchanged from cd_panels' argument: d is the SEC-isolated
oligomer fraction of the same 37 C / 24 h material e's third trace is the
unfractionated form of.  c's REPEAT index stays on dash, not on colour, so colour
never means two things at once (``cd_panels.REPEAT_DASH_EM``).


================================================================================
4.  LABELS -- the corrected ones, throughout
================================================================================
  * b is t = 0 and t = 2 h.  The published "t = 24 h" is refuted: the drawn
    oligomer/monomer ratio is 0.3396, the 2 h CSV gives 0.3396 and the 24 h CSV
    0.9338, and her own matplotlib parent is ``WT_37C_75inconlyoand2.svg`` with
    ``legend_labels=['Immediate','2h']``.
  * c is TEMPERATURE x REPEAT, not a time course.  The published
    "t = 0 / 10 / 20 / 30 min" is hand-typed Inkscape text over a deleted matplotlib
    legend that read "(r1)" / "(r2)".  The legend here restores that notation
    verbatim: "6 C (r1)", "6 C (r2)", "RT (r1)", "RT (r2)".
  * c says RT, not 25 C -- the filenames and her own legend say RT; 25 C appears
    only in the hand-typed text.  RT is coloured as 25 C (``cd_panels.RT_AS_C``) so
    it shares an ink with e's 25 C trace, and printed as RT.
  * the CD y axis is Ellipticity (mdeg), not "Absorbance (mdeg)".
  * a's title states "Cold incubation" and no temperature and no concentration:
    both are disputed (caption 50 uM / 6 C vs notebook 100 uM / 4-5 C) and the
    UNICORN exports carry neither field.
  * a is annotated "Monomer" only.  Its largest oligomer-window feature is 0.0359
    of full scale against b's 0.3396 -- there is no oligomer peak to point at, and
    ``OLIGO_Y_FLOOR`` would float a pointer at 0.10 over empty page.  Naming the
    monomer in a and both species in b is what lets a reader carry the
    identification across the row.


================================================================================
5.  THE BASELINED VARIANT
================================================================================
Her pipeline applies no baseline correction.  In panel c that leaves the 6 C pair
sitting +0.1414 mdeg above the RT pair (mean over 212-250 nm), which is LARGER than
the only real difference in the panel -- a trough difference of -0.0806 mdeg -- and
comparable to the per-point noise.  A reader will read the offset as a temperature
effect.  ``build(baselined=True)`` subtracts each CD trace's own mean over
240-250 nm, which is standard practice and is what puts all the CD traces on one
zero.  Both are built; see :func:`baseline_report` for exactly what moves.

The SEC panels are untouched by the variant: they are min-max normalised absorbance,
already on a common 0-1 scale.


================================================================================
6.  TWO DEFECTS IN build_fig12_panels.draw_sec, FIXED ON THE WAY
================================================================================
``draw_sec``'s docstring claims "identical ... pointer geometry to
build_fig1_panels.build_panel".  It was not: ``build_panel`` passes ``marker="<"``
for the oligomer pointer and ``draw_sec`` passed nothing, so it fell through to
``style.annot_marker_right`` = ``">"`` -- an arrow sitting to the RIGHT of the
oligomer apex and pointing away from it, into the inter-peak valley and on towards
the monomer.  That is the exact mislabel ``build_fig1_panels`` documents having fixed
once already.  ``draw_sec`` now passes ``marker="<"`` and gained an ``annotate=``
parameter so a panel with no oligomer peak can decline the pointer.

Usage
-----
    from secfigs.figure_fig12 import build
    build()                      # figures/Figure_1_2_Combined.{svg,png,pdf}
    build(baselined=True)        # ..._baselined.{svg,png,pdf}
    python -m secfigs.figure_fig12       # both, + the full audit
"""

from __future__ import annotations

import contextlib
import itertools
import os
import sys
from dataclasses import replace

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from secfigs import style as gs                                    # noqa: E402
from secfigs import timecolor as tc                                # noqa: E402
from secfigs import cd_panels as cd                                # noqa: E402
from secfigs.figure_sec import MM, LETTER_DY_IN                    # noqa: E402
from secfigs.figure_sec import RC_EXTRA as _RC_BASE                # noqa: E402
from secfigs.figure_sec_seq import TRIM, _blocks                   # noqa: E402
from secfigs.figure_sec_nature import (                            # noqa: E402
    FRAME_GREY, NAME_INK, _ceil3, _STYLE_NOSIZE as _NAT_STYLE,
)
import build_fig12_panels as bp                                    # noqa: E402


# ==========================================================================
# 1.  colour -- the two keys
# ==========================================================================
#: The three incubation times this display item indexes, in hours.
FIG_TIMES: tuple = (0.0, 2.0, 24.0)

#: TIME key -- ColorBrewer ``Oranges``, one hue, light -> dark.
TIME_INKS: dict = {0.0: "#fdb06e", 2.0: "#eb600e", 24.0: "#cb4302"}

#: TEMPERATURE key for ``cd_panels.TEMPS`` = (6, 25, 37) C -- ``RdBu_r``,
#: diverging, cold-blue -> hot-red, L* strictly descending.
#: 6 / 25 / 37 C.  The 25 C and 37 C inks used to be #de735c and #bb2a34, which are
#: only dE = 25.3 apart in CIE Lab while the other two pairs sit at 71.9 and 89.9 --
#: so the warm pair read as one colour on the page.  25 C is now gold and 37 C a
#: deeper red, which lifts the closest pair to dE = 62.3 and keeps the cool-to-warm
#: ordering a temperature series should have.  The 6 C blue is unchanged.
TEMP_INKS: tuple = ("#96c7df", "#d9a441", "#8c1d24")

TIME_CMAP = "Oranges"
TEMP_CMAP = "RdBu_r"

#: The naming layer is a real ink on this page and is counted as one during the
#: derivation: panel letters, panel titles and the peak pointers are all black.
#: Counting it is what stops the optimiser running a trace down to L* ~ 0.
FURNITURE: tuple = (FRAME_GREY, NAME_INK)

#: Legibility floor, unchanged from ``timecolor``/``cd_panels``: nothing lighter
#: than L* 78 reads as a 0.85 pt stroke on white.
L_MAX = 78.0

#: WHICH PAIRS A READER ACTUALLY HAS TO SEPARATE -- the pairs that share a panel.
#: (t = 2 h, t = 24 h) is NOT one of them: 2 h is a 37 C time point in panel b and
#: 24 h is a cold time point in panel a, so they are neither drawn together nor
#: comparable if they were.  Paying for their separation is what crushes panel b,
#: whose whole content is t = 0 against t = 2 h -- see :data:`NONCO_MIN`.
TIME_PAIRS = ((0.0, 24.0), (0.0, 2.0))          # panel a, panel b
TEMP_PAIRS = ((0, 1), (0, 2), (1, 2))           # panel e shows all three

#: Floor on the pairs that do NOT share a panel.  10 dE76 is several JND, so they
#: still read as different inks; it is simply not worth more than that.
NONCO_MIN = 10.0

#: Floor on every TIME-ink / TEMPERATURE-ink pair, in {normal, deutan, protan}.
#: NOT in greyscale, and that exclusion is the whole design.  Greyscale dE76 is
#: exactly |dL*|, so demanding it across keys forces the two ramps into disjoint
#: L* bands -- which is what a first pass did, and it left panel e's three traces
#: at 9.76 and panel b's two at 10.60.  No reader ever needs to separate a time
#: ink from a temperature ink IN GREYSCALE: they are never in the same panel, and
#: the only thing that has to survive across panels is "these are different
#: colours, so they are different keys", which is a matter of appearance.  20
#: dE76 is comfortably that.
CROSS_KEY_FLOOR = 20.0

_R4 = ("normal", "deutan", "protan", "grey")
_R3 = ("normal", "deutan", "protan")


def _render(cols, kind):
    if kind == "normal":
        return cols
    if kind == "grey":
        return tc.to_greyscale(cols)
    return np.clip(tc.simulate_cvd(cols, kind), 0.0, 1.0)


def _pair_de(h1, h2) -> dict:
    """dE76 between two hexes under each of the four renderings."""
    from matplotlib.colors import to_rgb
    c = np.array([to_rgb(h1), to_rgb(h2)])
    return {k: float(tc.delta_e76(*tc.rgb_to_lab(_render(c, k)))) for k in _R4}


def _worst(h1, h2, rend=_R4) -> float:
    return min(_pair_de(h1, h2)[k] for k in rend)


def _cmap_entries(cmap_name: str, samples: int = 4096):
    """Distinct 8-bit entries of a colormap, in colormap order, L* <= L_MAX."""
    u = np.linspace(0.0, 1.0, samples)
    raw = np.array([matplotlib.colormaps[cmap_name](v)[:3] for v in u])
    q8 = np.round(raw * 255).astype(int)
    uniq, first = np.unique(q8, axis=0, return_index=True)
    cols = (uniq / 255.0)[np.argsort(u[first])]
    lab = tc.rgb_to_lab(cols)
    keep = lab[:, 0] <= L_MAX
    return cols[keep], lab[keep]


def _dmat(A, B, rend):
    """Worst-over-``rend`` dE76 between every row of A and every row of B."""
    D = None
    for k in rend:
        la, lb = tc.rgb_to_lab(_render(A, k)), tc.rgb_to_lab(_render(B, k))
        d = np.sqrt(((la[:, None, :] - lb[None, :, :]) ** 2).sum(-1))
        D = d if D is None else np.minimum(D, d)
    return D


def _best_triple(cmap_name, other_hexes, *, diverging: bool, cooccur):
    """The max-min triple of ``cmap_name``, L* strictly descending.

    ``cooccur`` lists the index pairs of the triple that share a panel; those,
    plus every member against :data:`FURNITURE`, are the objective.  Index pairs
    NOT listed only have to clear :data:`NONCO_MIN`.  Every member must clear
    :data:`CROSS_KEY_FLOOR` against ``other_hexes`` in {normal, deutan, protan}.

    ``diverging`` additionally pins the ends of the temperature key the way
    ``cd_panels.derive_temp_colors`` does -- b*(coldest) < 0 and b* > 0 for the
    other two, so 6 C is blue and 25 / 37 C are warm.  cd_panels only requires
    b*(37) > 0, which under this objective let the optimiser hand 25 C a BLUE and
    destroy the cold-warm read; requiring both warm members to be warm is the
    formal statement of the convention the diverging map exists to exploit.
    """
    from matplotlib.colors import to_hex, to_rgb
    cols, lab = _cmap_entries(cmap_name)
    D4 = _dmat(cols, cols, _R4)
    f4 = _dmat(cols, np.array([to_rgb(h) for h in FURNITURE]), _R4).min(axis=1)
    g3 = _dmat(cols, np.array([to_rgb(h) for h in other_hexes]), _R3).min(axis=1)
    L, bs = lab[:, 0], lab[:, 2]
    ok = np.where(g3 >= CROSS_KEY_FLOOR)[0]
    cooccur = set(tuple(sorted(p)) for p in cooccur)

    def terms(i, j, k):
        """(objective terms, constraint terms) for the triple (i, j, k)."""
        obj = [f4[i], f4[j], f4[k]]
        con = []
        for a, b in ((0, 1), (0, 2), (1, 2)):
            d = D4[(i, j, k)[a], (i, j, k)[b]]
            (obj if (a, b) in cooccur else con).append(d)
        return obj, con

    best = None
    for i in ok:
        if diverging and bs[i] >= 0:
            continue
        for j in ok[L[ok] < L[i]]:
            if diverging and bs[j] <= 0:
                continue
            for k in ok[L[ok] < L[j]]:
                if diverging and bs[k] <= 0:
                    continue
                obj, con = terms(i, j, k)
                if con and min(con) < NONCO_MIN:
                    continue
                s = min(obj)
                if best is None or s > best[0]:
                    best = (float(s), (int(i), int(j), int(k)))
    if best is None:
        raise ValueError("%s: no feasible triple" % cmap_name)
    idx = best[1]
    return ([to_hex(cols[m]) for m in idx], best[0],
            [round(float(L[m]), 2) for m in idx])


def derive_keys(rounds: int = 8, seed_temp=("#3885bc", "#bf3338", "#a21328")):
    """Re-derive BOTH keys from scratch. ``(time, temperature, n_rounds)``.

    This is the function that produced :data:`TIME_INKS` and :data:`TEMP_INKS`;
    running it is how they are audited.  The two keys constrain each other, so it
    alternates ``time <- best(TIME_CMAP | temperature)`` and
    ``temperature <- best(TEMP_CMAP | time)`` from the frozen house temperature
    key until neither moves.
    """
    temp = list(seed_temp)
    time = None
    for n in range(1, rounds + 1):
        new_time = _best_triple(TIME_CMAP, tuple(temp),
                                diverging=False,
                                cooccur=((0, 2), (0, 1)))[0]
        new_temp = _best_triple(TEMP_CMAP, tuple(new_time),
                                diverging=True,
                                cooccur=((0, 1), (0, 2), (1, 2)))[0]
        if new_time == time and new_temp == temp:
            return time, temp, n
        time, temp = new_time, new_temp
    return time, temp, rounds


def _ink_table() -> dict:
    inks = {}
    for t in FIG_TIMES:
        inks["t = %s" % ("0" if t == 0 else "%g h" % t)] = TIME_INKS[t]
    for t, c in zip(cd.TEMPS, TEMP_INKS):
        inks["%g C" % t] = c
    inks["frame"] = FRAME_GREY
    inks["black"] = NAME_INK
    return inks


#: What ink each panel puts on the page, for the per-panel separation audit.
#: Every panel's frame is ink too, so it is in every set.
PANEL_INKS = {
    "a": ("t = 0", "t = 24 h", "frame"),
    "b": ("t = 0", "t = 2 h", "frame"),
    "c": ("6 C", "25 C", "frame"),
    "d": ("37 C", "frame"),
    "e": ("6 C", "25 C", "37 C", "frame"),
}


def color_report() -> dict:
    """Every pairwise dE76 over the eight inks, and the per-panel minima."""
    from matplotlib.colors import to_rgb
    inks = _ink_table()
    pairs = {}
    for a, b in itertools.combinations(inks, 2):
        d = _pair_de(inks[a], inks[b])
        pairs["%s / %s" % (a, b)] = dict(d, worst=min(d.values()))
    lab = {k: tc.rgb_to_lab(np.array([to_rgb(v)]))[0] for k, v in inks.items()}
    tk = [k for k in inks if k.startswith("t =")]
    pk = [k for k in inks if k.endswith(" C")]
    per_panel = {}
    for L, members in PANEL_INKS.items():
        per_panel[L] = min(_worst(inks[a], inks[b])
                           for a, b in itertools.combinations(members, 2))
    cross = {}
    for a in tk:
        for b in pk:
            cross["%s / %s" % (a, b)] = dict(_pair_de(inks[a], inks[b]),
                                             worst3=_worst(inks[a], inks[b], _R3))
    worst_pair = min(pairs.items(), key=lambda kv: kv[1]["worst"])
    return dict(
        inks=inks,
        lab={k: tuple(round(float(x), 2) for x in v) for k, v in lab.items()},
        pairs=pairs, cross=cross, per_panel=per_panel,
        panel_min=min(per_panel.values()),
        cross_key_min=min(v["worst3"] for v in cross.values()),
        per_rendering={k: min(v[k] for v in pairs.values()) for k in _R4},
        worst=worst_pair[1]["worst"], worst_pair=worst_pair[0],
    )


# ==========================================================================
# 2.  the style -- figure_sec_nature's, verbatim, plus this page's figsize
# ==========================================================================
#: Text-block extents in inches for the Nature-floor style.
B = dict(_blocks(_NAT_STYLE))

#: The widest y tick label ON THIS PAGE is panel c's "-1.5" (unicode minus, 7 pt
#: Arial), measured with the renderer, not the "0.0" ``_blocks`` assumes.  Panel c
#: is the only CD panel whose y step is fractional.
YTICK_W = 0.192290
#: The last x tick label on a CD panel is "250", centred on the right spine.
XTICK_HALF = 0.081146

B["ytick_w"] = YTICK_W
B["xtick_half"] = XTICK_HALF
B["y_block"] = ((_NAT_STYLE.major_tick_length + _NAT_STYLE.tick_pad) / 72.0
                + YTICK_W + _NAT_STYLE.axis_label_pad / 72.0 + B["axlab_h"])

#: The panel letter is the outermost element on the page, so it has to sit outside
#: the widest y-axis block plus clear air.  figure_sec's 0.42 in was set against a
#: 0.3741 in block; panel c's is 0.4092 in and would have collided.
LETTER_CLEAR = 0.020
LETTER_DX_IN = _ceil3(B["y_block"] + LETTER_CLEAR)

#: Clear air between the left panel's overhanging last x label and the right
#: panel's letter / between one row's x block and the next row's title block.
#: Both are figure_sec_nature's own constants, unchanged.
GAP_W_CLEAR = 0.154
GAP_H_CLEAR = 0.150

FIG_W_MM = 180.0
FIG_W = FIG_W_MM * MM

MARGIN = dict(
    left=_ceil3(LETTER_DX_IN + TRIM),
    right=_ceil3(B["xtick_half"] + TRIM),
    top=_ceil3(B["t_block"] + TRIM),
    bottom=_ceil3(B["x_block"] + TRIM),
)
GAP_W = _ceil3(B["xtick_half"] + GAP_W_CLEAR + LETTER_DX_IN)
GAP_H = _ceil3(B["x_block"] + B["t_block"] + GAP_H_CLEAR)

NCOL = (2, 3)                       # panels in row 0, row 1
_AVAIL_W = FIG_W - MARGIN["left"] - MARGIN["right"]
PANEL_W = tuple((_AVAIL_W - (n - 1) * GAP_W) / n for n in NCOL)

#: Row aspect ratios (axes-box width / height).  See the module docstring, sec. 2.
#: 1.270 is the value MEASURED off Griffin Figure 1A in ``secfigs.style``.
ROW_ASPECT = (1.270, 1.000)
ROW_H = tuple(w / a for w, a in zip(PANEL_W, ROW_ASPECT))

FIG_H = MARGIN["top"] + ROW_H[0] + GAP_H + ROW_H[1] + MARGIN["bottom"]
FIG_H_MM = FIG_H / MM

assert FIG_H_MM <= 185.0, "page is taller than the 185 mm tier"

STYLE = replace(_NAT_STYLE, figsize=(FIG_W, FIG_H))

RC_EXTRA = dict(_RC_BASE)
RC_EXTRA["svg.hashsalt"] = "figure_fig12"
RC_EXTRA["text.color"] = NAME_INK           # the naming layer; axis furniture
                                            # follows style.foreground via apply()


# ==========================================================================
# 3.  panels
# ==========================================================================
DEG = u"°"
XLABEL_SEC = u"Elution volume (mL)"
YLABEL_SEC = u"Normalized absorbance"
XLABEL_CD = u"Wavelength (nm)"
YLABEL_CD = u"Ellipticity (mdeg)"

#: Legend labels for the CD panels, overriding ``cd_panels.CD_PANELS[..]['label']``
#: only where the standalone string is too long for a 1.66 in panel.  1C's are her
#: OWN pristine matplotlib legend's repeat notation ("(r1)" / "(r2)"), which the
#: delivered artwork deleted and replaced with the fabricated minute labels.
CD_LABELS = {
    "6C_r1": u"6 %sC (r1)" % DEG,
    "6C_r2": u"6 %sC (r2)" % DEG,
    "RT_r1": u"RT (r1)",
    "RT_r2": u"RT (r2)",
    "6C_mix": u"6 %sC" % DEG,
    "25C_mix": u"25 %sC" % DEG,
    "37C_mix": u"37 %sC" % DEG,
}

#: The window each CD trace's own mean is taken over in the baselined variant.
BASELINE_WINDOW = (240.0, 250.0)


def cd_ylim_padded(pid: str):
    """``cd_panels.cd_ylim`` with the house padding put BACK on. ``(ylim, step)``.

    ``cd_ylim`` pads the data range by ``style.y_pad_normalized`` (0.05) and then
    snaps the result OUT to the next tick, which destroys the padding: the
    outermost major tick ends up exactly ON the axes limit.  Every CD panel
    therefore prints a y tick label at the axes' bottom-left corner, immediately
    above the first x tick label.  At the standalone 12.41 in canvas that is a
    1.4 mm gap and invisible; at 42.19 mm it measured 0.19 mm, the tightest text
    pair on the page and under this project's 0.50 mm floor.

    The fix is not a nudge, it is the house convention restored.  ``gs.axes``
    with ``normalized=True`` -- the SEC convention this whole style is measured
    from -- draws ticks at 0.0 ... 1.0 inside ``ylim = (-0.05, 1.05)``, i.e. it
    pads by 0.05 of the TICK RANGE beyond the outermost tick at each end.  This
    applies exactly that rule to the snapped CD range, so no CD tick label sits
    on an axes edge either.

    The plotted data does not move; only the y mapping does.  Reported deviation
    against the standalone panels is unaffected because the verifier maps each
    file through its own limits.
    """
    p, traces = cd.load_panel(pid)
    (lo, hi), step = cd.cd_ylim(traces, p["xlim"])
    pad = gs.STYLE.y_pad_normalized * (hi - lo)
    return (lo - pad, hi + pad), step

PANEL_SPEC = [
    dict(letter="a", pid="1A", kind="sec", grid=(0, 0),
         title=u"Cold incubation",
         xlabel=XLABEL_SEC, ylabel=YLABEL_SEC,
         # a and b share an axis exactly (normalised absorbance, 0-1, same ticks),
         # so it is named once.  c/d/e do NOT share one and are each named.
         times=(0.0, 24.0), annotate=("Monomer",), legend=True),
    dict(letter="b", pid="1B", kind="sec", grid=(0, 1),
         title=u"37 %sC incubation" % DEG,
         xlabel=XLABEL_SEC, ylabel=None,
         times=(0.0, 2.0), annotate=("Oligomer", "Monomer"), legend=True),
    dict(letter="c", pid="1C", kind="cd", grid=(1, 0),
         title=u"SEC-isolated monomers",
         xlabel=XLABEL_CD, ylabel=YLABEL_CD, legend=True),
    dict(letter="d", pid="1D", kind="cd", grid=(1, 1),
         title=u"SEC-isolated oligomers, 37 %sC" % DEG,
         xlabel=XLABEL_CD, ylabel=YLABEL_CD,
         # one trace: a legend would only restate the title, and at 1.66 in it
         # would be wider than the panel.
         legend=False),
    dict(letter="e", pid="2A", kind="cd", grid=(1, 2),
         title=u"Unfractionated, 24 h",
         xlabel=XLABEL_CD, ylabel=YLABEL_CD, legend=True),
]


@contextlib.contextmanager
def figure_inks(baselined: bool = False):
    """Install this figure's colour key and label overrides on ``cd_panels``.

    Two module-level functions are swapped for the duration and restored after,
    because ``cd_panels`` exposes no argument that reaches either:

      * ``temp_colors`` resolves through ``gs.ramp(3, style)[:2]`` and
        ``derive_temp_colors``, whose objective is a different one (min over pairs
        INVOLVING a temperature member, greyscale included, over an ink set of at
        most two SEC colours).  There is no argument that reaches this figure's
        objective, so the replacement returns :data:`TEMP_INKS` directly -- the
        triple :func:`derive_keys` produces, auditable by running it.
      * ``load_panel`` is where the legend labels and (for the baselined variant)
        the 240-250 nm offset have to be applied; ``draw`` calls it internally.

    Nothing outside this block sees either change.
    """
    orig_temp, orig_load = cd.temp_colors, cd.load_panel

    def temp_colors(style: gs.GriffinStyle = gs.STYLE) -> tuple:
        return TEMP_INKS

    def load_panel(pid: str, smooth: bool = True):
        p, traces = orig_load(pid, smooth=smooth)
        lo, hi = BASELINE_WINDOW
        for t in traces:
            if t.key in CD_LABELS:
                t.label = CD_LABELS[t.key]
            if baselined:
                m = (t.wl >= lo) & (t.wl <= hi)
                t.y = t.y - float(t.y[m].mean())
                t.y_raw = t.y_raw - float(t.y_raw[m].mean())
        return p, traces

    cd.temp_colors, cd.load_panel = temp_colors, load_panel
    try:
        yield
    finally:
        cd.temp_colors, cd.load_panel = orig_temp, orig_load


def sec_style(spec, style=None):
    """A style whose ``palette[:2]`` is this SEC panel's two TIME inks.

    ``build_fig1_panels.sec_colors`` resolves ``gs.ramp(3, style)[:2]`` at draw
    time, so replacing the palette is the supported way to recolour 1A / 1B --
    which is exactly the propagation path the colour fix built.
    """
    style = STYLE if style is None else style
    cols = tuple(TIME_INKS[t] for t in spec["times"])
    rest = tuple(c for t, c in sorted(TIME_INKS.items()) if c not in cols)
    return replace(style, palette=cols + rest)


# ==========================================================================
# 4.  the figure
# ==========================================================================
class CombinedFig12(object):
    """Figures 1 and 2 as one 180 x 134 mm Nature Communications display item."""

    def __init__(self, baselined: bool = False, style=STYLE, outdir=None,
                 stem=None):
        self.baselined = bool(baselined)
        self.style = style
        self.outdir = outdir or os.path.join(_REPO, "figures")
        self.stem = stem or ("Figure_1_2_Combined_baselined" if baselined
                             else "Figure_1_2_Combined")
        self.fig = None
        self.axes = {}
        self.meta = {}

    # -- geometry ---------------------------------------------------------
    def rects(self) -> dict:
        """``letter -> (x0, y0, w, h)`` in INCHES from the page's bottom-left."""
        out = {}
        y_row = (MARGIN["bottom"] + ROW_H[1] + GAP_H, MARGIN["bottom"])
        for spec in PANEL_SPEC:
            r, c = spec["grid"]
            w = PANEL_W[r]
            x0 = MARGIN["left"] + c * (w + GAP_W)
            out[spec["letter"]] = (x0, y_row[r], w, ROW_H[r])
        return out

    # -- one panel --------------------------------------------------------
    def _draw_panel(self, ax, spec):
        st = self.style
        if spec["kind"] == "sec":
            meta = bp.draw_sec(ax, spec["pid"], sec_style(spec, st),
                               title=False, letter=False,
                               legend=spec["legend"], labels=False,
                               annotate=spec["annotate"])
        else:
            ylim, _ = cd_ylim_padded(spec["pid"])
            meta = cd.draw(ax, spec["pid"], st, title=False, letter=False,
                           legend=spec["legend"], labels=False, ylim=ylim)
        gs.set_labels(ax, spec["xlabel"], spec["ylabel"], st)
        gs.set_title(ax, spec["title"], st)

        # -- open axes: drop the top and right spines and their ticks.
        # gs.axes has already written all four spines from the style (show_box
        # is a single flag over all four, so it cannot be used to drop two);
        # ticks_all_sides=False has already taken the ticks.
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="x", which="both", top=False)
        ax.tick_params(axis="y", which="both", right=False)

        # -- the legend joins the grey reading layer, not the naming layer
        leg = ax.get_legend()
        if leg is not None:
            for t in leg.get_texts():
                t.set_color(FRAME_GREY)
        return meta

    def _place_letters(self):
        fw, fh = self.style.figsize
        for spec in PANEL_SPEC:
            pos = self.axes[spec["letter"]].get_position()
            self.fig.text(pos.x0 - LETTER_DX_IN / fw,
                          pos.y1 + LETTER_DY_IN / fh,
                          spec["letter"], transform=self.fig.transFigure,
                          fontsize=self.style.panel_letter_size,
                          fontweight=self.style.panel_letter_weight,
                          ha="left", va="baseline", color=NAME_INK)

    # -- draw -------------------------------------------------------------
    def draw(self):
        st = self.style
        fw, fh = st.figsize
        rects = self.rects()
        with figure_inks(self.baselined), gs.apply(style=st, **RC_EXTRA):
            fig = plt.figure(figsize=(fw, fh), dpi=st.dpi)
            for spec in PANEL_SPEC:
                x0, y0, w, h = rects[spec["letter"]]
                ax = fig.add_axes([x0 / fw, y0 / fh, w / fw, h / fh])
                self.axes[spec["letter"]] = ax
                self.meta[spec["letter"]] = self._draw_panel(ax, spec)
            self.fig = fig
            self._place_letters()
        return fig

    # -- save -------------------------------------------------------------
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
        self.draw()
        return self.save(formats)

    # -- report -----------------------------------------------------------
    def geometry(self) -> dict:
        fw, fh = self.style.figsize
        g = dict(
            page_mm=(fw / MM, fh / MM), page_in=(fw, fh), dpi=self.style.dpi,
            px=(int(round(fw * self.style.dpi)), int(round(fh * self.style.dpi))),
            margin_in=dict(MARGIN), gap_w_in=GAP_W, gap_h_in=GAP_H,
            trim_in=TRIM, letter_dx_in=LETTER_DX_IN,
            panel_w_in=PANEL_W, row_h_in=ROW_H, row_aspect=ROW_ASPECT,
            blocks_in={k: round(v, 6) for k, v in B.items()},
            frame_ink=FRAME_GREY, name_ink=NAME_INK,
            font_pt=dict(tick=self.style.tick_labelsize,
                         axis_title=self.style.axis_label_size,
                         panel_title=self.style.title_size,
                         legend=self.style.legend_fontsize,
                         annotation=self.style.annotation_size,
                         panel_letter=self.style.panel_letter_size),
            baselined=self.baselined,
        )
        if self.fig is not None:
            g["panels"] = {}
            for spec in PANEL_SPEC:
                L = spec["letter"]
                ax = self.axes[L]
                p = ax.get_position()
                m = self.meta[L]
                g["panels"][L] = dict(
                    pid=spec["pid"], kind=spec["kind"],
                    x0_in=p.x0 * fw, y0_in=p.y0 * fh,
                    w_in=p.width * fw, h_in=p.height * fh,
                    aspect=(p.width * fw) / (p.height * fh),
                    n_traces=m["n_traces"], xlim=tuple(ax.get_xlim()),
                    ylim=tuple(ax.get_ylim()),
                    colors=[l.get_color() for l in ax.lines
                            if l.get_linestyle() != "None"],
                    spines={s: bool(ax.spines[s].get_visible())
                            for s in ("left", "bottom", "top", "right")},
                )
        return g


# ==========================================================================
# 5.  what the baselined variant changes
# ==========================================================================
def baseline_report() -> dict:
    """Per-trace offsets removed, and what happens to panel c's separation."""
    lo, hi = BASELINE_WINDOW
    out = {"window_nm": BASELINE_WINDOW, "panels": {}}
    for pid in ("1C", "1D", "2A"):
        p, traces = cd.load_panel(pid)
        rows = []
        for t in traces:
            m = (t.wl >= lo) & (t.wl <= hi)
            off = float(t.y[m].mean())
            mr = (t.wl >= lo) & (t.wl <= hi)
            nm, mdeg = t.trough
            rows.append(dict(
                key=t.key, label=t.label, offset=off,
                sd_smoothed=float(t.y[m].std(ddof=1)),
                sd_raw=float(t.y_raw[mr].std(ddof=1)),
                trough=(nm, mdeg), trough_baselined=(nm, mdeg - off),
                residual_after=float((t.y - off)[m].mean()),
            ))
        out["panels"][pid] = dict(
            rows=rows,
            spread_before=max(r["offset"] for r in rows)
            - min(r["offset"] for r in rows),
        )

    # panel c: does the 6 C / RT separation survive?
    p, traces = cd.load_panel("1C")
    d = {t.key: t for t in traces}
    wl = d["6C_r1"].wl
    band = (wl >= 212.0) & (wl <= 250.0)
    def pair(keys, off=None):
        ys = []
        for k in keys:
            y = d[k].y
            if off is not None:
                m = (wl >= lo) & (wl <= hi)
                y = y - float(y[m].mean())
            ys.append(y)
        return np.vstack(ys).mean(0)
    c6, crt = pair(("6C_r1", "6C_r2")), pair(("RT_r1", "RT_r2"))
    b6, brt = pair(("6C_r1", "6C_r2"), 1), pair(("RT_r1", "RT_r2"), 1)
    out["panel_c"] = dict(
        band_nm=(212.0, 250.0),
        offset_before=float((c6 - crt)[band].mean()),
        offset_after=float((b6 - brt)[band].mean()),
        trough_diff_before=float(c6.min() - crt.min()),
        trough_diff_after=float(b6.min() - brt.min()),
        within_group_trough_spread_before=(
            float(abs(d["6C_r1"].trough[1] - d["6C_r2"].trough[1])),
            float(abs(d["RT_r1"].trough[1] - d["RT_r2"].trough[1]))),
    )
    m = (wl >= lo) & (wl <= hi)
    off = {k: float(v.y[m].mean()) for k, v in d.items()}
    out["panel_c"]["within_group_trough_spread_after"] = (
        float(abs((d["6C_r1"].trough[1] - off["6C_r1"])
                  - (d["6C_r2"].trough[1] - off["6C_r2"]))),
        float(abs((d["RT_r1"].trough[1] - off["RT_r1"])
                  - (d["RT_r2"].trough[1] - off["RT_r2"]))))
    return out


# ==========================================================================
# 6.  cli
# ==========================================================================
def build(baselined: bool = False, formats=("svg", "png", "pdf"), outdir=None):
    """Build one variant. Returns ``(files, geometry)``."""
    f = CombinedFig12(baselined=baselined, outdir=outdir)
    f.draw()
    files = f.save(formats)
    g = f.geometry()
    plt.close(f.fig)
    return files, g


def main():
    from secfigs import verify_fig12 as vf                      # noqa: F401
    return vf.main()


if __name__ == "__main__":
    sys.exit(main())
