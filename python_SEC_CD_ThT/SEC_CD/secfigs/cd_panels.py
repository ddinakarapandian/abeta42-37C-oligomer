# -*- coding: utf-8 -*-
r"""
secfigs.cd_panels -- the three CIRCULAR DICHROISM panels of Figures 1 and 2 of
the 37 C Abeta42 oligomer manuscript (doc 10): 1C, 1D and 2A.

    1C   SEC-isolated WT monomers, 4 traces = 2 incubation temperatures x 2 repeats
    1D   SEC-isolated WT oligomers after 24 h at 37 C, 1 trace
    2A   unfractionated WT solutions after 24 h at 6 / 25 / 37 C, 3 traces

Everything here is loading + drawing.  Every size, weight, stroke width and ink is
taken from :mod:`secfigs.style`; nothing typographic is hardcoded, so the house-style
revision that is in flight (7 pt floor, spines, magma) propagates without editing
this file.  The two things this module DOES own are the CD axis conventions (a
wavelength axis is not an elution-volume axis) and the TEMPERATURE colour key --
both are derived below, with the measurements that chose them.


================================================================================
1.  PROVENANCE  (established by trace; see PROVENANCE_CD_FIGURES_1-2.md)
================================================================================
Raw data are headerless Chirascan CSV exports, comma separated, ``wavelength,value``,
with blank-line-delimited channel blocks.  **Block 0 is CD ellipticity in mdeg and is
the only block ever plotted.**

    1C   ...\Data\AbetaWT\052025data\   7-channel export, 91 pts, 280 -> 190 nm @ 1 nm
         WTmon6C052225.csv      6 C, repeat 1     single scan
         WTmonr26C052225.csv    6 C, repeat 2     single scan
         WTRTmonr1052225.csv    RT,  repeat 1     single scan
         WRrTmonr2052225.csv    RT,  repeat 2     single scan

    1D   ...\Data\AbetaWT\100325\       1-channel export, 61 pts, 250 -> 190 nm @ 1 nm
         37c oligo avg-autosub.csv      SEC-isolated oligomer fraction, 37 C / 24 h

    2A   ...\Data\AbetaWT\100325\       1-channel export, 61 pts, 250 -> 190 nm @ 1 nm
         6c mix avg-autosub.csv         unfractionated, 6 C / 24 h
         25c mix avg-autosub.csv        unfractionated, 25 C / 24 h
         37c mix avg-autosub.csv        unfractionated, 37 C / 24 h

PROCESSING CONTRACT, reimplemented here so that nothing under ``/path/to/Dropbox``
is imported (importing her module would drop a ``__pycache__`` into the read-only tree):

    1. read block 0 of the CSV -> (wavelength, mdeg)
    2. sort ascending by wavelength
    3. scipy.signal.savgol_filter(window_length=11, polyorder=2)
    4. nothing else.  No baseline subtraction, no normalisation, NO MRE conversion.

Step 3 is not optional in the source: ``load_cd_spectrum`` in both
``CD_Analysis_052025.ipynb`` (cell 4) and ``CD_Analysis_100325.ipynb`` (cell 4)
takes ``smooth=`` and **every call site passes ``smooth=True``**.  Reproducing her
published curves therefore requires it, and :data:`SAVGOL` records it as a fact of
the pipeline rather than a choice made here.  ``smooth=False`` is available on
:func:`load_trace` for auditing.

Verified: this loader reproduces every trough decoded from her delivered vector
artwork to <= 0.001 mdeg and 0.0 nm.  See :func:`verify`.

For 1C the smoothing runs over the FULL 91-point 190-280 nm array and the panel is
then cropped to 190-250 nm by ``xlim`` -- which is the order her artwork was made in
(matplotlib plotted 189-280, Inkscape deleted the 30 nodes above 250 nm by hand).
Smoothing after cropping would move the 250 nm end of those four curves.


================================================================================
2.  WHAT THE LABELS SAY, AND WHY THEY DIFFER FROM THE PUBLISHED ONES
================================================================================
Labels are derived from the FILENAMES and the notebook's own legend strings, never
from the hand-typed Inkscape text.  Three published labels do not survive that test;
each is carried in the panel's ``flags`` so the discrepancy travels with the figure.

FIG 1C IS NOT A TIME COURSE.  Doc 10's caption and the artwork label the four traces
``t = 0 / 10 / 20 / 30 min``.  That is fabricated.  The four files are two
temperatures x two repeats -- the filenames carry the repeat index
(``WTmon6C`` / ``WTmon_r2_6C``, ``WTRTmon_r1_`` / ``WRrTmon_r2_``), and the pristine
matplotlib export (``052025data\selected_overlay_WT052025.svg``, zero Inkscape
attributes) still carries the original legend ``6C WT monomers (r1)`` /
``(r2)`` / ``RT WT monomers (r1)`` / ``(r2)``.  The Inkscape-edited copy of the same
figure is that file with the matplotlib legend DELETED and the minute labels typed
in.  No archived record supports a 10-minute spacing, and the only temperature
channel in those CSVs (channel 6) reads 22.43-22.73 C for all four, i.e. they were
all measured at ambient after incubation -- so the panel cannot be a time course AT
temperature even in principle.  These panels therefore label by TEMPERATURE AND
REPEAT, which is what the files establish.

1C SAYS "RT", NOT "25 C".  The filenames say ``RT`` and so does her own notebook
legend.  ``25 C`` appears only in the hand-typed Inkscape text and in Doc 10.  The
label here is ``RT``; :data:`RT_AS_C` records that it is treated as 25 C *for
colouring only*, so 1C's RT traces and 2A's 25 C trace share an ink.

THE Y AXIS IS ELLIPTICITY, NOT ABSORBANCE.  Her ``plt.ylabel`` reads
``Absorbance (mdeg)`` in all three panels.  mdeg is a unit of ellipticity; absorbance
is dimensionless.  The label here is :data:`Y_LABEL`.

Two further defects are recorded in ``flags`` but NOT silently corrected, because
correcting them needs a number nobody archived:

  * no MRE conversion is applied anywhere, and 1C's amplitudes are 9.6x smaller than
    2A's for the same nominal species at the same nominal concentration.  The two
    panels' y axes are therefore not comparable.  Fixing it needs the cuvette
    pathlength and measured concentration for the 2025-05-22 and 2025-10-03 sessions.
  * 2A's three curves are not on a common zero: residual 240-250 nm offsets are
    +0.782 / +0.615 / +0.174 mdeg, a spread of 1.21 mdeg.


================================================================================
3.  COLOUR -- A TEMPERATURE KEY, DERIVED, NOT BORROWED
================================================================================
:mod:`secfigs.timecolor` maps INCUBATION TIME to colour and is deliberately not
reused here.  Its ``rank_for`` is a lookup over the eleven measured *hours* and
raises on anything else, so feeding it 6 / 25 / 37 is not even possible -- but the
real reason is semantic: reusing the time ramp for a temperature axis would put the
same ink on two different meanings, which is exactly the defect this project already
fixed once for time.

Is temperature ordinal here?  Yes.  6 < 25 < 37 C, and the physical claim the figure
makes IS the ordering -- more heat, more beta-strand.  So an ordered key, not a
categorical set, even at three levels.  1D gets ``color_for_temp(37)`` rather than a
colour of its own: its sample is the SEC-isolated oligomer fraction of the same 37 C
/ 24 h material 2A's 37 C curve is the unfractionated form of, and drawing them in
one ink is a true statement.

WHICH ramp.  Temperature, unlike time, has a colour convention a reader decodes with
no legend at all: cold is blue, hot is red.  That convention is DIVERGING, where the
time ramp is sequential -- so the two keys are unconfusable by construction and not
merely by hue accident.

The members were then chosen by measurement, reusing timecolor's own colour science
(``rgb_to_lab``, ``delta_e76``, ``simulate_cvd``, ``to_greyscale``).  See
:func:`derive_temp_colors`, which is the function that produced the frozen values and
is how they are audited.

    constraints  every member L* <= 78            (reads as a thin stroke on white)
                 L* strictly decreasing 6 -> 25 -> 37, step >= 8
                     greyscale dE76 is exactly |dL*|, so this is what makes the
                     ORDER survive a greyscale print
                 b*(6 C) < 0 and b*(37 C) > 0     (the cold -> warm read)
    objective    maximise the minimum dE76 of any pair INVOLVING a temperature
                 member, over the whole Figure-1 ink set, worst of
                 {normal, deuteranopia, protanopia, greyscale}

The ink set matters and is the reason the first pass was thrown away.  Figure 1 puts
SEC panels 1A/1B on the same page as CD panels 1C/1D, so the inks a reader must
separate are the two SEC trace colours ``style.ramp(3, style)[:2]`` PLUS the three
temperature members PLUS the frame.  Closed against the frame alone, the optimiser
returns ``BrBG_r`` and hands 6 C a teal ``#81cec2`` that collides with the SEC teal
``#66C2A4`` -- the exact defect being guarded against.  With the SEC inks in the set:

    cmap          WORST   hex 6 / 25 / 37              L* 6/25/37    normal/deutan/protan/grey
    RdBu_r        13.70   #2f79b5 #a51429 #6d0220      49.0/35.3/21.6   22.96 21.36 14.33  9.51
    BrBG_r        12.89   #21847c #804a09 #543005      49.8/36.9/23.7   18.75 18.58 16.38  9.51
    turbo         12.31   #476ee6 #af1801 #7a0403      49.7/37.4/24.4   23.09 18.33 16.77  9.51
    seismic       12.12   #6161ff #b60000 #800000      49.9/37.7/25.5   23.40 17.47 13.75  9.51
    Spectral_r     9.70   #3585bb #c1274a #9e0142      53.1/43.2/33.3
    PuOr_r         9.61   #8477af #a44f07 #813c08      53.2/43.5/33.8
    RdYlBu_r       9.23   #4f81ba #c82227 #a50026      52.8/43.6/34.3
    coolwarm       8.01   #5d7ce6 #c83836 #b40426      54.5/46.1/37.7
    bwr            4.60   #b8b8ff #ff7c7c #ff3e3e      77.1/67.4/57.1
    RdGy_r, PiYG_r          infeasible -- no cold->warm triple clears the L* ladder
    cividis, magma, viridis, plasma, inferno   infeasible -- sequential, no cold end

``RdBu_r`` wins outright at 13.70 and is also the canonical scientific cold-warm map,
so the measurement and the convention agree.  Note the ``grey`` column is 9.51 for
all four leaders: that is the SEC coral/teal pair's own greyscale separation, a
pre-existing property of Figure 1 that no choice made here can improve.

THE FRAME INK IS NOT FIXED, so neither is the key.  The house revision may move
``style.foreground`` from ``#000000`` to timecolor's ``FRAME_INK_GREY`` ``#3b3b3b``,
and that changes the answer: ``#6d0220`` sits at L* 21.6, only 3.29 from ``#3b3b3b``
at L* 24.9, so under a grey frame the 37 C trace would nearly vanish into the axes in
a greyscale print.  Re-run under that ink, RdBu_r returns
``#3885bc / #bf3338 / #a21328`` at worst 9.24 -- second by 1.4% to seismic's 9.37,
which is a tie, and holding one colormap family across both regimes is worth far
more than 1.4% of dE.  Both answers are frozen in :data:`TEMP_COLORS_BY_FRAME`, and
:func:`temp_colors` resolves through ``style.foreground`` at call time, so the panels
follow the style revision without being edited.  An unrecognised frame ink is
derived on demand and cached.
"""

from __future__ import annotations

import functools
import itertools
import os
from dataclasses import dataclass, field

import numpy as np
from scipy.signal import savgol_filter

from . import style as gs
from . import timecolor as tc

__all__ = [
    "CD_PANELS", "CDTrace", "SAVGOL", "RT_AS_C", "TEMPS", "Y_LABEL", "X_LABEL",
    "TEMP_COLORS_BY_FRAME", "HOUSE_SEC_INKS", "ARTWORK_REFERENCE",
    "load_cd_csv", "load_trace", "load_panel",
    "temp_colors", "color_for_temp", "derive_temp_colors",
    "cd_axes", "draw", "trough", "verify",
]


# ==========================================================================
# 1.  Raw data -- READ ONLY.  Nothing in this tree is ever written to.
# ==========================================================================
_DROPBOX = r"/path/to/Dropbox"
_AB = os.path.join(_DROPBOX, "For", "ForTarunya", "MyPapers",
                   "Abeta42WT37COligomers", "Data", "AbetaWT")
ROOT_052025 = os.path.join(_AB, "052025data")
ROOT_100325 = os.path.join(_AB, "100325")

#: Savitzky-Golay parameters, as hardcoded in her ``load_cd_spectrum``.
SAVGOL = dict(window_length=11, polyorder=2)

#: 1C's filenames say ``RT``.  It is treated as 25 C FOR COLOURING ONLY, so that
#: 1C's RT traces and 2A's 25 C trace share an ink.  The printed label stays ``RT``.
RT_AS_C: float = 25.0

#: The three incubation temperatures the key is defined over, in C. Ordered.
TEMPS: tuple = (6.0, 25.0, 37.0)

#: The quantity actually plotted.  Her ``plt.ylabel`` says "Absorbance (mdeg)";
#: mdeg is a unit of ellipticity and absorbance is dimensionless.
Y_LABEL: str = "Ellipticity (mdeg)"
X_LABEL: str = "Wavelength (nm)"


def load_cd_csv(path: str):
    """Block 0 of a Chirascan CSV export -> ``(wavelength_asc, mdeg)``.

    The export is headerless ``wavelength,value`` with blank-line-delimited channel
    blocks.  Block 0 is CD ellipticity in mdeg; every later block (HT volts,
    absorbance, error, temperature, ...) is ignored, exactly as her loader does --
    the ``.xlsx`` twin she actually reads contains block 0 and nothing else, and
    equals it to ``max|diff| = 0.0``.
    """
    text = open(path, "rb").read().decode("utf-8-sig", errors="strict")
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    wl, y = [], []
    for line in text.split("\n"):
        line = line.strip()
        if not line:
            if wl:
                break            # end of block 0
            continue
        a, b = line.split(",")[:2]
        wl.append(float(a))
        y.append(float(b))
    wl = np.asarray(wl, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(wl)       # her loader sorts ascending
    return wl[order], y[order]


@dataclass
class CDTrace:
    """One CD spectrum, loaded and processed."""
    key: str                     # stable id, e.g. "6C_r1"
    path: str
    label: str                   # what the legend prints
    temp_c: float                # incubation temperature used for COLOUR
    temp_label: str              # what the temperature is called in print
    repeat: int | None           # 1, 2 or None
    fraction: str                # "isolated monomer" / "isolated oligomer" / "unfractionated"
    wl: np.ndarray = field(repr=False)
    y: np.ndarray = field(repr=False)          # plotted (smoothed)
    y_raw: np.ndarray = field(repr=False)      # block 0, unsmoothed

    @property
    def trough(self):
        """``(wavelength_nm, mdeg)`` of the plotted minimum."""
        i = int(np.argmin(self.y))
        return float(self.wl[i]), float(self.y[i])

    @property
    def trough_raw(self):
        i = int(np.argmin(self.y_raw))
        return float(self.wl[i]), float(self.y_raw[i])


# ==========================================================================
# 2.  The manifest.  Labels derive from FILENAMES, never from the artwork.
# ==========================================================================
def _t(key, fname, root, temp_c, temp_label, repeat, fraction, label):
    return dict(key=key, path=os.path.join(root, fname), temp_c=temp_c,
                temp_label=temp_label, repeat=repeat, fraction=fraction,
                label=label)


CD_PANELS: dict = {

    # ------------------------------------------------------------------ 1C
    "1C": {
        "letter": "C",
        "title": "CD of SEC-isolated A\u03b242 monomers",
        "xlim": (190.0, 250.0),
        "style_by": "temperature + repeat",
        "traces": [
            _t("6C_r1", "WTmon6C052225.csv", ROOT_052025, 6.0, "6 \u00b0C", 1,
               "SEC-isolated monomer", "6 \u00b0C, repeat 1"),
            _t("6C_r2", "WTmonr26C052225.csv", ROOT_052025, 6.0, "6 \u00b0C", 2,
               "SEC-isolated monomer", "6 \u00b0C, repeat 2"),
            _t("RT_r1", "WTRTmonr1052225.csv", ROOT_052025, RT_AS_C, "RT", 1,
               "SEC-isolated monomer", "RT, repeat 1"),
            _t("RT_r2", "WRrTmonr2052225.csv", ROOT_052025, RT_AS_C, "RT", 2,
               "SEC-isolated monomer", "RT, repeat 2"),
        ],
        "flags": [
            "PUBLISHED LABEL IS WRONG. Doc 10 and the artwork read 't = 0 / 10 / 20 / "
            "30 min'. These four files are 2 temperatures x 2 repeats; the repeat index "
            "is in the filenames and in the pristine matplotlib legend "
            "(052025data\\selected_overlay_WT052025.svg: '6C WT monomers (r1)/(r2)', "
            "'RT WT monomers (r1)/(r2)'). Nothing archived supports a 10-minute grid.",
            "PUBLISHED LABEL IS UNSUPPORTED. '25 C' for the RT traces appears only in "
            "hand-typed Inkscape text and in Doc 10; the filenames and her own legend say "
            "'RT'. Printed as 'RT'; coloured as 25 C (RT_AS_C).",
            "The 6 C / RT figures are INCUBATION temperatures. Channel 6 of these CSVs "
            "reads 22.43 / 22.49 / 22.58 / 22.68 C for the four files -- all were MEASURED "
            "at ambient. That monotonic rise does independently fix the acquisition ORDER, "
            "which is the only time information that survives.",
            "These four are SINGLE SCANS, not averaged and not auto-subtracted -- unlike "
            "1D and 2A. Processing is not uniform across the CD panels.",
            "AMPLITUDES ARE NOT COMPARABLE WITH 2A. 1C's 6 C monomer troughs at "
            "-1.670 mdeg, 2A's 6 C monomer at -16.049 mdeg: same species, same nominal "
            "50 uM, same 197 nm trough, shape correlation 0.983 after peak-normalisation, "
            "but 9.61x apart. No MRE conversion is applied anywhere, so the y axes of 1C "
            "and 2A cannot be read against each other. Needs cuvette pathlength and "
            "measured concentration for the two sessions.",
            "Source data span 190-280 nm; the panel is cropped to 190-250 nm to match the "
            "artwork and 1D/2A. Above 250 nm the traces are baseline noise with residual "
            "offsets +0.141 / +0.094 / -0.037 / +0.072 mdeg against a noise sd of 0.07.",
        ],
    },

    # ------------------------------------------------------------------ 1D
    "1D": {
        "letter": "D",
        "title": "CD of SEC-isolated A\u03b242 oligomers, 37 \u00b0C / 24 h",
        "xlim": (190.0, 250.0),
        "style_by": "temperature",
        "traces": [
            _t("37C_oligo", "37c oligo avg-autosub.csv", ROOT_100325, 37.0,
               "37 \u00b0C", None, "SEC-isolated oligomer",
               "SEC-isolated oligomers, 37 \u00b0C"),
        ],
        "flags": [
            "Drawn in the 37 C ink, the same one 2A's 37 C curve uses: this is the "
            "SEC-isolated oligomer fraction of the same 37 C / 24 h material.",
            "'avg' averages an UNRECORDED number of repeats and 'autosub' is Pro-Data "
            "automatic reference subtraction against an UNARCHIVED blank. No .dsx session "
            "file exists for this project, so bandwidth, scan speed, pathlength and the "
            "identity of the subtracted reference are all unrecoverable.",
            "Doc 10 par.25 calls the beta-strand minimum '~220 nm'. Measured: 216 nm "
            "smoothed, 217 nm raw.",
            "Doc 10 par.19 ('Figure 1C and Figure 2') and par.25 ('Figure 2') both cite "
            "the wrong panel for the SEC-isolated oligomer spectrum; it is this one, 1D.",
            "In Figure1R0.svg this curve appears TWICE (path99-5 and path99-5-1) -- the "
            "assembled Figure 1 carries two 'D)' panels.",
        ],
    },

    # ------------------------------------------------------------------ 2A
    "2A": {
        "letter": "A",
        "title": "CD of A\u03b242 solutions after 24 h at three temperatures",
        "xlim": (190.0, 250.0),
        "style_by": "temperature",
        "traces": [
            _t("6C_mix", "6c mix avg-autosub.csv", ROOT_100325, 6.0, "6 \u00b0C",
               None, "unfractionated", "6 \u00b0C"),
            _t("25C_mix", "25c mix avg-autosub.csv", ROOT_100325, 25.0, "25 \u00b0C",
               None, "unfractionated", "25 \u00b0C"),
            _t("37C_mix", "37c mix avg-autosub.csv", ROOT_100325, 37.0, "37 \u00b0C",
               None, "unfractionated", "37 \u00b0C"),
        ],
        "flags": [
            "'mix' = the unfractionated incubated solution, as opposed to the "
            "SEC-separated 'mono' / 'oligo' fractions. All three are 24 h.",
            "THIS PANEL, NOT 1C, HOLDS THE 24 h / 6 C MEASUREMENT. Doc 10 par.19 "
            "('remains stable at 6 C for 24 h ... (Figure 1C)') cites the wrong figure; "
            "the 6 C / 24 h random-coil spectrum is this panel's 6 C trace, trough "
            "-16.05 mdeg @ 197 nm. The same sentence's 'at 25 C for 30 min' is supported "
            "by nothing at all -- see 1C's first flag.",
            "DEFECT, NOT CORRECTED HERE: the three curves are not on a common zero. "
            "Residual 240-250 nm offsets are +0.782 / +0.615 / +0.174 mdeg, a spread of "
            "1.21 mdeg, and the pipeline applies no baseline correction. Correcting it "
            "would depart from the established processing contract, so it is flagged "
            "instead.",
            "'37c mono avg-autosub.csv' (SEC-isolated monomer after 37 C, trough "
            "-10.283 mdeg @ 198 nm) was loaded by cell[12] and plotted by cell[16], but "
            "appears in no manuscript figure. Its omission may be unintentional.",
            "No Figure2R0.svg exists -- only Figure2R0.png, plus a Figures_Final "
            "Figure2_CD.pdf/.png dated 2025-07-24 that PREDATES the 2025-11-13 clean CD "
            "exports. Which Figure 2 artwork is canonical is unresolved.",
        ],
    },
}

#: Troughs decoded from her delivered vector artwork, ``(nm, mdeg)`` per trace key.
#: These are the numbers :func:`verify` checks a rebuilt panel against.  Method and
#: per-curve RMSE are in PROVENANCE_CD_FIGURES_1-2.md sec. 4; every match was exact
#: (RMSE <= 0.0009 mdeg) with margins of 56x to 1.9e6x over the next-best file.
ARTWORK_REFERENCE: dict = {
    "1C": {"6C_r1": (197.0, -1.670), "6C_r2": (196.0, -1.721),
           "RT_r1": (197.0, -1.592), "RT_r2": (197.0, -1.636)},
    "1D": {"37C_oligo": (216.0, -7.496)},
    "2A": {"6C_mix": (197.0, -16.049), "25C_mix": (202.0, -5.426),
           "37C_mix": (216.0, -7.617)},
}


# ==========================================================================
# 3.  Loading
# ==========================================================================
def load_trace(spec: dict, smooth: bool = True) -> CDTrace:
    """Load and process one manifest entry.  See the module docstring, sec. 1."""
    wl, y_raw = load_cd_csv(spec["path"])
    y = savgol_filter(y_raw, **SAVGOL) if smooth else y_raw.copy()
    return CDTrace(key=spec["key"], path=spec["path"], label=spec["label"],
                   temp_c=spec["temp_c"], temp_label=spec["temp_label"],
                   repeat=spec["repeat"], fraction=spec["fraction"],
                   wl=wl, y=y, y_raw=y_raw)


def load_panel(pid: str, smooth: bool = True):
    """``(panel_spec, [CDTrace, ...])`` for ``'1C'``, ``'1D'`` or ``'2A'``."""
    p = CD_PANELS[pid]
    return p, [load_trace(s, smooth=smooth) for s in p["traces"]]


def trough(t: CDTrace):
    """``(nm, mdeg)`` of the plotted minimum. Convenience for callers."""
    return t.trough


# ==========================================================================
# 4.  The temperature colour key
# ==========================================================================
#: Frozen results of :func:`derive_temp_colors`, keyed by frame ink.  See the
#: module docstring, sec. 3, for the selection table and the tie-break.
#: THE HOUSE TEMPERATURE TRIPLE, shared with ``figure_fig12.TEMP_INKS``.
#:
#: The derived values below were chosen by maximising the WORST pairwise dE against
#: the SEC inks, the frame ink and a greyscale L* ladder, over diverging colormaps.
#: That optimisation has a blind spot: sampling any cold-to-warm diverging map at
#: three points puts BOTH 25 C and 37 C on the warm arm, so those two are close by
#: construction no matter which map wins. Measured on the frozen answers:
#:
#:     frame #000000   25 vs 37   dE = 23.0
#:     frame #3b3b3b   25 vs 37   dE =  9.6      <- indistinguishable in print
#:
#: A.K.P. flagged exactly this on Figure 1 panel c. It was also an internal
#: inconsistency: panel b drew the same three temperatures in TEMP_INKS while panel c
#: drew them in these, so one legend had two colour schemes.
#:
#: Both are now the house triple: light blue / gold / deep red, 25 vs 37 at dE = 62.3,
#: and still clear of the SEC inks (nearest pair dE = 34.4, against 23.0 before). It
#: leaves the diverging-colormap family, which is what buys the separation.
HOUSE_TEMP_INKS: tuple = ("#96c7df", "#d9a441", "#8c1d24")

TEMP_COLORS_BY_FRAME: dict = {
    "#000000": HOUSE_TEMP_INKS,
    "#3b3b3b": HOUSE_TEMP_INKS,
}

#: the superseded derived answers, kept so the reasoning above can be checked
_DERIVED_TEMP_COLORS: dict = {
    "#000000": ("#2f79b5", "#a51429", "#6d0220"),   # worst dE 13.70, 25/37 at 23.0
    "#3b3b3b": ("#3885bc", "#bf3338", "#a21328"),   # worst dE  9.24, 25/37 at  9.6
}

#: The SEC ink pair :data:`TEMP_COLORS_BY_FRAME` was closed against.
#:
#: :func:`derive_temp_colors` is a function of BOTH the frame ink and the SEC
#: inks -- the whole point of sec. 3 above is that the two SEC trace colours are
#: in the ink set the temperature members are separated from.  Keying the frozen
#: table on the frame ink alone was therefore a latent bug: a style carrying a
#: different ``palette`` got the table's answer back unchanged, so the SEC panels
#: recoloured and the CD panels did not.  :func:`temp_colors` now checks this
#: pair as part of the key and re-derives when it does not hold.
HOUSE_SEC_INKS: tuple = tuple(gs.RAMP3[:2])

_RENDERINGS = ("normal", "deutan", "protan", "grey")
TEMP_CMAP = "RdBu_r"
TEMP_L_MAX = 78.0          # same legibility floor style.RAMP9 was closed against
TEMP_MIN_DL = 8.0          # minimum L* step, i.e. the greyscale ordering budget


def _render(cols, kind):
    if kind == "normal":
        return cols
    if kind == "grey":
        return tc.to_greyscale(cols)
    return np.clip(tc.simulate_cvd(cols, kind), 0.0, 1.0)


def _cmap_entries(cmap_name: str, samples: int = 4096):
    """Distinct 8-bit entries of a colormap, in colormap order."""
    import matplotlib as mpl
    u = np.linspace(0.0, 1.0, samples)
    raw = np.array([mpl.colormaps[cmap_name](v)[:3] for v in u])
    q8 = np.round(raw * 255).astype(int)
    uniq, first = np.unique(q8, axis=0, return_index=True)
    order = np.argsort(u[first])
    return (uniq / 255.0)[order], u[first][order]


@functools.lru_cache(maxsize=8)
def derive_temp_colors(frame_ink: str,
                       cmap_name: str = TEMP_CMAP,
                       sec_inks: tuple = None):
    """Re-derive the three temperature inks from scratch. ``(hexes, worst_dE, L*)``.

    This is the function that produced :data:`TEMP_COLORS_BY_FRAME`; running it is
    how those hexes are audited.  Constraints and objective are in the module
    docstring, sec. 3.  ``sec_inks`` defaults to ``style.ramp(3, style)[:2]``, the two
    colours SEC panels 1A and 1B draw with -- they share Figure 1 with the CD
    panels and so belong in the ink set.
    """
    from matplotlib.colors import to_hex, to_rgb

    if sec_inks is None:
        sec_inks = tuple(gs.ramp(3, gs.STYLE)[:2])
    fixed = np.array([to_rgb(h) for h in tuple(sec_inks) + (frame_ink,)])

    cols, us = _cmap_entries(cmap_name)
    lab = tc.rgb_to_lab(cols)
    keep = lab[:, 0] <= TEMP_L_MAX
    cols, us, lab = cols[keep], us[keep], lab[keep]
    if len(cols) < 3:
        raise ValueError(f"{cmap_name}: fewer than 3 entries below L* {TEMP_L_MAX}")

    # worst-rendering dE, entry-to-entry and entry-to-nearest-fixed-ink
    D = f = None
    for k in _RENDERINGS:
        la = tc.rgb_to_lab(_render(cols, k))
        lf = tc.rgb_to_lab(_render(fixed, k))
        d = np.sqrt(((la[:, None, :] - la[None, :, :]) ** 2).sum(-1))
        df = np.sqrt(((la[:, None, :] - lf[None, :, :]) ** 2).sum(-1)).min(axis=1)
        D = d if D is None else np.minimum(D, d)
        f = df if f is None else np.minimum(f, df)

    L, bstar = lab[:, 0], lab[:, 2]
    cold = np.where(bstar < 0)[0]      # 6 C must be blue
    warm = np.where(bstar > 0)[0]      # 37 C must be warm
    best = None
    for a in cold:
        okb = np.where(L <= L[a] - TEMP_MIN_DL)[0]
        okb = okb[okb > a]
        for b in okb:
            okc = warm[(warm > b) & (L[warm] <= L[b] - TEMP_MIN_DL)]
            if not len(okc):
                continue
            s = np.minimum.reduce([np.full(len(okc), min(D[a, b], f[a], f[b])),
                                   D[a, okc], D[b, okc], f[okc]])
            k = int(np.argmax(s))
            if best is None or s[k] > best[0]:
                c = int(okc[k])
                best = (float(s[k]),
                        (to_hex(cols[a]), to_hex(cols[b]), to_hex(cols[c])),
                        (float(L[a]), float(L[b]), float(L[c])))
    if best is None:
        raise ValueError(f"{cmap_name}: no cold->warm triple clears the L* ladder")
    return best[1], best[0], best[2]


def temp_colors(style: gs.GriffinStyle = gs.STYLE) -> tuple:
    """The three inks for :data:`TEMPS`, resolved through the STYLE.

    The key is ``(style.foreground, ramp(3, style)[:2])`` -- both of the things
    :func:`derive_temp_colors` closes against.  Frozen for the two frame inks the
    house style uses AT the house palette; derived and cached for anything else.
    This indirection is the whole point -- when the style revision moves
    ``foreground`` or ``palette``, the temperature key follows without this file
    changing.
    """
    ink = str(style.foreground).lower()
    sec_inks = tuple(gs.ramp(3, style)[:2])
    if ink in TEMP_COLORS_BY_FRAME and sec_inks == HOUSE_SEC_INKS:
        return TEMP_COLORS_BY_FRAME[ink]
    return derive_temp_colors(ink, sec_inks=sec_inks)[0]


def color_for_temp(temp_c: float, style: gs.GriffinStyle = gs.STYLE) -> str:
    """Ink for an incubation temperature in C. Identical in 1C, 1D and 2A."""
    key = float(temp_c)
    diffs = [abs(key - float(t)) for t in TEMPS]
    i = int(np.argmin(diffs))
    if diffs[i] > 1e-6:
        raise KeyError(
            f"T = {temp_c!r} C is not one of the measured incubation temperatures "
            f"{TEMPS}. cd_panels.color_for_temp is a lookup, not an interpolation: "
            "if a new temperature was collected, add it to TEMPS and re-run "
            "derive_temp_colors()."
        )
    return temp_colors(style)[i]


def temp_color_report(style: gs.GriffinStyle = gs.STYLE) -> dict:
    """Audit the resolved key: every pairwise dE over the Figure-1 ink set."""
    from matplotlib.colors import to_rgb
    inks = dict(zip(("6 C", "25 C", "37 C"), temp_colors(style)))
    inks["SEC a"], inks["SEC b"] = gs.ramp(3, style)[:2]
    inks["frame"] = style.foreground
    names = list(inks)
    cols = np.array([to_rgb(inks[n]) for n in names])
    out = {"inks": dict(inks), "pairs": {}}
    for i, j in itertools.combinations(range(len(names)), 2):
        per = {}
        for k in _RENDERINGS:
            la = tc.rgb_to_lab(_render(cols[[i, j]], k))
            per[k] = float(tc.delta_e76(la[0], la[1]))
        out["pairs"][f"{names[i]} / {names[j]}"] = dict(per, worst=min(per.values()))
    out["worst_involving_temperature"] = min(
        v["worst"] for k, v in out["pairs"].items()
        if any(p.strip() in ("6 C", "25 C", "37 C") for p in k.split("/")))
    return out


# ==========================================================================
# 5.  CD axis conventions
# ==========================================================================
# A wavelength axis is not an elution-volume axis, so style.x_major (2.0 mL) does
# not carry over.  Over the 190-250 nm window a 10 nm major with 2 nm minors gives
# seven labelled majors -- the same label density style.x_major gives over the SEC
# panels' 6-20 mL window -- and puts a minor on every even nanometre, matching the
# 1 nm sampling.  These are the only geometry numbers this module owns; every
# size, weight, width and ink still comes from the style object.
CD_X_MAJOR: float = 10.0
CD_X_MINOR: float = 2.0

#: y padding as a fraction of the data span, taken from the style so it tracks the
#: house look (style.y_pad_normalized is 0.05 of the normalised 0-1 SEC range).
_Y_PAD_FRAC = gs.STYLE.y_pad_normalized

#: Her artwork has NO zero line -- it has ``plt.grid(True)``, which the house style
#: forbids (``axes.grid: False``).  So no zero line is drawn by default.  A CD reader
#: who wants the sign reference back can pass ``zeroline=True``.
ZEROLINE_DEFAULT: bool = False

#: Dash pattern for the REPEAT dimension in 1C, in multiples of the trace linewidth
#: so it tracks the house stroke weight rather than fixing points.
REPEAT_DASH_EM: tuple = (2.6, 1.4)

#: Corners :func:`place_legend` will consider, in preference order for ties.
LEGEND_CORNERS: tuple = ("lower right", "upper right", "upper left", "lower left")

#: Clearance demanded between the legend box and the nearest trace vertex, as a
#: fraction of the axes box diagonal.  Not a look: it is the margin that stops the
#: 1D legend handle from being drawn across the 216 nm trough.
LEGEND_CLEARANCE: float = 0.02


def _nice_step(span: float, target: int = 6) -> float:
    """A 1/2/5-decade tick step giving roughly ``target`` majors across ``span``."""
    if span <= 0:
        return 1.0
    raw = span / float(target)
    mag = 10.0 ** np.floor(np.log10(raw))
    for m in (1.0, 2.0, 2.5, 5.0, 10.0):
        if raw <= m * mag:
            return m * mag
    return 10.0 * mag


def cd_ylim(traces, xlim, pad_frac: float = _Y_PAD_FRAC):
    """Data-driven y limits over the drawn x window, snapped out to a tick step."""
    lo, hi = np.inf, -np.inf
    for t in traces:
        m = (t.wl >= xlim[0]) & (t.wl <= xlim[1])
        if m.any():
            lo = min(lo, float(t.y[m].min()))
            hi = max(hi, float(t.y[m].max()))
    span = hi - lo
    lo -= pad_frac * span
    hi += pad_frac * span
    step = _nice_step(hi - lo)
    return (float(np.floor(lo / step) * step), float(np.ceil(hi / step) * step)), step


def cd_axes(ax, style: gs.GriffinStyle = gs.STYLE, *, xlim=(190.0, 250.0),
            ylim=None, y_major=None, zeroline: bool = ZEROLINE_DEFAULT):
    """The house box frame with CD (wavelength / ellipticity) conventions.

    Delegates every stroke, size and ink to :func:`secfigs.style.axes`; only the
    locators and limits are CD-specific.
    """
    y_minor = None if y_major is None else y_major / 2.0
    gs.axes(ax, style, normalized=False,
            x_major=CD_X_MAJOR, x_minor=CD_X_MINOR,
            y_major=y_major, y_minor=y_minor,
            x_range=xlim, y_range=ylim)
    if zeroline:
        ax.axhline(0.0, color=style.foreground,
                   lw=style.spine_linewidth, zorder=0.5)
    return ax


def place_legend(ax, traces, style: gs.GriffinStyle = gs.STYLE,
                 corners=LEGEND_CORNERS, clearance: float = LEGEND_CLEARANCE):
    """Put the legend in the emptiest corner, MEASURED, not chosen by hand.

    Draws the legend once to learn its true size, then scores each corner by the
    number of trace vertices its box would cover and, among corners that cover
    none, by how far the nearest vertex is.  A fixed corner would have to be
    revisited every time the house style changed the legend font size; this does
    not.  It is also what catches 1D, whose only sensible corner is upper right --
    at lower right the legend handle lands on the 216 nm trough.
    """
    leg = gs.legend(ax, style, loc=corners[0])
    ax.figure.canvas.draw()
    bb = leg.get_window_extent().transformed(ax.transAxes.inverted())
    w, h = bb.width, bb.height

    pts = []
    for t in traces:
        xy = ax.transLimits.transform(np.column_stack([t.wl, t.y]))
        pts.append(xy)
    pts = np.vstack(pts) if pts else np.zeros((0, 2))
    inside = (pts[:, 0] >= -0.02) & (pts[:, 0] <= 1.02)
    pts = pts[inside]

    best = None
    for loc in corners:
        vert, horiz = loc.split()
        x0 = 0.0 if horiz == "left" else 1.0 - w
        y0 = 0.0 if vert == "lower" else 1.0 - h
        x1, y1 = x0 + w, y0 + h
        pad = clearance
        cover = int(((pts[:, 0] >= x0 - pad) & (pts[:, 0] <= x1 + pad) &
                     (pts[:, 1] >= y0 - pad) & (pts[:, 1] <= y1 + pad)).sum())
        dx = np.maximum.reduce([x0 - pts[:, 0], np.zeros(len(pts)), pts[:, 0] - x1])
        dy = np.maximum.reduce([y0 - pts[:, 1], np.zeros(len(pts)), pts[:, 1] - y1])
        near = float(np.min(np.hypot(dx, dy))) if len(pts) else 1.0
        score = (cover, -near)
        if best is None or score < best[0]:
            best = (score, loc, cover, near)

    _, loc, cover, near = best
    if loc != corners[0]:
        leg.remove()
        leg = gs.legend(ax, style, loc=loc)
    return leg, dict(loc=loc, covered_vertices=cover, nearest_vertex=near,
                     size_axes_frac=(float(w), float(h)))


# ==========================================================================
# 6.  draw(ax) -- standalone or into a combined figure
# ==========================================================================
def draw(ax, pid: str, style: gs.GriffinStyle = gs.STYLE, *,
         smooth: bool = True, title: bool = True, letter: bool = True,
         legend: bool = True, labels: bool = True,
         zeroline: bool = ZEROLINE_DEFAULT, ylim=None):
    """Draw CD panel ``pid`` into an existing Axes. Returns a metadata dict.

    Every panel of Figures 1 and 2 exposes this signature, so the combined figure
    the user will ask for next is a grid of ``draw(ax, pid)`` calls.  Call it inside
    ``secfigs.style.apply()`` if you want the rcParams too; nothing here depends on
    them, because :func:`cd_axes` sets every property explicitly.
    """
    p, traces = load_panel(pid, smooth=smooth)
    xlim = p["xlim"]

    auto_ylim, step = cd_ylim(traces, xlim)
    if ylim is None:
        ylim = auto_ylim

    dash_on, dash_off = (v * style.trace_linewidth for v in REPEAT_DASH_EM)

    drawn = []
    for t in traces:
        kw = dict(color=color_for_temp(t.temp_c, style),
                  lw=style.trace_linewidth,
                  label=t.label,
                  solid_capstyle="round")
        # 1C carries a second, NESTED dimension -- the repeat.  Colour stays a pure
        # function of temperature; the repeat is carried by dash, so the key never
        # means two things at once.
        if t.repeat == 2:
            kw["dashes"] = (dash_on, dash_off)
            kw["dash_capstyle"] = "round"
        (line,) = ax.plot(t.wl, t.y, **kw)
        drawn.append((t, line))

    cd_axes(ax, style, xlim=xlim, ylim=ylim, y_major=step, zeroline=zeroline)
    if labels:
        gs.set_labels(ax, X_LABEL, Y_LABEL, style)
    if title:
        gs.set_title(ax, p["title"], style)
    leg_info = None
    if legend:
        _, leg_info = place_legend(ax, traces, style)
    if letter:
        gs.panel_letter(ax, p["letter"], style)

    return {
        "pid": pid,
        "panel": p,
        "traces": traces,
        "n_traces": len(traces),
        "xlim": xlim,
        "ylim": ylim,
        "y_major": step,
        "smoothed": smooth,
        "savgol": dict(SAVGOL) if smooth else None,
        "y_quantity": Y_LABEL,
        "legend": leg_info,
        "colors": {t.key: color_for_temp(t.temp_c, style) for t in traces},
        "labels": {t.key: t.label for t in traces},
        "troughs": {t.key: t.trough for t in traces},
        "troughs_raw": {t.key: t.trough_raw for t in traces},
        "flags": list(p["flags"]),
    }


# ==========================================================================
# 7.  Verification against the decoded artwork
# ==========================================================================
def verify(pid: str, smooth: bool = True) -> dict:
    """Rebuilt troughs vs the ones decoded from her vector artwork.

    Returns per-trace ``(nm, mdeg)`` for both, the deviations, and the panel-level
    maxima.  A rebuild that reproduces her curve gives 0.0 nm and <= 0.001 mdeg.
    """
    p, traces = load_panel(pid, smooth=smooth)
    ref = ARTWORK_REFERENCE[pid]
    rows = []
    for t in traces:
        nm, mdeg = t.trough
        rnm, rmdeg = ref[t.key]
        rows.append(dict(key=t.key, label=t.label,
                         built=(nm, mdeg), artwork=(rnm, rmdeg),
                         d_nm=abs(nm - rnm), d_mdeg=abs(mdeg - rmdeg),
                         file=os.path.basename(t.path)))
    return dict(pid=pid, n_traces=len(traces), n_reference=len(ref), rows=rows,
                max_d_nm=max(r["d_nm"] for r in rows),
                max_d_mdeg=max(r["d_mdeg"] for r in rows),
                verdict=("MATCH"
                         if max(r["d_nm"] for r in rows) == 0.0
                         and max(r["d_mdeg"] for r in rows) <= 1e-3
                         else "MISMATCH"))
