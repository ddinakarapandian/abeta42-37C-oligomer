# -*- coding: utf-8 -*-
r"""
build_fig1_panels.py -- rebuild SEC panels 1A and 1B of the 37 C Abeta42 oligomer
manuscript (doc 10) from RAW UNICORN exports, in the Griffin Figure 1A house style.

Reads   : the CSV manifest in secfigs/panel_data.py  (labels derived from FILENAMES)
Style   : secfigs/style.py                            (measured off Griffin Fig 1A)
Writes  : figures/panels/1A.{svg,png}, figures/panels/1B.{svg,png}   (PNG at 600 dpi)

Processing contract, reimplemented here so nothing under /path/to/Dropbox is
imported (importing her sec1.py would drop __pycache__ into the read-only tree):

    1. pd.read_csv(f, encoding='utf-16', sep='\t', skiprows=3)
    2. first two NUMERIC columns -> Volume_ml, Absorbance_mAU
    3. truncate 0 <= Volume_ml <= 25
    4. A = A - min(A)  only if min(A) < 0
    5. normalized: (A - min) / (max - min) over that same 0-25 mL window

No smoothing, no fitted baseline, no volume-axis transform.

Peak pointers are placed at DETECTED peaks (see find_peaks_for_panel), not at
hardcoded coordinates.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from secfigs import style as gs                      # noqa: E402
from secfigs.panel_data import PANELS                # noqa: E402

import matplotlib                                    # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402

# Live text in the SVG (Nature requires editable vector type) and no vertex
# decimation -- both were missing, so these panels shipped with their type baked
# to glyph outlines. path.simplify collapses in 72-dpi units, which visibly
# flattens the sharp oligomer peak flanks at this panel size.
import matplotlib as _mpl
_mpl.rcParams["svg.fonttype"] = "none"
_mpl.rcParams["path.simplify"] = False


OUTDIR = os.path.join(HERE, "figures", "panels")


def sec_colors(style: gs.GriffinStyle = gs.STYLE) -> tuple:
    """The two SEC trace inks, resolved through the STYLE at DRAW time.

    These used to be frozen into :data:`SPEC` at import
    (``colors=tuple(_ramp(3)[:2])``), evaluated once against the module-level
    constants.  A replaced palette could therefore never reach 1A or 1B no
    matter what ``ramp()`` did.  Resolving here means it does.
    """
    return tuple(gs.ramp(3, style)[:2])

# ---------------------------------------------------------------- window definitions
VMIN, VMAX = 0.0, 25.0        # processing window (step 3 of the contract)
OLIGO_WIN = (7.0, 11.0)       # elution window searched for the oligomer peak, mL
MONO_WIN = (12.0, 18.0)       # elution window searched for the monomer peak, mL

# ---------------------------------------------------------------- per-panel drawing spec
# Titles state only what the DATA / the filenames support.  1A's temperature and
# concentration are disputed (caption 50 uM / 6 C vs notebook 100 uM / 4-5 C) and the
# UNICORN exports carry neither field, so the title says "cold incubation" and nothing
# more.  1B's 37 C comes from the filenames themselves ("...reinj 37C 021925...").
#
# Trace colour is deliberately NOT in this dict -- see :func:`sec_colors`.
SPEC = {
    "1A": dict(
        letter="A",
        title="Normalized SEC Profiles of A\u03b242 after Cold Incubation",
        xlim=(6.0, 20.0),
    ),
    "1B": dict(
        letter="B",
        title="Normalized SEC Profiles of A\u03b242 at 37 \u00b0C",
        xlim=(6.0, 20.0),
    ),
}

# Pointer geometry, expressed as fractions so it tracks the axis limits.
# Measured on Griffin Fig 1A (xlim span 12, oligomer peak 8.58 / monomer peak 14.45):
#   oligomer triangle  +0.75 mL from the peak  = +0.0625 of the x span, at 0.48 * peak h
#   monomer  triangle  -0.83 mL from the peak  = -0.0692 of the x span, at 0.77 * peak h
OLIGO_DX_FRAC = +0.0625
MONO_DX_FRAC = -0.0692
OLIGO_Y_FRAC = 0.48
MONO_Y_FRAC = 0.77
OLIGO_Y_FLOOR = 0.10   # keep the pointer legible when there is no oligomer peak


# ============================================================== data
def load_trace(path: str, normalized: bool):
    """The documented three-step pipeline.  Returns (V, A_plotted, A_adjusted_raw)."""
    df = pd.read_csv(path, encoding="utf-16", sep="\t", skiprows=3)
    num = df.select_dtypes(include=[np.number]).columns[:2]
    if len(num) < 2:
        raise ValueError("%s: fewer than two numeric columns" % path)
    df = df[num].copy()
    df.columns = ["Volume_ml", "Absorbance_mAU"]
    df = df[(df.Volume_ml >= VMIN) & (df.Volume_ml <= VMAX)].dropna()
    df = df.sort_values("Volume_ml")
    V = df.Volume_ml.values.astype(float)
    A = df.Absorbance_mAU.values.astype(float)
    if A.min() < 0.0:                      # step 4 -- conditional offset only
        A = A - A.min()
    raw = A.copy()
    if normalized:                         # step 5 -- min-max over the SAME window
        A = (A - A.min()) / (A.max() - A.min())
    return V, A, raw


def load_panel(pid: str):
    p = PANELS[pid]
    out = []
    for path, label in zip(p["csvs"], p["labels"]):
        V, A, raw = load_trace(path, p["normalized"])
        out.append(dict(path=path, label=label, V=V, A=A, raw=raw))
    return p, out


# ============================================================== peak detection
def _window_peak(V, A, win):
    """Highest-prominence local maximum inside `win`; falls back to the window argmax."""
    m = (V >= win[0]) & (V <= win[1])
    if not m.any():
        return None
    Vw, Aw = V[m], A[m]
    idx, props = find_peaks(Aw, prominence=1e-9)
    if len(idx):
        j = idx[int(np.argmax(props["prominences"]))]
    else:
        j = int(np.argmax(Aw))
    return float(Vw[j]), float(Aw[j])


def find_peaks_for_panel(traces):
    """Detected oligomer/monomer peak per trace, plus the panel-level pointer anchors."""
    per = []
    for t in traces:
        per.append(dict(label=t["label"],
                        oligomer=_window_peak(t["V"], t["A"], OLIGO_WIN),
                        monomer=_window_peak(t["V"], t["A"], MONO_WIN)))
    ol = [d["oligomer"] for d in per if d["oligomer"]]
    mo = [d["monomer"] for d in per if d["monomer"]]
    # anchor the pointer on the trace with the STRONGEST feature -- that is the curve
    # the reader's eye goes to, and it is what the reference panel does.
    anchor_ol = max(ol, key=lambda p: p[1]) if ol else None
    anchor_mo = max(mo, key=lambda p: p[1]) if mo else None
    return per, anchor_ol, anchor_mo


# ============================================================== drawing
def build_panel(pid: str, style: gs.GriffinStyle = gs.STYLE):
    spec = SPEC[pid]
    p, traces = load_panel(pid)
    per, anchor_ol, anchor_mo = find_peaks_for_panel(traces)

    os.makedirs(OUTDIR, exist_ok=True)
    xlim = spec["xlim"]
    xspan = xlim[1] - xlim[0]

    with gs.apply(style):
        fig, ax = gs.figure(style)
        for t, c in zip(traces, sec_colors(style)):
            ax.plot(t["V"], t["A"], color=c, lw=style.trace_linewidth,
                    label=t["label"], solid_capstyle="round")

        gs.axes(ax, style, normalized=p["normalized"], x_major=2.0, x_range=xlim)
        gs.set_labels(ax, "Elution Volume (ml)", "Normalized Absorbance", style)
        gs.set_title(ax, spec["title"], style)
        gs.legend(ax, style)

        if anchor_ol is not None:
            xo, ho = anchor_ol
            # marker="<": the pointer sits OLIGO_DX_FRAC to the RIGHT of the
            # oligomer apex, so it has to point back left at the peak it names.
            # Defaulting to ">" (as the Griffin reference does) aims it into the
            # inter-peak valley and, followed literally, lands on the monomer --
            # the same defect already fixed once in secfigs/figure_sec.py.
            gs.annotate_peak(ax, xo + OLIGO_DX_FRAC * xspan,
                             max(OLIGO_Y_FRAC * ho, OLIGO_Y_FLOOR),
                             "Oligomer", side="right", style=style, marker="<")
        if anchor_mo is not None:
            xm, hm = anchor_mo
            gs.annotate_peak(ax, xm + MONO_DX_FRAC * xspan, MONO_Y_FRAC * hm,
                             "Monomer", side="left", style=style, marker=">")

        gs.panel_letter(ax, spec["letter"], style)

        svg = os.path.join(OUTDIR, "%s.svg" % pid)
        png = os.path.join(OUTDIR, "%s.png" % pid)
        fig.savefig(svg, facecolor="white")
        fig.savefig(png, dpi=600, facecolor="white")
        plt.close(fig)

    return dict(pid=pid, panel=p, traces=traces, peaks=per,
                anchor_oligomer=anchor_ol, anchor_monomer=anchor_mo,
                svg=svg, png=png, xlim=xlim)


def main():
    for pid in ("1A", "1B"):
        r = build_panel(pid)
        print("=" * 78)
        print("%s  %s" % (pid, r["panel"]["title"]))
        print("   xlim %s   traces %d" % (r["xlim"], len(r["traces"])))
        for t, pk in zip(r["traces"], r["peaks"]):
            o = pk["oligomer"]
            m = pk["monomer"]
            print("   %-9s oligomer %s   monomer %s   %s"
                  % (t["label"],
                     "V=%6.3f h=%.4f" % o if o else "none",
                     "V=%6.3f h=%.4f" % m if m else "none",
                     os.path.basename(t["path"])))
        print("   wrote %s" % r["svg"])
        print("   wrote %s" % r["png"])


if __name__ == "__main__":
    main()
