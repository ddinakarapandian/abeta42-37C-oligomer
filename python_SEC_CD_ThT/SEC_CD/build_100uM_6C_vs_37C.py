"""build_100uM_6C_vs_37C.py — 100 uM WT Abeta42 SEC at 6 C and at 37 C,
non-normalised and normalised, as four separate panels.

WHY THIS EXISTS.  Figure 2 is mostly 100 uM, so the useful comparison is the same
concentration at the two incubation temperatures.  The non-normalised panels are
the point: normalising rescales every trace to its own maximum and therefore
destroys exactly the quantity that separates the two temperatures.

WHAT THE DATA SAY (monomer AUC over 12-18 mL, relative to that series' t = 0):
    37 C   1.00  0.91  0.81  0.81  0.73  0.71  0.63  0.31  0.20
           0h    2h    4h    6h    8h    12h   24h   48h   72h
     6 C   1.00        0.95  0.97  0.49
           0     [1 d] 5 d   12 d  22 d
At 37 C four fifths of the monomer is gone by 72 h.  At 6 C the monomer is intact
to within 5 % at twelve DAYS.  Hours against weeks; the normalised panels cannot
show it.

THE ONE TRAP, and it is why the 1 d cold point is absent from the raw panel.
The cold series was not run on one loop:
    t = 0     500 uL      AUC 12.604
    t = 1 d   100 uL      AUC  2.276     <-- 0.18x, and 12.604/2.276 = 5.54 ~ 5x
    t = 5 d   500 uL      AUC 11.962
    t = 12 d  500 uL      AUC 12.256
    t = 22 d  500 uL      AUC  6.157
The 5.54 ratio tracks the 5x loop ratio, not the sample: the peptide concentration
was essentially unchanged.  Plotting that point in raw mAU beside the others would
manufacture a 5x drop at 24 h that did not happen.  It is therefore EXCLUDED from
the non-normalised cold panel and INCLUDED in the normalised one, where rescaling
removes the loop difference.  Every 37 C run used the 500 uL loop, so all nine are
raw-comparable and none is dropped.

TEMPERATURE LABEL.  The manuscript says 6 C; the source notebook says 4-5 C and two
of the cold files are named "4C".  The panels are labelled 6 C to match the
manuscript, and the discrepancy is unresolved in the bench record -- see
PROVENANCE_FIGURES_1-4.md.  Set COLD_LABEL below to change it in one place.

    python build_100uM_6C_vs_37C.py                 # all four panels
    python build_100uM_6C_vs_37C.py --only 37C_raw
    python build_100uM_6C_vs_37C.py --combined      # also a 2x2 overview page
"""
from __future__ import annotations

import argparse
import os
import sys

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from secfigs import style as gs                      # the project's Griffin style

RES = r"/path/to/Ab42oligomers37C/data/SEC\Run_Results"      # READ ONLY
E22K = os.path.join(RES, "E22K")
WT = os.path.join(RES, "WT")
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

MM = 1.0 / 25.4
COLD_LABEL = "6 °C"
AUC_WINDOW = (12.0, 18.0)     # the monomer peak, same window panel_data uses

#: (label, path, loop_uL).  Loop size is read from the filename and carried
#: explicitly, because it decides whether a trace may be drawn in raw mAU.
COLD = [
    ("t = 0",    f"{E22K}/WT immediate reinj 021325 500ul Loop no fractions 001.csv", 500),
    ("t = 1 d",  f"{E22K}/WT 0213 monomer reinj 021425 100ul Loop no fractions  001.csv", 100),
    ("t = 5 d",  f"{E22K}/500ul Loop no fractions WT 5 day reinjection 001.csv", 500),
    ("t = 12 d", f"{E22K}/500ul Loop no fractions WT 12d reinjection 4C 22525 001.csv", 500),
    ("t = 22 d", f"{E22K}/030725 WT noesy 22day 4C reinj  001.csv", 500),
]
HOT = [
    ("t = 0 h",  f"{WT}/WT 0 hr reinj 37C 021925 500ul Loop no fractions.csv", 500),
    ("t = 2 h",  f"{WT}/WT 2hr reinj 37C 021925 500ul Loop no fractions.csv", 500),
    ("t = 4 h",  f"{WT}/WT 4hr reinj 37C 021925 500ul Loop no fractions.csv", 500),
    ("t = 6 h",  f"{WT}/WT 6hr reinj 37C 021925 500ul Loop All fractions.csv", 500),
    ("t = 8 h",  f"{WT}/WT 8hr reinj 37C 021925 500ul Loop no fractions 001.csv", 500),
    ("t = 12 h", f"{WT}/WT 12hr reinj 37C 021925 500ul Loop no fractions  001 001.csv", 500),
    ("t = 24 h", f"{WT}/WT 24hr reinj 37C 022025 500ul Loop no fractions 001.csv", 500),
    ("t = 48 h", f"{WT}/WT 48hr reinj 37C 022125 500ul Loop no fractions  001 001.csv", 500),
    ("t = 72 h", f"{WT}/WT 72hr reinj 37 C 022225 500ul Loop no fractions  001 001.csv", 500),
]

trapz = getattr(np, "trapezoid", np.trapz)


# --------------------------------------------------------------------------- io

def load(path, vmax=25.0):
    """UNICORN export -> (volume mL, absorbance mAU), exactly as the project does:
    utf-16, tab-separated, three header rows; first two numeric columns; truncate
    to 0 <= V <= vmax; subtract the minimum ONLY if it is negative.  No smoothing,
    no fitted baseline."""
    df = pd.read_csv(path, encoding="utf-16", sep="\t", skiprows=3)
    num = [c for c in df.columns
           if pd.to_numeric(df[c], errors="coerce").notna().sum() > 10]
    if len(num) < 2:
        raise ValueError(f"{os.path.basename(path)}: fewer than two numeric columns")
    v = pd.to_numeric(df[num[0]], errors="coerce").to_numpy(float)
    a = pd.to_numeric(df[num[1]], errors="coerce").to_numpy(float)
    m = np.isfinite(v) & np.isfinite(a) & (v >= 0) & (v <= vmax)
    v, a = v[m], a[m]
    if a.size and a.min() < 0:
        a = a - a.min()
    return v, a


def auc(v, a, lo=AUC_WINDOW[0], hi=AUC_WINDOW[1]):
    s = (v >= lo) & (v <= hi)
    return float(trapz(a[s], v[s])) if s.sum() > 1 else float("nan")


def series(entries, normalized, drop_mismatched_loops):
    """Load a series. In raw mode, traces whose loop differs from the majority are
    dropped and returned separately -- a raw mAU comparison across loops is a
    volume artefact, not a measurement."""
    loops = [e[2] for e in entries]
    ref = max(set(loops), key=loops.count)
    out, dropped = [], []
    for label, path, loop in entries:
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        v, a = load(path)
        rec = dict(label=label, path=path, loop=loop, v=v, a=a, auc=auc(v, a))
        if (not normalized) and drop_mismatched_loops and loop != ref:
            dropped.append(rec)
            continue
        if normalized:
            s = (v >= AUC_WINDOW[0] - 4) & (v <= AUC_WINDOW[1] + 4)
            span = a[s].max() - a[s].min() if s.sum() else a.max() - a.min()
            rec["y"] = (a - a[s].min()) / span if span > 0 else a
        else:
            rec["y"] = a
        out.append(rec)
    return out, dropped, ref


# ------------------------------------------------------------------------ panel

#: Every size in GriffinStyle is calibrated against her 12 x 9.53 in reference
#: figure (tick labels 21.5 pt, axis labels 24 pt, panel letters 54 pt).  Dropped
#: onto an 88 mm panel unscaled, the type is roughly three times the page.  Scale
#: the whole style by one factor so the proportions she calibrated are preserved
#: and only the absolute size changes.
_SCALE_FIELDS = (
    "tick_labelsize", "tick_pad", "axis_label_size", "axis_label_pad",
    "title_size", "title_pad", "legend_fontsize", "annotation_size",
    "panel_letter_size", "trace_linewidth",
    "major_tick_length", "minor_tick_length",
)


def scaled_style(width_mm, base: "gs.GriffinStyle" = gs.STYLE):
    from dataclasses import replace
    k = (width_mm * MM) / base.figsize[0]
    kw = {f: getattr(base, f) * k for f in _SCALE_FIELDS if hasattr(base, f)}
    kw["figsize"] = (width_mm * MM, width_mm * MM * 0.773)
    return replace(base, **kw)


def panel(recs, dropped, ref_loop, title, normalized, stem,
          width_mm=88.0, height_mm=68.0, xlim=(5.0, 22.0),
          formats=("svg", "pdf", "png")):
    st = scaled_style(width_mm)
    # gs.apply is a @contextmanager -- calling it bare never enters the generator,
    # so no Griffin rcParams are installed and matplotlib falls back to DejaVu Sans.
    # That shipped: four delivered PDFs here were typeset in DejaVu while every
    # other figure in the project is Arial.
    with gs.apply(st):
        fig = plt.figure(figsize=(width_mm * MM, height_mm * MM))
        ax = fig.add_axes([0.165, 0.185, 0.805, 0.735])

        colors = gs.ramp(len(recs), st)
        for rec, c in zip(recs, colors):
            ax.plot(rec["v"], rec["y"], color=c, lw=st.trace_linewidth,
                    solid_capstyle="round", label=rec["label"], zorder=3)

        # normalized=False is what stops axes() forcing ylim (-0.05, 1.05); the raw
        # panels are in mAU and their range is set by the data.
        ymax = max(r["y"].max() for r in recs)
        gs.axes(ax, style=st, normalized=normalized, x_major=5.0, x_minor=1.0,
                x_range=xlim,
                y_range=None if normalized else (0.0, ymax * 1.06),
                y_major=None if normalized else 2.0,
                y_minor=None if normalized else 0.5)
        gs.set_labels(ax, "Elution volume (mL)",
                      "Normalised absorbance" if normalized else "Absorbance (mAU)",
                      style=st)
        gs.set_title(ax, title, style=st)
        gs.legend(ax, style=st, loc="upper right")

        # savefig MUST be inside the context.  matplotlib resolves the font family
        # when the file is written, not when the artist is created, so saving after
        # the contextmanager has restored the old rcParams silently typesets the
        # whole figure in the fallback face.  That is the actual mechanism behind
        # the four DejaVu Sans PDFs: entering the context was necessary and not
        # sufficient.  pair_24h() below always saved inside, which is why it alone
        # came out in Arial.
        paths = []
        for ext in formats:
            p = os.path.join(OUT, f"{stem}.{ext}")
            fig.savefig(p, dpi=600 if ext == "png" else None, facecolor="white")
            paths.append(p)
    plt.close(fig)
    return paths


# ------------------------------------------- the head-to-head 24 h comparison

def bold_style(base=None):
    """Delegates to secfigs.style.bold_variant so there is ONE definition of
    what 'bold' means across the project -- text weights AND the frame."""
    return gs.bold_variant(base)


def pair_24h(stem="SEC_100uM_24h_6C_vs_37C", xlim=(6.0, 20.0), panel_letter=None,
             bold=False, formats=("svg", "pdf", "png")):
    """One panel, two traces: 100 uM at 24 h, 6 C against 37 C, min-max normalised.

    Normalisation is the manuscript's, not a variant of it: this calls
    ``build_fig1_panels.load_trace(path, normalized=True)`` directly rather than
    reimplementing it, so the processing is by construction the same one behind
    Figures 1-4 -- utf-16 / tab / skiprows=3, first two numeric columns, sort
    ascending, truncate 0 <= V <= 25 mL, subtract the minimum ONLY if negative,
    then (A - A.min()) / (A.max() - A.min()) over that same window.

    Because each trace is scaled to its own maximum, the panel compares PEAK
    SHAPE AND POSITION -- where the material is -- and not how much of it there
    is.  The amplitude comparison lives in the non-normalised panels.
    """
    from build_fig1_panels import load_trace       # the manuscript's own pipeline

    picks = [
        ("6 °C, 24 h",  COLD[1]),      # 021425
        ("37 °C, 24 h", HOT[6]),       # 022025
    ]
    recs = []
    for label, (_l, path, loop) in picks:
        v, y, raw = load_trace(path, normalized=True)
        s = (v >= AUC_WINDOW[0]) & (v <= AUC_WINDOW[1])
        o = (v >= 7.0) & (v <= 10.0)
        recs.append(dict(label=label, v=v, y=y, loop=loop,
                         mono_peak_ml=float(v[s][np.argmax(y[s])]),
                         olig_frac=float(y[o].max()) if o.sum() else 0.0))

    # Draw exactly as build_fig1_panels does: the reference GriffinStyle at its own
    # figsize, x_major 2.0 with the style's inherited AutoMinorLocator, the
    # manuscript's own axis wording, a panel letter and no scaling.  An earlier
    # version of this function scaled the style to 88 mm and forced x_major 5 /
    # x_minor 1, which is why its ticks did not match the other SEC panels.
    style = bold_style() if bold else gs.STYLE
    with gs.apply(style):
        fig, ax = gs.figure(style)
        for rec, c in zip(recs, ("#2171B5", "#CB181D")):   # cold blue, hot red
            ax.plot(rec["v"], rec["y"], color=c, lw=style.trace_linewidth,
                    label=rec["label"], solid_capstyle="round")

        gs.axes(ax, style, normalized=True, x_major=2.0, x_range=xlim)
        gs.set_labels(ax, "Elution Volume (ml)", "Normalized Absorbance", style)
        gs.legend(ax, style)
        if panel_letter:
            gs.panel_letter(ax, panel_letter, style)

        paths = []
        for ext in formats:
            p = os.path.join(OUT, f"{stem}.{ext}")
            fig.savefig(p, dpi=600 if ext == "png" else None, facecolor="white")
            paths.append(p)
    plt.close(fig)

    print("\n=== 24 h, 100 µM, min-max normalised (build_fig1_panels.load_trace)")
    for r in recs:
        print(f"   {r['label']:12} monomer peak at {r['mono_peak_ml']:5.2f} mL   "
              f"oligomer max (7–10 mL) {r['olig_frac']:.3f} of its own maximum")
    for p in paths:
        print("   wrote", p)
    return paths


# ------------------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", choices=("6C_raw", "6C_norm", "37C_raw", "37C_norm", "pair24h"))
    ap.add_argument("--keep-loop-mismatch", action="store_true",
                    help="draw mismatched-loop traces in the raw panels too "
                         "(this is a volume artefact -- off by default)")
    ap.add_argument("--bold", action="store_true",
                    help="bold tick/axis/legend text (a DEPARTURE from the "
                         "reference, whose ticks and axis labels are Regular)")
    a = ap.parse_args(argv)

    jobs = [
        ("6C_raw",   COLD, False, f"100 µM, {COLD_LABEL}",  "SEC_100uM_6C_nonnormalized"),
        ("6C_norm",  COLD, True,  f"100 µM, {COLD_LABEL}",  "SEC_100uM_6C_normalized"),
        ("37C_raw",  HOT,  False, "100 µM, 37 °C",          "SEC_100uM_37C_nonnormalized"),
        ("37C_norm", HOT,  True,  "100 µM, 37 °C",          "SEC_100uM_37C_normalized"),
    ]
    if a.only == "pair24h":
        pair_24h(bold=a.bold)
        return 0
    if a.only:
        jobs = [j for j in jobs if j[0] == a.only]

    for key, entries, norm, title, stem in jobs:
        recs, dropped, ref = series(entries, norm,
                                    drop_mismatched_loops=not a.keep_loop_mismatch)
        print(f"\n=== {key}  ({'normalised' if norm else 'raw mAU'})  "
              f"majority loop {ref} µL")
        base = recs[0]["auc"] if recs else float("nan")
        for r in recs:
            print(f"   {r['label']:9} loop {r['loop']:3d} µL  max {r['a'].max():7.2f} mAU  "
                  f"AUC(12–18) {r['auc']:8.3f}  = {r['auc']/base:5.2f}× t0")
        for r in dropped:
            print(f"   {r['label']:9} loop {r['loop']:3d} µL  DROPPED from the raw panel: "
                  f"loop differs from the {ref} µL majority; its AUC {r['auc']:.3f} is "
                  f"{base/r['auc']:.2f}× lower, tracking the "
                  f"{ref/r['loop']:.0f}× loop ratio, not the sample.")
        paths = panel(recs, dropped, ref, title, norm, stem)
        for p in paths:
            print("   wrote", p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
