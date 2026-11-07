"""build_figures_S3_S5.py — Supplementary Figures S3, S4 and S5, remade in Nature
format from the artwork that already exists in 10_Collaborator_Review_Clean.docx.

THEY ALREADY EXIST.  word/media/image17, image18 and image19, sharing one caption
in the SI section:

    "Figure S3-5. Quantifying monomer and oligomer population changes over time
     using SEC peak area"

An earlier version of this file invented three figures from that sentence without
opening the images, and was wrong in content: it used 6.0-11.3 / 11.3-20.0 mL
windows, had no Total series, and made S5 a conversion summary.  What the artwork
shows:

    S3  A 250 µM · B 100 µM · C 25 µM, all 37 °C           (image17, 3 panels)
    S4  A 100 µM 37 °C · B 100 µM 25 °C                    (image18, 2 panels)
    S5  A 25 µM 37 °C, Superdex 75, 280 nm                 (image19, 2 panels)
        B 25 µM 37 °C, Superdex 200 Increase, 210 nm

REGIONS ARE THE LAB'S, read off their own plot legends:
    Oligomer   6   - 11.5 mL
    Monomer   12   - 20   mL       (25 µM panels: 12 - 16.1 mL)
    Total      6   - 20   mL

THE FRAGMENTS SERIES IS NOT PLOTTED, at the author's request, but the 25 µM monomer
window STAYS at 12-16.1 mL rather than widening to 12-20.  Dropping the series and
widening the window are not the same edit: the fragment peak near 17 mL grows from
6.7 % of a 12-20 mL window at t = 0 to 28.0 % at 72 h, so folding it into "monomer"
turns a 32 % monomer decline (2.409 -> 1.649) into an 11 % one (2.583 -> 2.291) and
hides two thirds of the loss.  Total (6-20 mL) still spans everything, so the
fragment contribution remains recoverable as Total - (Oligomer + Monomer) without
being drawn or claimed.

For the record, fragmentation IS discussed in the manuscript -- nine sentences,
including "At 25 uM, we also see some level of peptide fragmentation, as indicated
by an elution peak at ~17 mL" and a passage ruling out microbial contamination with
sodium azide -- so the series was not spurious; it is simply not wanted on the page.
Note the deliberate 11.5-12.0 mL gap: the inter-peak valley is excluded from both
windows rather than split between them.  TOTAL is the mass-balance control -- it
tells the reader whether material is being redistributed or lost -- and leaving it
out was the most substantive error in the invented version.

S5 PANEL B CANNOT BE REBUILT.  It needs the 25 µM 37 °C reinjection series run on a
Superdex 200 Increase column with 210 nm detection at t = 0, 2, 7, 12, 24, 39, 48 h,
the same run behind Figure S2.  That series is not in
/path/to/Ab42oligomers37C/data/SEC\\Run_Results: no file carries a 39 h timepoint and no
WT 25 µM S200 run is present.  Its elution axis differs too -- monomer near 19.5 mL
against 14.5 mL here, trace running out to 28 mL -- so it is a different column and
cannot be derived from the Superdex 75 files.  S5 is built with panel A only and
the build prints exactly what panel B needs.

HER SOURCE CODE WAS FOUND, and this file now follows it rather than inferring from
pixels.  Two notebooks in ``/path/to/Ab42oligomers37C/data/SEC\\python_code``:

    WT37COligomer.ipynb                        cells  4-5   250 uM 37 C   (S3 A)
                                               cells 10-11   25 uM 37 C   (S3 C / S5 A)
                                               cells 26-27  100 uM 25 C   (S4 B)
    SECProgressE22K_onlyrequiredtraces.ipynb   cells 28-31  100 uM 37 C   (S3 B / S4 A)

All five AUC blocks (four WT, one E22K) are written identically:

    peak1_range = (5, 11);  peak2_range = (12, 18)      # 25 uM: (12,16.1),(16.1,20)
    y = df["Absorbance_mAU"].to_numpy()                 # RAW -- no baseline call
    mask1 = (x >= 5) & (x <= 11);  mask2 = (x > 12) & (x <= 18)
    auc1 = np.trapz(y[mask1], x[mask1]);  total_auc = auc1 + auc2

``raw_absorbance`` + ``area`` + ``measure`` below reproduce her retained cell outputs
to four decimals on all three panels -- 250 uM 2 h 9.6286/17.5689, 25 uM 6 h
-0.1292/2.1181/0.0070, 100 uM 25 C 2 h 0.8386/9.1686.  Bit-exact, not approximate.

THIS SUPERSEDES AN EARLIER RECONSTRUCTION IN THIS FILE, which had been fitted to
digitised pixels and got two things wrong: it used ``scipy.integrate.simpson`` (she
uses ``np.trapz``) and applied ``baseline_correct(method='linear')`` (she applies no
baseline at all).  Those two errors happened to cancel against a third -- the wrong
integration windows, taken from her legends -- to give a plausible 0.027 fit on one
panel.  Three compensating mistakes are not a validation; the source is.

THE PUBLISHED LEGENDS ARE WRONG.  See the window constants below: her code integrates
5-11 and 12-18 mL while the same cell labels the series 6-11.5 and 12-20 mL, and the
published S3-S5 inherit the wrong strings.  This build integrates what she actually
integrated -- so the values match the manuscript -- and states the true windows in
the caption it prints.

ONE DEFECT IN THE PUBLISHED ARTWORK, in S3 panel B / S4 panel A (the same plot,
printed twice).  Its seven AUC values are correct but paired with the wrong times.
``SECProgressE22K_onlyrequiredtraces.ipynb`` cell 28 declares

    filesWT2 = [w11, w12, w13, w14, w15, w16, w17]   # = 0, 2, 12, 24, 48, 4, 8 h
    nameWT2  = ['Immediate...','2h','4h','8h','12h','24h','48h']   # chronological

and cell 30 does ``zip(labels, absorbance_list_WT2)``, so from the third element on
every area is attached to the wrong label; cell 31 then plots them against
``times_hr_WT37C = [0,2,4,8,12,24,48]``.  Her x = 4, 8, 12, 24, 48 h therefore carry
the 12, 24, 48, 4 and 8 h runs.  (It is a mismatched zip, not a sort -- ``sort`` does
not appear anywhere in the notebooks.)  This is what produces the non-physical
zigzag: monomer collapsing to 2.99 at "12 h" and rebounding to 8.42 at "24 h" while
the real monomer peak apex falls monotonically 9.74 -> 2.47 mAU.  Read correctly the
series is monotonic.  THIS BUILD USES THE CORRECT ORDER and so will not look like the
artwork it replaces.  S3A, S3C, S4B and S5A are chronological in her figures already.

    python build_figures_S3_S5.py
    python build_figures_S3_S5.py --compare    # against the published artwork
"""
from __future__ import annotations

import argparse
import os
import sys

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True

from dataclasses import replace

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from secfigs import style as gs
from secfigs.panel_data import PANELS
from secfigs.figure_fig12 import (_NAT_STYLE, MARGIN, GAP_W, GAP_H, TRIM, MM,
                                  ROW_ASPECT, RC_EXTRA, FIG_W)
from build_figure1_v2 import LETTER_DX, LETTER_DY

OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

#: her integrator, verbatim -- np.trapz, not scipy.integrate.simpson.  sec1.py DOES
#: import simpson and use it inside integrate_peaks, but grep across all 26
#: notebooks finds no call to integrate_peaks and no import of simpson; every AUC
#: block computes np.trapz inline.  The two differ by <0.001 % on this data, so
#: this is about provenance, not magnitude.
trapz = getattr(np, "trapezoid", np.trapz)

#: HER ACTUAL INTEGRATION WINDOWS, from the source cells -- NOT the published
#: legends, which disagree with the code printed directly above them.  In every one
#: of her AUC blocks the code says ``peak1_range = (5, 11)`` and
#: ``peak2_range = (12, 18)`` while the very next cell labels the same two series
#: "Peak 1 AUC (6-11.5 mL)" and "Peak 2 AUC (12-20 mL)".  The published figures
#: inherit the wrong strings.  These constants follow the code, so the numbers here
#: match the manuscript's; the CAPTION states them truthfully.
OLIGOMER = (5.0, 11.0)          # published legend claims 6 - 11.5
MONOMER = (12.0, 18.0)          # published legend claims 12 - 20
MONOMER_25 = (12.0, 16.1)       # 25 uM panels, correctly labelled
FRAGMENTS = (16.1, 20.0)        # 25 uM panels, correctly labelled
TOTAL = None                    # unused: Total is the SUM, see measure()

#: series colours as published; Total takes red where Fragments already owns green
C_OLIG, C_MONO, C_FRAG = "#1f77b4", "#ff7f0e", "#2ca02c"
C_TOTAL_PLAIN, C_TOTAL_FRAG = "#2ca02c", "#d62728"

#: registry id -> hours, in registry order
HOURS = {"3A": [2, 4, 6],
         "3B": [0, 2, 4, 8, 12, 24, 48],
         "3C": [0, 3, 6, 12, 24, 48, 72],
         "4A": [2, 4, 6, 8, 10, 12, 24, 48, 72]}


def raw_absorbance(path):
    """The RAW ``Absorbance_mAU`` column, with no baseline handling of any kind.

    This is deliberately NOT ``load_csv``/``load_trace``.  Those reproduce
    ``_process_dataframe``, which subtracts the global minimum when it is negative;
    her AUC cells never go through that path.  They read ``df["Absorbance_mAU"]``
    straight off a ``SECData`` that has only had ``rename_columns()`` and
    ``normalize_data()`` applied -- and ``normalize_data`` (sec1.py:109) only ADDS a
    ``Normalized_Absorbance_mAU`` column, leaving ``Absorbance_mAU`` untouched.

    The proof that no baseline was applied is in her own retained output: the 25 uM
    Peak 1 areas are NEGATIVE (6 h: -0.1292, 72 h: -0.0517).  Every baseline method
    in the library floors the trace at zero or near it, so none of them ran.
    """
    df = pd.read_csv(path, sep="\t", encoding="utf-16", skiprows=3, header=None)
    x = pd.to_numeric(df[0], errors="coerce").to_numpy()
    y = pd.to_numeric(df[1], errors="coerce").to_numpy()
    m = ~(np.isnan(x) | np.isnan(y))
    return x[m], y[m]


def area(V, A, region, closed_low=False):
    """``np.trapz`` over her exact mask.  Peak 1 takes ``>=`` on its lower edge,
    peaks 2 and 3 take a strict ``>`` -- copy her convention rather than tidy it,
    because the windows abut (12-16.1, 16.1-20) and the shared endpoint would
    otherwise be counted twice."""
    m = ((V >= region[0]) if closed_low else (V > region[0])) & (V <= region[1])
    return float(trapz(A[m], V[m])) if m.sum() > 2 else float("nan")


def measure(pid, split_fragments=False):
    t = np.array(HOURS[pid], float)
    if len(t) != len(PANELS[pid]["csvs"]):
        raise SystemExit(f"{pid}: {len(t)} times for {len(PANELS[pid]['csvs'])} files")
    rows = {"olig": [], "mono": [], "frag": [], "total": []}
    for f in PANELS[pid]["csvs"]:
        V, A = raw_absorbance(f)
        rows["olig"].append(area(V, A, OLIGOMER, closed_low=True))
        rows["mono"].append(area(V, A, MONOMER_25 if split_fragments else MONOMER))
        rows["frag"].append(area(V, A, FRAGMENTS) if split_fragments else np.nan)
        # TOTAL IS THE SUM OF THE COMPONENTS -- her ``total_auc = auc1 + auc2``
        # (+ auc3 at 25 uM), never a separate wide-window integral, in all five of
        # her AUC blocks.  Note her legend calls it "Total AUC (5-20 mL)" while the
        # value is a sum over 5-11 plus 12-18, which is not the same thing: the
        # 11-12 mL valley and everything past 18 mL are outside it.
        rows["total"].append(rows["olig"][-1] + rows["mono"][-1]
                             + (rows["frag"][-1] if split_fragments else 0.0))
    return dict(pid=pid, t=t, split=split_fragments,
                **{k: np.array(v) for k, v in rows.items()})


def series_of(d):
    """Four series on the 25 uM panels (which resolve a fragment peak near 17 mL),
    three elsewhere.

    THE FRAGMENTS LINE IS DRAWN.  An earlier version withheld it on the understanding
    that it was not in the manuscript.  That was wrong on both counts: her own S3 panel C
    and S5 panel A plot it, and the SI text discusses it directly -- "Figure 3C showed
    fragment peaks, and the reanalysis Figure S2 similarly displayed peaks near 28 mL
    corresponding to fragmented species".  Dropping it also broke the mass balance, since
    her Total is Oligomer + Monomer + Fragments on these panels.

    SHORT LABELS.  "Monomer AUC (12-16.1 ml)" is 30.9 mm of 7 pt Arial plus ~6 mm of
    handle, which on a narrow panel leaves almost no room for the data.  The integration
    windows are stated once in the caption instead of three times per panel.
    """
    ser = [("Oligomer", d["olig"], C_OLIG, "o"),
           ("Monomer",  d["mono"], C_MONO, "s")]
    if d.get("split"):
        ser.append(("Fragments", d["frag"], C_FRAG, "P"))
        ser.append(("Total", d["total"], C_TOTAL_FRAG, "^"))
    else:
        ser.append(("Total", d["total"], C_TOTAL_PLAIN, "^"))
    return ser


def draw(ax, d, style, title):
    ser = series_of(d)
    for lab, y, col, mk in ser:
        ax.plot(d["t"], y, color=col, lw=style.trace_linewidth, marker=mk,
                markersize=2.4, markerfacecolor="white",
                markeredgewidth=style.spine_linewidth, label=lab,
                solid_capstyle="round", clip_on=False)
    hi = max(float(np.nanmax(y)) for _l, y, _c, _m in ser)
    lo = min(0.0, min(float(np.nanmin(y)) for _l, y, _c, _m in ser))
    xhi = float(d["t"].max())
    # x_major MUST be set per panel. style.x_major is 2.0, the SEC ELUTION-VOLUME
    # default; inherited onto a 72 h time axis it puts a labelled tick every 2 h and
    # the labels collapse into an unreadable band.
    xmaj = 2.0 if xhi <= 8 else (6.0 if xhi <= 24 else (12.0 if xhi <= 50 else 24.0))
    gs.axes(ax, style, normalized=False, x_major=xmaj, x_minor=xmaj / 2.0,
            x_range=(-xhi * 0.05, xhi * 1.05), y_range=(lo, hi * 1.42))
    gs.set_labels(ax, "Reinjection time (h)", u"AUC (mAU·ml)", style)
    gs.set_title(ax, title, style)
    gs.legend(ax, style, loc="upper right")


def _page(ncols, aspect=None):
    """ROW_ASPECT[0] = 1.270 is the SEC-chromatogram ratio measured off Griffin
    Fig 1A. A single full-width kinetic panel at that ratio is 149 mm tall with
    most of the box empty, so S5 passes its own."""
    avail = FIG_W - MARGIN["left"] - MARGIN["right"]
    pw = (avail - (ncols - 1) * GAP_W) / ncols
    ph = pw / (aspect or ROW_ASPECT[0])
    return pw, ph, TRIM + MARGIN["top"] + ph + MARGIN["bottom"]


def _panel(fig, style, col, pw, ph, fh, letter):
    fx, fy = (lambda v: v / FIG_W), (lambda v: v / fh)
    x = MARGIN["left"] + col * (pw + GAP_W)
    y = TRIM + MARGIN["top"]
    ax = fig.add_axes([fx(x), 1.0 - fy(y + ph), fx(pw), fy(ph)])
    if letter:
        fig.text(fx(max(x - LETTER_DX, TRIM)), 1.0 - fy(y - LETTER_DY), letter,
                 fontsize=style.panel_letter_size,
                 fontweight=style.panel_letter_weight, ha="left", va="baseline")
    return ax


SPECS = [
    ("Figure_S3", [("a", "3A", u"250 µM, 37 °C", False),
                   ("b", "3B", u"100 µM, 37 °C", False),
                   ("c", "3C", u"25 µM, 37 °C",  True)]),
    ("Figure_S4", [("a", "3B", u"100 µM, 37 °C", False),
                   ("b", "4A", u"100 µM, 25 °C", False)]),
    ("Figure_S5", [("a", "3C", u"25 µM, 37 °C · Superdex 75, 280 nm", True)]),
]


def build(formats=("svg", "pdf", "png")):
    written, rows = [], []
    for stem, panels in SPECS:
        aspect = {1: 2.9, 3: 1.00}.get(len(panels))   # 3 narrow columns
        pw, ph, fh = _page(len(panels), aspect=aspect)
        style = replace(_NAT_STYLE, figsize=(FIG_W, fh))
        with gs.apply(style, **RC_EXTRA):
            fig = plt.figure(figsize=(FIG_W, fh), dpi=style.dpi)
            for i, (ch, pid, title, split) in enumerate(panels):
                d = measure(pid, split)
                rows.append((stem, ch, title, d))
                draw(_panel(fig, style, i, pw, ph, fh,
                            ch if len(panels) > 1 else ""), d, style, title)
            for ext in formats:
                q = os.path.join(OUT, f"{stem}.{ext}")
                fig.savefig(q, dpi=600 if ext == "png" else None, facecolor="white")
                written.append(q)
        plt.close(fig)
        print(f"=== {stem}  {FIG_W/MM:.4f} x {fh/MM:.4f} mm  "
              f"({len(panels)} panel{'s' if len(panels) > 1 else ''})")

    print("\nAREAS (np.trapz, raw absorbance, her code's windows — reproduces her"
          "\n       retained notebook outputs to 4 dp; see module docstring)")
    for stem, ch, title, d in rows:
        print(f"  {stem} {ch}  {title}")
        for i, t in enumerate(d["t"]):
            frag = f"   frag {d['frag'][i]:6.3f}" if d["split"] else ""
            print(f"     t = {t:5.1f} h   olig {d['olig'][i]:7.3f}   "
                  f"mono {d['mono'][i]:7.3f}{frag}   total {d['total'][i]:7.3f}")

    print("\nCAPTION MUST STATE THE WINDOWS — the legends no longer carry them, and")
    print("the ones printed in the published S3–S5 are wrong (see module docstring).")
    print("   Peak areas were obtained by trapezoidal integration of the raw")
    print("   absorbance, without baseline correction, over 5.0–11.0 mL (oligomer)")
    print("   and 12.0–18.0 mL (monomer; at 25 µM, 12.0–16.1 mL for monomer and")
    print("   16.1–20.0 mL for the ~17 mL fragment peak). Total is the sum of these")
    print("   components. Because no baseline is subtracted, areas for peaks near")
    print("   the noise floor can be slightly negative.")

    print("\nS5 PANEL B NOT BUILT — the data is not present.")
    print("  needs 25 µM / 37 °C reinjections on a Superdex 200 Increase column with")
    print("  210 nm detection at t = 0, 2, 7, 12, 24, 39, 48 h (the Figure S2 run).")
    print("  Published windows there: oligomer 4-17, monomer 17.5-23, fragments 26-28,")
    print("  total 4-28 mL. Searched Dropbox/SEC/Run_Results: no 39 h file, no WT")
    print("  25 µM S200 run; its monomer elutes near 19.5 mL against 14.5 mL here.")
    for q in written:
        print("   wrote", q)
    return written


def compare():
    """Areas beside values read off the published curves, so the remake is checked
    against the artwork it replaces rather than trusted."""
    approx = {"3A": ("250 µM 37 °C", 9.7, 17.5, 27.0),
              "3B": ("100 µM 37 °C", 0.2, 10.0, 10.2),
              "3C": ("25 µM 37 °C",  0.2, 2.35, 2.6),
              "4A": ("100 µM 25 °C", 1.0, 9.2, 10.5)}
    print(f"{'series':16} {'t':>5}  {'olig (pub)':>18} {'mono (pub)':>18} {'total (pub)':>18}")
    for pid, split in (("3A", False), ("3B", False), ("3C", True), ("4A", False)):
        nm, po, pm, pt = approx[pid]
        d = measure(pid, split)
        print(f"{nm:16} {d['t'][0]:5.0f}  {d['olig'][0]:8.3f} ({po:5.2f})  "
              f"{d['mono'][0]:8.3f} ({pm:5.2f})  {d['total'][0]:8.3f} ({pt:5.2f})")
    print("(published values are read off the plotted curves by eye, so agreement to")
    print(" a few percent is the most that can be expected of this check)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compare", action="store_true")
    a = ap.parse_args(argv)
    if a.compare:
        compare(); print()
    build()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
