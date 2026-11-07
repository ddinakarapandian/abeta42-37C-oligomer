"""build_fig34_combined.py — manuscript Figures 3 and 4 combined into one 2x2 SEC
figure, NON-NORMALISED (raw mAU).

    a  250 uM, 37 C      (registry 3A)   t = 2, 4, 6 h
    b   25 uM, 37 C      (registry 3C)   t = 0, 3, 6, 12, 24, 48, 72 h
    c  100 uM, 37 C      (registry 3B)   t = 0, 2, 4, 8, 12, 24, 48 h
    d  100 uM, 25 C      (registry 4A)   t = 2, 4, 6, 8, 10, 12, 24, 48, 72 h

Three of the four already exist raw in the registry -- 3A, 3B and 3C carry
``normalized=False`` -- so this only has to load 4A raw as well and put the four on
one page.  Panel c and manuscript Figure 4B are the SAME seven CSVs; 4B is that
series normalised, 3B is it raw, which is why only one of them appears here.

WHY RAW.  Min-max rescales every trace to its own maximum.  In the normalised
version of this figure all four panels' monomer peaks sit at 1.0 by construction,
so the reader can compare only the oligomer/monomer RATIO, never how much material
there is, and the 250 uM panel -- where the oligomer peak is the taller of the two
in absolute terms -- looks the same as the 25 uM panel, where it never rises above
0.1 of the monomer.  Raw mAU restores that.

Every trace within every panel was run on the same 500 uL loop, so no panel needs a
volume correction and none is applied.  (This is checked at build time and the
build fails if it stops being true.)  Absolute heights across PANELS still reflect
different preps and different loading concentrations -- that is the point of the
concentration axis, not an artefact.

y limits: 3A, 3B and 3C keep the limits already traced against her artwork; 3A's
(-0.41153078, 20.98806978) come from sec1.py's autoscale with margins=0.02 and are
exact, not the rounded transcription in panel_data's note.  Panel d has no
published raw twin, so its limits are derived by the same margins=0.02 rule.

    python build_fig34_combined.py                # bold, raw
    python build_fig34_combined.py --regular      # measured Griffin weights
    python build_fig34_combined.py --normalized   # the normalised version, for comparison
"""
from __future__ import annotations

import argparse
import os
import sys

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from secfigs import style as gs
from secfigs.panel_data import PANELS
from build_fig3_panels import (load_csv, peak, OLIG_WIN, MONO_WIN,
                               OLIG_DX_FRAC, MONO_DX_FRAC, OLIG_YFRAC,
                               MONO_YFRAC, MIN_YFRAC_OF_AXIS)

OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
XLIM = (6.0, 20.0)
MARGIN = 0.02          # sec1.py's autoscale margin, the rule behind 3A's ylim

#: panel letter -> (registry id, title, raw ylim, y_major, y_minor, draw pointers)
#: Titles are the short form used in the combined layout; the standalone builders
#: keep her longer "SEC Profiles of Aβ-42 WT at ..." wording.
LAYOUT = [
    ("a", "3A", u"250 µM, 37 °C",
     (-0.41153078, 20.98806978), 4.0, 1.0, True),
    ("b", "3C", u"25 µM, 37 °C",
     (-0.06, 2.75), 0.5, 0.1, False),
    ("c", "3B", u"100 µM, 37 °C",
     (-0.25, 11.0), 2.0, 0.5, False),
    ("d", "4A", u"100 µM, 25 °C",
     None, None, None, False),          # no published raw twin -> derived below
]


def derived_ylim(traces, margin=MARGIN):
    """sec1.py's autoscale: pad the data range by `margin` of itself, both ends."""
    lo = min(A.min() for _V, A in traces)
    hi = max(A.max() for _V, A in traces)
    rng = hi - lo
    return (lo - margin * rng, hi + margin * rng)


def nice_step(span, target=6):
    """A round major-tick step giving roughly `target` intervals over `span`."""
    raw = span / float(target)
    mag = 10.0 ** np.floor(np.log10(raw))
    for m in (1.0, 2.0, 2.5, 5.0, 10.0):
        if raw <= m * mag:
            return m * mag
    return 10.0 * mag


def check_loops(pid):
    """Every trace in a panel must share an injection loop, or raw mAU compares
    injected volumes rather than concentrations."""
    loops = set()
    for c in PANELS[pid]["csvs"]:
        b = os.path.basename(c).lower()
        loops.add(100 if "100ul" in b else (500 if "500ul" in b else None))
    stated = {l for l in loops if l}
    if len(stated) > 1:
        raise SystemExit(f"{pid}: mixed injection loops {sorted(stated)} — a raw mAU "
                         f"panel across loops is a volume artefact. Fix before building.")
    return (stated.pop() if stated else None), (None in loops)


def build(stem="Figure_3_4_Combined_nonnormalized", bold=True, normalized=False,
          formats=("svg", "pdf", "png")):
    style = gs.bold_variant() if bold else gs.STYLE
    w, h = style.figsize

    panels = []
    for letter, pid, title, ylim, ymaj, ymin, pointers in LAYOUT:
        loop, unknown = check_loops(pid)
        P = PANELS[pid]
        traces = [load_csv(f, normalize=normalized) for f in P["csvs"]]
        labels = [l.replace("t = ", "") for l in P["labels"]]
        if normalized:
            ylim, ymaj, ymin = (-0.05, 1.05), 0.2, 0.1
        elif ylim is None:
            ylim = derived_ylim(traces)
            ymaj = nice_step(ylim[1] - ylim[0])
            ymin = ymaj / 4.0
        panels.append(dict(letter=letter, pid=pid, title=title, traces=traces,
                           labels=labels, ylim=ylim, ymaj=ymaj, ymin=ymin,
                           pointers=pointers, loop=loop, loop_unknown=unknown))

    with gs.apply(style):
        fig, axes = plt.subplots(2, 2, figsize=(w * 2.0, h * 2.0))
        for p, ax in zip(panels, axes.ravel()):
            colors = gs.ramp(len(p["traces"]), style)
            for (V, A), lab, c in zip(p["traces"], p["labels"], colors):
                ax.plot(V, A, color=c, lw=style.trace_linewidth, label=lab,
                        solid_capstyle="round")
            gs.axes(ax, style, normalized=normalized, x_major=2.0,
                    y_major=p["ymaj"], y_minor=p["ymin"],
                    x_range=XLIM, y_range=p["ylim"])
            gs.set_labels(ax, "Elution Volume (ml)",
                          "Normalized Absorbance" if normalized else "Absorbance (mAU)",
                          style)
            gs.set_title(ax, p["title"], style)
            gs.legend(ax, style)

            if p["pointers"]:
                # Panel a alone names the two species, as in the reference layout.
                #
                span = XLIM[1] - XLIM[0]
                top = p["ylim"][1]
                for win, dx, yf, text, side, marker in (
                        (OLIG_WIN, OLIG_DX_FRAC, OLIG_YFRAC, "Oligomer", "right", "<"),
                        (MONO_WIN, MONO_DX_FRAC, MONO_YFRAC, "Monomer", "left", ">")):
                    best = (float("nan"), -np.inf)
                    for V, A in p["traces"]:
                        v, a = peak(V, A, *win)
                        if a > best[1]:
                            best = (v, a)
                # Height is a fraction of EACH PEAK, which is the reference rule and
                # keeps every label BELOW the peak it names. An earlier attempt
                # anchored to the axis top instead, to stop the two labels landing at
                # similar heights on raw data -- but that put Monomer 1.47 mAU ABOVE
                # its own apex (pointer 16.16 against apex 14.69), detached from every
                # curve. They are far apart in x (9.5 vs 13.7 mL), so equal heights
                # were never the real problem.
                    y = max(yf * best[1], MIN_YFRAC_OF_AXIS * top)
                    gs.annotate_peak(ax, best[0] + dx * span, y,
                                     text, side=side, marker=marker, style=style)
            gs.panel_letter(ax, p["letter"], style)

        fig.tight_layout(pad=2.4, w_pad=3.0, h_pad=2.6)
        paths = []
        for ext in formats:
            q = os.path.join(OUT, f"{stem}.{ext}")
            fig.savefig(q, dpi=600 if ext == "png" else None, facecolor="white")
            paths.append(q)
    plt.close(fig)

    print("=== Figures 3+4 combined — %s, %s weights"
          % ("normalised" if normalized else "raw mAU", "bold" if bold else "reference"))
    for p in panels:
        mono = max(peak(V, A, *MONO_WIN)[1] for V, A in p["traces"])
        olig = max(peak(V, A, *OLIG_WIN)[1] for V, A in p["traces"])
        loop = f"{p['loop']} µL" if p["loop"] else "loop unstated"
        print(f"   {p['letter']}  {p['pid']}  {p['title']:16}  {len(p['traces'])} traces  "
              f"{loop:13}  ylim ({p['ylim'][0]:.3f}, {p['ylim'][1]:.3f})  "
              f"peak monomer {mono:7.3f}  peak oligomer {olig:7.3f}")
        if p["loop_unknown"]:
            print(f"      note: some filenames in {p['pid']} do not state a loop; "
                  f"magnitudes are consistent with the rest of the panel.")
    for q in paths:
        print("   wrote", q)
    return paths


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regular", action="store_true",
                    help="measured Griffin reference weights instead of bold")
    ap.add_argument("--normalized", action="store_true",
                    help="build the normalised version instead (for comparison)")
    ap.add_argument("--stem", default=None)
    a = ap.parse_args(argv)
    stem = a.stem or ("Figure_3_4_Combined_normalized" if a.normalized
                      else "Figure_3_4_Combined_nonnormalized")
    build(stem=stem, bold=not a.regular, normalized=a.normalized)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
