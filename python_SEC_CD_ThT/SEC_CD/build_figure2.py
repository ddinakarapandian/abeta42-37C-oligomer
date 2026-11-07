"""build_figure2.py — Figure 2: oligomer formation kinetics depend on
concentration and temperature.  Consolidates manuscript Figures 3 and 4.

    a  250 µM, 37 °C     (registry 3A)   t = 2, 4, 6 h
    b  100 µM, 37 °C     (registry 3B)   t = 0, 2, 4, 8, 12, 24, 48 h
    c   25 µM, 37 °C     (registry 3C)   t = 0, 3, 6, 12, 24, 48, 72 h
    d  100 µM, 25 °C     (registry 4A)   t = 2, 4, 6, 8, 10, 12, 24, 48, 72 h

PANEL ORDER IS THE ARGUMENT.  a -> b -> c is a descending concentration series at
one temperature; d sits directly under b as its temperature control, same 100 µM.
The published pair ordered these 250 / 25 / 100 / 100@25, which splits the
concentration series across the diagonal and leaves the temperature comparison
without a neighbour to be compared against.

TWO BUILDS FROM ONE CODE PATH.  ``Figure_2`` is raw mAU for the main text;
``Figure_S1_normalized`` is the same four panels min-max normalised, for the SI.
Raw is the main-text choice because normalisation rescales every trace to its own
maximum, and at 250 µM the oligomer peak is genuinely TALLER than the monomer peak
-- the single most striking thing in the figure, and the thing normalisation
deletes.  Every trace within a panel shares a 500 µL loop, so raw is honest at the
panel level; only absolute height ACROSS panels reflects different loading, which
is what the concentration axis means.

Processing is sec3.py's formula throughout, via ``build_fig3_panels.load_csv``;
``verify_sec3_equivalence.py`` proves that path matches sec3 exactly.

STYLE matches Figure 1: ``figure_fig12._NAT_STYLE`` (7 pt type, 8 pt letters,
0.85 pt traces) with that module's measured margins, and the same one-rule panel
letter placement, so the two pages sit together without a visible seam.

    python build_figure2.py                # raw, main text
    python build_figure2.py --normalized   # the SI variant
    python build_figure2.py --verify       # font census of the delivered PDF
"""
from __future__ import annotations

import argparse
import os
import sys

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True

from dataclasses import replace

import numpy as np
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
from build_fig3_panels import (load_csv, peak, OLIG_WIN, MONO_WIN,
                               OLIG_DX_FRAC, MONO_DX_FRAC, OLIG_YFRAC, MONO_YFRAC)
from build_figure1_v2 import LETTER_DX, LETTER_DY

OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
XLIM = (6.0, 20.0)
MARGIN_FRAC = 0.02          # sec1.py's autoscale margin, the rule behind 3A's ylim

#: letter, registry id, title, raw ylim, y_major, y_minor, draw the two pointers
LAYOUT = [
    ("a", "3A", u"250 µM, 37 °C", (-0.41153078, 20.98806978), 4.0, 1.0, True),
    ("b", "3B", u"100 µM, 37 °C", (-0.25, 11.0), 2.0, 0.5, False),
    ("c", "3C", u"25 µM, 37 °C",  (-0.06, 2.75), 0.5, 0.1, False),
    ("d", "4A", u"100 µM, 25 °C", None, None, None, False),   # no published raw twin
]


def derived_ylim(traces, margin=MARGIN_FRAC):
    lo = min(A.min() for _V, A in traces)
    hi = max(A.max() for _V, A in traces)
    rng = hi - lo
    return (lo - margin * rng, hi + margin * rng)


def nice_step(span, target=6):
    raw = span / float(target)
    mag = 10.0 ** np.floor(np.log10(raw))
    for m in (1.0, 2.0, 2.5, 5.0, 10.0):
        if raw <= m * mag:
            return m * mag
    return 10.0 * mag


def check_loops(pid):
    """Within one panel every trace must share a loop, or raw mAU compares
    injected volumes rather than concentrations."""
    stated = set()
    for c in PANELS[pid]["csvs"]:
        b = os.path.basename(c).lower()
        if "100ul" in b:
            stated.add(100)
        elif "500ul" in b:
            stated.add(500)
    if len(stated) > 1:
        raise SystemExit(f"{pid}: mixed loops {sorted(stated)} — a raw panel across "
                         f"loops is a volume artefact.")
    return (sorted(stated)[0] if stated else None)


def build(stem=None, normalized=False, formats=("svg", "pdf", "png")):
    stem = stem or ("Figure_S1_normalized" if normalized else "Figure_2")

    avail_w = FIG_W - MARGIN["left"] - MARGIN["right"]
    panel_w = (avail_w - GAP_W) / 2.0
    row_h = panel_w / ROW_ASPECT[0]
    fig_h = (TRIM + MARGIN["top"] + row_h + GAP_H + row_h + MARGIN["bottom"])

    style = replace(_NAT_STYLE, figsize=(FIG_W, fig_h))
    fx = lambda v: v / FIG_W
    fy = lambda v: v / fig_h

    panels = []
    for ch, pid, title, ylim, ymaj, ymin, pointers in LAYOUT:
        loop = check_loops(pid)
        traces = [load_csv(f, normalize=normalized) for f in PANELS[pid]["csvs"]]
        labels = [l.replace("t = ", "") for l in PANELS[pid]["labels"]]
        if normalized:
            ylim, ymaj, ymin = (-0.05, 1.05), 0.2, 0.1
        elif ylim is None:
            ylim = derived_ylim(traces)
            ymaj = nice_step(ylim[1] - ylim[0]); ymin = ymaj / 4.0
        panels.append(dict(ch=ch, pid=pid, title=title, traces=traces, labels=labels,
                           ylim=ylim, ymaj=ymaj, ymin=ymin, pointers=pointers,
                           loop=loop))

    with gs.apply(style, **RC_EXTRA):
        fig = plt.figure(figsize=(FIG_W, fig_h), dpi=style.dpi)
        for i, p in enumerate(panels):
            r, c = divmod(i, 2)
            x = MARGIN["left"] + c * (panel_w + GAP_W)
            y = TRIM + MARGIN["top"] + r * (row_h + GAP_H)
            ax = fig.add_axes([fx(x), 1.0 - fy(y + row_h), fx(panel_w), fy(row_h)])

            colors = gs.ramp(len(p["traces"]), style)
            for (V, A), lab, col in zip(p["traces"], p["labels"], colors):
                ax.plot(V, A, color=col, lw=style.trace_linewidth, label=lab,
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
                    yy = max(yf * best[1], 0.0)
                    gs.annotate_peak(ax, best[0] + dx * span, yy, text,
                                     side=side, marker=marker, style=style)

            fig.text(fx(max(x - LETTER_DX, TRIM)), 1.0 - fy(y - LETTER_DY), p["ch"],
                     fontsize=style.panel_letter_size,
                     fontweight=style.panel_letter_weight, ha="left", va="baseline")

        paths = []
        for ext in formats:
            q = os.path.join(OUT, f"{stem}.{ext}")
            fig.savefig(q, dpi=600 if ext == "png" else None, facecolor="white")
            paths.append(q)
    plt.close(fig)

    print(f"=== {stem}   {FIG_W/MM:.4f} x {fig_h/MM:.4f} mm   "
          f"({'normalised' if normalized else 'raw mAU'})")
    for p in panels:
        mono = max(peak(V, A, *MONO_WIN)[1] for V, A in p["traces"])
        olig = max(peak(V, A, *OLIG_WIN)[1] for V, A in p["traces"])
        print(f"   {p['ch']}  {p['pid']}  {p['title']:16} {len(p['traces'])} traces  "
              f"loop {str(p['loop']) + ' µL' if p['loop'] else 'unstated':9} "
              f"ylim ({p['ylim'][0]:7.3f}, {p['ylim'][1]:7.3f})  "
              f"peak mono {mono:7.3f}  olig {olig:7.3f}")
    tier = 185.0 if fig_h/MM <= 185.0 else (210.0 if fig_h/MM <= 210.0 else 225.0)
    print(f"   page {fig_h/MM:.2f} mm -> {tier:.0f} mm tier, "
          f"{230.0 - fig_h/MM:.1f} mm of caption room")
    for q in paths:
        print("   wrote", q)
    return paths


def verify(stem="Figure_2"):
    import collections, fitz
    p = fitz.open(os.path.join(OUT, f"{stem}.pdf"))[0]
    sizes, letters = collections.Counter(), []
    for b in p.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                sizes[round(s["size"], 2)] += 1
                if s["text"].strip() in list("abcd"):
                    letters.append((s["text"].strip(), round(s["size"], 2),
                                    round(s["bbox"][0]/72*25.4, 2)))
    pix = p.get_pixmap(dpi=1200)
    a = np.frombuffer(pix.samples, np.uint8).reshape(pix.h, pix.w, pix.n)
    ink = a[:, :, :3].min(2) < 245
    r, c = np.where(ink); k = 25.4/1200
    print(f"\n{stem}: {p.rect.width/72*25.4:.4f} x {p.rect.height/72*25.4:.4f} mm")
    print("  sizes:", dict(sorted(sizes.items())),
          "| outside 5–7 pt (8 pt letters ok):",
          [s for s in sizes if (s > 7.01 and s != 8.0) or s < 5.0] or "none")
    print("  letters:", sorted(letters), )
    print(f"  margins L {c.min()*k:.3f} R {(pix.w-1-c.max())*k:.3f} "
          f"T {r.min()*k:.3f} B {(pix.h-1-r.max())*k:.3f} mm")
    print("  rasters:", len(p.get_images()))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--normalized", action="store_true", help="build the SI variant")
    ap.add_argument("--both", action="store_true", help="build main text and SI")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--stem", default=None)
    a = ap.parse_args(argv)
    todo = [False, True] if a.both else [a.normalized]
    for nrm in todo:
        stem = a.stem if (a.stem and not a.both) else None
        build(stem=stem, normalized=nrm)
        if a.verify:
            verify(stem or ("Figure_S1_normalized" if nrm else "Figure_2"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
