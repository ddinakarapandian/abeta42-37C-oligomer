# -*- coding: utf-8 -*-
"""build_figure_S2_tht.py — Supplementary Figure S2: ThT kinetics for 25 uM WT Abeta42
MONOMERS ALONE.

Requested by A.K.P.: the SI figure shows only the monomer trace.

WHY THIS IS THE RIGHT SPLIT. Figure 1e now carries the quantitative comparison as three
bars (fibril / oligomer / monomer, each read where its species is defined). What that
panel cannot show is WHERE the fibril number comes from -- a bar is one number, and the
reader has to take on trust that the monomer well actually converted and then held
steady. This figure is that evidence: the sigmoidal rise, the three replicates, and the
plateau window the fibril bar is averaged over, drawn across the trace.

WHAT WAS REMOVED, and why it is not a loss:
  * the oligomer trace       -- its value is a single t = 0 reading, which a kinetic
                                trace communicates no better than the bar in Fig. 1e
  * the plateau-window zoom  -- the window is now marked directly on the one trace
  * the three-bar comparison -- that IS Fig. 1e; repeating it in the SI duplicates a
                                main-text display item

WIDTH. 88 mm, Nature's single-column measure. One trace stretched to 180 mm would be
mostly empty; the other SI figures stay at 180 mm because they are multi-panel.

    python build_figure_S2_tht.py
"""
from __future__ import annotations

import os
import sys

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import tht_plateau as tp        # noqa: E402
import tht_panels as pan        # noqa: E402

OUT = os.path.join(HERE, "figures")
STEM = "Figure_S2_ThT_25uM"
MM = 1.0 / 25.4


def main():
    df, t = tp.load()
    groups = tp.analyse_all(df, t)
    mono, mono_c = groups[0], tp.CONDITIONS[0]

    st = pan.Style(name="s2", fig_w_in=88.0 * MM, fig_h_in=68.0 * MM,
                   left_pad=0.155, bottom_pad=0.155, right_pad=0.035)
    pan.apply_rc(st)
    fig = plt.figure(figsize=(st.fig_w_in, st.fig_h_in), dpi=600)
    ax = fig.add_axes([st.left_pad, st.bottom_pad,
                       1.0 - st.left_pad - st.right_pad,
                       0.985 - st.bottom_pad])

    # one condition, so the list passed in has to be a one-element list as well --
    # panel_traces zips groups against conditions.
    #
    # The built-in window note is centred on the shaded band and is wider than it, so
    # on a single-column panel it runs off the right-hand frame. Drawn here instead,
    # right-anchored to the axes, where it cannot overflow whatever the panel width.
    pan.panel_traces(ax, t, [mono], st, conditions=[mono_c],
                     show_legend=True, show_window_note=False)
    ax.text(0.975, 0.045, f"plateau window\nlast {tp.PLAT_H:g} h",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=st.annot * 0.92, color="#555555", linespacing=1.25)

    paths = []
    for ext in ("pdf", "svg", "png"):
        p = os.path.join(OUT, f"{STEM}.{ext}")
        fig.savefig(p, dpi=600 if ext == "png" else None, facecolor="white")
        paths.append(p)
    plt.close(fig)

    b, bsd = mono.mean_sd("baseline")
    pl, plsd = mono.mean_sd("plateau")
    print(f"=== {STEM}   {st.fig_w_in*25.4:.4f} x {st.fig_h_in*25.4:.4f} mm   "
          f"25 µM monomer only, n = {len(mono.per_well)}")
    print(f"   wells    {', '.join(w.well for w in mono.per_well)}")
    print(f"   t = 0    {b:7.1f} ± {bsd:5.1f} mAU   (the monomer bar in Fig. 1e)")
    print(f"   plateau  {pl:7.1f} ± {plsd:5.1f} mAU   (the fibril bar in Fig. 1e)")
    print(f"   window   last {tp.PLAT_H:g} h of a {t[-1]:.1f} h run")
    for w in mono.per_well:
        flag = "  still rising at the end -> plateau is a lower bound" if w.still_rising else ""
        print(f"     {w.well}  baseline {w.baseline:6.1f}  plateau {w.plateau:7.1f}"
              f"  peak {w.peak:7.1f} @ {w.t_peak:5.1f} h{flag}")
    for p in paths:
        print("   wrote", p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
