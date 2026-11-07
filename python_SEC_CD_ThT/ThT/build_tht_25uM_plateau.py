"""build_tht_25uM_plateau.py — make the 25 uM ThT plateau panel.

    python build_tht_25uM_plateau.py              # both styles
    python build_tht_25uM_plateau.py --style talk # presentation only
    python build_tht_25uM_plateau.py --plateau-h 6 --base-h 1

Writes into ./figures and prints the numbers behind every mark.
Nothing under Dropbox is ever modified.
"""
from __future__ import annotations

import argparse
import os
import sys

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True

import matplotlib
matplotlib.use("Agg")

import tht_plateau as tp
import tht_panels as pan


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--style", choices=("nature", "talk", "both"), default="both")
    ap.add_argument("--base-h", type=float, default=tp.BASE_H,
                    help=f"baseline window in hours from the start (default {tp.BASE_H:g})")
    ap.add_argument("--plateau-h", type=float, default=tp.PLAT_H,
                    help=f"plateau window in hours back from the end (default {tp.PLAT_H:g})")
    ap.add_argument("--metric", choices=("amplitude", "plateau"), default="amplitude",
                    help="what panel b plots (default amplitude = plateau - baseline)")
    ap.add_argument("--stem", default="ThT_25uM_plateau")
    ap.add_argument("--panels", default=",".join(pan.DEFAULT_PANELS),
                    help="comma-separated, from " + ",".join(pan.PANEL_KINDS)
                         + f" (default {','.join(pan.DEFAULT_PANELS)})")
    ap.add_argument("--xy-start", type=float, default=None,
                    help="hour the x-y plateau panel starts at "
                         "(default: the earliest mean peak time across conditions)")
    a = ap.parse_args(argv)
    panels = tuple(p.strip() for p in a.panels.split(",") if p.strip())

    df, t = tp.load()
    print(f"loaded {tp.data_path()}")
    print(f"{len(t)} timepoints, 0 -> {t[-1]:.2f} h, "
          f"{(t[1]-t[0])*60:.0f} min spacing\n")

    groups = tp.analyse_all(df, t, base_h=a.base_h, plat_h=a.plateau_h)
    stats_ = tp.fold_difference(groups)
    print(tp.report(groups, stats_))

    rising = [w.well for g in groups for w in g.per_well if w.still_rising]
    drifting = [w.well for g in groups for w in g.per_well if w.drift_pct < -15.0]
    if rising:
        print(f"\nNOTE  {', '.join(rising)} peak inside the plateau window — still "
              f"rising at the end of the run, so the plateau is a LOWER BOUND there.")
    if drifting:
        print(f"NOTE  {', '.join(drifting)} fall by more than 15% across the plateau "
              f"window, which pulls that condition's mean down.")

    styles = {"nature": pan.NATURE, "talk": pan.TALK}
    chosen = list(styles) if a.style == "both" else [a.style]
    for name in chosen:
        st = styles[name]
        stem = a.stem if name == "nature" else f"{a.stem}_talk"
        title = (None if name == "nature" else
                 "ThT: WT Aβ42 25 µM — monomer vs oligomer (no seeding, 37 °C shaking)")
        fig, paths = pan.build(t, groups, stats_, style=st, stem=stem, title=title,
                               panels=panels, xy_t_start=a.xy_start, metric=a.metric)
        w, h = fig.get_size_inches()
        print(f"\n[{name}] {w*25.4:.4f} x {h*25.4:.4f} mm")
        for p in paths:
            print("   wrote", p)
        fig.clf()

    # Each panel also standalone, at Nature single-column width, so any one of them
    # can be dropped in on its own or reused in a talk.
    for kind in panels:
        fig, paths = pan.build(t, groups, stats_, style=pan.SOLO, panels=(kind,),
                               letters=(""), stem=f"{a.stem}_{kind}",
                               xy_t_start=a.xy_start, metric=a.metric)
        w, h = fig.get_size_inches()
        print(f"\n[solo:{kind}] {w*25.4:.4f} x {h*25.4:.4f} mm")
        for p in paths:
            print("   wrote", p)
        fig.clf()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
