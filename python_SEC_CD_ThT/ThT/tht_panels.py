"""tht_panels.py — drawing for the 25 uM ThT plateau comparison.

Two panels, one question:
    a  the kinetics, with the plateau window shaded and each condition's plateau
       level drawn across it, so the reader sees where the number comes from
    b  the plateau amplitude itself, every replicate shown, mean +- 1 SD

Two styles from the same drawing code:
    "talk"    matches the existing wt_0713.ipynb look (large bold Arial, 12 x 7.2 in)
    "nature"  180.0 mm wide, nothing above 7 pt except 8 pt bold panel letters,
              live vector text, exact page size (never bbox_inches="tight")

Every panel-scoped size lives on the Style object, so a caller can change any of
it without editing this file.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from tht_plateau import BASE_H, PLAT_H, CONDITIONS

MM = 1.0 / 25.4


@dataclass
class Style:
    name: str = "nature"
    fig_w_in: float = 180.0 * MM
    fig_h_in: float = 78.0 * MM
    margin_mm: float = 3.5
    axis_label: float = 7.0
    tick_label: float = 6.0
    legend: float = 6.0
    annot: float = 6.0
    panel_letter: float = 8.0
    letter_weight: str = "bold"
    label_weight: str = "normal"
    title: float = 0.0            # 0 disables titles (Nature captions carry the text)
    mean_lw: float = 1.4
    rep_lw: float = 0.5
    rep_alpha: float = 0.55
    band_alpha: float = 0.18
    spine_lw: float = 0.7
    tick_len: float = 2.4
    tick_w: float = 0.7
    marker_size: float = 16.0
    err_lw: float = 1.0
    cap_size: float = 2.2
    plateau_shade: str = "#000000"
    plateau_shade_alpha: float = 0.055
    grid: bool = False
    font: str = "Arial"
    width_ratios: tuple = (2.15, 1.0)
    left_pad: float = 0.055     # extra room for the y-axis label, figure fraction
    bottom_pad: float = 0.13    # extra room for the x-axis label
    right_pad: float = 0.0      # room for the plateau value labels, which sit OUTSIDE
                                # the axes; needed when the xy panel is the rightmost
    wspace: float = 0.30


TALK = Style(
    name="talk", fig_w_in=12.0, fig_h_in=6.4, margin_mm=6.0,
    axis_label=17.0, tick_label=14.0, legend=13.0, annot=13.0,
    panel_letter=20.0, label_weight="bold", title=18.0,
    mean_lw=3.0, rep_lw=1.3, rep_alpha=0.45, band_alpha=0.15,
    spine_lw=2.0, tick_len=8.0, tick_w=2.0,
    marker_size=70.0, err_lw=2.0, cap_size=5.0, width_ratios=(2.0, 1.0),
    left_pad=0.075, bottom_pad=0.135, wspace=0.24,
)
NATURE = Style()

#: single-column standalone, for using one panel on its own (88 mm = Nature single column)
SOLO = Style(name="solo", fig_w_in=88.0 * MM, fig_h_in=70.0 * MM,
             left_pad=0.125, bottom_pad=0.165, right_pad=0.025, wspace=0.0)


def apply_rc(st: Style):
    mpl.rcParams.update({
        "font.family": st.font,
        "font.size": st.tick_label,
        "svg.fonttype": "none",          # keep text live and editable
        "pdf.fonttype": 42,
        "axes.linewidth": st.spine_lw,
        "axes.unicode_minus": False,
    })


def _frame(ax, st: Style):
    for s in ax.spines.values():
        s.set_linewidth(st.spine_lw)
    ax.tick_params(which="major", length=st.tick_len, width=st.tick_w,
                   labelsize=st.tick_label)
    ax.tick_params(which="minor", length=st.tick_len * 0.55, width=st.tick_w * 0.7)


def panel_traces(ax, t, groups, st: Style, conditions=CONDITIONS,
                 show_plateau_lines=True, shade_window=True,
                 show_legend=True, show_window_note=True):
    """Panel a — smoothed kinetics, replicates, +-1 SD band, plateau level marked."""
    from tht_plateau import smooth, POLY_DISP

    t0 = t[-1] - PLAT_H
    if shade_window:
        ax.axvspan(t0, t[-1], color=st.plateau_shade,
                   alpha=st.plateau_shade_alpha, lw=0, zorder=0)

    for g, c in zip(groups, conditions):
        reps = g.traces
        mean_s = reps.mean(0)
        sd = reps.std(0, ddof=1) if reps.shape[0] > 1 else np.zeros_like(mean_s)
        mean_d = smooth(mean_s, poly=POLY_DISP)
        up = smooth(mean_s + sd, poly=POLY_DISP)
        lo = smooth(mean_s - sd, poly=POLY_DISP)

        ax.fill_between(t, lo, up, color=c.light, alpha=st.band_alpha, lw=0, zorder=1)
        for j in range(reps.shape[0]):
            ax.plot(t, reps[j], color=c.light, alpha=st.rep_alpha,
                    lw=st.rep_lw, zorder=2, solid_capstyle="butt")
        ax.plot(t, mean_d, color=c.dark, lw=st.mean_lw, zorder=4,
                label=f"{c.label} (n={reps.shape[0]})", solid_capstyle="round")

        if show_plateau_lines:
            p, _ = g.mean_sd("plateau")
            ax.plot([t0, t[-1]], [p, p], color=c.dark, lw=st.mean_lw * 0.85,
                    ls=(0, (2.6, 1.6)), zorder=5, solid_capstyle="butt")

    ax.set_xlabel("Time (h)", fontsize=st.axis_label, fontweight=st.label_weight)
    ax.set_ylabel("ThT fluorescence (mAU)", fontsize=st.axis_label,
                  fontweight=st.label_weight)
    ax.set_xlim(0, t[-1])
    ax.xaxis.set_major_locator(MultipleLocator(12))
    ax.xaxis.set_minor_locator(MultipleLocator(3))
    ax.set_ylim(bottom=0)
    _frame(ax, st)

    if show_legend:
        leg = ax.legend(fontsize=st.legend, frameon=False, loc="upper left",
                        handlelength=1.5, borderpad=0.2, labelspacing=0.35,
                        handletextpad=0.6)
        if st.label_weight == "bold":
            for tx in leg.get_texts():
                tx.set_weight("bold")

    if show_window_note:
        ax.annotate(f"plateau window\nlast {PLAT_H:g} h",
                    xy=((t0 + t[-1]) / 2, ax.get_ylim()[1] * 0.045),
                    ha="center", va="bottom", fontsize=st.annot * 0.92,
                    color="#555555", linespacing=1.25)
    return ax


def panel_plateau_xy(ax, t, groups, st: Style, conditions=CONDITIONS,
                     t_start=None, show_mean_level=True, annotate_drift=True,
                     show_values=True):
    """Panel — the plateau itself as an x-y trace: intensity vs time after the maximum.

    A bar chart collapses this region to one number and hides whether it is actually
    flat. Drawn out, the two wells that misbehave become visible rather than a
    footnote: A10 is still climbing at the end of the run, G12 decays through the
    window. Shared y with panel a, so the ~4x separation is read directly off the axis.

    t_start defaults to the earliest mean peak time across the conditions, i.e. the
    point from which every condition has reached its maximum.
    """
    if t_start is None:
        t_start = min(np.mean([w.t_peak for w in g.per_well]) for g in groups)
        t_start = float(np.floor(t_start))
    m = t >= t_start
    tp_ = t[m]

    for g, c in zip(groups, conditions):
        reps = g.traces[:, m]
        mean_s = reps.mean(0)
        sd = reps.std(0, ddof=1) if reps.shape[0] > 1 else np.zeros_like(mean_s)
        ax.fill_between(tp_, mean_s - sd, mean_s + sd, color=c.light,
                        alpha=st.band_alpha, lw=0, zorder=1)
        for j in range(reps.shape[0]):
            ax.plot(tp_, reps[j], color=c.light, alpha=st.rep_alpha,
                    lw=st.rep_lw, zorder=2, solid_capstyle="butt")
        ax.plot(tp_, mean_s, color=c.dark, lw=st.mean_lw, zorder=4,
                solid_capstyle="round")
        if show_mean_level:
            p, _ = g.mean_sd("plateau")
            ax.plot([tp_[0], tp_[-1]], [p, p], color=c.dark, lw=st.mean_lw * 0.8,
                    ls=(0, (2.6, 1.6)), zorder=5, solid_capstyle="butt")
            if show_values:
                ax.annotate(f"{p:.0f}", xy=(tp_[-1], p), xytext=(2.0, 0),
                            textcoords="offset points", ha="left", va="center",
                            fontsize=st.annot * 0.92, color=c.dark,
                            fontweight=st.label_weight)

    if annotate_drift:
        rising = [w.well for g in groups for w in g.per_well if w.still_rising]
        falling = [w.well for g in groups for w in g.per_well if w.drift_pct < -15.0]
        notes = []
        if rising:
            notes.append(f"{', '.join(rising)} still rising")
        if falling:
            notes.append(f"{', '.join(falling)} decaying")
        if notes:
            ax.text(0.02, 0.02, "\n".join(notes), transform=ax.transAxes,
                    ha="left", va="bottom", fontsize=st.annot * 0.9,
                    color="#555555", linespacing=1.25)

    ax.set_xlabel("Time (h)", fontsize=st.axis_label, fontweight=st.label_weight)
    ax.set_ylabel("ThT fluorescence (mAU)", fontsize=st.axis_label,
                  fontweight=st.label_weight)
    ax.set_xlim(tp_[0], tp_[-1])
    ax.set_ylim(bottom=0)
    ax.xaxis.set_major_locator(MultipleLocator(4))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    _frame(ax, st)
    return ax


def panel_plateau(ax, groups, st: Style, conditions=CONDITIONS, stats_=None,
                  metric="amplitude", show_points=True, brackets=None):
    """Panel b — the plateau number itself: every replicate, mean +- 1 SD."""
    xs = np.arange(len(groups), dtype=float)
    rng = np.random.default_rng(0)          # deterministic jitter
    top = 0.0

    for x, g, c in zip(xs, groups, conditions):
        # metric="as_measured" reads each species where it is defined: the plateau
        # for the monomer, t = 0 for the oligomer.  See Condition.measure_at.
        key = ("baseline" if getattr(g, "measure_at", "plateau") == "baseline"
               else "plateau") if metric == "as_measured" else metric
        vals = g._arr(key)
        m, s = g.mean_sd(key)
        top = max(top, float(vals.max()), m + s)

        ax.bar(x, m, width=0.56, color=c.light, edgecolor=c.dark,
               lw=st.spine_lw, zorder=2)
        ax.errorbar(x, m, yerr=s, fmt="none", ecolor=c.dark, elinewidth=st.err_lw,
                    capsize=st.cap_size, capthick=st.err_lw, zorder=3)
        if show_points:
            jitter = (rng.random(vals.size) - 0.5) * 0.20
            ax.scatter(x + jitter, vals, s=st.marker_size, facecolor="white",
                       edgecolor=c.dark, linewidth=st.spine_lw, zorder=4, clip_on=False)

    ax.set_xticks(xs)
    ax.set_xticklabels([(c.measured_label or c.label.replace("25 µM ", ""))
                        for c in conditions],
                       fontsize=st.tick_label, fontweight=st.label_weight)
    ax.set_xlim(-0.62, len(groups) - 0.38)
    ax.set_ylim(0, top * 1.22)
    ax.set_ylabel("ThT fluorescence (mAU)" if metric == "as_measured"
                  else ("Plateau amplitude (mAU)" if metric == "amplitude"
                        else f"Plateau {metric} (mAU)"),
                  fontsize=st.axis_label, fontweight=st.label_weight)
    _frame(ax, st)

    if brackets:
        # Three-state panel (fibril / oligomer / monomer): each bracket is computed
        # from the bars it spans, so the number printed always describes the two bars
        # under it. Stacked innermost first so they do not cross.
        import tht_bars as tb
        for k, (i, j) in enumerate(brackets):
            r = tb.pairwise(groups, i, j)
            y = top * (1.03 + 0.13 * k)
            lo, hi = sorted((xs[i], xs[j]))
            ax.plot([lo, lo, hi, hi], [y * 0.975, y, y, y * 0.975],
                    color="black", lw=st.spine_lw, zorder=5)
            ptxt = f"p = {r['p']:.3f}" if r["p"] >= 0.001 else "p < 0.001"
            ax.text((lo + hi) / 2.0, y * 1.01, f"{r['fold']:.1f}×  {ptxt}",
                    ha="center", va="bottom", fontsize=st.annot,
                    fontweight=st.label_weight)
        ax.set_ylim(0.0, top * (1.03 + 0.13 * len(brackets) + 0.09))
        return ax

    if stats_ is not None and len(groups) == 2:
        y = top * 1.06
        ax.plot([xs[0], xs[0], xs[1], xs[1]],
                [y * 0.965, y, y, y * 0.965], color="black", lw=st.spine_lw, zorder=5)
        # the p must belong to the quantity actually plotted: metric="as_measured"
        # compares fibril plateau against oligomer t = 0 on the log scale, which is
        # a different test from the amplitude comparison and a different p.
        p = {"as_measured": stats_.get("p_as_measured"),
             "plateau":     stats_.get("fold_plateau_p"),
             "amplitude":   stats_.get("p")}.get(metric) or stats_["p"]
        ptxt = f"p = {p:.3f}" if p >= 0.001 else "p < 0.001"
        # THE BRACKET MUST DESCRIBE THE BARS. A two-way test on "as_measured" let
        # metric="plateau" fall through to the amplitude statistics: bars showing
        # plateau (3.38x, Welch p = 0.057) captioned 4.6x / p = 0.046 -- a different
        # comparison, and across the 0.05 line. Select on the metric explicitly.
        fold = {"as_measured": stats_.get("fold_as_measured"),
                "plateau":     stats_.get("fold_plateau"),
                "amplitude":   stats_.get("fold_amplitude")}.get(
                    metric, stats_.get("fold_amplitude"))
        ax.text(xs.mean(), y * 1.02,
                f"{fold:.1f}×   {ptxt}",
                ha="center", va="bottom", fontsize=st.annot,
                fontweight=st.label_weight)
    return ax


#: what each panel key draws, and its natural width share
PANEL_KINDS = {
    "traces":  ("full kinetics, plateau window shaded",      2.15),
    "xy":      ("the plateau region as an x-y trace",        1.35),
    "bar":     ("plateau amplitude, replicates + mean ± SD", 1.00),
}
DEFAULT_PANELS = ("traces", "xy", "bar")


def build(t, groups, stats_, style: Style = NATURE, conditions=CONDITIONS,
          stem="ThT_25uM_plateau", outdir=None, formats=("svg", "pdf", "png"),
          letters="abcdefg", title=None, panels=DEFAULT_PANELS,
          width_ratios=None, xy_t_start=None, metric="amplitude", clean=False):
    """Draw the requested panels and write the files. Returns (fig, paths).

    panels selects and orders the panels from PANEL_KINDS; width_ratios overrides
    their default width shares. Passing panels=("traces", "bar") reproduces the
    original two-panel layout.
    """
    from tht_plateau import FIGURES
    outdir = FIGURES if outdir is None else outdir
    st = style
    apply_rc(st)

    bad = [p for p in panels if p not in PANEL_KINDS]
    if bad:
        raise ValueError(f"unknown panel(s) {bad}; choose from {list(PANEL_KINDS)}")
    if width_ratios is None:
        width_ratios = tuple(PANEL_KINDS[p][1] for p in panels)

    # The xy panel writes its plateau values just outside the right spine. If it is the
    # rightmost panel those labels are the last ink on the page, so reserve room for
    # them; otherwise the neighbouring panel's own gutter already covers it.
    auto_right_pad = st.right_pad
    if panels and panels[-1] == "xy":
        digits = 4  # widest plateau label, e.g. "1234"
        auto_right_pad = max(auto_right_pad,
                             (digits * st.annot * 0.62 / 72.0) / st.fig_w_in + 0.012)

    fig = plt.figure(figsize=(st.fig_w_in, st.fig_h_in))
    m = st.margin_mm * MM
    left = m / st.fig_w_in
    right = 1.0 - m / st.fig_w_in
    bottom = m / st.fig_h_in
    topfrac = 1.0 - (m + (4.0 * MM if title else 0.0)) / st.fig_h_in

    # The panel letter is the topmost ink on the page, so its BASELINE is placed from
    # the trim inward in millimetres: margin + cap height. Anchoring it to the axes
    # instead lets an 8 pt ascender overhang the margin (measured 1.29 mm against a
    # 3.3 mm budget before this was fixed).
    cap_mm = st.panel_letter * 0.716 / 72.0 * 25.4
    letter_y = 1.0 - (st.margin_mm + cap_mm) / (st.fig_h_in / MM)
    axes_top = letter_y - (2.0 * MM) / st.fig_h_in      # 2 mm clear beneath the letter

    gs = fig.add_gridspec(
        1, len(panels), width_ratios=width_ratios,
        left=left + st.left_pad, right=right - auto_right_pad,
        bottom=bottom + st.bottom_pad,
        top=min(axes_top, topfrac - 0.035), wspace=st.wspace,
    )
    axes = [fig.add_subplot(gs[0, i]) for i in range(len(panels))]

    # clean=True draws the data and the axes only: no legend, no in-panel notes, no
    # value or significance labels. For figures whose caption carries the text.
    for ax, kind in zip(axes, panels):
        if kind == "traces":
            panel_traces(ax, t, groups, st, conditions,
                         show_legend=not clean, show_window_note=not clean)
        elif kind == "xy":
            panel_plateau_xy(ax, t, groups, st, conditions, t_start=xy_t_start,
                             annotate_drift=not clean, show_values=not clean)
        elif kind == "bar":
            panel_plateau(ax, groups, st, conditions,
                          None if clean else stats_, metric=metric)

    # Letters anchored in FIGURE coordinates so they cannot drift with panel width.
    # The offset is the y-label gutter, which is wider for the first panel only if
    # it is the one carrying the tick labels with the most digits.
    for ax, ltr in zip(axes, letters):
        bb = ax.get_position()
        dx = 0.052 if ax is axes[0] else 0.046
        fig.text(bb.x0 - dx, letter_y, ltr,
                 fontsize=st.panel_letter, fontweight=st.letter_weight,
                 ha="left", va="baseline")

    if title and st.title:
        fig.text(0.5, 1.0 - m / st.fig_h_in, title, ha="center", va="top",
                 fontsize=st.title, fontweight="bold")

    paths = []
    for ext in formats:
        p = f"{outdir}/{stem}.{ext}"
        # No bbox_inches="tight": it re-crops and destroys the exact page width.
        fig.savefig(p, dpi=600 if ext == "png" else None,
                    facecolor="white", transparent=False)
        paths.append(p)
    return fig, paths
