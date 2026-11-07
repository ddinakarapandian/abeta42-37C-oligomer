# -*- coding: utf-8 -*-
"""tht_bars.py — the three-state ThT comparison: fibril, oligomer, monomer.

Requested by A.K.P.: show the monomer alongside the fibril/oligomer pair, with each
species read where it is defined.

    fibril     the monomer well at its PLATEAU -- by then the well contains fibril
    oligomer   the oligomer well at t = 0     -- SEC-isolated oligomer is already
                                                 assembled when the plate starts, and
                                                 its signal only decays afterwards
    monomer    the monomer well at t = 0      -- the same well as the fibril bar,
                                                 read before it converts

Bars one and three are therefore THE SAME WELLS at opposite ends of the run. That is
the point: it makes the fibril bar an internal comparison against its own starting
material, so the difference cannot be attributed to well-to-well variation.

WHY THIS IS A SEPARATE CONDITION LIST from ``tht_plateau.CONDITIONS``. The kinetic
traces panel must draw the monomer well ONCE; the bar panel draws it twice. Sharing
one list would either duplicate a trace or lose a bar.

No data is re-read to build the third bar: every ``GroupStats`` already carries both
``baseline`` and ``plateau`` for each well, so the monomer-at-t=0 bar is the existing
monomer group read at its other end.
"""
from __future__ import annotations

import copy

import numpy as np

from tht_plateau import Condition, CONDITIONS

#: ColorBrewer single-hue pairs (light fill / dark edge), one hue per state. Greys for
#: the monomer because it is the inactive reference the other two are measured against.
BAR_CONDITIONS = (
    Condition("25 µM monomer", ("A10", "A11", "A12"), "#FCAE91", "#DE2D26",
              measure_at="plateau", measured_label="fibril"),
    Condition("25 µM oligomer", ("G10", "G11", "G12"), "#9ECAE1", "#3182BD",
              measure_at="baseline", measured_label="oligomer"),
    Condition("25 µM monomer", ("A10", "A11", "A12"), "#BDBDBD", "#636363",
              measure_at="baseline", measured_label="monomer"),
)


def read_key(g):
    """The attribute this group's species is defined by."""
    return "baseline" if getattr(g, "measure_at", "plateau") == "baseline" else "plateau"


def bar_groups(groups):
    """[fibril, oligomer, monomer] from the two analysed kinetic groups."""
    mono, olig = groups[0], groups[1]
    m0 = copy.deepcopy(mono)
    m0.measure_at = "baseline"
    m0.label = "25 µM monomer (t = 0)"
    return [mono, olig, m0]


def pairwise(groups, i, j):
    """Fold change and Welch p for two bars, each read where its species is defined.

    The test is on the LOG of the replicate values. A fold change is multiplicative,
    so the log scale is the one on which the difference is additive and the two
    variances are comparable; a t-test on the raw values would be testing a
    difference of means, which is not the quantity being quoted.
    """
    from scipy import stats
    a = groups[i]._arr(read_key(groups[i]))
    b = groups[j]._arr(read_key(groups[j]))
    la, lb = np.log(np.clip(a, 1e-9, None)), np.log(np.clip(b, 1e-9, None))
    tt = stats.ttest_ind(la, lb, equal_var=False)
    return dict(fold=float(a.mean() / b.mean()), p=float(tt.pvalue),
                a_mean=float(a.mean()), a_sd=float(a.std(ddof=1)),
                b_mean=float(b.mean()), b_sd=float(b.std(ddof=1)))


#: which comparisons to bracket, innermost first so the brackets stack without crossing
BRACKETS = ((1, 2), (0, 1))


def report(bg):
    out = []
    for i, (g, c) in enumerate(zip(bg, BAR_CONDITIONS)):
        k = read_key(g)
        m, s = g.mean_sd(k)
        out.append("   %-9s read at %-8s  %8.1f +- %6.1f mAU   (%s)"
                   % (c.measured_label, k, m, s, ", ".join(w.well for w in g.per_well)))
    for i, j in BRACKETS:
        r = pairwise(bg, i, j)
        out.append("   %-9s / %-9s = %5.1fx   p = %.4f (Welch on log)"
                   % (BAR_CONDITIONS[i].measured_label,
                      BAR_CONDITIONS[j].measured_label, r["fold"], r["p"]))
    return "\n".join(out)


if __name__ == "__main__":
    import tht_plateau as tp
    df, t = tp.load()
    print(report(bar_groups(tp.analyse_all(df, t))))
