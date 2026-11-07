"""tht_plateau.py — 25 uM ThT plateau comparison, WT Abeta42 monomer vs oligomer.

The question this answers: after the ThT signal has reached its maximum, how much
fluorescence is there?  Monomer converts to fibril and ends ThT-bright; oligomer
starts ThT-positive but never gets there.  The plateau level is the observable.

Data
----
/path/to/Ab42oligomers37C/data/ThT/ForTerry/Abeta0708_0711.xlsx  (READ ONLY)
    Sheet1, merged 24 h + 40 h runs, 250 points over 62.48 h at 15 min spacing.
    Column Elapsed_Time_s holds seconds; well columns are plate coordinates.

Plate map, taken verbatim from ThT_plot0708_37C_Shaking.ipynb (WT, 37 C, shaking,
20 mM NaPi).  Only the 25 uM no-seed arms are used here:
    25 uM monomer  no seed   A10 A11 A12
    25 uM oligomer no seed   G10 G11 G12
(For reference, not plotted: monomer 10% seed A1-A3; oligomer 10% seed G1-G3;
 monomer 10% seed with NaCl-prepped seeds G4-G6; 50 uM arms B10-B12 / H10-H12.)

Processing matches the existing wt_0713.ipynb so the traces are the same curves:
Savitzky-Golay, window 15, polyorder 3 per replicate, polyorder 5 for the display
mean and band.  No baseline subtraction is applied to the traces; the plateau
AMPLITUDE is baseline-subtracted, and both are reported.

Definitions, fixed and identical for both conditions so the comparison is
like-for-like:
    baseline  = median of the smoothed trace over the first BASE_H hours
    plateau   = mean of the smoothed trace over the last PLAT_H hours
    amplitude = plateau - baseline

A fixed late window is used rather than "everything after each well's own
maximum" because two wells make that definition misbehave: A10's maximum falls on
the final timepoint (it is still rising), and G12 decays 45% after its maximum.
A common window keeps one rule for every well and lets those two behaviours show
up as drift, which is reported rather than smoothed away.

Colours follow the 25 uM entries of the existing scheme (monomer #FCAE91,
oligomer #9ECAE1) for replicates and bands, with the next darker step of the same
ColorBrewer ramp for the mean lines so they stay legible.
"""
from __future__ import annotations

import os
import pathlib
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")

# --------------------------------------------------------------------------- data

DATA_CANDIDATES = (
    r"/path/to/Ab42oligomers37C/data/ThT/ForTerry/Abeta0708_0711.xlsx",
)
HERE = pathlib.Path(__file__).resolve().parent
FIGURES = HERE / "figures"
FIGURES.mkdir(exist_ok=True)

SMOOTH_WINDOW, POLY_REP, POLY_DISP = 15, 3, 5
BASE_H = 2.0     # baseline window, hours from the start
PLAT_H = 10.0    # plateau window, hours back from the end


@dataclass(frozen=True)
class Condition:
    label: str
    wells: tuple
    light: str          # replicates + SD band
    dark: str           # mean line
    measure_at: str = "plateau"
    measured_label: str = ""
    #: what the sample IS at the point it is measured.  "25 uM monomer" is what was
    #: added; by the plateau it is fibril, so that is what the plateau bar is named.
    #: The kinetic legend keeps the sample name, because at t = 0 the red trace
    #: really is monomer -- calling the whole trace "fibril" would be wrong at the
    #: left-hand edge, which is exactly where it is ThT-dark.
    #: WHERE THIS SPECIES' ThT INTENSITY IS READ, and the two are not the same
    #: question.  A monomer sample is ThT-dark until it converts, so its intensity
    #: is the PLATEAU -- what the fibril reaches.  An isolated oligomer is already
    #: ThT-positive when the run starts and then slowly converts too, so its
    #: plateau describes whatever it has become by 62 h, not the oligomer.  The
    #: oligomer as prepared is the t = 0 reading.  Comparing oligomer plateau to
    #: monomer plateau (3.38x; 4.6x is the AMPLITUDE ratio, a different quantity)
#: compares two end states; comparing oligomer t = 0 to
    #: monomer plateau (11.9x) compares the oligomer to the fibril, which is the
    #: question the figure is asking.


CONDITIONS = (
    Condition("25 µM monomer",  ("A10", "A11", "A12"), "#FCAE91", "#DE2D26",
              measure_at="plateau",  measured_label="fibril"),
    Condition("25 µM oligomer", ("G10", "G11", "G12"), "#9ECAE1", "#3182BD",
              measure_at="baseline", measured_label="oligomer"),
)

#: the 10 % self-seeded arms of the same plate, for the seeding question
SEEDED = (
    Condition("25 µM monomer + seed",  ("A1", "A2", "A3"), "#FCAE91", "#DE2D26",
              measure_at="plateau",  measured_label="fibril + seed"),
    Condition("25 µM oligomer + seed", ("G1", "G2", "G3"), "#9ECAE1", "#3182BD",
              measure_at="baseline", measured_label="oligomer + seed"),
)


def data_path() -> pathlib.Path:
    override = os.environ.get("THT_XLSX")
    if override:
        return pathlib.Path(override)
    for c in DATA_CANDIDATES:
        if pathlib.Path(c).is_file():
            return pathlib.Path(c)
    raise FileNotFoundError(
        "Abeta0708_0711.xlsx not found. Tried: " + "; ".join(DATA_CANDIDATES)
        + ". Set THT_XLSX to point at it."
    )


def load():
    """Return (dataframe, time in hours). The source workbook is never written to."""
    df = pd.read_excel(data_path(), sheet_name="Sheet1", engine="openpyxl")
    t = df["Elapsed_Time_s"].to_numpy(float) / 3600.0
    return df, t


def _odd(n: int) -> int:
    return n if n % 2 else max(1, n - 1)


def smooth(y, window: int = SMOOTH_WINDOW, poly: int = POLY_REP):
    """Savitzky-Golay with clamped window/poly; returns the input if too short."""
    y = np.asarray(y, float)
    n = y.size
    if n < 3:
        return y
    w = min(_odd(window), _odd(n))
    return savgol_filter(y, w, min(poly, w - 1)) if w >= 3 else y


# ----------------------------------------------------------------------- analysis

@dataclass
class WellStats:
    well: str
    baseline: float
    plateau: float
    amplitude: float
    peak: float
    t_peak: float
    drift_pct: float      # change across the plateau window, % of its first value
    still_rising: bool    # peak lands in the plateau window -> plateau is a lower bound


@dataclass
class GroupStats:
    label: str
    wells: list = field(default_factory=list)
    traces: np.ndarray = None      # (n_reps, n_time) smoothed
    per_well: list = field(default_factory=list)

    def _arr(self, attr):
        return np.array([getattr(w, attr) for w in self.per_well], float)

    def mean_sd(self, attr):
        a = self._arr(attr)
        return a.mean(), (a.std(ddof=1) if a.size > 1 else 0.0)

    #: the reading this condition's species is defined by -- see Condition.measure_at
    measure_at: str = "plateau"

    def as_measured(self):
        """(mean, sd) of the reading appropriate to THIS species: the plateau for a
        monomer that has to convert first, the t = 0 baseline for an oligomer that
        is already ThT-positive when the run starts."""
        return self.mean_sd("baseline" if self.measure_at == "baseline" else "plateau")


def analyse(df, t, cond: Condition, base_h=BASE_H, plat_h=PLAT_H) -> GroupStats:
    wells = [w for w in cond.wells if w in df.columns]
    missing = [w for w in cond.wells if w not in df.columns]
    if missing:
        raise KeyError(f"{cond.label}: wells not in the workbook: {missing}")

    base_m = t <= base_h
    plat_m = t >= t[-1] - plat_h
    traces, per_well = [], []
    for w in wells:
        y = smooth(df[w].to_numpy(float))
        traces.append(y)
        b = float(np.median(y[base_m]))
        p = float(y[plat_m].mean())
        i = int(np.argmax(y))
        seg = y[plat_m]
        per_well.append(WellStats(
            well=w, baseline=b, plateau=p, amplitude=p - b,
            peak=float(y[i]), t_peak=float(t[i]),
            drift_pct=float((seg[-1] - seg[0]) / seg[0] * 100.0),
            still_rising=bool(t[i] >= t[-1] - plat_h),
        ))
    return GroupStats(label=cond.label, wells=wells, measure_at=cond.measure_at,
                      traces=np.vstack(traces), per_well=per_well)


def analyse_all(df, t, conditions=CONDITIONS, **kw):
    return [analyse(df, t, c, **kw) for c in conditions]


def fold_difference(groups):
    """Ratio of the first group's mean amplitude to the second's, plus a Welch test."""
    from scipy import stats
    a = groups[0]._arr("amplitude")
    b = groups[1]._arr("amplitude")
    tt = stats.ttest_ind(a, b, equal_var=False)
    # each species read where it is defined: monomer plateau vs oligomer t = 0
    am = groups[0]._arr("plateau" if groups[0].measure_at != "baseline" else "baseline")
    bm = groups[1]._arr("plateau" if groups[1].measure_at != "baseline" else "baseline")
    ttm = stats.ttest_ind(am, bm, equal_var=False)
    # ON THE LOG SCALE.  ThT intensities are positive and multiplicative, and the
    # raw variances here differ by ~3600:1 (fibril sd 190 against oligomer sd 3.2),
    # which drives Welch's dof to 2.00 and makes the raw p fragile.  Logs bring the
    # variance ratio to ~31:1 and turn the quantity being tested into the ratio
    # itself, which is what the figure claims.  A rank test cannot help: at n = 3
    # vs 3 the two-sided Mann-Whitney floor is p = 0.10.
    la, lb = np.log(am), np.log(bm)
    ltt = stats.ttest_ind(la, lb, equal_var=False)
    va, vb = la.var(ddof=1) / la.size, lb.var(ddof=1) / lb.size
    se = float(np.sqrt(va + vb))
    # THE INTERVAL MUST USE THE t QUANTILE AT WELCH'S OWN dof, NOT 1.96.  With
    # n = 3 per group the Satterthwaite dof is 2.13 and t(0.975, 2.13) = 4.063,
    # more than twice the normal quantile.  Using 1.96 here reported the ratio as
    # 7.8-18.2x when the honest interval is 5.0-28.6x -- an interval about half as
    # wide as the data support.
    dof = (va + vb) ** 2 / (va ** 2 / (la.size - 1) + vb ** 2 / (lb.size - 1))
    tcrit = float(stats.t.ppf(0.975, dof))
    ratio = float(am.mean() / bm.mean())
    # a Welch p for EVERY metric a panel can plot, so the bracket can never quote a
    # test of a different comparison from the bars it sits over
    pl_a, pl_b = groups[0]._arr("plateau"), groups[1]._arr("plateau")
    tt_pl = stats.ttest_ind(pl_a, pl_b, equal_var=False)
    return dict(
        fold_plateau_p=float(tt_pl.pvalue),
        fold_as_measured=ratio,
        p_as_measured=float(ltt.pvalue),               # the log-scale Welch p
        t_as_measured=float(ltt.statistic),
        p_as_measured_raw=float(ttm.pvalue),           # raw scale, for reference
        ci_as_measured=(float(np.exp(np.log(ratio) - tcrit * se)),
                        float(np.exp(np.log(ratio) + tcrit * se))),
        ci_dof=float(dof), ci_t=tcrit,
        measured_at=(groups[0].measure_at, groups[1].measure_at),
        fold_amplitude=float(a.mean() / b.mean()),
        fold_plateau=float(groups[0]._arr("plateau").mean() / groups[1]._arr("plateau").mean()),
        baseline_ratio=float(groups[1]._arr("baseline").mean() / groups[0]._arr("baseline").mean()),
        t=float(tt.statistic), p=float(tt.pvalue), n=(a.size, b.size),
    )


def report(groups, stats_=None) -> str:
    """Plain-text summary; printed by the builder and embedded in the notebook."""
    L = []
    L.append(f"baseline = median over the first {BASE_H:g} h; "
             f"plateau = mean over the last {PLAT_H:g} h; amplitude = plateau - baseline")
    for g in groups:
        L.append("")
        L.append(f"{g.label}  (n={len(g.per_well)})")
        for w in g.per_well:
            flag = "  <- peak is in the plateau window, still rising" if w.still_rising else ""
            L.append(f"   {w.well}: baseline {w.baseline:7.1f}  plateau {w.plateau:8.1f}  "
                     f"amplitude {w.amplitude:8.1f}  peak {w.peak:7.1f} @ {w.t_peak:5.1f} h  "
                     f"drift {w.drift_pct:+6.1f}%{flag}")
        for name in ("baseline", "plateau", "amplitude"):
            m, s = g.mean_sd(name)
            L.append(f"   mean {name:<9} {m:8.1f} +- {s:6.1f}")
    if stats_:
        L.append("")
        L.append(f"amplitude ratio  {stats_['fold_amplitude']:.2f}x   "
                 f"raw plateau ratio {stats_['fold_plateau']:.2f}x")
        L.append(f"oligomer baseline is {stats_['baseline_ratio']:.2f}x the monomer baseline")
        L.append(f"Welch t-test on amplitude: t = {stats_['t']:.3f}, p = {stats_['p']:.4f} "
                 f"(n = {stats_['n'][0]} vs {stats_['n'][1]})")
    return "\n".join(L)
