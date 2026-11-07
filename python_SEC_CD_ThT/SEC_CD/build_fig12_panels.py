# -*- coding: utf-8 -*-
r"""
build_fig12_panels.py -- build ALL FIVE panels of Figures 1 and 2 of the 37 C
Abeta42 oligomer manuscript (doc 10) from raw instrument exports, as individually
checkable files.

    1A.svg/.png   SEC   cold incubation, t = 0 and 24 h
    1B.svg/.png   SEC   37 C, t = 0 and t = 2 h   (published as "24 h" -- see below)
    1C.svg/.png   CD    SEC-isolated monomers, 2 temperatures x 2 repeats
    1D.svg/.png   CD    SEC-isolated oligomers, 37 C / 24 h
    2A.svg/.png   CD    unfractionated solutions, 24 h at 6 / 25 / 37 C

Style   : secfigs/style.py      -- every size, weight, width and ink
Colour  : secfigs/cd_panels.py  -- an ordered TEMPERATURE key, derived, not the
                                   time ramp (see that module's sec. 3)
SEC data: secfigs/panel_data.py
CD data : secfigs/cd_panels.py

Output contract, per the Nature vector requirements this project already adopted:
``svg.fonttype = "none"`` (live, editable text), ``path.simplify = False`` (no vertex
decimation of the sharp peak flanks), no ``bbox_inches="tight"``, PNG at 600 dpi.

Nothing under ``/path/to/Dropbox`` is written, moved or imported.  The raw
trees are opened read-only; the processing contracts are reimplemented in
``build_fig1_panels`` (SEC) and ``secfigs.cd_panels`` (CD) rather than imported from
her notebooks, because importing a module that lives under Dropbox would drop a
``__pycache__`` into it.

Every panel exposes ``draw(ax, pid)``, so the combined Figure 1 + Figure 2 layout is
a grid of calls into :func:`draw_panel` and needs none of this file's saving code.
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import matplotlib                                          # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                            # noqa: E402

from secfigs import style as gs                            # noqa: E402
from secfigs import cd_panels as cd                        # noqa: E402
import build_fig1_panels as sec                            # noqa: E402

OUTDIR = os.path.join(HERE, "figures", "panels")

#: Live text and un-decimated paths in the SVG.  Set globally rather than inside
#: ``gs.apply()`` because ``savefig`` reads them at write time.
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["path.simplify"] = False

PNG_DPI = 600

SEC_PANELS = ("1A", "1B")
CD_PANEL_IDS = ("1C", "1D", "2A")
ALL_PANELS = SEC_PANELS + CD_PANEL_IDS


# ==========================================================================
# SEC draw(ax) -- the same pipeline build_fig1_panels already verified, lifted
# out of its figure-owning build_panel() so a combined figure can call it.
# ==========================================================================
#: Which peak pointers :func:`draw_sec` draws by default.  A panel whose oligomer
#: window holds no peak should pass a shorter tuple -- ``OLIGO_Y_FLOOR`` would
#: otherwise float an "Oligomer" arrow at 0.10 of full scale over empty page.
ANNOTATE_BOTH = ("Oligomer", "Monomer")


def draw_sec(ax, pid: str, style: gs.GriffinStyle = gs.STYLE, *,
             title: bool = True, letter: bool = True, legend: bool = True,
             labels: bool = True, annotate=ANNOTATE_BOTH):
    """Draw SEC panel ``1A`` or ``1B`` into an existing Axes. Metadata dict back.

    Identical processing, colours, peak detection and pointer geometry to
    ``build_fig1_panels.build_panel`` -- this only removes that function's
    ``fig, ax`` ownership, so 1A/1B render byte-comparably here and are reusable in
    the combined figure.

    ``annotate`` selects which peak pointers are drawn; see :data:`ANNOTATE_BOTH`.

    DEFECT FIXED 2026-08-08.  The claim of identical pointer geometry was false.
    ``build_fig1_panels.build_panel`` passes ``marker="<"`` for the oligomer
    pointer; this function passed nothing and fell through to
    ``style.annot_marker_right`` = ``">"``.  Since the pointer is placed
    ``OLIGO_DX_FRAC`` to the RIGHT of the oligomer apex, ">" aims it away from the
    peak it names, into the inter-peak valley and on towards the monomer -- the
    exact mislabel ``build_fig1_panels`` documents having fixed once already, and
    the one ``secfigs/figure_sec.py`` fixed before that.
    """
    spec = sec.SPEC[pid]
    p, traces = sec.load_panel(pid)
    per, anchor_ol, anchor_mo = sec.find_peaks_for_panel(traces)
    xlim = spec["xlim"]
    xspan = xlim[1] - xlim[0]

    for t, c in zip(traces, sec.sec_colors(style)):
        ax.plot(t["V"], t["A"], color=c, lw=style.trace_linewidth,
                label=t["label"], solid_capstyle="round")

    gs.axes(ax, style, normalized=p["normalized"], x_major=2.0, x_range=xlim)
    if labels:
        gs.set_labels(ax, "Elution Volume (ml)", "Normalized Absorbance", style)
    if title:
        gs.set_title(ax, spec["title"], style)
    if legend:
        gs.legend(ax, style)

    if anchor_ol is not None and "Oligomer" in annotate:
        xo, ho = anchor_ol
        gs.annotate_peak(ax, xo + sec.OLIGO_DX_FRAC * xspan,
                         max(sec.OLIGO_Y_FRAC * ho, sec.OLIGO_Y_FLOOR),
                         "Oligomer", side="right", style=style, marker="<")
    if anchor_mo is not None and "Monomer" in annotate:
        xm, hm = anchor_mo
        gs.annotate_peak(ax, xm + sec.MONO_DX_FRAC * xspan,
                         sec.MONO_Y_FRAC * hm,
                         "Monomer", side="left", style=style, marker=">")
    if letter:
        gs.panel_letter(ax, spec["letter"], style)

    return dict(pid=pid, kind="SEC", panel=p, traces=traces, peaks=per,
                n_traces=len(traces), xlim=xlim,
                labels={t["label"]: os.path.basename(t["path"]) for t in traces},
                anchor_oligomer=anchor_ol, anchor_monomer=anchor_mo)


# ==========================================================================
# Canvas fitting
# ==========================================================================
def fit_canvas(fig, style: gs.GriffinStyle = gs.STYLE, pad_pt: float = None,
               max_iter: int = 4):
    """Grow the canvas until nothing is clipped, WITHOUT resizing the axes box.

    ``style.panel_letter_xy`` is ``(-0.20, 1.13)`` in axes fraction, i.e. the panel
    letter is deliberately drawn OUTSIDE the axes -- and the style's own usage
    example saves with ``bbox_inches="tight"``, which is what used to absorb it.
    The Nature output contract forbids ``bbox_inches``, so the letter fell off the
    left edge of the canvas instead; the existing 1A and 1B both shipped with a
    clipped letter.

    The fix has to leave the axes box alone.  Every point value in
    :class:`~secfigs.style.GriffinStyle` is calibrated against a 9.30 x 7.34 in axes
    box (aspect 1.266 vs the reference 1.270), which the default subplot params give
    at ``figsize=(12, 9.53)``; shrinking the axes to make room would silently
    decalibrate the whole style.  So this measures the real overflow with the
    renderer, grows the FIGURE by exactly that much, and repositions the axes so its
    absolute size in inches is unchanged.  The only free parameter is the padding,
    and that is taken from ``style.tick_pad``.
    """
    pad = (style.tick_pad if pad_pt is None else pad_pt) / 72.0   # pt -> in
    for _ in range(max_iter):
        fig.canvas.draw()
        fw, fh = fig.get_size_inches()
        tb = fig.get_tightbbox(fig.canvas.get_renderer())
        dl = max(0.0, pad - tb.x0)
        db = max(0.0, pad - tb.y0)
        dr = max(0.0, tb.x1 + pad - fw)
        dt = max(0.0, tb.y1 + pad - fh)
        if max(dl, db, dr, dt) < 1e-6:
            return fig
        nfw, nfh = fw + dl + dr, fh + db + dt
        boxes = [(ax, ax.get_position()) for ax in fig.axes]
        fig.set_size_inches(nfw, nfh, forward=False)
        for ax, p in boxes:
            ax.set_position([(p.x0 * fw + dl) / nfw, (p.y0 * fh + db) / nfh,
                             p.width * fw / nfw, p.height * fh / nfh])
    return fig


def axes_box_inches(fig, ax):
    """The axes box in inches -- the quantity the style is calibrated against."""
    p = ax.get_position()
    fw, fh = fig.get_size_inches()
    return (p.width * fw, p.height * fh)


# ==========================================================================
# One entry point for all five
# ==========================================================================
def draw_panel(ax, pid: str, style: gs.GriffinStyle = gs.STYLE, **kw):
    """``draw(ax, pid)`` for any of 1A, 1B, 1C, 1D, 2A."""
    if pid in SEC_PANELS:
        return draw_sec(ax, pid, style, **kw)
    if pid in CD_PANEL_IDS:
        return cd.draw(ax, pid, style, **kw)
    raise KeyError(f"unknown panel {pid!r}; have {ALL_PANELS}")


def build_panel(pid: str, style: gs.GriffinStyle = gs.STYLE,
                outdir: str = None, stem: str = None, **kw):
    """Render one panel standalone to ``<stem>.svg`` and ``<stem>.png``."""
    outdir = outdir or OUTDIR
    stem = stem or pid
    os.makedirs(outdir, exist_ok=True)
    with gs.apply(style):
        fig, ax = gs.figure(style)
        meta = draw_panel(ax, pid, style, **kw)
        box_before = axes_box_inches(fig, ax)
        fit_canvas(fig, style)
        box_after = axes_box_inches(fig, ax)
        svg = os.path.join(outdir, f"{stem}.svg")
        png = os.path.join(outdir, f"{stem}.png")
        fig.savefig(svg, facecolor="white")            # no bbox_inches
        fig.savefig(png, dpi=PNG_DPI, facecolor="white")
        meta["figsize_in"] = tuple(round(float(v), 4) for v in fig.get_size_inches())
        meta["axes_box_in"] = tuple(round(v, 4) for v in box_after)
        meta["axes_box_preserved"] = all(
            abs(a - b) < 1e-6 for a, b in zip(box_before, box_after))
        meta["axes_box_aspect"] = round(box_after[0] / box_after[1], 4)
        plt.close(fig)
    meta["svg"], meta["png"] = svg, png
    return meta


# ==========================================================================
# Reporting
# ==========================================================================
def _report_sec(meta):
    print("   traces %d   xlim %s" % (meta["n_traces"], meta["xlim"]))
    for t, pk in zip(meta["traces"], meta["peaks"]):
        o, m = pk["oligomer"], pk["monomer"]
        print("     %-9s oligomer %-22s monomer %-22s %s"
              % (t["label"],
                 ("V=%6.3f h=%.4f" % o) if o else "none",
                 ("V=%6.3f h=%.4f" % m) if m else "none",
                 os.path.basename(t["path"])))


def _report_cd(meta):
    v = cd.verify(meta["pid"])
    print("   traces %d   xlim %s   ylim (%.2f, %.2f) step %.2f   savgol %s"
          % (meta["n_traces"], meta["xlim"], meta["ylim"][0], meta["ylim"][1],
             meta["y_major"], meta["savgol"]))
    lg = meta["legend"]
    if lg:
        print("   legend '%s' (measured: covers %d trace vertices, nearest %.3f)"
              % (lg["loc"], lg["covered_vertices"], lg["nearest_vertex"]))
    print("     %-22s %-9s %-20s %-20s %s"
          % ("label", "colour", "trough built", "trough artwork", "dev"))
    for r in v["rows"]:
        print("     %-22s %-9s %8.1f nm %8.3f  %8.1f nm %8.3f   %.1f nm %.4f mdeg"
              % (r["label"], meta["colors"][r["key"]],
                 r["built"][0], r["built"][1],
                 r["artwork"][0], r["artwork"][1],
                 r["d_nm"], r["d_mdeg"]))
    print("   VERDICT %s   max dev %.1f nm / %.4f mdeg"
          % (v["verdict"], v["max_d_nm"], v["max_d_mdeg"]))
    return v


def main(panels=ALL_PANELS):
    rep = cd.temp_color_report()
    print("=" * 86)
    print("TEMPERATURE COLOUR KEY (resolved through style.foreground = %s)"
          % gs.STYLE.foreground)
    for t, c in zip(cd.TEMPS, cd.temp_colors()):
        print("   %5.0f C  ->  %s" % (t, c))
    print("   worst pairwise dE76 involving a temperature member, over the whole")
    print("   Figure-1 ink set, worst of {normal, deutan, protan, greyscale}: %.2f"
          % rep["worst_involving_temperature"])

    results = {}
    for pid in panels:
        meta = build_panel(pid)
        results[pid] = meta
        print("=" * 86)
        print("%s   %s" % (pid, meta["panel"]["title"]))
        if pid in SEC_PANELS:
            _report_sec(meta)
        else:
            meta["verification"] = _report_cd(meta)
            for f in meta["flags"]:
                print("   FLAG %s" % f)
        print("   canvas %.3f x %.3f in   axes box %.3f x %.3f in "
              "(aspect %.3f, preserved=%s)"
              % (meta["figsize_in"] + meta["axes_box_in"]
                 + (meta["axes_box_aspect"], meta["axes_box_preserved"])))
        print("   wrote %s" % meta["svg"])
        print("   wrote %s" % meta["png"])
    return results


if __name__ == "__main__":
    main(sys.argv[1:] or ALL_PANELS)
