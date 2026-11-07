# -*- coding: utf-8 -*-
r"""
secfigs.verify_fig12 -- build and audit figures/Figure_1_2_Combined{,_baselined}.

Everything is MEASURED off the delivered files, or off the renderer that wrote
them, not read back out of :mod:`secfigs.figure_fig12`.

  1  page size in mm, from the SVG root AND the PDF MediaBox
  2  PNG pixel size and embedded dpi
  3  font census -- every <text> in the SVG.  FAILS on anything above 7.0 pt that
     is not an 8.0 pt bold panel letter, or anything below 5.0 pt.
  4  live text -- <text> element count, and no rasterised <image>
  5  ink at the trim edge -- outermost 1, 2 and 3 px rings of the PNG
  6  panel boxes and every gap, decoded from the SVG clip rects
  7  open axes -- which spines each panel actually drew
  8  renderer-level collision check over every text bounding box
  9  data identity -- every plotted vertex of the combined figure against the
     standalone panel SVG it came from (main variant), and the exact 240-250 nm
     offset (baselined variant)

    python -m secfigs.verify_fig12
"""

from __future__ import annotations

import os
import re
import sys

import numpy as np
from lxml import etree

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

import matplotlib                                              # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402
from PIL import Image                                          # noqa: E402

Image.MAX_IMAGE_PIXELS = None

from secfigs import style as gs                                # noqa: E402
from secfigs import cd_panels as cd                            # noqa: E402
from secfigs import figure_fig12 as ff                         # noqa: E402
import build_fig12_panels as bp                                # noqa: E402

SVGNS = "{http://www.w3.org/2000/svg}"
FIGDIR = os.path.join(_REPO, "figures")
PANELDIR = os.path.join(FIGDIR, "panels")
PT_MM = 25.4 / 72.0
FONT_MAX, FONT_MIN, LETTER_PT = 7.0, 5.0, 8.0

_fails = []


def check(ok, msg):
    print(("  PASS  " if ok else "  FAIL  ") + msg)
    if not ok:
        _fails.append(msg)
    return ok


# ==========================================================================
# SVG decoding
# ==========================================================================
def _num(s):
    return float(re.sub(r"[^0-9.eE+-]", "", s))


def load_svg(path):
    root = etree.parse(path).getroot()
    return root, _num(root.get("width")), _num(root.get("height"))


def clip_rects(root):
    out = {}
    for cp in root.iter(SVGNS + "clipPath"):
        r = cp.find(SVGNS + "rect")
        if r is not None:
            out[cp.get("id")] = tuple(_num(r.get(k))
                                      for k in ("x", "y", "width", "height"))
    return out


def parse_d(d):
    return np.asarray(re.findall(r"[-+0-9.eE]+", d), dtype=float).reshape(-1, 2)


def trace_paths(root, min_vertices=40):
    """``[(axes_id, [(clip_rect, verts, stroke_hex), ...]), ...]``, doc order."""
    rects = clip_rects(root)
    out = []
    for g in root.iter(SVGNS + "g"):
        gid = g.get("id") or ""
        if not gid.startswith("axes_"):
            continue
        found = []
        for p in g.iter(SVGNS + "path"):
            d = p.get("d") or ""
            if d.count("L") < min_vertices:
                continue
            m = re.search(r"url\(#([^)]+)\)", p.get("clip-path") or "")
            sm = re.search(r"stroke:\s*(#[0-9a-fA-F]{6})", p.get("style") or "")
            found.append((rects.get(m.group(1)) if m else None, parse_d(d),
                          sm.group(1).lower() if sm else None))
        if found:
            out.append((gid, found))
    return out


def to_data(verts, rect, xlim, ylim):
    rx, ry, rw, rh = rect
    x = xlim[0] + (verts[:, 0] - rx) / rw * (xlim[1] - xlim[0])
    y = ylim[1] - (verts[:, 1] - ry) / rh * (ylim[1] - ylim[0])
    return np.column_stack([x, y])


# ==========================================================================
def check_page(stem):
    print("\n[1] PAGE SIZE")
    root, w_pt, h_pt = load_svg(os.path.join(FIGDIR, stem + ".svg"))
    w_mm, h_mm = w_pt * PT_MM, h_pt * PT_MM
    print("  SVG root      %.4f x %.4f pt = %.4f x %.4f mm"
          % (w_pt, h_pt, w_mm, h_mm))
    check(abs(w_mm - 180.0) < 5e-3, "width is 180.000 mm exactly")
    check(h_mm <= 185.0, "height %.3f mm is inside the 185 mm tier "
                         "(%.2f mm to spare)" % (h_mm, 185.0 - h_mm))
    print("  computed      %.4f x %.4f mm" % (ff.FIG_W_MM, ff.FIG_H_MM))
    check(abs(h_mm - ff.FIG_H_MM) < 5e-3,
          "delivered height matches the derived height to 0.005 mm")

    import fitz
    doc = fitz.open(os.path.join(FIGDIR, stem + ".pdf"))
    mb = doc[0].mediabox
    print("  PDF MediaBox  %.4f %.4f %.4f %.4f pt = %.4f x %.4f mm"
          % (mb.x0, mb.y0, mb.x1, mb.y1, mb.width * PT_MM, mb.height * PT_MM))
    check(abs(mb.width - w_pt) < 0.02 and abs(mb.height - h_pt) < 0.02,
          "PDF MediaBox equals the SVG page to 0.02 pt")
    nf = len(doc[0].get_fonts())
    check(nf > 0, "PDF carries %d embedded font(s) -- text is live" % nf)
    doc.close()

    print("\n[2] PNG")
    im = Image.open(os.path.join(FIGDIR, stem + ".png"))
    dpi = im.info.get("dpi")
    print("  %d x %d px, embedded dpi %s" % (im.width, im.height, dpi))
    exp = (int(round(w_pt / 72.0 * 600)), int(round(h_pt / 72.0 * 600)))
    check(abs(im.width - exp[0]) <= 1 and abs(im.height - exp[1]) <= 1,
          "pixel size matches %d x %d at 600 dpi" % exp)
    check(dpi is not None and abs(dpi[0] - 600) < 1,
          "embedded dpi is 600")
    return root, w_pt, h_pt


_FONT_RE = re.compile(r"font:\s*(?:(\d{3})\s+)?([0-9.]+)px")


def check_fonts(root):
    print("\n[3] FONT CENSUS -- every <text> in the SVG")
    texts = root.findall(".//" + SVGNS + "text")
    census = {}
    for t in texts:
        m = _FONT_RE.search(t.get("style") or "")
        if m is None:
            census.setdefault(("UNPARSED", (t.get("style") or "")[:40]),
                              []).append(t)
            continue
        census.setdefault((float(m.group(2)),
                           "bold" if m.group(1) == "700" else "regular"),
                          []).append(t)
    for k in sorted(census, key=lambda k: (isinstance(k[0], str), k[0])):
        v = census[k]
        s = "".join(x.text or "" for x in v[0].iter())[:30]
        print("   %-6s %-8s  n=%-4d  e.g. %r" % (k[0], k[1], len(v), s))
    check(not [k for k in census if isinstance(k[0], str)],
          "every <text> declares a parseable font size")
    sizes = [k for k in census if not isinstance(k[0], str)]
    over = [k for k in sizes if k[0] > FONT_MAX + 1e-9
            and not (abs(k[0] - LETTER_PT) < 1e-9 and k[1] == "bold")]
    under = [k for k in sizes if k[0] < FONT_MIN - 1e-9]
    check(not over, "nothing above the %.1f pt ceiling except 8 pt BOLD panel "
                    "letters (offenders: %s)" % (FONT_MAX, over or "none"))
    check(not under, "nothing below the %.1f pt floor (offenders: %s)"
                     % (FONT_MIN, under or "none"))
    letters = [k for k in sizes if abs(k[0] - LETTER_PT) < 1e-9]
    check(all(k[1] == "bold" for k in letters) and len(census[letters[0]]) == 5
          if letters else False,
          "the only 8 pt type is 5 bold panel letters")

    print("\n[4] LIVE TEXT")
    print("  %d <text> elements" % len(texts))
    check(len(texts) > 0, "SVG text is live (svg.fonttype='none'), not paths")
    check(not root.findall(".//" + SVGNS + "image"),
          "no rasterised <image> in the SVG -- fully vector")
    return len(texts)


def check_edges(stem):
    print("\n[5] INK AT THE TRIM EDGE")
    a = np.asarray(Image.open(os.path.join(FIGDIR, stem + ".png")).convert("L"))
    for k in (1, 2, 3):
        ring = np.concatenate([a[:k, :].ravel(), a[-k:, :].ravel(),
                               a[:, :k].ravel(), a[:, -k:].ravel()])
        n = int((ring < 250).sum())
        print("   outermost %d px: %d of %d below L 250 (min L %d)"
              % (k, n, ring.size, ring.min()))
        check(n == 0, "outermost %d px ring is pure white" % k)
    ink = a < 250
    rows = np.where(ink.any(axis=1))[0]
    cols = np.where(ink.any(axis=0))[0]
    px_mm = 25.4 / 600.0
    d = dict(top=rows[0], bottom=a.shape[0] - 1 - rows[-1], left=cols[0],
             right=a.shape[1] - 1 - cols[-1])
    print("   first ink: " + "  ".join("%s %d px (%.2f mm)" % (k, v, v * px_mm)
                                       for k, v in d.items()))
    worst = min(d.values()) * px_mm
    check(worst >= 3.0, "closest ink to any trim edge is %.2f mm (>= 3.00 mm; "
                        "TRIM budget %.2f mm)" % (worst, ff.TRIM * 25.4))


def check_geometry(root, page_h):
    print("\n[6] PANEL BOXES, decoded from the SVG clip rects")
    tp = trace_paths(root)
    print("   %-3s %-9s %-9s %-9s %-9s %-7s %s"
          % ("", "x0 (in)", "y0 (in)", "w (in)", "h (in)", "aspect", "traces"))
    boxes = {}
    for L, (gid, found) in zip("abcde", tp):
        rx, ry, rw, rh = found[0][0]
        x0, w = rx / 72.0, rw / 72.0
        h = rh / 72.0
        y0 = (page_h - ry - rh) / 72.0
        boxes[L] = (x0, y0, w, h)
        print("   %-3s %-9.4f %-9.4f %-9.4f %-9.4f %-7.4f %d"
              % (L, x0, y0, w, h, w / h, len(found)))
    exp = ff.CombinedFig12().rects()
    ok = all(max(abs(np.array(boxes[L]) - np.array(exp[L]))) < 2e-3
             for L in "abcde")
    check(ok, "every panel box matches the derived rect to 0.002 in")

    print("\n   derived lengths (in):")
    for k, v in [("TRIM", ff.TRIM), ("letter dx", ff.LETTER_DX_IN),
                 ("margin left", ff.MARGIN["left"]),
                 ("margin right", ff.MARGIN["right"]),
                 ("margin top", ff.MARGIN["top"]),
                 ("margin bottom", ff.MARGIN["bottom"]),
                 ("GAP_W", ff.GAP_W), ("GAP_H", ff.GAP_H),
                 ("panel w row 1", ff.PANEL_W[0]),
                 ("panel w row 2", ff.PANEL_W[1]),
                 ("row h 1", ff.ROW_H[0]), ("row h 2", ff.ROW_H[1])]:
        print("      %-14s %.4f  (%.3f mm)" % (k, v, v / ff.MM))

    print("\n   measured gaps between adjacent panel boxes (in):")
    gw1 = boxes["b"][0] - (boxes["a"][0] + boxes["a"][2])
    gw2 = boxes["d"][0] - (boxes["c"][0] + boxes["c"][2])
    gw3 = boxes["e"][0] - (boxes["d"][0] + boxes["d"][2])
    gh = boxes["a"][1] - (boxes["c"][1] + boxes["c"][3])
    for nm, v in [("a|b", gw1), ("c|d", gw2), ("d|e", gw3),
                  ("row1|row2", gh)]:
        print("      %-10s %.4f" % (nm, v))
    check(max(abs(g - ff.GAP_W) for g in (gw1, gw2, gw3)) < 2e-3,
          "all three column gaps equal GAP_W %.4f in to 0.002 in" % ff.GAP_W)
    check(abs(gh - ff.GAP_H) < 2e-3,
          "the row gap equals GAP_H %.4f in to 0.002 in" % ff.GAP_H)
    check(abs(boxes["a"][0] - boxes["c"][0]) < 2e-3
          and abs((boxes["b"][0] + boxes["b"][2])
                  - (boxes["e"][0] + boxes["e"][2])) < 2e-3,
          "the two rows are flush left and flush right")
    return tp, boxes


def check_open_axes(f):
    print("\n[7] OPEN AXES")
    ok = True
    for L in "abcde":
        sp = {s: bool(f.axes[L].spines[s].get_visible())
              for s in ("left", "bottom", "top", "right")}
        print("   %s  left %-5s bottom %-5s top %-5s right %-5s"
              % (L, sp["left"], sp["bottom"], sp["top"], sp["right"]))
        ok &= sp["left"] and sp["bottom"] and not sp["top"] and not sp["right"]
    check(ok, "all five panels keep left+bottom and drop top+right")


def check_collisions(f):
    """Renderer bbox audit. MUST run inside gs.apply so Arial metrics are live."""
    print("\n[8] TEXT COLLISIONS -- renderer bboxes, every string on the page")
    from matplotlib import font_manager as _fm
    print("   font resolved for metrics: %s"
          % os.path.basename(_fm.findfont(
              _fm.FontProperties(family=["Arial", "DejaVu Sans"]))))
    fig = f.fig
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    items = []

    def add(t, owner):
        s = t.get_text()
        if not s or not t.get_visible():
            return
        bb = t.get_window_extent(r).transformed(fig.dpi_scale_trans.inverted())
        if bb.width <= 0 or bb.height <= 0:
            return
        items.append((owner, s, bb))

    def in_view(axis, lo, hi):
        return [lab for loc, lab in zip(axis.get_majorticklocs(),
                                        axis.get_ticklabels())
                if lo - 1e-9 <= loc <= hi + 1e-9]

    for t in fig.texts:
        add(t, "figure")
    for L, ax in f.axes.items():
        add(ax.title, L)
        add(ax.xaxis.label, L)
        add(ax.yaxis.label, L)
        for t in (in_view(ax.xaxis, *sorted(ax.get_xlim()))
                  + in_view(ax.yaxis, *sorted(ax.get_ylim()))):
            add(t, L)
        for t in ax.texts:
            add(t, L)
        leg = ax.get_legend()
        if leg is not None:
            for t in leg.get_texts():
                add(t, L + ":legend")
    print("   %d text bounding boxes tested, %d pairs"
          % (len(items), len(items) * (len(items) - 1) // 2))

    hits, worst = [], (1e9, None)
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = items[i][2], items[j][2]
            dx = min(a.x1, b.x1) - max(a.x0, b.x0)
            dy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if dx > 0 and dy > 0:
                hits.append((items[i], items[j], dx, dy))
            else:
                sep = max(-dx, -dy)
                if sep < worst[0]:
                    worst = (sep, (items[i], items[j]))
    for a, b, dx, dy in hits:
        print("   OVERLAP  %r (%s) x %r (%s) by %.4f x %.4f in"
              % (a[1], a[0], b[1], b[0], dx, dy))
    check(not hits, "no two text bounding boxes intersect")
    sep, pair = worst
    print("   tightest non-touching pair: %r (%s) / %r (%s), %.4f in = %.2f mm"
          % (pair[0][1], pair[0][0], pair[1][1], pair[1][0], sep, sep * 25.4))
    check(sep * 25.4 >= 0.5,
          "tightest text pair is %.2f mm apart (>= 0.50 mm)" % (sep * 25.4))

    out = [(o, s, bb) for o, s, bb in items
           if bb.x0 < ff.TRIM - 1e-6 or bb.y0 < ff.TRIM - 1e-6
           or bb.x1 > ff.FIG_W - ff.TRIM + 1e-6
           or bb.y1 > ff.FIG_H - ff.TRIM + 1e-6]
    for o, s, bb in out:
        print("   OUTSIDE TRIM  %r (%s) at [%.4f %.4f %.4f %.4f]"
              % (s, o, bb.x0, bb.y0, bb.x1, bb.y1))
    check(not out, "every string sits inside the %.2f mm trim box"
                   % (ff.TRIM * 25.4))


def panel_lims(combined: bool):
    """``pid -> (xlim, ylim)``.

    SEC is the normalised convention in both files.  CD differs BY DESIGN: the
    combined figure restores the house 5%-of-tick-range padding that
    ``cd_panels.cd_ylim``'s snap-to-tick destroys (see
    ``figure_fig12.cd_ylim_padded``), so each file has to be mapped through its
    own limits.  The DATA is what is being compared, not the mapping.
    """
    lims = {"1A": ((6.0, 20.0), (-0.05, 1.05)),
            "1B": ((6.0, 20.0), (-0.05, 1.05))}
    for pid in ("1C", "1D", "2A"):
        p, traces = cd.load_panel(pid)
        if combined:
            yl, _ = ff.cd_ylim_padded(pid)
        else:
            yl, _ = cd.cd_ylim(traces, p["xlim"])
        lims[pid] = (tuple(p["xlim"]), tuple(yl))
    return lims


def _compare_curves(dc, ds, xlim):
    """max |dy| between two decoded polylines, on their shared x support.

    An elementwise diff is NOT available and the reason is a property of the
    SVG backend, not of the data: ``RendererSVG.draw_path`` clips an unfilled
    path to the FIGURE rectangle before serialising it, so a trace whose data
    runs past the axes (panel c's spectra reach 280 nm against a 250 nm window)
    is written with a different number of vertices on a 12.41 in standalone
    canvas than on this 7.09 in page -- 69 against 91 for panel c.  Both
    polylines still describe the same single-valued curve, so the honest
    comparison is to evaluate one on the other's abscissae.
    """
    oc = np.argsort(dc[:, 0])
    os_ = np.argsort(ds[:, 0])
    xc, yc = dc[oc, 0], dc[oc, 1]
    xs, ys = ds[os_, 0], ds[os_, 1]
    lo = max(xc.min(), xs.min(), xlim[0])
    hi = min(xc.max(), xs.max(), xlim[1])
    m = (xc >= lo) & (xc <= hi)
    if not m.any():
        return np.nan, 0
    return float(np.abs(yc[m] - np.interp(xc[m], xs, ys)).max()), int(m.sum())


def check_data(tp, page_h):
    print("\n[9] DATA IDENTITY -- combined figure vs the standalone panels")
    lc, ls = panel_lims(True), panel_lims(False)
    worst_all, total = 0.0, 0
    print("   panel pid  traces  vertices here/standalone  compared   max |dy|")
    for L, (gid, found) in zip("abcde", tp):
        pid = dict((s["letter"], s["pid"]) for s in ff.PANEL_SPEC)[L]
        sroot, sw, sh = load_svg(os.path.join(PANELDIR, pid + ".svg"))
        stp = trace_paths(sroot)
        if len(stp) != 1 or len(stp[0][1]) != len(found):
            check(False, "%s: %d traces here vs %d standalone"
                         % (L, len(found), len(stp[0][1]) if stp else 0))
            continue
        xlim, ylim_c = lc[pid]
        _, ylim_s = ls[pid]
        worst, nc, ns, ncmp = 0.0, 0, 0, 0
        for (rc, vc, _), (rs, vs, _) in zip(found, stp[0][1]):
            d, n = _compare_curves(to_data(vc, rc, xlim, ylim_c),
                                   to_data(vs, rs, xlim, ylim_s), xlim)
            worst = max(worst, d)
            nc += len(vc)
            ns += len(vs)
            ncmp += n
        total += ncmp
        print("     %s    %-4s   %2d     %6d / %-6d          %7d   %.3e"
              % (L, pid, len(found), nc, ns, ncmp, worst))
        worst_all = max(worst_all, worst)
    print("   %d vertices compared; worst deviation %.3e" % (total, worst_all))
    check(worst_all < 5e-3, "every plotted vertex reproduces its standalone "
                            "panel to better than 5e-3 in data units "
                            "(the SVG carries 6 significant figures at two very "
                            "different axes sizes, so this is quantisation)")

    print("\n   CD troughs vs the decoded artwork (cd_panels.verify):")
    for pid in ("1C", "1D", "2A"):
        v = cd.verify(pid)
        print("     %-4s %s   max dev %.1f nm / %.4f mdeg"
              % (pid, v["verdict"], v["max_d_nm"], v["max_d_mdeg"]))
        check(v["verdict"] == "MATCH", "%s reproduces her artwork" % pid)


def check_baselined(stem_bl, stem_main):
    print("\n[9b] BASELINED VARIANT -- the offset actually applied")
    lims = panel_lims(True)
    with ff.figure_inks(True):
        lims_b = panel_lims(True)
    same = all(np.allclose(lims[k][1], lims_b[k][1]) for k in lims)
    print("   y limits identical in both variants: %s  %s"
          % (same, {k: tuple(round(v, 3) for v in lims[k][1]) for k in lims}))
    check(same, "baselining does not move any panel's y limits, so the two "
                "variants are panel-for-panel comparable")
    broot, _, bh = load_svg(os.path.join(FIGDIR, stem_bl + ".svg"))
    mroot, _, mh = load_svg(os.path.join(FIGDIR, stem_main + ".svg"))
    btp, mtp = trace_paths(broot), trace_paths(mroot)
    rep = ff.baseline_report()
    off_by_pid = {pid: [r["offset"] for r in rep["panels"][pid]["rows"]]
                  for pid in ("1C", "1D", "2A")}
    worst = 0.0
    for L, (gb, fb), (gm, fm) in zip("abcde", btp, mtp):
        pid = dict((s["letter"], s["pid"]) for s in ff.PANEL_SPEC)[L]
        xlim, ylim = lims[pid]
        offs = off_by_pid.get(pid, [0.0] * len(fb))
        for (rb, vb, _), (rm, vm, _), off in zip(fb, fm, offs):
            db = to_data(vb, rb, xlim, ylim)
            dm = to_data(vm, rm, xlim, ylim)
            dm[:, 1] -= off
            d, _n = _compare_curves(db, dm, xlim)
            worst = max(worst, d)
        print("     %s  %-4s  %d traces, offsets %s"
              % (L, pid, len(fb), np.round(offs, 4)))
    print("   worst |(main - offset) - baselined| = %.3e mdeg" % worst)
    check(worst < 5e-3, "every baselined CD vertex is its main-variant vertex "
                        "minus that trace's own 240-250 nm mean")


# ==========================================================================
def audit(stem, f, data=True):
    root, w, h = check_page(stem)
    check_fonts(root)
    check_edges(stem)
    tp, boxes = check_geometry(root, h)
    check_open_axes(f)
    with gs.apply(style=f.style, **ff.RC_EXTRA):
        check_collisions(f)
    if data:
        check_data(tp, h)
    return root, boxes


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    print("=" * 78)
    print("BUILD -- standalone panels first, so [9] compares like with like")
    print("=" * 78)
    bp.main()

    print("\n" + "=" * 78)
    print("COLOUR KEY")
    print("=" * 78)
    time, temp, rounds = ff.derive_keys()
    print("  re-derived in %d round(s): time %s   temperature %s"
          % (rounds, time, temp))
    check(list(time) == [ff.TIME_INKS[t] for t in ff.FIG_TIMES]
          and list(temp) == list(ff.TEMP_INKS),
          "the frozen TIME_INKS / TEMP_INKS are what derive_keys() produces")
    rep = ff.color_report()
    for k, v in rep["inks"].items():
        print("   %-9s %-9s L* %6.2f  a* %7.2f  b* %7.2f"
              % (k, v, rep["lab"][k][0], rep["lab"][k][1], rep["lab"][k][2]))
    print("\n   PER-PANEL minimum dE76 (traces + frame ink, worst of the four "
          "renderings):")
    for L in "abcde":
        print("     %s  %6.2f   [%s]"
              % (L, rep["per_panel"][L], ", ".join(ff.PANEL_INKS[L])))
    check(rep["panel_min"] >= 17.0,
          "the weakest WITHIN-PANEL separation on the page is %.2f dE76"
          % rep["panel_min"])
    print("\n   min over ALL %d pairs of the eight inks, per rendering: "
          "normal %.2f  deutan %.2f  protan %.2f  grey %.2f"
          % (len(rep["pairs"]), rep["per_rendering"]["normal"],
             rep["per_rendering"]["deutan"], rep["per_rendering"]["protan"],
             rep["per_rendering"]["grey"]))
    print("   worst pair overall %.2f (%s) -- a cross-key, cross-panel pair, "
          "see the module docstring sec. 3" % (rep["worst"], rep["worst_pair"]))
    print("\n   CROSS-KEY pairs (time ink vs temperature ink):")
    for k, v in sorted(rep["cross"].items(), key=lambda kv: kv[1]["worst3"]):
        print("     %-22s worst3 %6.2f   normal %6.2f deutan %6.2f protan %6.2f"
              "   (grey %6.2f, not constrained)"
              % (k, v["worst3"], v["normal"], v["deutan"], v["protan"],
                 v["grey"]))
    check(rep["cross_key_min"] >= ff.CROSS_KEY_FLOOR - 1e-9,
          "every cross-key pair clears the %.0f dE76 floor (min %.2f)"
          % (ff.CROSS_KEY_FLOOR, rep["cross_key_min"]))
    print("\n   tightest six pairs overall:")
    for k, v in sorted(rep["pairs"].items(), key=lambda kv: kv[1]["worst"])[:6]:
        print("     %-24s worst %5.2f  normal %6.2f deutan %6.2f protan %6.2f "
              "grey %6.2f" % (k, v["worst"], v["normal"], v["deutan"],
                              v["protan"], v["grey"]))

    figs = {}
    for baselined in (False, True):
        f = ff.CombinedFig12(baselined=baselined)
        f.draw()
        files = f.save()
        figs[baselined] = (f, files)
        print("\n" + "=" * 78)
        print("AUDIT -- %s%s" % (f.stem, "  (baselined)" if baselined else ""))
        print("=" * 78)
        for k, v in files.items():
            print("   wrote %s" % v)
        audit(f.stem, f, data=not baselined)

    print("\n" + "=" * 78)
    print("BASELINE REPORT")
    print("=" * 78)
    rep = ff.baseline_report()
    for pid in ("1C", "1D", "2A"):
        d = rep["panels"][pid]
        print("  %s   offsets removed (240-250 nm mean), spread %.4f mdeg"
              % (pid, d["spread_before"]))
        for r in d["rows"]:
            print("     %-12s %-11s offset %+8.4f -> %+.2e   trough %8.4f -> "
                  "%8.4f   sd(240-250) smoothed %.4f raw %.4f"
                  % (r["key"], r["label"], r["offset"], r["residual_after"],
                     r["trough"][1], r["trough_baselined"][1],
                     r["sd_smoothed"], r["sd_raw"]))
    c = rep["panel_c"]
    print("\n  PANEL c, 6 C vs RT (means of the two repeats):")
    print("     offset over 212-250 nm   %+.4f  ->  %+.4f mdeg"
          % (c["offset_before"], c["offset_after"]))
    print("     trough difference        %+.4f  ->  %+.4f mdeg"
          % (c["trough_diff_before"], c["trough_diff_after"]))
    print("     within-group repeat spread (6 C, RT)  %s -> %s"
          % (tuple(round(v, 4) for v in c["within_group_trough_spread_before"]),
             tuple(round(v, 4) for v in c["within_group_trough_spread_after"])))

    check_baselined("Figure_1_2_Combined_baselined", "Figure_1_2_Combined")

    for f, _ in figs.values():
        plt.close(f.fig)

    print("\n" + "=" * 78)
    if _fails:
        print("%d CHECK(S) FAILED:" % len(_fails))
        for m in _fails:
            print("   - " + m)
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
