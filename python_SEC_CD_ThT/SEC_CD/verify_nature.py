# -*- coding: utf-8 -*-
r"""
verify_nature.py -- audit figures/Figure_SEC_natstyle.{svg,png,pdf}.

Everything here is MEASURED off the delivered files (or off a re-render of the
same figure through the renderer), not read back out of the module that wrote
them.  Nine checks:

  1  page size, in mm, from the SVG root and from the PDF MediaBox
  2  PNG pixel size and embedded dpi
  3  font census -- every <text> in the SVG, by point size; FAILS if anything
     is under 7.0 pt
  4  live text -- <text> element count (0 would mean the type got pathed)
  5  ink at the trim edge -- outermost 1, 2 and 3 px rings of the PNG
  6  panel boxes and every derived gap, from the SVG clip rects
  7  open axes -- which spines each panel actually drew
  8  renderer-level collision check over every text bounding box
  9  DATA IDENTITY -- decode the 26 trace paths out of BOTH
     Figure_SEC_natstyle.svg and Figure_SEC_normalized.svg, map each back to
     data coordinates through its own axes clip rect, and difference them.
     This is a restyle: the numbers must not move.

    python verify_nature.py
"""

from __future__ import annotations

import os
import re
import sys

import numpy as np
from lxml import etree

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import matplotlib                                              # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402
from PIL import Image                                          # noqa: E402

Image.MAX_IMAGE_PIXELS = None

from secfigs import style as gs                                # noqa: E402
from secfigs import figure_sec_nature as fsn                   # noqa: E402
from secfigs.figure_sec import XLIM, YLIM                      # noqa: E402

SVGNS = "{http://www.w3.org/2000/svg}"
FIGDIR = os.path.join(_HERE, "figures")
NAT = os.path.join(FIGDIR, "Figure_SEC_natstyle")
REF = os.path.join(FIGDIR, "Figure_SEC_normalized")
PT_MM = 25.4 / 72.0
FONT_FLOOR = 7.0

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
    """``clipPath id -> (x, y, w, h)`` in SVG user units (= points here)."""
    out = {}
    for cp in root.iter(SVGNS + "clipPath"):
        r = cp.find(SVGNS + "rect")
        if r is not None:
            out[cp.get("id")] = tuple(_num(r.get(k))
                                      for k in ("x", "y", "width", "height"))
    return out


def parse_d(d):
    """'M x y L x y ...' -> (N, 2) array.  Trace paths are M/L only."""
    toks = re.findall(r"[-+0-9.eE]+", d)
    v = np.asarray(toks, dtype=float)
    return v.reshape(-1, 2)


def trace_paths(root, min_vertices=500):
    """``[(axes_id, [(clip_id, verts), ...]), ...]`` in document order."""
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
            cid = (p.get("clip-path") or "")
            m = re.search(r"url\(#([^)]+)\)", cid)
            found.append((rects.get(m.group(1)) if m else None,
                          parse_d(d),
                          re.search(r"stroke:\s*(#[0-9a-fA-F]{6})",
                                    p.get("style") or "")))
        if found:
            out.append((gid, found))
    return out


def to_data(verts, rect, page_h):
    """SVG user units -> (elution volume mL, normalized absorbance)."""
    rx, ry, rw, rh = rect
    x = XLIM[0] + (verts[:, 0] - rx) / rw * (XLIM[1] - XLIM[0])
    # SVG y grows downward; the rect's TOP edge is YLIM[1]
    y = YLIM[1] - (verts[:, 1] - ry) / rh * (YLIM[1] - YLIM[0])
    return np.column_stack([x, y])


# ==========================================================================
# 1-2  page size, dpi
# ==========================================================================
def check_page():
    print("\n[1] PAGE SIZE")
    root, w_pt, h_pt = load_svg(NAT + ".svg")
    w_mm, h_mm = w_pt * PT_MM, h_pt * PT_MM
    print("  SVG root      %.4f x %.4f pt = %.4f x %.4f mm"
          % (w_pt, h_pt, w_mm, h_mm))
    check(abs(w_mm - 180.0) < 5e-3, "width is 180.000 mm (Nat Comms double column)")
    check(h_mm <= 230.0, "height %.3f mm is inside the 230 mm cap "
                         "(%.2f mm to spare)" % (h_mm, 230.0 - h_mm))
    print("  computed      %.4f x %.4f mm" % (fsn.FIG_W_MM, fsn.FIG_H_MM))
    check(abs(h_mm - fsn.FIG_H_MM) < 5e-3, "delivered height matches the "
          "derived height to 0.005 mm")

    import fitz
    doc = fitz.open(NAT + ".pdf")
    mb = doc[0].mediabox
    print("  PDF MediaBox  %.4f %.4f %.4f %.4f pt = %.4f x %.4f mm"
          % (mb.x0, mb.y0, mb.x1, mb.y1,
             mb.width * PT_MM, mb.height * PT_MM))
    check(abs(mb.width - w_pt) < 0.02 and abs(mb.height - h_pt) < 0.02,
          "PDF MediaBox equals the SVG page to 0.02 pt")
    nfonts = len(doc[0].get_fonts())
    check(nfonts > 0, "PDF carries %d embedded font(s) -- text is live, "
                      "not outlined" % nfonts)
    doc.close()

    print("\n[2] PNG")
    im = Image.open(NAT + ".png")
    dpi = im.info.get("dpi")
    print("  %d x %d px, embedded dpi %s" % (im.width, im.height, dpi))
    exp = (int(round(w_pt / 72.0 * 600)), int(round(h_pt / 72.0 * 600)))
    check(abs(im.width - exp[0]) <= 1 and abs(im.height - exp[1]) <= 1,
          "pixel size matches %d x %d at 600 dpi (+/-1 px rounding)" % exp)
    check(dpi is not None and abs(dpi[0] - 600) < 1 and abs(dpi[1] - 600) < 1,
          "embedded dpi is 600 (Nature asks >= 300)")
    return root, w_pt, h_pt


# ==========================================================================
# 3-4  fonts
# ==========================================================================
_FONT_RE = re.compile(r"font:\s*(?:(\d{3})\s+)?([0-9.]+)px")


def check_fonts(root):
    print("\n[3] FONT CENSUS -- every <text> in the SVG")
    texts = root.findall(".//" + SVGNS + "text")
    census = {}
    for t in texts:
        st = t.get("style") or ""
        m = _FONT_RE.search(st)
        if m is None:
            census.setdefault(("UNPARSED", st[:40]), []).append(t)
            continue
        weight = "bold" if m.group(1) == "700" else "regular"
        census.setdefault((float(m.group(2)), weight), []).append(t)

    def sample(v):
        s = "".join(x.text or "" for x in v[0].iter())
        return s[:26]

    for k in sorted(census, key=lambda k: (isinstance(k[0], str), k[0])):
        v = census[k]
        print("   %-6s %-8s  n=%-4d  e.g. %r"
              % (k[0], k[1], len(v), sample(v)))
    sizes = [k[0] for k in census if not isinstance(k[0], str)]
    check(not [k for k in census if isinstance(k[0], str)],
          "every <text> declares a parseable font size")
    check(sizes and min(sizes) >= FONT_FLOOR,
          "smallest type on the page is %.1f pt (floor %.1f pt)"
          % (min(sizes), FONT_FLOOR))

    print("\n[4] LIVE TEXT")
    print("  %d <text> elements" % len(texts))
    check(len(texts) > 0, "SVG text is live (svg.fonttype='none'), not paths")
    check(not root.findall(".//" + SVGNS + "image"),
          "no rasterised <image> in the SVG -- fully vector")


# ==========================================================================
# 5  trim-edge ink
# ==========================================================================
def check_edges():
    print("\n[5] INK AT THE TRIM EDGE")
    a = np.asarray(Image.open(NAT + ".png").convert("L"))
    for k in (1, 2, 3):
        ring = np.concatenate([a[:k, :].ravel(), a[-k:, :].ravel(),
                               a[:, :k].ravel(), a[:, -k:].ravel()])
        n = int((ring < 250).sum())
        print("   outermost %d px: %d of %d pixels below L 250 (min L %d)"
              % (k, n, ring.size, ring.min()))
        check(n == 0, "outermost %d px ring is pure white" % k)
    # how far in does the ink actually start?
    ink = a < 250
    rows = np.where(ink.any(axis=1))[0]
    cols = np.where(ink.any(axis=0))[0]
    px_mm = 25.4 / 600.0
    print("   first ink: top %d px (%.2f mm), bottom %d px (%.2f mm), "
          "left %d px (%.2f mm), right %d px (%.2f mm)"
          % (rows[0], rows[0] * px_mm, a.shape[0] - 1 - rows[-1],
             (a.shape[0] - 1 - rows[-1]) * px_mm, cols[0], cols[0] * px_mm,
             a.shape[1] - 1 - cols[-1], (a.shape[1] - 1 - cols[-1]) * px_mm))
    worst = min(rows[0], a.shape[0] - 1 - rows[-1], cols[0],
                a.shape[1] - 1 - cols[-1]) * px_mm
    check(worst >= 2.0, "closest ink to any trim edge is %.2f mm "
                        "(TRIM budget %.2f mm)" % (worst, fsn.TRIM * 25.4))


# ==========================================================================
# 6-7  geometry and open axes, from the SVG itself
# ==========================================================================
def check_geometry(root, page_h):
    print("\n[6] PANEL BOXES, decoded from the SVG clip rects")
    tp = trace_paths(root)
    boxes = []
    for gid, found in tp:
        rect = found[0][0]
        boxes.append((gid, rect, len(found)))
    letters = ["a", "b", "c", "d"]
    print("   panel  n   x0 (in)   y0 (in)   w (in)    h (in)    aspect")
    for L, (gid, r, n) in zip(letters, boxes):
        x0, y0 = r[0] / 72.0, (page_h - r[1] - r[3]) / 72.0
        w, h = r[2] / 72.0, r[3] / 72.0
        print("     %s   %2d  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f"
              % (L, n, x0, y0, w, h, w / h))
    check([b[2] for b in boxes] == [3, 7, 7, 9],
          "trace counts are 3 / 7 / 7 / 9 = 26")
    ws = {round(b[1][2] / 72.0, 6) for b in boxes}
    hs = {round(b[1][3] / 72.0, 6) for b in boxes}
    check(len(ws) == 1 and len(hs) == 1,
          "all four panels are the same size, %.4f x %.4f in"
          % (list(ws)[0], list(hs)[0]))
    check(abs(list(ws)[0] - fsn.PANEL_W) < 1e-4,
          "panel width %.4f in matches the derived PANEL_W %.4f in"
          % (list(ws)[0], fsn.PANEL_W))
    check(abs(list(hs)[0] - fsn.ROW_H) < 1e-4,
          "panel height %.4f in is ROW_H, unchanged from the parent figure"
          % list(hs)[0])

    xa, xb = boxes[0][1][0] / 72.0, boxes[1][1][0] / 72.0
    gap_w = xb - xa - list(ws)[0]
    ytop = (page_h - boxes[0][1][1] - boxes[0][1][3]) / 72.0
    ybot = (page_h - boxes[2][1][1] - boxes[2][1][3]) / 72.0
    gap_h = ytop - ybot - list(hs)[0]
    print("\n   measured gaps and margins (in), and what sets each:")
    B = fsn.B
    rows = [
        ("left",   xa, fsn.MARGIN["left"],
         "panel letter 0.420 + trim %.3f" % fsn.TRIM),
        ("right",  fsn.FIG_W - (boxes[1][1][0] + boxes[1][1][2]) / 72.0,
         fsn.MARGIN["right"],
         "half of \"20\" %.4f + trim %.3f" % (B["xtick_half"], fsn.TRIM)),
        ("top",    fsn.FIG_H - (ytop + list(hs)[0]), fsn.MARGIN["top"],
         "T_BLOCK %.4f + trim %.3f" % (B["t_block"], fsn.TRIM)),
        ("bottom", ybot, fsn.MARGIN["bottom"],
         "X_BLOCK %.4f + cbar block %.4f" % (B["x_block"], fsn.CB_BLOCK)),
        ("GAP_W",  gap_w, fsn.GAP_W,
         "letter 0.420 + \"20\" overhang %.4f + clear %.3f"
         % (B["xtick_half"], fsn.GAP_W_CLEAR)),
        ("GAP_H",  gap_h, fsn.GAP_H,
         "X_BLOCK %.4f + T_BLOCK %.4f + clear %.3f"
         % (B["x_block"], B["t_block"], fsn.GAP_H_CLEAR)),
    ]
    for name, meas, derived, why in rows:
        print("     %-7s measured %.4f   derived %.4f   = %s"
              % (name, meas, derived, why))
        check(abs(meas - derived) < 2e-3,
              "%s matches its derivation to 0.002 in" % name)
    print("     %-7s measured %s" % ("blocks", ""))
    for k in ("tick_h", "axlab_h", "title_h", "ytick_w", "xtick_half",
              "x_block", "y_block", "t_block"):
        print("        %-11s %.6f in" % (k, B[k]))
    return tp


def check_open_axes():
    print("\n[7] OPEN AXES")
    f = fsn.NatureSECFigure()
    f.load()
    f.draw()
    ok = True
    for L, ax in f.axes.items():
        vis = {s: ax.spines[s].get_visible() for s in
               ("left", "bottom", "top", "right")}
        nx = sum(1 for t in ax.xaxis.get_major_ticks() if t.tick2line.get_visible())
        ny = sum(1 for t in ax.yaxis.get_major_ticks() if t.tick2line.get_visible())
        print("   panel %s  spines %s   top-side x ticks %d   "
              "right-side y ticks %d" % (L, vis, nx, ny))
        ok &= (vis["left"] and vis["bottom"]
               and not vis["top"] and not vis["right"]
               and nx == 0 and ny == 0)
    check(ok, "all four panels: left+bottom drawn, top+right gone, "
              "no ticks on the two dropped sides")
    cb = f.cbar
    n_out = cb.outline.get_visible()
    print("   colourbar outline visible: %s, edge %s, linewidth %.2f pt"
          % (n_out, cb.outline.get_edgecolor(), cb.outline.get_linewidth()))
    check(bool(n_out), "the colourbar keeps its own closed rectangle -- it is "
                       "a swatch strip, not a plot")
    p = f.cax.get_position()
    fw, fh = f.style.figsize
    span = 2 * fsn.PANEL_W + fsn.GAP_W
    print("   colourbar x %.4f .. %.4f in, span %.4f in; axes extent %.4f in"
          % (p.x0 * fw, p.x1 * fw, p.width * fw, span))
    check(abs(p.width * fw - span) < 2e-3 and abs(p.x0 * fw - fsn.MARGIN["left"]) < 2e-3,
          "colourbar spans exactly the axes extent, so it stays visually "
          "attached with the right spines gone")
    return f


# ==========================================================================
# 8  renderer-level text collisions
# ==========================================================================
def check_collisions(f):
    """Renderer-level bbox audit.

    MUST run inside ``gs.apply``.  A Text resolves its font file at DRAW time
    from the live rcParams, and ``draw()`` has already exited the style
    context by the time we get here -- measured outside it, matplotlib falls
    back to DejaVu Sans, whose digits are 0.636 em against Arial's 0.556, and
    every tick label comes out 14% too wide.  ``save()`` re-enters the
    context, so Arial is what the delivered files were laid out with.
    """
    print("\n[8] TEXT COLLISIONS -- renderer bboxes, every string on the page")
    fig = f.fig
    from matplotlib import font_manager as _fm
    print("   font resolved for metrics: %s"
          % os.path.basename(_fm.findfont(
              _fm.FontProperties(family=["Arial", "DejaVu Sans"]))))
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

    for t in fig.texts:
        add(t, "figure")
    axes_all = list(f.axes.items()) + [("cbar", f.cax)]
    def in_view(axis, lo, hi):
        """Tick labels matplotlib actually DRAWS.

        A MultipleLocator hands back ticks past the view limits (x = 4 and 22
        for an xlim of (6, 20)); the Tick objects exist and report
        get_visible() True, but ``Axis._update_ticks`` filters them out before
        anything is rendered, so they are not on the page and must not be
        tested against the trim box.
        """
        out = []
        for loc, lab in zip(axis.get_majorticklocs(), axis.get_ticklabels()):
            if lo - 1e-9 <= loc <= hi + 1e-9:
                out.append(lab)
        return out

    for L, ax in axes_all:
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

    def overlap(a, b):
        dx = min(a.x1, b.x1) - max(a.x0, b.x0)
        dy = min(a.y1, b.y1) - max(a.y0, b.y0)
        return dx, dy

    hits = []
    worst = (1e9, None)
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            dx, dy = overlap(items[i][2], items[j][2])
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
    (sep, pair) = worst
    print("   tightest non-touching pair: %r (%s) / %r (%s), %.4f in = %.2f mm"
          % (pair[0][1], pair[0][0], pair[1][1], pair[1][0], sep, sep * 25.4))
    check(sep * 25.4 >= 0.5,
          "tightest text pair is %.2f mm apart (>= 0.50 mm)" % (sep * 25.4))

    # everything must also stay inside the trim box
    out = [(o, s, bb) for o, s, bb in items
           if bb.x0 < fsn.TRIM - 1e-6 or bb.y0 < fsn.TRIM - 1e-6
           or bb.x1 > fsn.FIG_W - fsn.TRIM + 1e-6
           or bb.y1 > fsn.FIG_H - fsn.TRIM + 1e-6]
    for o, s, bb in out:
        print("   OUTSIDE TRIM  %r (%s) at [%.4f %.4f %.4f %.4f]"
              % (s, o, bb.x0, bb.y0, bb.x1, bb.y1))
    check(not out, "every string sits inside the %.2f mm trim box"
                   % (fsn.TRIM * 25.4))


# ==========================================================================
# 9  data identity
# ==========================================================================
def check_data(nat_tp, nat_h):
    print("\n[9] DATA IDENTITY -- natstyle vs Figure_SEC_normalized")
    ref_root, _, ref_h = load_svg(REF + ".svg")
    ref_tp = trace_paths(ref_root)
    check(len(nat_tp) == len(ref_tp) == 4,
          "both SVGs carry 4 trace-bearing axes")
    worst = 0.0
    worst_at = None
    total = 0
    print("   panel  traces  vertices    max |dV| (mL)   max |dA| (norm)")
    for L, (gn, fn_), (gr, fr) in zip("abcd", nat_tp, ref_tp):
        if len(fn_) != len(fr):
            check(False, "panel %s: %d traces vs %d" % (L, len(fn_), len(fr)))
            continue
        mv = ma = 0.0
        nv = 0
        for (rn, vn, _), (rr, vr, _) in zip(fn_, fr):
            if len(vn) != len(vr):
                check(False, "panel %s: vertex count %d vs %d"
                             % (L, len(vn), len(vr)))
                continue
            dn = to_data(vn, rn, nat_h)
            dr = to_data(vr, rr, ref_h)
            d = np.abs(dn - dr)
            mv = max(mv, float(d[:, 0].max()))
            ma = max(ma, float(d[:, 1].max()))
            nv += len(vn)
            total += len(vn)
        print("     %s      %2d    %7d      %.3e       %.3e"
              % (L, len(fn_), nv, mv, ma))
        if max(mv, ma) > worst:
            worst, worst_at = max(mv, ma), L
    print("   %d vertices compared in total; worst deviation %.3e "
          "(panel %s)" % (total, worst, worst_at))
    # 6 decimal places of SVG coordinate precision on a ~206 pt axis is
    # 14 mL * 1e-6 / 206 = 6.8e-8 mL, so anything at 1e-7 is quantisation.
    check(worst < 1e-6, "every plotted vertex agrees to better than 1e-6 "
                        "(SVG coordinate quantisation alone is ~7e-8); "
                        "the restyle moved no data")

    # and the colours DID change, which is the point
    nat_c = [m.group(1).lower() for _, f_ in nat_tp for _, _, m in f_ if m]
    ref_c = [m.group(1).lower() for _, f_ in ref_tp for _, _, m in f_ if m]
    from secfigs import timecolor as tc
    check(set(nat_c) <= {c.lower() for c in tc.COLORS_MAGMA},
          "all 26 natstyle strokes come from COLORS_MAGMA")
    check(set(ref_c) <= {c.lower() for c in tc.COLORS},
          "all 26 reference strokes still come from the cividis COLORS")
    print("   natstyle distinct stroke colours: %d, reference: %d"
          % (len(set(nat_c)), len(set(ref_c))))


# ==========================================================================
def main():
    print("=" * 74)
    print("verify_nature.py -- %s" % NAT)
    print("=" * 74)
    root, w_pt, h_pt = check_page()
    check_fonts(root)
    check_edges()
    tp = check_geometry(root, h_pt)
    f = check_open_axes()
    with gs.apply(style=f.style, **f._rc()):
        check_collisions(f)
    plt.close(f.fig)
    check_data(tp, h_pt)
    print("\n" + "=" * 74)
    if _fails:
        print("%d CHECK(S) FAILED:" % len(_fails))
        for m in _fails:
            print("   - " + m)
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
