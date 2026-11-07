# -*- coding: utf-8 -*-
"""Extract data traces from a matplotlib-exported SEC SVG (Inkscape-saved, relative paths)."""
import re, sys, io
import numpy as np
from lxml import etree
SVGNS='http://www.w3.org/2000/svg'
def q(t): return '{%s}%s'%(SVGNS,t)

def _use_xy(group):
    u=group.find('.//'+q('use'))
    if u is None: return None
    return float(u.get('x')), float(u.get('y'))

def _tick_value(group):
    # first numeric comment inside the group
    for el in group.iter():
        prev=el.getprevious()
    for node in group.iter():
        pass
    # comments: iterate and collect
    for c in group.iter(etree.Comment):
        s=(c.text or '').strip()
        if re.fullmatch(r'-?\d+(\.\d+)?', s): return float(s)
    return None

def parse_ticks(root, prefix):
    ticks=[]
    for g in root.iter(q('g')):
        gid=g.get('id') or ''
        if gid.startswith(prefix):
            xy=_use_xy(g); val=_tick_value(g)
            if xy is not None and val is not None:
                ticks.append((xy, val))
    return ticks

def parse_rel_path(d):
    """Parse an Inkscape relative path 'm x,y dx,dy dx,dy ...' into absolute pixel points."""
    toks=re.findall(r'[-+]?\d*\.?\d+(?:e[-+]?\d+)?', d)
    nums=list(map(float,toks))
    # command letters
    cmds=re.findall(r'[MmLl]', d)
    pts=[]; x=y=0.0; i=0
    # handle possibly mixed; assume starts with m/M then pairs
    # We'll re-walk using regex over commands+coords
    it=re.finditer(r'([MmLl])?\s*([-+]?\d*\.?\d+(?:e[-+]?\d+)?)[ ,]+([-+]?\d*\.?\d+(?:e[-+]?\d+)?)', d)
    rel=True; started=False
    for m in it:
        cmd=m.group(1); dx=float(m.group(2)); dy=float(m.group(3))
        if cmd in ('M','L'): rel=False
        elif cmd in ('m','l'): rel=True
        if not started:
            x,y=dx,dy; started=True  # first pair after m/M is the start (relative m first pair is absolute in SVG)
        else:
            if rel: x+=dx; y+=dy
            else: x,y=dx,dy
        pts.append((x,y))
    return np.array(pts)

def calib(ticks):
    P=np.array([t[0] for t in ticks]); V=np.array([t[1] for t in ticks])
    A=np.vstack([P,np.ones_like(P)]).T
    m,b=np.linalg.lstsq(A,V,rcond=None)[0]
    return m,b

def _yticks_fallback(root):
    """Find y-axis tick (pixel_y, value) when there are no ytick_ groups:
    scan matplotlib.axis_2 (the y axis) for <use> marks and pair with numeric comments 0.0/1.0..."""
    ax=None
    for g in root.iter(q('g')):
        if (g.get('id') or '')=='matplotlib.axis_2': ax=g; break
    if ax is None: return []
    uses=[u for u in ax.iter(q('use'))]
    ys=sorted(set(round(float(u.get('y')),2) for u in uses if u.get('y')))
    vals=sorted(set(float((c.text or '').strip()) for c in ax.iter(etree.Comment)
                    if re.fullmatch(r'-?\d+(\.\d+)?',(c.text or '').strip())))
    if len(ys)>=2 and len(vals)>=2:
        # ys ascending in pixels = descending in value (y flipped): pair max-pixel->min-val
        ys_desc=sorted(ys, reverse=True)
        return list(zip([(0,y) for y in ys_desc][:len(vals)], vals))
    return []

def extract(svgpath, min_verts=40):
    root=etree.parse(svgpath).getroot()
    xt=parse_ticks(root,'xtick_'); yt=parse_ticks(root,'ytick_')
    xtx=[( t[0][0], t[1]) for t in xt]; yty=[(t[0][1], t[1]) for t in yt]
    if len(yty)<2:
        yfb=_yticks_fallback(root)
        yty=[(xy[1], v) for xy,v in yfb]
    mx,bx=calib(xtx); my,by=calib(yty)
    traces=[]
    for p in root.iter(q('path')):
        st=p.get('style') or ''
        col=re.search(r'stroke:\s*(#[0-9a-fA-F]{6})', st)
        if not col:
            par=p.getparent()
            if par is not None:
                col=re.search(r'stroke:\s*(#[0-9a-fA-F]{6})', par.get('style') or '')
        if not col: continue
        color=col.group(1).lower()
        d=p.get('d') or ''
        if not d: continue
        pts=parse_rel_path(d)
        # black data traces are long (>400 verts); black frame/ticks/grid are short -> keep only long black
        thresh = 400 if color=='#000000' else min_verts
        if len(pts)<thresh: continue
        X=mx*pts[:,0]+bx; Y=my*pts[:,1]+by
        traces.append((color, X, Y, len(pts)))
    return traces, (mx,bx,my,by)

def extract_pix(svgpath, min_verts=40, black_thresh=400):
    """Like extract() but returns RAW pixel-Y (not y-calibrated), for normalized plots that
    have no y-axis calibration. Returns list of (color, X_calibrated, Ypix_raw, n)."""
    root=etree.parse(svgpath).getroot()
    xt=parse_ticks(root,'xtick_')
    xtx=[(t[0][0], t[1]) for t in xt]
    mx,bx=calib(xtx)
    out=[]
    for p in root.iter(q('path')):
        st=p.get('style') or ''
        col=re.search(r'stroke:\s*(#[0-9a-fA-F]{6})', st)
        if not col:
            par=p.getparent()
            if par is not None:
                col=re.search(r'stroke:\s*(#[0-9a-fA-F]{6})', par.get('style') or '')
        if not col: continue
        color=col.group(1).lower()
        d=p.get('d') or ''
        if not d: continue
        pts=parse_rel_path(d)
        thresh = black_thresh if color=='#000000' else min_verts
        if len(pts)<thresh: continue
        X=mx*pts[:,0]+bx
        out.append((color, X, pts[:,1], len(pts)))
    return out, (mx,bx)

def norm_trace(Ypix):
    """Flip SVG pixel-Y (down=+) and normalize to [0,1]: baseline->0, peak->1."""
    import numpy as _np
    y = Ypix.max() - Ypix          # flip so peak (small pixel) becomes large
    rng = y.max() - y.min()
    return y/ rng if rng>0 else y*0.0

if __name__=='__main__':
    sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding='utf-8')
    f=r"/path/to/Ab42oligomers37C/lab_main_folder/Writing/Inkscape/WT_37c_100uM_normalized.svg"
    tr,cal=extract(f)
    print("calibration (mx,bx,my,by):", [round(c,4) for c in cal])
    print("traces extracted:", len(tr))
    for color,X,Y,n in tr:
        print(f"  {color}: {n} pts, x[{X.min():.1f},{X.max():.1f}] y[{Y.min():.2f},{Y.max():.2f}]")
    # verification plot
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    plt.figure(figsize=(6,4))
    for color,X,Y,n in tr: plt.plot(X,Y,color=color,lw=1)
    plt.xlabel("Elution Volume (ml)"); plt.ylabel("Normalized Absorbance"); plt.xlim(6,20)
    plt.title("EXTRACTED WT 100uM 37C")
    plt.savefig(r"/path/to/output.png",dpi=130,bbox_inches='tight')
    print("wrote verify_sec.png")
