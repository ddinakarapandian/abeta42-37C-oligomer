"""figstyle.py — one shared typography for all main figures (ROX321-style: native panels, one type set,
panel letters in the top-left MARGIN, never over content, sized proportionate to the axis labels)."""
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib import font_manager

# pick Arial if present
for _f in ("Arial", "Helvetica", "DejaVu Sans"):
    try: font_manager.findfont(_f, fallback_to_default=False); FONT=_f; break
    except Exception: FONT="DejaVu Sans"

# font sizes are tuned so that, embedded at 180 mm double-column, everything lands ~7-9 pt (Nature ok)
LABEL_FS = 19      # axis titles
TICK_FS  = 16      # tick labels
LEG_FS   = 15      # legend
LETTER_FS= 24      # panel letter (bold) — ~1.25x the axis title, like ROX321
LW       = 2.6     # data line width
SPINE_LW = 1.6
TICK_LW  = 1.4
LETTER = dict(fontsize=LETTER_FS, fontweight="bold", fontfamily=FONT, va="top", ha="left")

def apply_axes(ax):
    """House style on one axes: Arial, bold-ish ticks all four sides pointing in, thick spines."""
    for s in ax.spines.values(): s.set_linewidth(SPINE_LW)
    ax.tick_params(which="major", direction="in", length=6, width=TICK_LW,
                   labelsize=TICK_FS, top=True, right=True)
    ax.tick_params(which="minor", direction="in", length=3, width=TICK_LW*0.8,
                   top=True, right=True)
    for lbl in ax.get_xticklabels()+ax.get_yticklabels():
        lbl.set_fontfamily(FONT)
    ax.xaxis.label.set_fontsize(LABEL_FS); ax.yaxis.label.set_fontsize(LABEL_FS)
    ax.xaxis.label.set_fontfamily(FONT);  ax.yaxis.label.set_fontfamily(FONT)

def panel_letter(ax, tag, x=-0.16, y=1.04):
    ax.text(x, y, tag, transform=ax.transAxes, **LETTER)

def new_fig(w, h):
    return plt.figure(figsize=(w, h))
