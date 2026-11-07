# *** SUPERSEDED — DO NOT BUILD FROM THIS FILE. ***
# Figure 3 (morphology) is now built by build_figure3_morphology.py, which is on the
# Nature contract: exactly 180 mm, Arial, live vector text, no bbox_inches="tight"
# (which silently changes the width), and the gel annotation redrawn as 7 pt vector.
# This file predates that and produces a non-conforming figure; kept only for history.
"""fig3_morphology.py — Figure 3 (size & morphology): SDS-PAGE + cryo-EM.
These are micrographs/photos (not SEC/CD data), so they are composed, not replotted.
(a) SDS-PAGE gel  (b) cryo-EM 100 nm field  (c) cryo-EM 10 nm zoom.
The embedded caption band on the 'Finalized' image is cropped off (it belongs in the legend).
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

SRC = r"/path/to/Ab42oligomers37C/lab_main_folder"
OUT = r"/path/to/Ab42oligomers37C/code/37C_SEC_CD_Python/figures"
GEL   = f"{SRC}/Data/SDS/SDS_gel.png"
CRYO1 = f"{SRC}/Writing/Inkscape/FromShirin/image2.png"                 # 100 nm field
CRYO2 = f"{SRC}/Writing/Inkscape/FromShirin/Finalized cryo-EM Image.png"  # 10 nm zoom (+caption)

def load(p, crop_bottom=0):
    im = Image.open(p).convert("RGB")
    if crop_bottom:
        w, h = im.size
        im = im.crop((0, 0, w, h - crop_bottom))
    return np.asarray(im)

def panel(ax, arr, s):
    ax.imshow(arr); ax.axis("off")
    ax.text(0.01, 0.99, s, transform=ax.transAxes, fontsize=15, fontweight="bold",
            va="top", ha="left", color="black",
            bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.85))

MM = 1/25.4
gel  = load(GEL)
c100 = load(CRYO1)
c10  = load(CRYO2, crop_bottom=22)   # drop the "Single particles ~6-7 nm" caption band

fig = plt.figure(figsize=(180*MM, 92*MM))
gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0], height_ratios=[1, 1],
                      wspace=0.04, hspace=0.06)
axg = fig.add_subplot(gs[:, 0])   # gel spans both rows
panel(axg, gel, "a")
panel(fig.add_subplot(gs[0, 1]), c100, "b")
panel(fig.add_subplot(gs[1, 1]), c10,  "c")

fig.savefig(f"{OUT}/Figure3_morphology.png", dpi=450, bbox_inches="tight", facecolor="white")
fig.savefig(f"{OUT}/Figure3_morphology.pdf", bbox_inches="tight", facecolor="white")
plt.close(fig); print("wrote Figure3_morphology")
