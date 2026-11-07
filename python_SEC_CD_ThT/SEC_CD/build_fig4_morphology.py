# *** SUPERSEDED — DO NOT BUILD FROM THIS FILE. ***
# Figure 3 (morphology) is now built by build_figure3_morphology.py, which is on the
# Nature contract: exactly 180 mm, Arial, live vector text, no bbox_inches="tight"
# (which silently changes the width), and the gel annotation redrawn as 7 pt vector.
# This file predates that and produces a non-conforming figure; kept only for history.
"""build_fig4_morphology.py — Figure 4 (size & morphology), revised per Shirin Kamalaldinezabadi
(email 2026-08-05).

She asked for the gridded cryo-EM zoom to be dropped ("does not add sufficient value and may
potentially confuse readers") and replaced with a negative-stain image, keeping the cryo-EM
image WITHOUT gridlines:
    a  silver-stained SDS-PAGE
    b  negative-stain TEM   <- her `Asset 1.png` (the version she prepared for the paper)
    c  cryo-EM, no gridlines <- `FromShirin/image2.png`
Both micrographs carry their own 100 nm scale bars.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from PIL import Image
import figstyle as FS

WRITING = r"/path/to/Ab42oligomers37C/code"
SRC = r"/path/to/Ab42oligomers37C/lab_main_folder"
OUT = WRITING + r"/37C_Oligomers_NatComms_Package/Figures_Final"

GEL   = f"{SRC}/Data/SDS/SDS_gel.png"
NEGST = f"{WRITING}/Asset 1.png"                       # Shirin's negative stain
CRYO  = f"{SRC}/Writing/Inkscape/FromShirin/image2.png"  # cryo-EM, no gridlines

def load(p):
    return np.asarray(Image.open(p).convert("RGB"))

def img_panel(ax, arr, tag, lx=-0.02):
    ax.imshow(arr); ax.axis("off")
    ax.text(lx, 1.02, tag, transform=ax.transAxes, **FS.LETTER)

gel, neg, cryo = load(GEL), load(NEGST), load(CRYO)

fig = FS.new_fig(13.0, 6.4)
gs = gridspec.GridSpec(1, 3, width_ratios=[1.05, 1.0, 1.0], wspace=0.04,
                       left=0.02, right=0.99, top=0.93, bottom=0.02)
img_panel(fig.add_subplot(gs[0]), gel,  "a", lx=-0.01)
img_panel(fig.add_subplot(gs[1]), neg,  "b")
img_panel(fig.add_subplot(gs[2]), cryo, "c")

for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/Figure4_morphology.{ext}", dpi=600, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("wrote Figure4_morphology (gel | negative stain | cryo-EM without gridlines)")
