"""build_fig8_cryo.py — original-numbering Figure 8 (cryo-EM), revised per Shirin
Kamalaldinezabadi's email of 5 Aug 2026.

She asked to drop the gridded cryo-EM zoom ("does not add sufficient value and may potentially
confuse readers") and instead show:
    A  negative-stain TEM   (left)   <- her `Asset 1.png`
    B  cryo-EM, no gridlines (right) <- `FromShirin/image2.png`
Both carry their own 100 nm scale bars. In this (pre-Nature-formatting) version the
silver-stained gel remains its own Figure 7, so Figure 8 is cryo/negative-stain only.
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

NEGST = f"{WRITING}/Asset 1.png"
CRYO  = f"{SRC}/Writing/Inkscape/FromShirin/image2.png"

def load(p): return np.asarray(Image.open(p).convert("RGB"))

def panel(ax, arr, tag):
    ax.imshow(arr); ax.axis("off")
    ax.text(-0.02, 1.02, tag, transform=ax.transAxes, **FS.LETTER)

fig = FS.new_fig(11.0, 5.6)
gs = gridspec.GridSpec(1, 2, wspace=0.05, left=0.02, right=0.99, top=0.93, bottom=0.02)
panel(fig.add_subplot(gs[0]), load(NEGST), "A")
panel(fig.add_subplot(gs[1]), load(CRYO),  "B")
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/Figure8_cryoEM_revised.{ext}", dpi=600, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("wrote Figure8_cryoEM_revised (A negative stain | B cryo-EM without gridlines)")
