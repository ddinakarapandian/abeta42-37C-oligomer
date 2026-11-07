"""compose.py — assemble Figure 1 (CD) and Figure 2 (SEC) from the panels made with Daniel's code."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
F = "figures"
def img(p): return np.asarray(Image.open(p).convert("RGB"))
def lab(ax, s): ax.text(0.005, 0.995, s, transform=ax.transAxes, fontsize=15, fontweight="bold", va="top", ha="left")
def panel(ax, p, s): ax.imshow(img(p)); ax.axis("off"); lab(ax, s)
MM = 1/25.4

# Figure 1 — CD (a monomer/oligomer, b temperature overlay, c oligomer)
fig = plt.figure(figsize=(180*MM, 64*MM))
gs = fig.add_gridspec(1, 3, wspace=0.02)
panel(fig.add_subplot(gs[0]), f"{F}/Fig1a_CD_monomer_vs_oligomer.png", "a")
panel(fig.add_subplot(gs[1]), f"{F}/Fig1b_CD_temperature_overlay.png", "b")
panel(fig.add_subplot(gs[2]), f"{F}/Fig1c_CD_oligomer.png", "c")
fig.savefig(f"{F}/Figure1_CD.png", dpi=450, bbox_inches="tight", facecolor="white")
fig.savefig(f"{F}/Figure1_CD.pdf", bbox_inches="tight", facecolor="white")
plt.close(fig); print("wrote Figure1_CD")

# Figure 2 — SEC (a 250 uM, b 100 uM at 37 C)  [25 uM / 25 C pending raw AKTA data]
fig = plt.figure(figsize=(180*MM, 74*MM))
gs = fig.add_gridspec(1, 2, wspace=0.04)
panel(fig.add_subplot(gs[0]), f"{F}/Fig2a_SEC_250uM_37C.png", "a")
panel(fig.add_subplot(gs[1]), f"{F}/Fig2b_SEC_100uM_37C.png", "b")
fig.savefig(f"{F}/Figure2_SEC.png", dpi=450, bbox_inches="tight", facecolor="white")
fig.savefig(f"{F}/Figure2_SEC.pdf", bbox_inches="tight", facecolor="white")
plt.close(fig); print("wrote Figure2_SEC")
