"""fig2_sec.py — Figure 2 SEC panels via DANIEL'S SEC code (sec1.SECData.list_plot).

NOTE: the raw ÄKTA/UNICORN exports for the 37 C 250/100/25 µM concentration series are NOT in the
Tarunya folder (only the finished SVGs exist). So the SEC traces are recovered from the SVGs with the
verified digitizer (37C_NMR_Python/digitize.py) and fed to Daniel's SECData as DataFrames.
The 250 µM and 100 µM SVGs digitize cleanly; the 25 µM / 25 C SVGs lack axis calibration in the SVG
(single/no y-tick) and need the raw ÄKTA data to be plotted natively.
"""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

DAN = r"/path/to/Ab42oligomers37C/code/shared/griffin_051926"
NMRPY = r"/path/to/Ab42oligomers37C/code/37C_NMR_Python"
sys.path.insert(0, DAN); sys.path.insert(0, NMRPY)
from sec1 import SECData
from digitize import extract, extract_pix, norm_trace

SVG = r"/path/to/Ab42oligomers37C/lab_main_folder/Writing/Inkscape/"
OUT = r"/path/to/Ab42oligomers37C/code/37C_SEC_CD_Python/figures"

def _oligo(T):
    X, Y = T[1], T[2]
    m = (X >= 7.5) & (X <= 10.5)
    return Y[m].max() if m.any() else -1

def dfs_from_svg(svg, normalized=False):
    """normalized=True -> plot has no y-calibration (0-1 axis); use pixel-Y flipped+normalized per trace."""
    if normalized:
        tr, _cal = extract_pix(SVG+svg)
        tr = [(c, X, norm_trace(Yp), n) for c, X, Yp, n in tr]
    else:
        tr, _cal = extract(SVG+svg)
    tr = sorted(tr, key=_oligo)          # ascending oligomer-peak height ~ incubation time
    out = []
    for c, X, Y, n in tr:
        o = np.argsort(X)
        out.append(pd.DataFrame({"Volume_ml": X[o], "Absorbance_mAU": Y[o]}))
    return out

VIR = lambda n: [plt.cm.viridis(v) for v in np.linspace(0.05, 0.92, n)]

def panel(svg, times, stem, title, normalized=False):
    dfs = dfs_from_svg(svg, normalized=normalized)
    n = len(dfs); labels = times if len(times) == n else [f"t{i+1}" for i in range(n)]
    sd = SECData(dfs); sd.rename_columns(); sd.normalize_data()
    filt = sd.filter_by_volume(min_volume=6, max_volume=20); filt.normalize_data()
    filt.list_plot(normalized=True, xlabel="Elution Volume (ml)", ylabel="Normalized Absorbance",
                   legend_labels=labels, colors=[VIR(n)[i] for i in range(n)], linewidth=2.2,
                   x_range=(6, 20), figsize=(9, 7), dpi=600, output_file=f"{OUT}/{stem}.png")
    print("wrote", stem, f"({n} traces)")

# 250 / 100 / 25 uM at 37 C — concentration series (old Fig 3 A,B,C):
panel("WT_4c_250uM.svg", ["2 h","4 h","6 h"], "Fig2a_SEC_250uM_37C", "250 µM, 37 °C")
panel("WT_37c_100uM_normalized.svg", ["0 h","2 h","4 h","8 h","12 h","24 h","48 h"], "Fig2b_SEC_100uM_37C", "100 µM, 37 °C")
# 25 uM SVG is a normalized plot (no y-calibration) -> digitize via pixel-Y flip+normalize.
panel("WT_37C_25uM.svg", ["0 h","2 h","4 h","8 h","12 h","24 h","48 h"], "Fig2c_SEC_25uM_37C", "25 µM, 37 °C",
      normalized=True)
print("DONE 250/100/25 uM via Daniel's sec1. (Time labels for 25 uM inferred by oligomer-peak order.)")
