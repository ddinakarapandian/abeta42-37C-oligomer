"""Figure 1 (SEC) — rebuilt NATIVELY: concentration series (250/100/25 uM, 37 C) as three panels that
share one typography and plot box; viridis time coloring; panel letters in the margin. Data = digitized
from Tarunya's Inkscape SVGs (raw AKTA exports absent from the folder)."""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
sys.path.insert(0, ".")
from digitize import extract, extract_pix, norm_trace
import figstyle as FS

SVG=r"/path/to/Ab42oligomers37C/lab_main_folder/Writing/Inkscape/"
OUT=r"/path/to/Ab42oligomers37C/manuscript/37C_Oligomers_NatComms_Package/Figures_Final"

def _oligo(T):
    X,Y=T[1],T[2]; m=(X>=7.5)&(X<=10.5); return Y[m].max() if m.any() else -1
def traces(svg, normalized=False):
    if normalized:
        tr,_=extract_pix(SVG+svg); tr=[(c,X,norm_trace(Yp),n) for c,X,Yp,n in tr]
    else:
        tr,_=extract(SVG+svg)
    tr=sorted(tr,key=_oligo)
    out=[]
    for c,X,Y,n in tr:
        o=np.argsort(X); x,y=X[o],Y[o]
        mx=np.nanmax(y);  y=y/mx if mx>0 else y     # per-trace normalize (peak=1), like Daniel's list_plot
        out.append((x,y))
    return out

def sec_panel(ax, tr, times, ylab=True):
    n=len(tr); cols=[plt.cm.viridis(v) for v in np.linspace(0.05,0.92,n)]
    labs=times if len(times)==n else [f"t{i+1}" for i in range(n)]
    for (x,y),c,lab in zip(tr,cols,labs):
        ax.plot(x,y,color=c,lw=FS.LW,label=lab)
    ax.set_xlim(6,20); ax.set_ylim(-0.03,1.05)
    ax.set_xticks(range(6,21,2)); ax.set_yticks(np.arange(0,1.01,0.2))
    ax.set_xlabel("Elution Volume (ml)")
    if ylab: ax.set_ylabel("Normalized Absorbance")
    FS.apply_axes(ax)
    leg=ax.legend(fontsize=FS.LEG_FS, frameon=False, handlelength=1.3, labelspacing=0.3,
                  borderaxespad=0.3, loc="upper right")
    for t in leg.get_texts(): t.set_fontfamily(FS.FONT)

fig=FS.new_fig(16.5,5.0)
gs=gridspec.GridSpec(1,3,wspace=0.26,left=0.055,right=0.99,top=0.90,bottom=0.17)
T7=["0 h","2 h","4 h","8 h","12 h","24 h","48 h"]
for i,(svg,times,norm,tag) in enumerate([
        ("WT_4c_250uM.svg",["2 h","4 h","6 h"],False,"a"),
        ("WT_37c_100uM_normalized.svg",T7,False,"b"),
        ("WT_37C_25uM.svg",T7,True,"c")]):
    ax=fig.add_subplot(gs[i]); sec_panel(ax, traces(svg,norm), times, ylab=(i==0)); FS.panel_letter(ax,tag)

fig.savefig(OUT+"/Figure1_SEC.png",dpi=600,bbox_inches="tight",facecolor="white")
fig.savefig(OUT+"/Figure1_SEC.pdf",bbox_inches="tight",facecolor="white")
print("wrote Figure1_SEC (native, 250/100/25 uM)")
