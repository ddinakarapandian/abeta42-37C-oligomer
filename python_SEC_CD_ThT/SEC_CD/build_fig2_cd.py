"""Figure 2 (CD) — rebuilt NATIVELY into one figure (ROX321 style): three CD panels share one
typography and one plot box; panel letters sit in the margin. Data = Tarunya's 100325 CD, savgol-smoothed."""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.signal import savgol_filter
import figstyle as FS

CD = r"/path/to/Ab42oligomers37C/lab_main_folder/Data/AbetaWT/100325/"
OUT = r"/path/to/Ab42oligomers37C/manuscript/37C_Oligomers_NatComms_Package/Figures_Final"
WL=(190,250)
def load(fn):
    df=pd.read_csv(CD+fn, header=None, names=["wl","cd"]).apply(pd.to_numeric, errors="coerce").dropna()
    df=df.sort_values("wl"); df["cd"]=savgol_filter(df["cd"].values, 11, 2)
    return df[(df.wl>=WL[0])&(df.wl<=WL[1])]
S={k:load(v) for k,v in {"mono":"37c mono avg-autosub.csv","oligo":"37c oligo avg-autosub.csv",
    "6c":"6c mix avg-autosub.csv","25c":"25c mix avg-autosub.csv","37c":"37c mix avg-autosub.csv"}.items()}

def cd_panel(ax, traces, ylab=True):
    for x,y,c,ls,lab in traces:
        ax.plot(x,y,color=c,lw=FS.LW,ls=ls,label=lab)
    ax.axhline(0,color="k",lw=1.0,ls=":")
    ax.set_xlim(*WL); ax.set_xticks(range(190,251,10)); ax.set_xticks(range(190,251,2),minor=True)
    ax.set_xlabel("Wavelength (nm)")
    if ylab: ax.set_ylabel("CD (mdeg)")
    FS.apply_axes(ax)
    leg=ax.legend(fontsize=FS.LEG_FS, frameon=False, handlelength=1.4, borderaxespad=0.4)
    for t in leg.get_texts(): t.set_fontfamily(FS.FONT)

fig=FS.new_fig(16.5, 5.0)
gs=gridspec.GridSpec(1,3,wspace=0.28,left=0.06,right=0.99,top=0.90,bottom=0.17)
BLUE,RED,BLACK="#1f4e79","#c0392b","#000000"
ax=fig.add_subplot(gs[0]); cd_panel(ax,[
    (S["6c"].wl,S["6c"].cd,BLACK,"-","6 °C"),
    (S["25c"].wl,S["25c"].cd,BLUE,"--","25 °C"),
    (S["37c"].wl,S["37c"].cd,RED,"-","37 °C")]); FS.panel_letter(ax,"a")
ax=fig.add_subplot(gs[1]); cd_panel(ax,[
    (S["mono"].wl,S["mono"].cd,BLUE,"-","Monomer"),
    (S["oligo"].wl,S["oligo"].cd,RED,"-","37 °C oligomer")]); FS.panel_letter(ax,"b")
ax=fig.add_subplot(gs[2]); cd_panel(ax,[
    (S["oligo"].wl,S["oligo"].cd,RED,"-","37 °C oligomer")]); FS.panel_letter(ax,"c")

fig.savefig(OUT+"/Figure2_CD.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(OUT+"/Figure2_CD.pdf", bbox_inches="tight", facecolor="white")
print("wrote Figure2_CD (native)")
