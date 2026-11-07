"""fig1_cd.py — Figure 1 CD panels for the 37 C oligomer paper, generated with DANIEL'S CD code
(cd_spectra_workflow.CDSpectraWorkflow.plot_publication_overlay) applied to Tarunya's 100325 CD data.
"""
import sys, os
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- import Daniel's CD plotting code ---
DAN = r"/path/to/Ab42oligomers37C/code/shared/griffin_051926"
sys.path.insert(0, DAN)
from cd_spectra_workflow import CDSpectraWorkflow

CD = r"/path/to/Ab42oligomers37C/lab_main_folder/Data/AbetaWT/100325/"
OUT = r"/path/to/Ab42oligomers37C/code/37C_SEC_CD_Python/figures"

from scipy.signal import savgol_filter
def load(fn):
    # Match Tarunya's CD_Analysis_100325.ipynb: read, coerce, sort, Savitzky-Golay smooth (win=11, poly=2)
    df = pd.read_csv(CD+fn, header=None, names=["wavelength_nm", "cd_mdeg"])
    df = df.apply(pd.to_numeric, errors="coerce").dropna().sort_values("wavelength_nm").reset_index(drop=True)
    df["cd_mdeg"] = savgol_filter(df["cd_mdeg"].values, window_length=11, polyorder=2)
    return df

spec = {k: load(v) for k, v in {
    "mono": "37c mono avg-autosub.csv", "oligo": "37c oligo avg-autosub.csv",
    "6c": "6c mix avg-autosub.csv", "25c": "25c mix avg-autosub.csv", "37c": "37c mix avg-autosub.csv",
}.items()}

wf = CDSpectraWorkflow(root_dir=CD)   # empty auto-discovery (no WT_/Deam_ subfolders) — fine

def save(fig, stem):
    fig.savefig(f"{OUT}/{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT}/{stem}.svg", bbox_inches="tight")
    plt.close(fig); print("wrote", stem)

# (a) monomer vs 37 C oligomer
fig, ax = wf.plot_publication_overlay(
    {"Monomer": spec["mono"], "37 °C oligomer": spec["oligo"]},
    colors={"Monomer": "#1f4e79", "37 °C oligomer": "#c0392b"},
    wavelength_range=(190, 250), ylabel="CD (mdeg)", zero_line=True)
save(fig, "Fig1a_CD_monomer_vs_oligomer")

# (b) temperature overlay 6/25/37 C
fig, ax = wf.plot_publication_overlay(
    {"6 °C": spec["6c"], "25 °C": spec["25c"], "37 °C": spec["37c"]},
    colors={"6 °C": "#000000", "25 °C": "#1f4e79", "37 °C": "#c0392b"},
    linestyles={"25 °C": "--"},
    wavelength_range=(190, 250), ylabel="CD (mdeg)", zero_line=True)
save(fig, "Fig1b_CD_temperature_overlay")

# (c) 37 C oligomer alone (beta-strand)
fig, ax = wf.plot_publication_overlay(
    {"37 °C oligomer": spec["oligo"]},
    colors={"37 °C oligomer": "#c0392b"},
    wavelength_range=(190, 250), ylabel="CD (mdeg)", zero_line=True)
save(fig, "Fig1c_CD_oligomer")
print("DONE CD panels")
