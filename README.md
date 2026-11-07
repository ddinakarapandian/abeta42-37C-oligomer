# Aβ42 37 °C oligomer: analysis code and NMR raw data

Code and raw NMR data for the manuscript "Stimulating microglial neuroinflammation with amyloid-beta oligomers compatible with solid-state NMR structural characterization" (Dinakarapandian, Rao Sudarshan, Griffin et al.).

## Contents
| Folder | What it is |
|---|---|
| `cytokines_R/` | R scripts for the multiplexed cytokine and principal component analyses, and the cytokine data table they read (`Monomer and Oligomer WT ABeta.xlsx`). Main script: `Mono Oligo WT Abeta Cytokine.R`. |
| `python_SEC_CD_ThT/` | Python scripts used to plot the circular dichroism (CD), size-exclusion chromatography (SEC) and thioflavin-T (ThT) data. See `python_SEC_CD_ThT/README.md`. |
| `mathematica_NMR/` | Mathematica notebooks for NMR peak fitting, secondary chemical-shift analysis and chemical-shift RMSD calculations (2D DARR, CHHC, CPMAS). Saved without their cached output, so they must be evaluated to regenerate plots. |

## NMR raw data
The Bruker raw data (CPMAS, DARR 20 ms and 500 ms, CHHC) are in `NMR_raw_data.zip`, attached to the [v1.0.0 release](../../releases/tag/v1.0.0). It contains a `README.txt` that describes each folder. The notebooks build their data paths from a `Dropbox` variable plus a folder name; set these to the location of the unzipped data.

## Licence
Code is released under the MIT licence (`LICENSE`). The NMR raw data are released under CC BY 4.0.

## Notes
- The scripts are provided as used in the study. R packages are loaded at the top of the R script (readxl, tidyverse, ropls and others).
- Processed values behind each figure are in the Source Data file of the paper.
