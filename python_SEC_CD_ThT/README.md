# Python plotting scripts (CD, SEC, ThT)

These are the scripts used to plot the CD, SEC and ThT data of the manuscript.

- `SEC_CD/` has the SEC and CD figure scripts and the `secfigs` helper package they import. Script names and some docstrings use an earlier figure numbering; the manuscript figures are Fig. 1b-e (SEC, CD, ThT), Fig. 3 (SEC) and the matching Supplementary Figures.
- `ThT/` has the ThT scripts (`tht_plateau.py`, `tht_panels.py`, `tht_bars.py` and the figure builders).
- `shared/` has `sec1.py` and `cd_spectra_workflow.py`, which the SEC and CD scripts import.

## Before running
- Paths are placeholders (`/path/to/Ab42oligomers37C/...`) and must be edited to the location of the instrument exports.
- The scripts read the original instrument exports (ÄKTA/UNICORN SEC results, CD CSV files, the ThT plate-reader workbook). These raw exports are not in this repository; the processed values that were plotted are in the paper's Source Data file.
- Python 3.12 with numpy, pandas, scipy, matplotlib, openpyxl, Pillow and lxml.
