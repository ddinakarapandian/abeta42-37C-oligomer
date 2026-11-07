# -*- coding: utf-8 -*-
r"""
panel_data.py -- authoritative CSV manifest for the seven SEC panels of the
37 C Abeta42 oligomer manuscript (doc 10).

Every time-point label here is DERIVED FROM THE CSV FILENAME, never from Tarunya's
legend text and never from the published caption.  Where a filename cannot establish
a fact (temperature, concentration), the fact is left in `notes` as a flag rather
than guessed.

Traces are listed CHRONOLOGICALLY, which for several panels is NOT the order in
her notebook's file list.  Because label and path travel together in one tuple,
the positional legend/file-order desynchronisation that produced the Fig 3B and
Fig 4B mislabels cannot recur.

Processing contract (see PROVENANCE_FIGURES_1-4.md sec. 3):
    1. pd.read_csv(f, encoding='utf-16', sep='\t', skiprows=3)
    2. take the first two NUMERIC columns -> Volume_ml, Absorbance_mAU
    3. truncate 0 <= Volume_ml <= 25
    4. A = A - min(A) only if min(A) < 0
    5. if normalized: (A - min) / (max - min) over that same 0-25 mL window
No smoothing, no fitted baseline, no volume-axis transform.

Evidence tier used below:
    CONFIRMED  -- the curve in Tarunya's delivered vector artwork was matched to this
                  exact CSV by decoding the SVG <path d="..."> node coordinates,
                  inverse-mapping through the axes rectangle, and comparing point by
                  point against the reprocessed raw data.  RMSE and the margin over
                  the next-best candidate are quoted per panel.
    INFERRED   -- established from notebook code only (no surviving vector to check).

Written 2026-08-07.
"""

import os

# ---------------------------------------------------------------- roots (READ-ONLY)
RES = r"/path/to/Ab42oligomers37C/data/SEC\Run_Results"
WT = os.path.join(RES, "WT")
WT250 = os.path.join(WT, "06-04-25 Results")
WT25 = os.path.join(WT, "WT 25 uM 37 Inc")   # NOT used by any published panel -- see notes
E22K = os.path.join(RES, "E22K")             # holds the two WT cold-incubation runs
RT = os.path.join(RES, "RT runs")

# Filename traps that have bitten before -- the on-disk spellings, verbatim:
#   * "WT prim run 122025 001.csv"      SINGLE space before 001 on disk;
#                                       WT37COligomer.ipynb cell 24 types a DOUBLE space.
#                                       It is not in any panel, so nothing breaks -- but do
#                                       not copy the notebook string.
#   * "wt 24h reinj 500ul Loop no fractions  001 25uM.csv"
#                                       ONE "001" token; its six siblings have "001 001".
#   * "WT 72hr reinj 37 C 022225 ..."   space inside "37 C" (siblings write "37C").
#   * Several names carry a DOUBLE space before "001". pathlib will not normalise it.
# WT\ also contains G33V and N27D runs -- never glob this folder.


def _f(folder, name):
    return os.path.join(folder, name)


PANELS = {

    # ------------------------------------------------------------------ FIGURE 1A
    "1A": {
        "title": "SEC, cold incubation, t = 0 and 24 h",
        "normalized": True,
        "csvs": [
            _f(E22K, "WT immediate reinj 021325 500ul Loop no fractions 001.csv"),
            _f(E22K, "WT 0213 monomer reinj 021425 100ul Loop no fractions  001.csv"),
        ],
        "labels": ["t = 0", "t = 24 h"],
        "confidence": "CONFIRMED",
        "notes": [
            "Both CSVs live in Run_Results\\E22K\\ despite being WT runs -- that folder is "
            "the Feb-2025 NOESY-prep tree, not an E22K-only tree. Do not go looking in WT\\.",
            "Labels derived from the filenames' run dates: monomers isolated 02/13/25 "
            "('WT immediate reinj 021325') and reinjected 02/14/25 ('WT 0213 monomer reinj "
            "021425') = +1 day = 24 h. The filename literally carries both the isolation "
            "date (0213) and the reinjection date (021425).",
            "Vector match against Figure1R0.svg: red path81-7 -> the 021325 immediate run, "
            "rmse = 0.00000 normalised units (exact). Blue path82-8 -> the 021425 run, "
            "rmse = 0.00243 vs 0.01004 for the next candidate in the cold series (3 d), "
            "a 4.1x margin; every other cold time point (5 d, 6 d, 12 d, 22 d) is worse still.",
            "Independent structural confirmation: the Inkscape id lineage runs "
            "line2d_99-3/path154-0 -> g82-7/path81 -> g82-7-6/path81-7 (red) and "
            "line2d_100-2/path155-4 -> g83-4/path82 -> g83-4-5/path82-8 (blue), i.e. both "
            "published curves are clone-descendants of one pasted 2-trace matplotlib group, "
            "which is SECProgressE22K_onlyrequiredtraces.ipynb cell 4 -> cells 5/6 "
            "(filesWTprep = [f11, f12], legend 'Immediate' / '1d').",
            "TIME POINTS ARE RIGHT; TEMPERATURE AND CONCENTRATION ARE DISPUTED AND THE DATA "
            "CANNOT SETTLE THEM. Doc-10 caption: '50 uM incubated at 6 C'. Notebook: cell-2 "
            "markdown 'WT incubations on 4C for NOSEY 2D experiments'; cell-4 comments "
            "'Monomers had a concentration of 100uM' and '1d reinjection after cold 5C "
            "incubation'. The UNICORN exports carry NO temperature and NO concentration "
            "field -- the only non-UV channels are conductivity, %B, injection marks and a "
            "Logbook whose entire content is Method Settings / Equilibration / Sample "
            "Application / Elution. The 'UV_CUT_TEMP@100' curve name is a UNICORN internal "
            "label, not a sample temperature. ONLY THE BENCH RECORD FOR RUNS 021325 / 021425 "
            "CAN DECIDE. Weight of written evidence favours the notebook (100 uM, 4-5 C).",
            "The two runs used DIFFERENT injection loops: 500 uL (t = 0) vs 100 uL (t = 24 h). "
            "Raw monomer AUC over 12-18 mL is 12.604 vs 2.276 mAU*mL, a ratio of 5.54 ~ the "
            "5x loop ratio -- so the peptide concentration was essentially unchanged between "
            "the two injections, but the two traces are NOT comparable in raw mAU. This panel "
            "MUST be normalised; plotting it in raw mAU would be a 5x artefact.",
        ],
    },

    # ------------------------------------------------------------------ FIGURE 1B
    "1B": {
        "title": "SEC, 37 C, t = 0 and t = 2 h  (caption says 24 h -- caption is wrong)",
        "normalized": True,
        "csvs": [
            _f(WT, "WT 0 hr reinj 37C 021925 500ul Loop no fractions.csv"),
            _f(WT, "WT 2hr reinj 37C 021925 500ul Loop no fractions.csv"),
        ],
        "labels": ["t = 0", "t = 2 h"],
        "confidence": "CONFIRMED",
        "notes": [
            "THE PUBLISHED CAPTION IS WRONG. Doc 10 Figure 1B reads 't = 0 and t = 24 h'. "
            "The blue curve in the artwork is the 2-HOUR reinjection.",
            "Vector match against Figure1R0.svg with the page-x -> elution-volume map "
            "calibrated on the red (t = 0) curve and then HELD FIXED, so peak volumes are "
            "absolute rather than fitted. Red path154-1-4-9-2 -> the 0 hr CSV, rmse = 0.00007. "
            "Blue path155-2-5-7-1 scored against the whole 100 uM / 37 C series:",
            "      2h  rmse 0.00019   <-- best, and best over the entire candidate pool",
            "      4h  rmse 0.01796   (95x worse)",
            "      6h  rmse 0.02745",
            "      8h  rmse 0.03323",
            "      0h  rmse 0.03782",
            "     12h  rmse 0.04293",
            "     24h  rmse 0.07133   (375x worse -- the captioned assignment)",
            "     48h  rmse 0.09509",
            "     72h  rmse 0.10058",
            "Margin best/next-in-series = 95x; margin best/captioned-24h = 375x. The only "
            "candidates anywhere in Run_Results that come close (0.00027, 0.00038) are two "
            "E22K monomer runs from a different peptide and a different month, and both need "
            "a volume shift railed at the +0.10 mL tolerance bound, the signature of a "
            "shape-degenerate near-match between two featureless monomer peaks rather than "
            "an identity. The true source needs zero shift.",
            "Independent structural confirmation: path154-1-4-9-2 and path155-2-5-7-1 are "
            "Inkscape clone-descendants (suffix chain -1-4-9-2 / -2-5-7-1) of path154-1 and "
            "path155-2, which are traces 1 and 2 of the Fig 3B seven-trace group inside "
            "37CcompE22KvsWT.svg. Traces 1 and 2 of that group are w11 (0 h) and w12 (2 h). "
            "The panel was literally built by copying the first two curves of Fig 3B.",
            "Code corroboration: SECProgressE22K_onlyrequiredtraces.ipynb cell 32 builds "
            "filesWT12 = [w11 (0 hr), w12 (2hr)] and cell 33 saves it as "
            "'WT_37C_75inconlyoand2.svg' -- 'only 0 and 2'. (Cell 32 has a bug: it passes "
            "SECData(filesWT2) instead of filesWT12, so re-running it as it stands yields 7 "
            "traces, not 2. Fix the argument if you re-run it; the intent is unambiguous.)",
            "The 24 h CSV exists and is used in Figs 3B and 4B, so this is an authoring "
            "choice, not missing data. DECISION NEEDED: fix the caption to t = 2 h, or swap "
            "in 'WT 24hr reinj 37C 022025 500ul Loop no fractions 001.csv' and keep the "
            "caption. The new panel must not silently perpetuate the current pairing.",
        ],
    },

    # ------------------------------------------------------------------ FIGURE 3A
    "3A": {
        "title": "SEC, 250 uM, 37 C, raw mAU",
        "normalized": False,
        "csvs": [
            _f(WT250, "WT 2 hr 37C inc reinj 001.csv"),
            _f(WT250, "WT 4hr 37C inc reinj 001.csv"),
            _f(WT250, "500ul Loop no fractions 6h WT 37C inc high conc 001.csv"),
        ],
        "labels": ["t = 2 h", "t = 4 h", "t = 6 h"],
        "confidence": "CONFIRMED",
        "notes": [
            "Labels read straight off the filenames: '2 hr', '4hr', '6h'. All three also say "
            "'37C inc' / '37C inc high conc', independently confirming the temperature.",
            "Vector match against WT_4c_250uM.svg (pure matplotlib, axes rect "
            "x=84.915 y=47.583 w=749.16 h=584.14, xlim (6,20), ylim (-0.41151, 20.98812)), "
            "no free parameters: red line2d_77 -> 2 hr, rmse 0.00002 mAU (margin 65917x); "
            "blue line2d_78 -> 4hr, rmse 0.00003 (72323x); black line2d_79 -> 6h high conc, "
            "rmse 0.00002 (74904x).",
            "Panel has NO t = 0. 'wt immediate reinj frac34 test 060325 500ul Loop no "
            "fractions  001 001.csv' sits in the same folder but was deliberately not plotted "
            "(WT37COligomer.ipynb cell 2, filesWT1 = [WT1, WT2, WT3]). Adding it is a "
            "scientific choice, not a correction.",
            "The parent SVG is misnamed 'WT_4c_250uM.svg' -- it is 37 C data, not 4 C. The "
            "AUC cell (cell 5) titles its plot '260uM' while everything else says 250 uM; "
            "the caption's 250 uM is the value to use.",
        ],
    },

    # ------------------------------------------------------------------ FIGURE 3B
    "3B": {
        "title": "SEC, 100 uM, 37 C, raw mAU",
        "normalized": False,
        "csvs": [
            _f(WT, "WT 0 hr reinj 37C 021925 500ul Loop no fractions.csv"),
            _f(WT, "WT 2hr reinj 37C 021925 500ul Loop no fractions.csv"),
            _f(WT, "WT 4hr reinj 37C 021925 500ul Loop no fractions.csv"),
            _f(WT, "WT 8hr reinj 37C 021925 500ul Loop no fractions 001.csv"),
            _f(WT, "WT 12hr reinj 37C 021925 500ul Loop no fractions  001 001.csv"),
            _f(WT, "WT 24hr reinj 37C 022025 500ul Loop no fractions 001.csv"),
            _f(WT, "WT 48hr reinj 37C 022125 500ul Loop no fractions  001 001.csv"),
        ],
        "labels": ["t = 0", "t = 2 h", "t = 4 h", "t = 8 h",
                   "t = 12 h", "t = 24 h", "t = 48 h"],
        "confidence": "CONFIRMED",
        "notes": [
            "HER PUBLISHED LEGEND IS WRONG FOR 5 OF THE 7 TRACES, and the bug is in the CODE, "
            "inherited into the artwork. SECProgressE22K_onlyrequiredtraces.ipynb cell 34 "
            "zips legend_labels ['Immediate','2h','4h','8h','12h','24h','48h'] positionally "
            "onto the cell-28/33 file list, whose order is [0h, 2h, 12h, 24h, 48h, 4h, 8h].",
            "Verified numerically on the surviving vector (the 7-trace flattened group "
            "path154-1 / 155-2 / 156-6 / 157-08 / 158-1 / 159-1 / 160-0 inside "
            "37CcompE22KvsWT.svg, whose stroke colours red/blue/black/grey/orange/green/"
            "yellow are exactly cell 34's `colors` list). Colour -> true time point:",
            "     red    #ff0000 ->  0h   rmse 0.00006  (margin 1040x)   legend said 'Immediate'  OK",
            "     blue   #0000ff ->  2h   rmse 0.00033  (51x)            legend said '2h'         OK",
            "     black  #000000 -> 12h   rmse 0.00048  (16x)            legend said '4h'         WRONG",
            "     grey   #808080 -> 24h   rmse 0.00058  (36x)            legend said '8h'         WRONG",
            "     orange #ffa500 -> 48h   rmse 0.00051  (41x)            legend said '12h'        WRONG",
            "     green  #008000 ->  4h   rmse 0.00039  (24x)            legend said '24h'        WRONG",
            "     yellow #ffff00 ->  8h   rmse 0.00042  (18x)            legend said '48h'        WRONG",
            "The same seven CSVs were also matched at rmse = 0.00000 (margins 5e6 - 2e7 x) "
            "through the normalised sibling plot WT_37c_100uM.svg, which is pure matplotlib "
            "and therefore has an exact axes transform. The CSV list is not in doubt.",
            "No standalone SVG for this panel exists. Cell 34 declares "
            "output_file='WT_37C_75inc.svg'; that file is nowhere in the SEC or For trees. "
            "Only Run_Results\\E22K\\WT_37C_75inc.png (an earlier 5-trace variant) survives. "
            "The panel is fully reproducible from raw, but her exact artwork is not "
            "re-exportable.",
            "6 h and 72 h CSVs for this series EXIST and were not plotted: "
            "'WT 6hr reinj 37C 021925 500ul Loop All fractions.csv' (note: 'All fractions', "
            "unlike its siblings) and 'WT 72hr reinj 37 C 022225 500ul Loop no fractions  "
            "001 001.csv' (note the space in '37 C'). Adding them would make 3B a 9-point "
            "series matching 3C's density. Scientific choice, flagged not taken.",
            "Fig 3B and Fig 4B are the SAME seven CSVs -- raw mAU here, min-max normalised "
            "there. Worth deciding whether the paper needs both.",
        ],
    },

    # ------------------------------------------------------------------ FIGURE 3C
    "3C": {
        "title": "SEC, 25 uM, 37 C, raw mAU",
        "normalized": False,
        "csvs": [
            _f(WT, "wt immed reinj 500ul Loop no fractions  001 001 25uM.csv"),
            _f(WT, "wt 3h reinj 500ul Loop no fractions  001 001 25uM.csv"),
            _f(WT, "wt 6h reinj 500ul Loop no fractions  001 001 25uM.csv"),
            _f(WT, "wt 12h reinj 500ul Loop no fractions  001 001 25uM.csv"),
            _f(WT, "wt 24h reinj 500ul Loop no fractions  001 25uM.csv"),
            _f(WT, "wt 48h reinj 500ul Loop no fractions  001 001 25uM.csv"),
            _f(WT, "wt 72h reinj 500ul Loop no fractions  001 001 25uM.csv"),
        ],
        "labels": ["t = 0", "t = 3 h", "t = 6 h", "t = 12 h",
                   "t = 24 h", "t = 48 h", "t = 72 h"],
        "confidence": "CONFIRMED",
        "notes": [
            "Labels read straight off the filenames ('immed', '3h', '6h', '12h', '24h', "
            "'48h', '72h'); every filename also carries its own '25uM' concentration token, "
            "so this is the one panel whose concentration is self-documenting.",
            "TRAP: the 24 h file is 'wt 24h reinj 500ul Loop no fractions  001 25uM.csv' -- "
            "ONE '001' token where all six siblings have '001 001'. Copy it verbatim.",
            "All seven names carry a DOUBLE space before the first '001'.",
            "Vector match against WT_37C_25uM.svg (pure matplotlib, axes rect x=71.56875 "
            "y=47.583125 w=762.48 h=584.14, xlim (6,20), ylim (0,3)), no free parameters: "
            "every one of the seven curves matched its CSV at rmse = 0.00000 mAU with "
            "next-best margins of 6.2e6 - 1.8e7 x. Plot order red/blue/black/orange/purple/"
            "grey/green = immed/3h/6h/12h/24h/48h/72h, i.e. her file order was already "
            "chronological here and her legend (cell 12) is CORRECT for this panel.",
            "DO NOT confuse this series with Run_Results\\WT\\WT 25 uM 37 Inc\\, a SEPARATE "
            "Oct-2025 25 uM run (0/2/4/6/8/12 h, six files) loaded by WT37COligomer.ipynb "
            "cell 20. It was never published in doc 10.",
        ],
    },

    # ------------------------------------------------------------------ FIGURE 4A
    #
    # ==================================================================================
    # DECISION RECORD -- WHY PANEL 4A (combined-figure panel d) HAS NO t = 0
    # Decided by the user, 2026-08-07, on the evidence below.  Option taken: leave the
    # series starting at t = 2 h and say so in the caption.  DO NOT "fix" this by
    # borrowing a zero point from another prep, and do not let a later editor assume
    # the 0 h file was simply forgotten.  It does not exist.
    #
    # 1. NO ZERO-TIME REINJECTION WAS COLLECTED FOR THE DEC-2025 25 C CAMPAIGN.
    #    Run_Results\RT runs\ holds exactly ten CSVs: the nine reinjections listed
    #    below (2, 4, 6, 8, 10, 12, 24, 48, 72 h) plus 'WT prim run 122025 001.csv'.
    #
    # 2. THE PRIMARY RUN IS NOT A ZERO POINT.  It is the crude-peptide isolation run,
    #    not a reinjection of isolated monomer.  Measured on the reprocessed raw data
    #    (areas by np.trapz over the 6-20 mL plotted range; "monomer" = 13.0-16.0 mL,
    #    "post-monomer" = 16.0-20.0 mL):
    #        WT prim run 122025 001.csv     monomer apex 15.68 mL, 39.3 % of area
    #                                       monomer, 30.9 % post-monomer salt/small
    #                                       molecule
    #        WT 0 hr reinj 37C 021925 ...   monomer apex 14.57 mL, 89.0 % of area
    #                                       monomer      <- what a genuine 0 h looks like
    #    A 1.1 mL later apex and less than half the monomer purity: this is crude
    #    material, not the t = 0 of the incubation series.
    #
    # 3. BORROWING A t = 0 FROM ANOTHER PREP WAS REJECTED.  The WT t = 0
    #    oligomer/monomer PEAK-HEIGHT ratio (apex in 7.0-11.5 mL over apex in
    #    12.0-18.0 mL, on min-subtracted raw traces) spans ~17x across preps:
    #        0.0099  E22K\WT immediate reinj 021325 500ul Loop no fractions 001.csv
    #        0.0185  WT\WT 0 hr reinj 37C 021925 500ul Loop no fractions.csv   (3B/4B t=0)
    #        0.0363  WT\wt immed reinj 500ul Loop no fractions  001 001 25uM.csv (3C t=0)
    #        0.0400  WT\06-04-25 Results\wt immediate reinj frac34 test 060325 ... .csv
    #        0.1677  WT\WT 25 uM 37 Inc\wt 25uM 0hr reinj 100125 500ul Loop no
    #                fractions  001 001.csv        (Oct-2025 series, never published)
    #    spread 0.1677 / 0.0099 = 16.9x.  The largest of them, 0.1677, is
    #    indistinguishable from THIS series' own 2 h point, 0.1601.  An imported t = 0
    #    could therefore sit on top of the 2 h trace and manufacture a "no growth in
    #    the first two hours" result that the data do not support.
    #
    # 4. CONSEQUENCE FOR THE CAPTION: state as fact that no zero-time reinjection was
    #    collected for the 25 C series and that d begins at 2 h.  Not an apology, not
    #    a limitation paragraph -- a statement of what was measured.
    # ==================================================================================
    "4A": {
        "title": "SEC, 100 uM, 25 C, normalized",
        "normalized": True,
        "csvs": [
            _f(RT, "WT 100uM RT 2hr 122025  001.csv"),
            _f(RT, "WT 100uM RT 4hr 122025  001.csv"),
            _f(RT, "WT 100uM RT 6hr 122025  001.csv"),
            _f(RT, "WT 100uM RT 8hr 122025  001.csv"),
            _f(RT, "WT 100uM RT 10hr 122025  001.csv"),
            _f(RT, "WT 100uM RT 12hr 122025  001.csv"),
            _f(RT, "24h RT WT reinjection 001.csv"),
            _f(RT, "WT 100uM RT 48hr 122225  001.csv"),
            _f(RT, "WT 100uM RT 72hr 122325  001.csv"),
        ],
        "labels": ["t = 2 h", "t = 4 h", "t = 6 h", "t = 8 h", "t = 10 h",
                   "t = 12 h", "t = 24 h", "t = 48 h", "t = 72 h"],
        "confidence": "CONFIRMED",
        "notes": [
            "Labels read off the filenames ('2hr' ... '72hr'; the 24 h file is named "
            "'24h RT WT reinjection 001.csv' and is the odd one out -- no '100uM', no date "
            "token). Eight of nine also carry '100uM' and 'RT', self-documenting the "
            "concentration and the nominal temperature.",
            "Vector match against WT_25c_100uM_normalized.svg (pure matplotlib, axes rect "
            "x=91.5825 y=47.241875 w=742.44 h=584.62, xlim (5,22), ylim (0,1.1)), no free "
            "parameters: all nine curves matched at rmse = 0.00000 normalised units, margins "
            "3.9e6 - 1.3e7 x. Colour -> time point in the delivered SVG: red 2h, blue 4h, "
            "green 6h, orange 8h, yellow 10h, black 12h, PINK 48h, GREY 72h, PURPLE 24h.",
            "NOTE THE ORDERING TRAP: the SVG's own plot order is 2,4,6,8,10,12,48,72,24 -- "
            "24 h is LAST in the drawn data, matching cell 31's legend_labels "
            "['2h','4h','6h','8h','10h','12h','48h','72h','24h']. But cell 24 as it stands "
            "today builds filesWT4 = [2,4,6,8,10,12,24,48,72] with 24 h in SEVENTH position. "
            "The file list was edited after the figure was made. RE-RUNNING CELL 31 AGAINST "
            "TODAY'S CELL 24 WOULD PRODUCE A 3-WAY LABEL SWAP (24h<->48h<->72h). Her "
            "published legend happens to be correct because the retype was done against the "
            "SVG, not against the notebook. Derive from filenames only.",
            "Panel has NO t = 0. 'WT prim run 122025 001.csv' is a primary crude run, not a "
            "reinjection time point, and is not in filesWT4. TRAP: on disk that name has a "
            "SINGLE space before '001'; cell 24 types a DOUBLE space. It is unused, so the "
            "typo is harmless -- but never copy filenames from the notebook.",
            "DECIDED 2026-08-07 (user's choice): leave the series starting at t = 2 h and "
            "say so in the caption. See the DECISION RECORD comment block immediately above "
            "this dict for the evidence -- primary run 39.3 % monomer at 15.68 mL vs 89.0 % "
            "at 14.57 mL for a genuine 0 h, and a 16.9x spread in the WT t = 0 "
            "oligomer/monomer ratio across preps (0.0099 / 0.0185 / 0.0363 / 0.0400 / "
            "0.1677) whose maximum equals this series' own 2 h value of 0.1601. "
            "Recorded in PROVENANCE_FIGURES_1-4.md sec. 9 and in "
            "reports\\FIGURE_SEC_CAPTION.md.",
            "'RT' in the filenames means room temperature; the caption calls it 25 C. The "
            "CSVs record no temperature. If Methods must state 25 C, that comes from the "
            "bench record, not from these files.",
            "The published artwork recoloured the 10 h trace from yellow to teal (yellow was "
            "invisible on white). That recolour is cosmetic and carries no label change.",
        ],
    },

    # ------------------------------------------------------------------ FIGURE 4B
    "4B": {
        "title": "SEC, 100 uM, 37 C, normalized  (same CSVs as 3B)",
        "normalized": True,
        "csvs": [
            _f(WT, "WT 0 hr reinj 37C 021925 500ul Loop no fractions.csv"),
            _f(WT, "WT 2hr reinj 37C 021925 500ul Loop no fractions.csv"),
            _f(WT, "WT 4hr reinj 37C 021925 500ul Loop no fractions.csv"),
            _f(WT, "WT 8hr reinj 37C 021925 500ul Loop no fractions 001.csv"),
            _f(WT, "WT 12hr reinj 37C 021925 500ul Loop no fractions  001 001.csv"),
            _f(WT, "WT 24hr reinj 37C 022025 500ul Loop no fractions 001.csv"),
            _f(WT, "WT 48hr reinj 37C 022125 500ul Loop no fractions  001 001.csv"),
        ],
        "labels": ["t = 0", "t = 2 h", "t = 4 h", "t = 8 h",
                   "t = 12 h", "t = 24 h", "t = 48 h"],
        "confidence": "CONFIRMED",
        "notes": [
            "Identical CSV set to Fig 3B; the only difference is min-max normalisation.",
            "Vector match against WT_37c_100uM.svg (pure matplotlib -- the file is misnamed, "
            "it IS the normalised plot; it was later hand-renamed to "
            "WT_37c_100uM_normalized.svg in the Inkscape folder). Axes rect x=91.5825 "
            "y=47.241875 w=742.44 h=584.62, xlim (5,22), ylim (0,1.1), no free parameters. "
            "All seven curves matched at rmse = 0.00000 normalised units, margins "
            "5.1e6 - 2.4e7 x. Drawn order red 0h, blue 2h, green 12h, orange 24h, yellow 48h, "
            "black 4h, red 8h.",
            "COLOUR-CYCLING BUG in the parent: WT37COligomer.ipynb cell 34 passes only SIX "
            "colours for SEVEN traces and sec1.py:408 does colors = colors * (n//len+1), so "
            "the 7th trace (8 h) is drawn in RED, identical to the 1st (0 h). Confirmed in "
            "the SVG: line2d_75 and line2d_81 are both #ff0000. Any faithful re-plot must "
            "supply seven distinct colours.",
            "HER PUBLISHED LEGEND IS WRONG FOR 3 OF 6 TRACES (green, black, red). Unlike 3B "
            "this one was introduced purely in Inkscape -- the matplotlib parent's own legend "
            "(cell 34 legend_labels ['0h','2h','12h','24h','48h','4h','8h']) is CORRECT for "
            "the file order. The orphaned original legend still sits off-artboard in "
            "WT25Cvs37CNormalized.svg reading '0h 2h 12h 24h 48h 4h 8h'. Correct published "
            "mapping: blue 2h, black 4h, red 8h, green 12h, orange 24h, lavender 48h.",
            "THE PUBLISHED PANEL SHOWS ONLY SIX TRACES: the t = 0 curve was DELETED by hand "
            "in WT25Cvs37CNormalized.svg. The matplotlib parent has seven. This manifest "
            "lists all seven; drop index 0 if the intent is to reproduce the published panel "
            "exactly, but keeping t = 0 is the better figure and makes 4A/4B consistent "
            "(4A also starts without t = 0 -- neither panel is currently self-consistent).",
        ],
    },
}


# --------------------------------------------------------------------------- checks
def _iter_rows():
    for pid in ("1A", "1B", "3A", "3B", "3C", "4A", "4B"):
        p = PANELS[pid]
        for path, label in zip(p["csvs"], p["labels"]):
            yield pid, p, path, label


def check():
    """Assert every listed file exists and every panel is internally consistent."""
    missing = []
    for pid, p in PANELS.items():
        assert len(p["csvs"]) == len(p["labels"]), \
            "%s: %d csvs but %d labels" % (pid, len(p["csvs"]), len(p["labels"]))
        assert isinstance(p["normalized"], bool), "%s: normalized must be bool" % pid
        for f in p["csvs"]:
            if not os.path.isfile(f):
                missing.append((pid, f))
    if missing:
        for pid, f in missing:
            print("MISSING  %-3s  %s" % (pid, f))
        raise AssertionError("%d listed CSV(s) not found on disk" % len(missing))
    return True


def main():
    check()
    order = ("1A", "1B", "3A", "3B", "3C", "4A", "4B")
    print("=" * 100)
    print("SEC PANEL RESOLUTION -- 37 C Abeta42 oligomer manuscript (doc 10)")
    print("=" * 100)
    print("%-4s %-6s %-7s %-11s %s" % ("PANEL", "TRACES", "Y-AXIS", "EVIDENCE", "DESCRIPTION"))
    print("-" * 100)
    for pid in order:
        p = PANELS[pid]
        print("%-4s %-6d %-7s %-11s %s" % (
            pid, len(p["csvs"]),
            "norm" if p["normalized"] else "raw mAU",
            p["confidence"], p["title"]))
    print("-" * 100)
    total = 0
    for pid in order:
        p = PANELS[pid]
        print()
        print("%s  --  %s  [%s, %s]" % (
            pid, p["title"], "normalized" if p["normalized"] else "raw mAU", p["confidence"]))
        for i, (f, lab) in enumerate(zip(p["csvs"], p["labels"])):
            total += 1
            print("   %d. %-10s %s" % (i + 1, lab, f))
    print()
    print("-" * 100)
    print("%d panels, %d trace/CSV pairs, all files present." % (len(order), total))
    print("Distinct CSVs referenced: %d" % len({f for _, _, f, _ in _iter_rows()}))
    print("-" * 100)
    print("FLAGS THAT NEED A HUMAN DECISION:")
    print("  1B  caption says t = 24 h; the plotted curve is the 2 h reinjection")
    print("      (rmse 0.00019 vs 0.07133, a 375x margin). Fix caption or swap CSV.")
    print("  1A  time points confirmed (0 / 24 h); TEMPERATURE and CONCENTRATION unresolved.")
    print("      Caption 50 uM / 6 C vs notebook 100 uM / 4-5 C. The CSVs carry neither field.")
    print("  3B  her legend mislabels 5 of 7 traces; corrected above.")
    print("  4B  her legend mislabels 3 of 6 traces, and the published panel silently drops")
    print("      t = 0. Corrected and restored above.")
    print("  4A  her published legend is correct, but re-running notebook cell 31 against")
    print("      today's cell 24 would introduce a 24h/48h/72h swap. Labels here come from")
    print("      filenames, so this manifest is immune.")
    print("  3B/4B are the same seven CSVs plotted twice (raw vs normalized).")
    print("-" * 100)


if __name__ == "__main__":
    main()
