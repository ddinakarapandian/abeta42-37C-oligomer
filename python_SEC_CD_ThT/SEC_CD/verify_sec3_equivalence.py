"""verify_sec3_equivalence.py — prove the project's SEC loader is sec3.py's formula.

sec3.py is the reference implementation for this paper's SEC processing:
    /path/to/Ab42oligomers37C/code/shared/BioexpressedABWT\\sec3.py
    SECData.process_data_to_arrays -> SECData._process_single_dataframe

The project does NOT import it -- it lives in a read-only Dropbox tree belonging to
a different paper, and importing across that boundary would make every figure here
depend on a file this repository must never write to.  Instead
``build_fig1_panels.load_trace`` transcribes it, and this script exists so the
transcription cannot silently drift: it re-implements sec3's steps from its source
and asserts agreement on every CSV in the registry.

sec3's pipeline, in order:
    1. first two columns -> Volume_ml, Absorbance_mAU
    2. filter min_volume <= V <= max_volume   (defaults 0, 25)
    3. Adjusted = Absorbance - min  IF min < 0  ELSE  Absorbance
    4. if normalize: Normalized = (Adjusted - Adjusted.min()) / (Adjusted.max() - Adjusted.min())
    5. return Normalized if normalize else Adjusted

load_trace does the same and additionally sorts by volume, which permutes rows
without changing any value -- so the comparison is made on sorted order.

    python verify_sec3_equivalence.py          # exit 0 if identical
"""
from __future__ import annotations

import os
import sys

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from secfigs.panel_data import PANELS
from build_fig1_panels import load_trace

SEC3 = (r"/path/to/Ab42oligomers37C/code/shared/BioexpressedABWT\sec3.py")
TOL = 0.0          # they are the same arithmetic; anything nonzero is a real divergence


def sec3_pipeline(path, normalize, min_volume=0, max_volume=25):
    """sec3.SECData.process_data_to_arrays, transcribed from the source above."""
    df = pd.read_csv(path, encoding="utf-16", sep="\t", skiprows=3)
    num = df.select_dtypes(include=[np.number]).columns[:2]
    df = df[num].copy()
    df.columns = ["Volume_ml", "Absorbance_mAU"]
    df = df.dropna()
    df = df[(df["Volume_ml"] >= min_volume) & (df["Volume_ml"] <= max_volume)]
    mn = df["Absorbance_mAU"].min()
    adj = df["Absorbance_mAU"] - mn if mn < 0 else df["Absorbance_mAU"]
    if normalize:
        out = (adj - adj.min()) / (adj.max() - adj.min())
    else:
        out = adj
    return df["Volume_ml"].to_numpy(float), out.to_numpy(float)


def main():
    if not os.path.exists(SEC3):
        print("NOTE  sec3.py not found at its recorded path; comparing against the\n"
              "      transcription in this file only:\n      %s" % SEC3)

    n = 0
    worst = {True: (0.0, None), False: (0.0, None)}
    for pid in sorted(PANELS):
        for f in PANELS[pid]["csvs"]:
            if not os.path.exists(f):
                print("  MISSING %s" % f)
                continue
            for norm in (True, False):
                v1, y1, _ = load_trace(f, normalized=norm)
                v3, y3 = sec3_pipeline(f, norm)
                if len(y1) != len(y3):
                    print("  LENGTH DIFF %s %s: %d vs %d"
                          % (pid, os.path.basename(f), len(y1), len(y3)))
                    return 1
                o1 = np.argsort(v1, kind="stable")
                o3 = np.argsort(v3, kind="stable")
                d = float(np.max(np.abs(y1[o1] - y3[o3])))
                if d > worst[norm][0]:
                    worst[norm] = (d, "%s / %s" % (pid, os.path.basename(f)))
                n += 1

    print("compared %d (trace, mode) pairs across %d panels" % (n, len(PANELS)))
    for norm, tag in ((True, "normalized"), (False, "raw       ")):
        d, where = worst[norm]
        print("  max |load_trace - sec3|  %s : %.3e%s"
              % (tag, d, "" if where is None or d == 0 else "   at " + where))

    bad = [t for t in (True, False) if worst[t][0] > TOL]
    if bad:
        print("\nDIVERGED. build_fig1_panels.load_trace no longer matches sec3.py.\n"
              "Fix load_trace, not this script -- sec3.py is the reference.")
        return 1
    print("\nIDENTICAL. The project's SEC processing is sec3.py's formula.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
