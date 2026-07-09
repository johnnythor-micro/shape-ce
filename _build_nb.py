import nbformat as nbf
nb = nbf.v4.new_notebook()
cells = []
def md(t): cells.append(nbf.v4.new_markdown_cell(t))
def code(t): cells.append(nbf.v4.new_code_cell(t))

md("""# SHAPE-CE Analysis Notebook (`shapece`)
### A modern, notebook-based pipeline for capillary-electrophoresis nucleic-acid probing

This notebook takes you from raw capillary-electrophoresis trace files (`.fsa`) to
**per-nucleotide reactivities** and a **SHAPE-directed secondary structure**, following the
analysis logic of **QuShape** (Karabiber et al., *RNA* 2013) reimplemented in modern Python 3,
with added internal size-standard alignment and ViennaRNA structure prediction.

It works for **any** probing experiment — SHAPE, DMS, CMCT, hydroxyl-radical, or enzymatic —
and makes **no assumptions** about your RNA, organism, or hypothesis.

**How to use this notebook (beginner-friendly):**
- Run cells top to bottom (Shift+Enter).
- You only edit **one cell** — the `CONFIG` cell in Step 1. Everything else runs as-is.
- After each step there is a **plot to check** — read the short "What to look for" note.
- No coding experience required beyond editing file paths and numbers in the CONFIG cell.

*Not affiliated with or endorsed by the original QuShape authors. MIT-licensed.*
""")

md("## Step 0 — Install and import\nRun this once. On Google Colab it installs everything; locally, it uses what you have.")
code("""# If running on Google Colab, this installs the package + dependencies.
try:
    import google.colab  # noqa
    !pip -q install ViennaRNA biopython openpyxl pandas scipy matplotlib
    !pip -q install git+https://github.com/YOURNAME/shape-ce.git   # <-- your repo URL after you publish
except ImportError:
    pass

import numpy as np, matplotlib.pyplot as plt
import shapece as sc
print("shapece version:", sc.__version__)""")

md("""## Step 1 — Configure your experiment  ✏️ **(the only cell you edit)**

Fill in your file paths, the reagent (+) / background (−) pairs, the DATA channel numbers,
your RNA sequence, and (optionally) your sequencing ladders.

**Finding your DATA channel numbers:** run Step 2 first with any guess — it prints and plots
all channels in each file so you can identify which holds your reagent trace and which holds
the size standard. Then come back and set them here.

**Separate-capillary design** (most common; e.g. + and − in different files): give each sample
its own file and the same `data_channel`. **One-capillary multi-dye design:** put +, −, and
ladders in one file with different `data_channel` numbers.""")
code('''CONFIG = {
    # --- your RNA, 5'->3' (U or T both fine) ---
    "rna_sequence": "GGGAUCGUCAC...",   # <-- EDIT

    # --- which DATA channels hold what (see Step 2 to identify) ---
    "data_channel": 9,     # the analyzed reagent/probe trace   <-- EDIT
    "size_channel": 205,   # the internal size standard (e.g. LIZ); set None if you have none

    # --- samples: label -> file path ---
    "samples": {           # <-- EDIT (add as many as you have)
        "plus_rep1":  "/path/to/plus_rep1.fsa",
        "minus_rep1": "/path/to/minus_rep1.fsa",
        "plus_rep2":  "/path/to/plus_rep2.fsa",
        "minus_rep2": "/path/to/minus_rep2.fsa",
    },
    # --- (+) reagent / (-) background pairs, by label ---
    "pairs": [("plus_rep1", "minus_rep1"),
              ("plus_rep2", "minus_rep2")],   # <-- EDIT

    # --- reference lane that defines the common axis (any good (+) lane) ---
    "reference": "plus_rep1",   # <-- EDIT

    # --- reactivity model: "area_difference" (QuShape-style, default) or "stop_fraction" ---
    "reactivity_model": "area_difference",

    # --- OPTIONAL sequencing ladders for base registration (set to None if none) ---
    "ladders": {
        # "S1_file": "/path/to/ddATP_ladder.fsa", "S1_base": "U",
        # "S2_file": "/path/to/ddGTP_ladder.fsa", "S2_base": "C",
    },

    # --- analysis window in sample units (set after viewing Step 2; None = auto) ---
    "read_window": (None, None),
}
print("Configured %d samples, %d pairs." % (len(CONFIG["samples"]), len(CONFIG["pairs"])))''')

md("""## Step 2 — Load files and inspect channels
This prints the DATA channels in each file and plots them, so you can confirm your
`data_channel` and `size_channel` choices in CONFIG.

**What to look for:** the reagent channel shows many sharp peaks (your RT stops); the size
standard shows an evenly spaced ladder of peaks. If they look swapped, fix the numbers in CONFIG.""")
code('''from shapece.io import read_abif, get_channel, list_data_channels

entries = {}
for label, path in CONFIG["samples"].items():
    entries[label] = read_abif(path)
    print(f'{label:14s} DATA channels: {list_data_channels(entries[label])}')

# quick look at the reference lane
e = entries[CONFIG["reference"]]
fig, ax = plt.subplots(2, 1, figsize=(13, 5))
ax[0].plot(get_channel(e, CONFIG["data_channel"]), lw=0.5)
ax[0].set_title(f'reagent trace (DATA{CONFIG["data_channel"]}) — {CONFIG["reference"]}')
if CONFIG["size_channel"]:
    ax[1].plot(get_channel(e, CONFIG["size_channel"]), lw=0.5, color="darkgreen")
    ax[1].set_title(f'size standard (DATA{CONFIG["size_channel"]})')
plt.tight_layout(); plt.show()''')

md("""## Step 3 — Preprocess (saturation → smooth → baseline)
Cleans each reagent trace. These corrections are standard and safe; defaults work for most data.

**What to look for:** peaks preserved, flat baseline at zero, no giant clipped plateaus.""")
code('''from shapece import preprocess as pp

traces = {}   # label -> {"RX": cleaned reagent trace}
size = {}     # label -> size-standard trace
for label, e in entries.items():
    d = get_channel(e, CONFIG["data_channel"])
    d = pp.baseline(pp.smooth(pp.correct_saturation(d)))
    traces[label] = {"RX": d}
    if CONFIG["size_channel"]:
        size[label] = get_channel(e, CONFIG["size_channel"])

ref = CONFIG["reference"]
plt.figure(figsize=(13, 3))
plt.plot(traces[ref]["RX"], lw=0.5); plt.title(f'cleaned reagent trace — {ref}')
plt.tight_layout(); plt.show()''')

md("""## Step 4 — Align lanes to a common axis
If you have an internal size standard, we warp every lane onto the reference lane using the
matched standard peaks (**recommended** — this is what makes lanes truly comparable and fixes
instrument calibration drift). If you have no size standard, set `size_channel=None` and this
falls back to DTW alignment on the reagent traces.

**What to look for:** after alignment, same-condition lanes should overlap almost perfectly.""")
code('''from shapece import align

if CONFIG["size_channel"]:
    aligned = align.size_standard_align(traces, size, reference=CONFIG["reference"], roles=("RX",))
else:
    aligned = align.dtw_align(traces, reference=CONFIG["reference"], align_on="RX", roles=("RX",))

# overlay the (+) lanes to check co-registration
plus_labels = [p for p, m in CONFIG["pairs"]]
plt.figure(figsize=(13, 3))
for lab in plus_labels:
    plt.plot(aligned[lab]["RX"], lw=0.4, alpha=0.7, label=lab)
lo, hi = CONFIG["read_window"]
plt.xlim(lo, hi); plt.legend(fontsize=8); plt.title("aligned (+) lanes (should overlap)")
plt.tight_layout(); plt.show()''')

md("""## Step 5 — Detect and quantify peaks
Peaks are detected once on the mean of the (+) lanes, then quantified in every lane at those
same positions (Gaussian deconvolution by default). Adjust `read_window` in CONFIG to the region
with clean signal if needed.

**What to look for:** red marks sit on peak apices across the read; a few hundred peaks is typical.""")
code('''from shapece import peaks

ref_trace = np.mean([aligned[lab]["RX"] for lab in plus_labels], axis=0)
lo, hi = CONFIG["read_window"]
pos = peaks.detect_peaks(ref_trace, lo=lo or 0, hi=hi, min_spacing=12, prominence_frac=0.03)
print("peaks detected:", len(pos))

quant = {lab: peaks.quantify(aligned[lab]["RX"], pos, mode="gaussian") for lab in entries}

plt.figure(figsize=(13, 3))
plt.plot(ref_trace, lw=0.5, color="0.4")
plt.plot(pos, ref_trace[pos], "rv", ms=3)
plt.xlim(lo, hi); plt.title(f"{len(pos)} consensus peaks"); plt.tight_layout(); plt.show()''')

md("""## Step 6 — Compute reactivity
For each (+)/(−) pair we compute reactivity, then apply model-free boxplot normalization
(reactive nucleotides ≈ 1.0). Default model is **area-difference** (QuShape-style); switch to
**stop-fraction** in CONFIG if your readout is stop-based.

**What to look for:** most values 0–2; replicates of the same condition tracking each other.""")
code('''from shapece import reactivity as rx

profiles = {}   # pair label -> normalized reactivity
for plus, minus in CONFIG["pairs"]:
    if CONFIG["reactivity_model"] == "area_difference":
        r = rx.area_difference(quant[plus]["area"], quant[minus]["area"], scale=True)
    else:
        r = (rx.stop_fraction(aligned[plus]["RX"], pos)
             - rx.stop_fraction(aligned[minus]["RX"], pos))
        r = np.clip(r, 0, None)
    norm, factor = rx.boxplot_normalize(r)
    profiles[f"{plus}"] = norm

plt.figure(figsize=(13, 3))
for lab, r in profiles.items():
    plt.plot(r, lw=0.7, alpha=0.8, label=lab)
plt.axhline(1, color="k", lw=0.4, ls="--"); plt.legend(fontsize=8)
plt.ylabel("normalized reactivity"); plt.title("reactivity profiles"); plt.tight_layout(); plt.show()''')

md("""## Step 7 — Register to nucleotide numbers
Assign each peak an absolute position in your RNA. **With ladders** (set in CONFIG) we align the
ladder-called bases to your sequence (accurate). **Without ladders** we skip this and number peaks
by order from the 3′ end (clearly approximate — a ladder is strongly recommended for real numbering).""")
code('''from shapece import sequence as seqmod

seq = CONFIG["rna_sequence"].replace("T", "U")
if CONFIG["ladders"]:
    lad = CONFIG["ladders"]
    ladder_area = {}; base_of = {}
    for role, base_key, file_key in [("S1", "S1_base", "S1_file"), ("S2", "S2_base", "S2_file")]:
        if file_key in lad:
            le = read_abif(lad[file_key])
            lt = pp.baseline(pp.smooth(pp.correct_saturation(get_channel(le, CONFIG["data_channel"]))))
            # warp ladder onto reference axis via its size standard if present
            if CONFIG["size_channel"]:
                la = align.size_standard_align({"L": {"RX": lt}},
                        {"L": get_channel(le, CONFIG["size_channel"])},
                        reference="L", roles=("RX",))  # ladder on its own axis; see docs to co-warp
                lt = la["L"]["RX"]
            ladder_area[role] = peaks.quantify(lt, pos, mode="gaussian")["area"]
            base_of[role] = lad[base_key]
    called = seqmod.call_ladder_bases(pos, ladder_area, base_of)
    nt_number = seqmod.register_to_sequence(pos, called, seq)
else:
    nt_number = np.arange(len(seq), len(seq) - len(pos), -1)   # approximate, from 3' end
    print("No ladder: numbering is APPROXIMATE (from the 3' end).")
print("nucleotide range:", int(np.nanmin(nt_number)), "-", int(np.nanmax(nt_number)))''')

md("""## Step 8 — Export (.shape, .map, and a GraphPad Prism–ready Excel file)
Writes reactivity files for RNAstructure/ViennaRNA and a tidy `.xlsx` you can import straight into
GraphPad Prism to make your own graphs and run statistics (see the **Statistics Guide** in the repo).""")
code('''from shapece import report

nt = np.asarray(nt_number)
base = [seq[int(n) - 1] if np.isfinite(n) and 1 <= int(n) <= len(seq) else "N" for n in nt]

df = report.reactivity_table(np.nan_to_num(nt).astype(int), base, profiles)
# group replicate columns per condition for mean/SD (edit to match your labels):
# df = report.add_group_stats(df, {"conditionA": ["plus_rep1", "plus_rep2"]})

report.write_excel_for_prism("reactivity_for_prism.xlsx", df)
first = list(profiles)[0]
report.write_shape(f"{first}.shape", nt, profiles[first], seq_len=len(seq))
report.write_map(f"{first}.map", nt, profiles[first], seq)
print("wrote reactivity_for_prism.xlsx,", f"{first}.shape,", f"{first}.map")
df.head()''')

md("""## Step 9 — SHAPE-directed structure prediction (ViennaRNA)
Fold your RNA using the reactivities as constraints (Deigan method, standard parameters). We show
one condition here; repeat for others. A **thermodynamic (no-SHAPE) fold** is also shown as a baseline.

**What to look for:** in a good SHAPE-directed model, high-reactivity nucleotides (red) fall in loops
and low-reactivity ones (dark) in stems.""")
code('''from shapece import structure

# map the chosen condition's reactivity onto the full sequence (NaN where no data)
r_full = np.full(len(seq), np.nan)
for n, v in zip(nt, profiles[first]):
    if np.isfinite(n) and 1 <= int(n) <= len(seq):
        r_full[int(n) - 1] = v

ss_shape, mfe_shape = structure.fold(seq, r_full, method="deigan")
ss_thermo, mfe_thermo = structure.fold(seq, None)
print(f"SHAPE-directed: {mfe_shape:.1f} kcal/mol")
print(f"thermodynamic : {mfe_thermo:.1f} kcal/mol")

fig, ax = plt.subplots(2, 1, figsize=(14, 6))
structure.plot_arc(seq, ss_shape, r_full, ax=ax[0], title="SHAPE-directed (colored by reactivity)")
structure.plot_arc(seq, ss_thermo, None, ax=ax[1], title="thermodynamic baseline")
plt.tight_layout(); plt.show()
structure.write_ct(f"{first}_shape.ct", seq, ss_shape)
print("wrote", f"{first}_shape.ct")''')

md("""## Appendix — parameters, tips, and statistics

**Key parameters you may tune** (all have sensible defaults):
- `smooth(window, method)`, `baseline(window)` — preprocessing strength.
- `detect_peaks(min_spacing, prominence_frac)` — peak sensitivity.
- `quantify(mode="gaussian"|"trapezoid")` — peak integration.
- `fold(method="deigan"|"zarringhalam", slope, intercept, temperature)` — folding.

**Statistics:** this notebook deliberately does **not** run hypothesis tests — it exports a clean
Prism-ready workbook so you control the analysis. See **`docs/STATISTICS_GUIDE.md`** in the repo for
recommended tests (replicate reproducibility, per-nucleotide comparisons with multiple-testing
correction, and how to set up the Prism tables).

**Citation:** if this pipeline is used in a publication, please cite QuShape (Karabiber et al.,
*RNA* 2013), ViennaRNA (Lorenz et al., 2011), the SHAPE normalization/folding references
(Deigan et al., 2009; Low & Weeks, 2010), and this repository (see README).""")

nb["cells"] = cells
nb.metadata = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}}
nbf.write(nb, "notebooks/SHAPE_CE_Analysis.ipynb")
print("notebook written with", len(cells), "cells")
