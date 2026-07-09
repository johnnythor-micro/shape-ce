# shapece — modern capillary-electrophoresis probing analysis
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/johnnythor-micro/shape-ce/blob/main/notebooks/SHAPE_CE_Analysis.ipynb)

**`shapece`** is a modern, notebook-based pipeline for analyzing capillary-electrophoresis (CE)
nucleic-acid chemical-probing experiments — **SHAPE, DMS, CMCT, hydroxyl-radical, and enzymatic** —
from raw `.fsa` traces to **per-nucleotide reactivities** and **SHAPE-directed secondary structure**.

It is a Python 3 reimplementation of the analysis workflow introduced by **QuShape**
(Karabiber, McGinnis, Favorov \& Weeks, *RNA* 2013), rebuilt without the legacy Python 2 / PyQt4 GUI,
and extended with **internal size-standard alignment** and **ViennaRNA-based structure prediction**.
It is intended as a *modern, scriptable alternative* to QuShape — not a replacement or an official
successor. It is **not affiliated with or endorsed by the original QuShape authors**.

* **No GUI** — a documented Jupyter/Google-Colab notebook drives the analysis.
* **Beginner-friendly** — you edit a single configuration cell; every step has a diagnostic plot.
* **Generalizable \& unbiased** — no assumptions about your RNA, organism, probe, or hypothesis.
* **Reproducible \& citable** — small, documented, MIT-licensed importable package.

\---

## Why this exists

QuShape remains the reference tool for CE-probing analysis, but its implementation
(Python 2.7, PyQt4, Matplotlib-Qt4Agg) no longer installs on modern systems. `shapece` ports the
scientific algorithms to maintained libraries (NumPy/SciPy/pandas/matplotlib/ViennaRNA), runs in a
notebook (including free Google Colab), and adds two capabilities QuShape lacks: alignment via an
internal size standard (e.g. GeneScan LIZ) and one-click SHAPE-directed folding.

## Pipeline

|Stage|Module|What it does|
|-|-|-|
|Read|`io`|Parse ABIF `.fsa`; map dye/DATA channels to logical roles (RX, BG, ladders, size std).|
|Preprocess|`preprocess`|Saturation repair, Savitzky–Golay/triangle smoothing, morphological baseline, exponential decay correction, dye mobility-shift.|
|Align|`align`|Co-register lanes by **internal size-standard warping** (recommended) or **banded dynamic time warping** (QuShape-style, when no standard).|
|Peaks|`peaks`|Derivative peak detection + iterative **Gaussian deconvolution** (area = amp × width); trapezoid option.|
|Reactivity|`reactivity`|**Area-difference** (QuShape-style, default) or **stop-fraction** (`-ln(1-f)`); background scaling on low-reactivity peaks; model-free **boxplot (2–8%) normalization**.|
|Sequence|`sequence`|Ladder-driven base registration by Needleman–Wunsch alignment to the reference sequence; approximate size-standard fallback.|
|Structure|`structure`|**ViennaRNA** SHAPE-directed folding (Deigan / Zarringhalam), base-pair probabilities, arc plots, `.ct` export.|
|Report|`report`|`.shape`, `.map` (RNAstructure), and **GraphPad Prism–ready `.xlsx`** exports.|

## Reactivity models

* **Area-difference (default).** `reactivity = area(RX) − factor·area(BG)`, where `factor` scales the
background to the reagent using the lowest-reactivity peaks (QuShape's approach). Chosen as the
default for continuity with the existing CE-SHAPE literature.
* **Stop-fraction.** Per-lane `-ln(1 − f)` where `f` is the fraction of molecules reaching a
nucleotide that stop there (Aviran/Weeks). Loading-independent and depletion-corrected; useful for
stop-based readouts.

Both are followed by the same **model-free boxplot normalization** (Low \& Weeks): exclude high
outliers (> Q3 + 1.5·IQR), divide by the mean of the top \~10% of the remainder.

## Alignment

* **`size\_standard\_align`** (recommended) — detects the size-standard peaks *directly and identically*
in every lane and warps each lane's data channels onto a reference lane. Because the standard runs in
the same capillary as the data, this co-registers the data exactly and is robust to instruments that
call inconsistent numbers of standard peaks per lane (a real failure mode we encountered).
* **`dtw\_align`** — QuShape-style banded dynamic time warping on a shared channel, for designs without
an internal size standard.

## Installation

```bash
pip install numpy scipy pandas matplotlib openpyxl ViennaRNA
pip install git+https://github.com/johnnythor-micro/shape-ce.git
```

All dependencies are open-source and install on Linux/macOS/Windows and Google Colab.

## Quickstart

**Notebook (recommended):** open `notebooks/SHAPE\_CE\_Analysis.ipynb` (locally or via the
"Open in Colab" badge), edit the `CONFIG` cell with your file paths / channels / sequence, and run
top to bottom.

**As a library:**

```python
import numpy as np, shapece as sc
from shapece.io import read\_abif, get\_channel

e   = read\_abif("plus\_rep1.fsa")
rx  = sc.preprocess.baseline(sc.preprocess.smooth(get\_channel(e, 9)))
# ... align, detect peaks, quantify, compute reactivity ...
ss, mfe = sc.structure.fold(sequence, reactivity, method="deigan")   # SHAPE-directed fold
```



## Between-condition comparison and visualization (v0.2)

`shapece.stats` and `shapece.plots` add differential analysis and publication plots.

### Statistics — `shapece.stats`

|Function|Question answered|
|-|-|
|`probe\_signal\_qc`|**Did the probe work at all?** Run first — see warning below.|
|`replicate\_stats`|Are replicates reproducible (leave-one-out r)?|
|`delta\_shape`|**Which regions changed?** Faithful deltaSHAPE port (Smola et al. 2015).|
|`ttest\_fdr`|Which individual nucleotides changed? (Welch t-test + BH-FDR + effect-size gate)|
|`permutation\_test`|Did anything change globally? (exact label permutation)|

`delta\_shape` is a verbatim port of the published algorithm — smoothing, Z-factor, standard
score, and the "≥3 significant nucleotides in a 5-nt window" site rule — validated
bit-for-bit against a transcription of the original Python-2 source (`tests/test\_stats.py`).
The one adaptation for CE: since capillary electrophoresis has no per-nucleotide error, the
**SEM across replicates** is used (≥2 replicates required).

> ⚠️ \*\*A reproducible background is not a measurement.\*\* When probe signal is weak, the (+) and
> (−) lanes are nearly identical, background subtraction cannot isolate reactivity, and the
> surviving natural-stop pattern is \*highly reproducible\* — so replicates agree beautifully and
> conditions look identical. `probe\_signal\_qc` detects exactly this. Run it before any
> comparison; if it fails, no downstream statistic can rescue the experiment.

### Visualization — `shapece.plots`

`skyline` (step plot; overlay conditions without occlusion) · `profile` (single condition,
SHAPE-colored, with error bars) · `delta\_shape\_plot` (smoothed Δ with significant sites shaded)
· `arc` (base pairs as arcs; two structures mirrored; accepts dot-bracket, pair list, or
probability matrix) · `circle` (chords on a circle — for long RNAs and long-range contacts) ·
`linreg` (replicate R²) · `heatmap` (spot a bad lane).

Plot vocabulary follows RNAvigate (Ehrhardt \& Weeks, *NAR* 2024). Every plotter returns a
matplotlib `Axes`.

### Point-and-click notebooks (v0.3)

`shapece.ui` turns configuration into widgets — no file paths, no code editing:
an **upload button**, a **sequence box** that cleans FASTA/whitespace/`T`→`U` and reports length
and GC, **channel dropdowns** populated from your actual file, and a **file table** where you set
each lane's role, condition and replicate (pairs are matched for you).
Requires `ipywidgets` (pre-installed on Colab); everything still works without it by writing a
plain `CONFIG` dict.

### Notebook and docs

* `notebooks/SHAPE\_CE\_Analysis.ipynb` — guided raw-trace → reactivity → structure.
* `notebooks/SHAPE\_CE\_Comparison.ipynb` — guided two-condition comparison.
* [**`docs/DIFFERENTIAL\_ANALYSIS.md`**](docs/DIFFERENTIAL_ANALYSIS.md) — the statistics, with power caveats.
* [**`docs/VISUALIZATION\_GUIDE.md`**](docs/VISUALIZATION_GUIDE.md) — what each plot demonstrates and when to use it.

## Statistics

`shapece` intentionally does **not** perform hypothesis testing; it exports a tidy Prism-ready
workbook so you control the analysis. See [**`docs/STATISTICS\_GUIDE.md`**](docs/STATISTICS_GUIDE.md)
for recommended tests (replicate reproducibility, per-nucleotide comparison with FDR correction, and
Prism table setup).

## Validation

`shapece` is validated against **QuShape's own published TPP-riboswitch reference analysis**
(the finished `test1001\_done.qushape` project shipped with QuShape), at four levels of strictness:

|Level|What is tested|Result|
|-|-|-|
|1|reactivity arithmetic (`areaDiff`)|**exact** (residual 5e−21)|
|2|boxplot normalization vs `normDiff`|**r = 1.0000000000** (see note below)|
|3|Gaussian peak areas|**r = 0.995 / 0.993** (RX / BG)|
|3|peaks → reactivity → normalization|**r = 0.9889**|
|4|end-to-end from raw `.fsa` (DTW alignment)|**r = 0.9061**|
|4|+ fully independent peak detection|**r = 0.8869**, 78/88 peaks matched|

!\[validation](docs/shapece\_TPP\_validation.png)

Level 2 is exact in shape but the normalization *divisor* differs by 6.5%: we traced this to the
shipped reference project having been generated by an **earlier QuShape build** whose boxplot outlier
count on these data was 1, whereas the currently published `findPOutlierBox` source computes 2.
`shapece` follows the current published source. Full details, including how to reproduce every number,
are in [**`docs/VALIDATION.md`**](docs/VALIDATION.md); the tests are in
[`tests/test\_validation\_tpp.py`](tests/test_validation_tpp.py).

## Citing

If you use this pipeline, please cite the repository and the underlying methods:

* Karabiber F, McGinnis JL, Favorov OV, Weeks KM. **QuShape**: rapid, accurate, and best-practices
quantification of nucleic acid probing information, resolved by capillary electrophoresis.
*RNA* 19:63–73 (2013).
* Lorenz R, et al. **ViennaRNA Package 2.0**. *Algorithms Mol Biol* 6:26 (2011).
* Deigan KE, Li TW, Mathews DH, Weeks KM. Accurate SHAPE-directed RNA structure determination.
*PNAS* 106:97–102 (2009).
* Low JT, Weeks KM. SHAPE-directed RNA secondary structure prediction. *Methods* 52:150–158 (2010).
* (stop-fraction) Aviran S, et al. *PNAS* 108:11069–11074 (2011).
* deltaSHAPE (differential SHAPE; ported in `shapece.stats`): Smola MJ, Calabrese JM, Weeks KM.
*Biochemistry* 54:6867–6875 (2015).
* RNAvigate (plot vocabulary): Ehrhardt JE, Weeks KM. *Nucleic Acids Research* (2024).

## License

MIT — see [`LICENSE`](LICENSE).

