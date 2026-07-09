"""
shapece -- a modern, notebook-friendly pipeline for capillary-electrophoresis
nucleic-acid probing analysis (SHAPE, DMS, CMCT, hydroxyl-radical, enzymatic).

A lightweight, dependency-modest reimplementation of the QuShape analysis
workflow (Karabiber et al., RNA 2013) in Python 3, extended with internal
size-standard alignment and ViennaRNA-based SHAPE-directed structure prediction.
Not affiliated with or endorsed by the original QuShape authors.

Modules
-------
io           ABIF (.fsa) reading and channel mapping
preprocess   saturation / smoothing / baseline / decay / mobility corrections
align        cross-lane co-registration (size-standard warp; DTW)
peaks        peak detection + Gaussian-deconvolution quantification
reactivity   area-difference (default) & stop-fraction models; boxplot norm
sequence     ladder-driven base registration to the reference sequence
structure    ViennaRNA SHAPE-directed folding + 2D plots
report       .shape / .map / Prism-ready .xlsx exports
stats        probe-signal QC, deltaSHAPE differential analysis, FDR + permutation
plots        skyline, profile, deltaSHAPE, arc, circle, linreg, heatmap
ui           notebook widgets: upload buttons, sequence box, channel pickers
"""
__version__ = "0.3.0"
from . import io, preprocess, align, peaks, reactivity, sequence, structure, report  # noqa
from . import stats, plots  # noqa
# `ui` imports ipywidgets lazily (inside functions), so importing it here is cheap
# and safe even when ipywidgets is not installed.
from . import ui  # noqa
