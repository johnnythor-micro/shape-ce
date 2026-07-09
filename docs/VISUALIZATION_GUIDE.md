# Visualization guide: `shapece.plots`

Every plotter returns a matplotlib `Axes`, so plots compose into multi-panel figures and stay
fully editable for publication. Plot vocabulary follows RNAvigate (Weeks lab), so figures made
here are legible to readers of the RNA-structure literature.

**Color convention** (shared with `shapece.structure`):

| Reactivity | Color | Interpretation |
|---|---|---|
| ≥ 0.85 | red | flexible / likely unpaired |
| 0.40 – 0.85 | orange | intermediate |
| < 0.40 | black | constrained / likely paired |
| no data | gray | not measured (e.g. primer site) |

---

## Skyline (step) plot — `plots.skyline`

```python
plots.skyline({"37 °C": mean37, "32 °C": mean32})
```

**What it shows.** Reactivity as a step function: each nucleotide is a full-width step rather
than a thin bar.

**Why use it.** Because the steps tile the axis, **several profiles can be superimposed without
occluding one another**. Bars from two conditions would overlap and hide each other; skylines
do not. This is the standard way to put two conditions — or all your replicates — on one axis
and see where they diverge.

**Read it as:** where the two traces separate, reactivity changed. Where they track together,
the local structure is unchanged (or unmeasured).

## Profile plot — `plots.profile`

```python
plots.profile(mean37, error=sem37)
```

**What it shows.** One condition, as bars colored by the SHAPE convention, with error bars
(typically the SEM across replicates).

**Why use it.** The clearest way to **report a single dataset**: red bars mark flexible,
likely-unpaired nucleotides; black bars mark constrained, likely-paired ones. The dashed line
at 1.0 is the normalization reference — a well-behaved profile has most nucleotides below it
with a minority of clearly reactive positions above.

## deltaSHAPE plot — `plots.delta_shape_plot`

```python
res = stats.delta_shape(None, None, replicates1=reps37, replicates2=reps32)
plots.delta_shape_plot(res, labels=("37 °C", "32 °C"))
```

**What it shows.** The smoothed difference between conditions, with **statistically significant
sites shaded**: orange where condition 1 is more reactive, blue where condition 2 is. Dots mark
nucleotides whose 95% confidence intervals don't overlap (Z-factor > 0).

**Why use it.** It is the one figure that answers "what changed, where, and is it real?" A
positive (orange) region means the RNA became *more* flexible/accessible in condition 1 — a
helix melting, a ligand leaving, a protein dissociating. Unshaded differences did not pass the
significance filters and should not be interpreted.

## Arc plot — `plots.arc`

```python
plots.arc(dot_bracket, reactivity=mean37, sequence=seq)
plots.arc(shape_directed, structure2=thermodynamic, sequence=seq)   # two models compared
```

**What it shows.** Base pairs as arcs spanning the sequence axis. Accepts a dot-bracket string,
a list of `(i, j)` pairs, or a **base-pair probability matrix** (arc opacity ∝ probability).
A second structure is drawn mirrored *below* the axis.

**Why use it.** Arcs preserve the linear sequence, so you can align structure directly against
reactivity: in a good model, red (reactive) nucleotides fall in loops and bulges — the gaps
between arcs — while black nucleotides sit under stacked arcs. **Nested arcs are helices;
crossing arcs indicate pseudoknots.** Mirroring two structures makes agreements appear as
reflections and disagreements jump out.

**Limitation.** Arc height grows with pairing distance, so on long RNAs (≳ 500 nt) long-range
pairs dominate the figure. Use a circle plot instead.

## Circle plot — `plots.circle`

```python
plots.circle(dot_bracket, reactivity=mean37, sequence_length=len(seq))
plots.circle(bpp_matrix, sequence_length=len(seq), threshold=0.3)    # any pair matrix
```

**What it shows.** The sequence wrapped around a circle, with pairs drawn as chords across the
interior. Rim nucleotides are colored by reactivity; chord opacity scales with pair weight.

**Why use it.** Every contact is rendered at the same visual scale regardless of how far apart
the partners are, so **long-range interactions stay legible** where they would overwhelm an arc
plot. This is the right choice for long RNAs, for comparing base-pair probability matrices, and
for any pair/correlation matrix (it accepts an arbitrary square matrix, not just structures).

**Read it as:** chords crossing the center are long-range; short chords hugging the rim are
local helices.

## Replicate regression — `plots.linreg`

```python
plots.linreg(rep1, rep2)
```

**What it shows.** One replicate against another, one point per nucleotide, with a fitted line
and R².

**Why use it.** The reproducibility check that should **precede any comparison**. Points should
hug the diagonal. If R² is low, the between-condition analysis that follows is not
interpretable — fix the experiment, not the statistics.

## Heatmap — `plots.heatmap`

```python
plots.heatmap(np.vstack([reps37, reps32]), row_labels=[...])
```

**What it shows.** Replicates (rows) × nucleotides (columns), colored by reactivity.

**Why use it.** The fastest way to **spot a bad lane**: a row that looks unlike its neighbors is
an outlier replicate. Also reveals systematic column-wise features (e.g. strong natural stops
present in every lane).

---

## A suggested figure order for a paper

1. `linreg` + `heatmap` — replicate quality (often supplementary).
2. `profile` — reactivity for each condition.
3. `skyline` — the two conditions overlaid.
4. `delta_shape_plot` — what changed and where.
5. `arc` (short RNA) or `circle` (long RNA) — reactivity mapped onto structure.

## References

- Ehrhardt JE, Weeks KM. RNAvigate: efficient exploration of RNA chemical probing datasets.
  *Nucleic Acids Research* (2024). — plot vocabulary
- Smola MJ, Calabrese JM, Weeks KM. *Biochemistry* 54:6867–6875 (2015). — deltaSHAPE plot
