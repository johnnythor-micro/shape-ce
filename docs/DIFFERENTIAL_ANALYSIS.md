# Comparing two conditions: `shapece.stats`

This guide covers the built-in statistics for asking **"did reactivity change between
condition A and condition B?"** — for example two temperatures, ± ligand, or ± protein.

It complements [`STATISTICS_GUIDE.md`](STATISTICS_GUIDE.md), which describes how to export
data and run your own tests in GraphPad Prism. Use whichever suits your workflow; the numbers
are the same.

---

## Step 0 (mandatory): did the probe work at all?

```python
from shapece import stats
qc = stats.probe_signal_qc(area_rx, area_bg)
print(qc["passed"], qc["reasons"])
```

| Diagnostic | Healthy | Meaning if it fails |
|---|---|---|
| `rx_bg_correlation` | **< ~0.95** | (+) and (−) lanes are nearly the same lane |
| `mean_area_ratio` | **> 1** (clearly) | the reagent added few stops |
| `background_carryover` | **≈ 0** | "reactivity" is still the background pattern, rescaled |

**Why this matters.** When probe signal is weak, background subtraction cannot isolate
reactivity, and what survives is the natural-stop pattern. That pattern is *highly
reproducible*, so replicates agree beautifully and conditions look nearly identical —
a result that is easy to mistake for "clean data showing no difference."

> **A reproducible background is not a measurement.** If `probe_signal_qc` fails, stop.
> No downstream statistic can rescue the experiment. Repeat it with fresh reagent.

## Step 1: check replicate reproducibility

```python
rs = stats.replicate_stats(reps_condition1)   # (n_reps, n_nt)
print(rs["loo_r"])        # each replicate vs. the mean of the others
print(rs["pairwise_r"])
```

Good CE-SHAPE replicates typically reach **r ≈ 0.8–0.9**. A single outlier replicate should be
investigated (and usually excluded, transparently) before any comparison. Visualize with
`plots.linreg` and `plots.heatmap`.

## Step 2: deltaSHAPE — which *regions* changed?

A faithful port of deltaSHAPE (Smola et al. 2015), adapted for CE by using the
**standard error of the mean across replicates** in place of MaP's per-nucleotide error.

```python
res = stats.delta_shape(None, None,
                        replicates1=reps_37C,    # (n_reps, n_nt)
                        replicates2=reps_32C)
print(res["n_sites"], res["sites"])
```

What it does, in order:

1. **Smooth** reactivities and errors over a centered `2*pad+1` window (`pad=1` → 3 nt).
2. **Difference**: `diff = condition1 − condition2`, then smooth it.
3. **Z-factor** = `1 − 1.96·(err₁+err₂)/|Δ|`. A value **> 0** means the 95% confidence
   intervals of the two conditions do **not** overlap at that nucleotide.
4. **Standard score** of the smoothed difference: how large is this change relative to
   *all other* changes in the molecule?
5. **Site calling**: a nucleotide passes if `Z-factor > 0` **and** `|standard score| ≥ 1`.
   A **site** requires **≥ 3 passing nucleotides within a 5-nt window**.

Sign convention: `diff = condition1 − condition2`, so a **positive** site means *more reactive
in condition 1* (more flexible/accessible — e.g. a helix melted, or a protein left).

**Why the windowing matters.** Requiring changes to cluster is not an arbitrary filter: real
structural rearrangements affect several adjacent nucleotides, while noise does not. This is
where deltaSHAPE gets its power, and it is why it can find real signal where a per-nucleotide
test corrected for multiple comparisons finds nothing.

Key parameters (all exposed): `pad`, `z_coeff`, `z_thresh`, `ss_thresh`, `site_pad`, `site_min`,
`mask5`, `mask3`. Use `mask5`/`mask3` to exclude primer-binding nucleotides.

Plot it with `plots.delta_shape_plot(res, labels=("37 °C", "32 °C"))`.

## Step 3: per-nucleotide testing with FDR control

```python
tt = stats.ttest_fdr(reps_37C, reps_32C, min_delta=0.3, q_thresh=0.10)
print(tt["n_significant"])
```

Welch t-test per nucleotide, Benjamini–Hochberg q-values, **plus an effect-size gate**
(`|Δ| ≥ min_delta`) so that statistically detectable but biologically trivial differences are
not called.

> ⚠️ **Power warning.** With n = 3 replicates across hundreds of nucleotides, this test is
> usually underpowered after FDR correction — **zero** significant nucleotides is a common
> outcome *even when a real global difference exists*. Read a null result as "insufficient
> power", not "no difference". Plan for n ≥ 4–6 if per-nucleotide resolution matters.

`delta_shape` and `ttest_fdr` answer different questions and are **complementary**:
the former asks "which regions changed?" and borrows strength across neighbors;
the latter asks "which individual nucleotides changed?" and controls the FDR across them.

## Step 4: is there *any* global difference?

```python
pt = stats.permutation_test(reps_37C, reps_32C)
print(pt["p_value"], pt["min_attainable_p"])
```

Exact label-permutation test on (within-condition correlation) − (between-condition correlation).

> ⚠️ With 3 vs 3 replicates there are only **10** distinct partitions, so the smallest possible
> p-value is **0.10**. A non-significant result under this design cannot establish that the
> conditions are identical. `min_attainable_p` is returned so you can state this honestly.

---

## Reporting checklist

- [ ] `probe_signal_qc` reported and passed.
- [ ] Replicate reproducibility (`loo_r`) reported per condition.
- [ ] Normalization method stated (model-free boxplot / 2–8%).
- [ ] For per-nucleotide claims: FDR-corrected q-values **and** an effect-size threshold.
- [ ] For regional claims: deltaSHAPE parameters (`pad`, thresholds, site rule) stated.
- [ ] Power limitations noted (n, `min_attainable_p`, probe strength).
- [ ] Tools cited: this pipeline; deltaSHAPE (Smola et al. 2015); QuShape (Karabiber et al. 2013).

## References

- Smola MJ, Calabrese JM, Weeks KM. Detection of RNA–protein interactions in living cells with
  SHAPE. *Biochemistry* 54:6867–6875 (2015). — deltaSHAPE
- Karabiber F, McGinnis JL, Favorov OV, Weeks KM. QuShape. *RNA* 19:63–73 (2013).
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *JRSS B* 57:289–300 (1995).
