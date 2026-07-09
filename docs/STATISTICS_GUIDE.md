# Statistics guide: comparing SHAPE-CE reactivities

This pipeline **does not** run hypothesis tests. It exports a clean, tidy
`reactivity_for_prism.xlsx` so you keep full control of the statistics in GraphPad
Prism (or any tool). This guide describes the correct, defensible way to analyze
CE probing reactivities and how to set up the Prism tables.

> Reactivities are already **normalized** (model-free boxplot / 2–8% rule), so they
> are comparable across replicates and conditions on a common ~0–2 scale.

---

## 1. Always check replicate reproducibility first

Before comparing conditions, confirm your replicates agree. Poor reproducibility
means any "difference" may be noise.

- **In Prism:** make an XY table with nucleotide on X and each replicate's
  reactivity as a Y column. Under **Analyze → Correlation**, compute the Pearson
  correlation between replicate columns. Good CE-SHAPE replicates typically give
  **r ≳ 0.8–0.9**. Flag and investigate any outlier replicate before pooling.
- A quick scatter of rep 1 vs rep 2 (XY, one point per nucleotide) should fall on
  the diagonal.

## 2. Choose the comparison that matches your question

**(a) Global: "are two conditions different overall?"**
The cleanest, best-powered test with few replicates. Summarize each replicate by
its whole-profile similarity to the other condition and compare within- vs
between-condition. With only n = 3 this is often a **permutation/randomization**
argument rather than a t-test (few possible label permutations → limited power;
report the exact permutation p-value and be honest about power).

**(b) Per-nucleotide: "which nucleotides change?"**
For each nucleotide, compare the replicate reactivities between conditions.

- **In Prism:** use a **Grouped** table — rows = nucleotides, sub-columns =
  replicates, grouped by condition. Run **multiple t-tests (one per row)** or,
  better, **Analyze → Multiple t-tests / Two-stage step-up (Benjamini, Krieger,
  Yekutieli) FDR**. Report **q-values**, not raw p-values, because you are testing
  hundreds of nucleotides. A common threshold is **q < 0.05** (or < 0.10) **and**
  a minimum effect size (e.g. |Δreactivity| ≥ 0.3 on the normalized scale) so you
  don't call tiny, meaningless differences.

**(c) Established tools for SHAPE differences.** For per-nucleotide differential
SHAPE specifically, **deltaSHAPE** (Smola et al.) implements a standard, published
Z-factor + magnitude criterion; it is a good complement to a generic t-test/FDR.

## 3. How much power do you have?

With **n = 3** replicates and hundreds of nucleotides, per-nucleotide tests are
**underpowered after multiple-testing correction** — it is common for *zero*
individual nucleotides to survive FDR even when a real global difference exists.
If per-nucleotide resolution matters, plan for **more replicates** (n ≥ 4–6),
prioritizing the noisier condition, and make sure the probing chemistry is strong
(a weak (+)/(−) separation limits everything downstream).

## 4. Setting up the Prism import

The exported workbook has a `reactivity` sheet: one **row per nucleotide**, columns
= `nucleotide`, `base`, then one column per replicate.

- **XY graph of a profile:** New table → XY → paste `nucleotide` into X and one
  reactivity column into Y. Repeat/overlay for replicates or conditions.
- **Grouped comparison + stats:** New table → Grouped → rows = nucleotides;
  create a sub-column group per condition and paste that condition's replicate
  columns as the sub-columns. Then run multiple t-tests with FDR (Section 2b).
- **Mean ± SD bar/line:** use the `add_group_stats()` export option (adds
  `<condition>_mean` and `<condition>_SD` columns) and plot those directly.

## 5. Reporting checklist

- [ ] Replicate correlations reported (within each condition).
- [ ] Normalization method stated (model-free boxplot / 2–8%).
- [ ] Multiple-testing correction applied for per-nucleotide claims (report q).
- [ ] Effect-size threshold stated alongside significance.
- [ ] Power/limitations noted honestly (n, chemistry strength, registration).
- [ ] Tools cited (this pipeline; QuShape; ViennaRNA; deltaSHAPE if used).
