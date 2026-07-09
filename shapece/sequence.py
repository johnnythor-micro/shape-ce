"""
shapece.sequence
================
Assign an absolute nucleotide number to every quantified peak, using sequencing
ladder channel(s) and the known RNA sequence -- QuShape's "Sequence" stage,
reimplemented as a Needleman-Wunsch alignment of the ladder-called base pattern
to the reference sequence.

Reverse transcription reads 3'->5', so peaks are matched to the sequence in
reverse. Ladders mark specific bases (you specify which): e.g. ddATP terminates
at template U, ddGTP at template C.

If no ladder is available, :func:`register_by_size_standard` gives an
approximate linear map from an internal size standard instead (clearly flagged
as approximate).
"""
from __future__ import annotations
import numpy as np


def call_ladder_bases(peak_pos: np.ndarray, ladder_area: dict, base_of: dict,
                      enrich: float = 1.5) -> list[str]:
    """Tentatively call each peak's base from ladder enrichment.

    Parameters
    ----------
    peak_pos : np.ndarray
        Consensus peak positions.
    ladder_area : dict[str, np.ndarray]
        Role -> per-peak area for each ladder channel (e.g. {"S1":..., "S2":...}),
        already background-subtracted or envelope-normalized.
    base_of : dict[str, str]
        Which base each ladder marks, e.g. {"S1": "U", "S2": "C"}.
    enrich : float
        A peak is 'called' for a ladder if its area exceeds ``enrich`` x the
        median area of that ladder.

    Returns
    -------
    list[str]
        Per-peak base call: the marked base, or 'N' if ambiguous/uncalled.
    """
    n = len(peak_pos)
    calls = ["N"] * n
    for role, areas in ladder_area.items():
        med = np.median(areas[areas > 0]) if np.any(areas > 0) else 0.0
        thr = enrich * med
        for i in range(n):
            if areas[i] > thr:
                calls[i] = base_of[role] if calls[i] == "N" else calls[i]
    return calls


def _nw_align(called: list[str], seq: str, match=2.0, mismatch=-1.0,
              n_neutral=0.0, gap=-2.0):
    """Needleman-Wunsch alignment of a base-call list to a sequence string.

    Returns the aligned index pairs (call_index, seq_index) for matched columns.
    """
    n, m = len(called), len(seq)
    S = np.zeros((n + 1, m + 1))
    S[0, :] = np.arange(m + 1) * gap
    S[:, 0] = np.arange(n + 1) * gap
    arrow = np.zeros((n + 1, m + 1), int)
    arrow[0, :] = 2; arrow[:, 0] = 1
    for i in range(1, n + 1):
        c = called[i - 1]
        for j in range(1, m + 1):
            if c == "N":
                sc = n_neutral
            elif c == seq[j - 1]:
                sc = match
            else:
                sc = mismatch
            diag = S[i - 1, j - 1] + sc
            up = S[i - 1, j] + gap
            left = S[i, j - 1] + gap
            best = max(diag, up, left)
            S[i, j] = best
            arrow[i, j] = 0 if best == diag else (1 if best == up else 2)
    # backtrace
    i, j = n, m
    pairs = []
    while i > 0 and j > 0:
        if arrow[i, j] == 0:
            pairs.append((i - 1, j - 1)); i -= 1; j -= 1
        elif arrow[i, j] == 1:
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


def register_to_sequence(peak_pos: np.ndarray, called: list[str], rna_seq: str) -> np.ndarray:
    """Return the 1-based nucleotide number for each consensus peak (NaN if unmapped).

    Peaks are matched to ``rna_seq`` read 3'->5'. The alignment is driven by the
    ladder base-calls; uncalled ('N') peaks are carried along between anchors.
    """
    rev = rna_seq[::-1]
    n = len(rna_seq)
    pairs = _nw_align(called, rev, )
    num = np.full(len(peak_pos), np.nan)
    for pi, sj in pairs:
        num[pi] = n - sj          # reverse index -> forward 1-based position
    # fill uncalled peaks by interpolation between anchored ones (monotone)
    anc = np.where(~np.isnan(num))[0]
    if len(anc) >= 2:
        num = np.interp(np.arange(len(peak_pos)), anc, num[anc])
        num = np.round(num)
    return num


def register_by_size_standard(peak_pos: np.ndarray, size_map, seq_len: int,
                              primer_len: int = 0):
    """Approximate nucleotide numbering from a size->position model (no ladder).

    ``size_map`` maps a peak's electrophoretic size to a nucleotide number; here
    we accept a callable pos(size) or a linear (alpha, beta) tuple. Flagged as
    approximate: prefer a ladder when available.
    """
    if callable(size_map):
        return np.round(size_map(peak_pos))
    alpha, beta = size_map
    return np.round(seq_len - (peak_pos - alpha) / beta)
