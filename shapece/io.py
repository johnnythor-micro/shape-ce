"""
shapece.io
==========
Reading capillary-electrophoresis (CE) trace files (ABIF ``.fsa``) and mapping
their dye channels onto the four logical roles used throughout the pipeline.

A modern, dependency-light ABIF reader is included (no Biopython required),
adapted in spirit from the ABIF specification used by Applied Biosystems
instruments. It returns the raw processed traces (``DATA`` tags) plus any
instrument-called size-standard peaks.

Logical channel roles
----------------------
* ``RX``   : the (+) reagent lane          (e.g. 1M7 / DMS / CMCT modified)
* ``BG``   : the (-) background lane        (vehicle / no-reagent control)
* ``S1``   : first sequencing ladder        (e.g. ddATP -> marks U, optional)
* ``S2``   : second sequencing ladder       (e.g. ddGTP -> marks C, optional)
* ``SIZE`` : the internal size standard     (e.g. GeneScan LIZ, optional)

Every role maps to a dye/``DATA`` channel number that you choose in the
notebook. Nothing about a specific probe, organism, or experiment is assumed.
"""
from __future__ import annotations
import struct
from dataclasses import dataclass, field
import numpy as np

# ABIF element type codes -> (struct format, byte size)
_ABIF_TYPES = {
    1: ("B", 1), 2: ("s", 1), 3: ("H", 2), 4: ("h", 2),
    5: ("i", 4), 6: ("I", 4), 7: ("f", 4), 8: ("d", 8),
    10: ("HBB", 4), 11: ("BBBB", 4), 12: ("iiHB", 10),
    18: ("s", 1), 19: ("s", 1),
}


def read_abif(path: str) -> dict:
    """Read an ABIF (``.fsa``/``.ab1``) file into a dict keyed by (tag, number).

    Parameters
    ----------
    path : str
        Path to the ABIF file.

    Returns
    -------
    dict
        Maps ``(tag_name, tag_number)`` -> value. Numeric ``DATA`` channels come
        back as lists of ints; use :func:`get_channel` to pull one as a float array.

    Notes
    -----
    Processed CE traces live in ``DATA`` tags. On 4-dye chemistries the analyzed
    dyes are usually ``DATA1..DATA4``; on 5/6-dye chemistries the analyzed traces
    are commonly ``DATA9..DATA12`` (+ ``DATA205`` for the size standard on some
    instruments). Inspect the file rather than assuming — the notebook plots all
    channels so you can pick the right ones.
    """
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw[:4] != b"ABIF":
        raise ValueError(f"{path} is not an ABIF file (bad magic bytes).")
    # Header: the directory's own entry starts at offset 6.
    num_elements = struct.unpack(">i", raw[18:22])[0]
    dir_offset = struct.unpack(">i", raw[26:30])[0]
    entries: dict = {}
    for i in range(num_elements):
        start = dir_offset + i * 28
        name = raw[start:start + 4].decode("latin-1")
        number = struct.unpack(">i", raw[start + 4:start + 8])[0]
        etype = struct.unpack(">h", raw[start + 8:start + 10])[0]
        count = struct.unpack(">i", raw[start + 12:start + 16])[0]
        dsize = struct.unpack(">i", raw[start + 16:start + 20])[0]
        doffset = struct.unpack(">i", raw[start + 20:start + 24])[0]
        data_start = doffset if dsize > 4 else start + 20
        entries[(name, number)] = _parse_value(raw, data_start, etype, count, dsize)
    return entries


def _parse_value(raw, offset, etype, count, dsize):
    if etype in (2, 18, 19):                # strings / char arrays
        return raw[offset:offset + dsize]
    if etype not in _ABIF_TYPES:
        return raw[offset:offset + dsize]
    fmt, size = _ABIF_TYPES[etype]
    if etype in (10, 11, 12):               # compound records
        out = []
        for k in range(count):
            chunk = raw[offset + k * size: offset + (k + 1) * size]
            out.append(struct.unpack(">" + fmt, chunk))
        return out
    values = struct.unpack(">" + fmt * count, raw[offset:offset + size * count])
    return list(values)


def get_channel(entries: dict, data_number: int) -> np.ndarray:
    """Return one ``DATA`` channel as a float array (e.g. ``get_channel(e, 9)``)."""
    key = ("DATA", data_number)
    if key not in entries:
        raise KeyError(f"DATA{data_number} not found. Available DATA channels: "
                       f"{sorted(n for (t, n) in entries if t == 'DATA')}")
    return np.asarray(entries[key], dtype=float)


def list_data_channels(entries: dict) -> list[int]:
    """List available ``DATA`` channel numbers, for interactive inspection."""
    return sorted(n for (t, n) in entries if t == "DATA")


def instrument_size_peaks(entries: dict):
    """Return (apex_scan, called_size) arrays from instrument size-calling, if present.

    These are convenient but *inconsistent across lanes* on some instruments, so
    the pipeline can instead detect the standard peaks directly (see
    :mod:`shapece.align`). Returns (None, None) if the file has no called peaks.
    """
    if ("Peak", 2) in entries and ("Peak", 12) in entries:
        return (np.asarray(entries[("Peak", 2)], float),
                np.asarray(entries[("Peak", 12)], float))
    return None, None


@dataclass
class Trace:
    """A single CE lane mapped to logical roles.

    Attributes
    ----------
    name : str
        User label (e.g. ``"sample_37C_rep1"``).
    channels : dict[str, np.ndarray]
        Role -> raw trace. Keys among {"RX","BG","S1","S2","SIZE"}.
    meta : dict
        Free-form metadata (file path, dye numbers, etc.).
    """
    name: str
    channels: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)


def load_trace(path: str, name: str, channel_map: dict) -> Trace:
    """Load one ABIF file and assign its DATA channels to logical roles.

    Parameters
    ----------
    path : str
        ABIF file path.
    name : str
        Label for this lane.
    channel_map : dict
        Role -> DATA channel number, e.g.
        ``{"RX": 9, "BG": 9, "S1": 9, "S2": 9, "SIZE": 205}``.
        (In separate-capillary designs RX and BG live in *different files* but the
        same DATA number; in one-capillary multi-dye designs they differ.)

    Returns
    -------
    Trace
    """
    e = read_abif(path)
    channels = {}
    for role, num in channel_map.items():
        if num is None:
            continue
        channels[role] = get_channel(e, num)
    return Trace(name=name, channels=channels,
                 meta={"path": path, "channel_map": dict(channel_map),
                       "entries": e})
