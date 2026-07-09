"""
shapece.ui
==========
Point-and-click notebook widgets, so users never have to type a file path or edit
code to configure an analysis.

Everything here is **optional sugar**. Each widget exposes a plain-Python value
(``.value`` / ``.to_config()``), so a user who prefers editing a dict can ignore
this module entirely and the rest of ``shapece`` behaves identically.

Components
----------
:func:`upload_files`      A "Choose files" button. Saves uploads to disk and
                          returns ``{filename: path}``. Uses Colab's native file
                          picker when available, otherwise ``ipywidgets``.
:class:`SequenceInput`    A text box to paste an RNA sequence or FASTA record.
                          Cleans and validates on the fly.
:class:`ChannelPicker`    Dropdowns to choose the reagent / size-standard
                          channels, populated from the actual file.
:class:`SamplePairer`     Assign each uploaded file a role, condition and
                          replicate number; builds the (+)/(-) pairs for you.
:class:`ConfigPanel`      Bundles the above into one form; ``.to_config()``
                          returns the dict the analysis notebook consumes.

Requires ``ipywidgets`` (pre-installed on Google Colab):  ``pip install ipywidgets``
"""
from __future__ import annotations

import os

__all__ = ["widgets_available", "upload_files", "SequenceInput", "ChannelPicker",
           "SamplePairer", "ConfigPanel", "clean_sequence"]

_BOX = ("border:1px solid #d0d7de; border-radius:8px; padding:10px 14px; "
        "margin:6px 0px; background:#fbfcfd;")
_TITLE = "font-weight:600; font-size:14px; color:#1f2328; margin-bottom:4px;"
_HINT = "color:#57606a; font-size:12px; margin-bottom:6px;"


def widgets_available() -> bool:
    """True if ``ipywidgets`` can be imported."""
    try:
        import ipywidgets  # noqa: F401
        return True
    except ImportError:
        return False


def _w():
    try:
        import ipywidgets as w
        return w
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError(
            "shapece.ui needs ipywidgets. Install it with:  pip install ipywidgets\n"
            "(On Google Colab it is already available.)") from exc


def _panel(title: str, hint: str, *children):
    w = _w()
    header = [w.HTML(f"<div style='{_TITLE}'>{title}</div>")]
    if hint:
        header.append(w.HTML(f"<div style='{_HINT}'>{hint}</div>"))
    return w.VBox(header + list(children),
                  layout=w.Layout(border="1px solid #d0d7de", border_radius="8px",
                                  padding="10px 14px", margin="6px 0px"))


# ---------------------------------------------------------------------------
# Sequence handling
# ---------------------------------------------------------------------------
def clean_sequence(text: str, to_rna: bool = True) -> str:
    """Normalize a pasted sequence: drop FASTA headers, whitespace, digits.

    Accepts a bare sequence or a FASTA record. Lowercase is upper-cased and, by
    default, ``T`` is converted to ``U``.

    Raises
    ------
    ValueError
        If characters outside ``ACGUTN`` remain after cleaning.
    """
    lines = [ln.strip() for ln in str(text).splitlines()]
    lines = [ln for ln in lines if ln and not ln.startswith(">")]
    seq = "".join(lines)
    seq = "".join(ch for ch in seq if not ch.isspace() and not ch.isdigit())
    seq = seq.upper()
    if to_rna:
        seq = seq.replace("T", "U")
    bad = sorted(set(seq) - set("ACGUTN"))
    if bad:
        raise ValueError(f"Sequence contains invalid characters: {bad}")
    return seq


class SequenceInput:
    """A text box for pasting an RNA sequence or FASTA record.

    Validates as you type and reports length / GC content.

    Examples
    --------
    >>> seq_box = SequenceInput()        # doctest: +SKIP
    >>> seq_box.display()                # doctest: +SKIP
    >>> seq_box.value                    # cleaned sequence  # doctest: +SKIP
    """

    def __init__(self, title: str = "RNA sequence",
                 hint: str = "Paste a sequence or FASTA record. "
                             "Headers, whitespace, digits and T→U are handled for you.",
                 placeholder: str = "GGGAUCGUCAC…  or  >my_rna\\nGGGAUC…"):
        w = _w()
        self._area = w.Textarea(
            placeholder=placeholder,
            layout=w.Layout(width="100%", height="120px"))
        self._status = w.HTML("")
        self._area.observe(self._on_change, names="value")
        self._panel = _panel(title, hint, self._area, self._status)
        self._clean = ""

    def _on_change(self, _=None):
        try:
            self._clean = clean_sequence(self._area.value)
        except ValueError as e:
            self._clean = ""
            self._status.value = f"<span style='color:#cf222e'>✗ {e}</span>"
            return
        if not self._clean:
            self._status.value = ""
            return
        gc = sum(self._clean.count(b) for b in "GC") / len(self._clean) * 100
        self._status.value = (
            f"<span style='color:#1a7f37'>✓ {len(self._clean)} nt, "
            f"{gc:.1f}% GC</span>")

    @property
    def value(self) -> str:
        """The cleaned, validated sequence ('' if empty or invalid)."""
        return self._clean

    def display(self):
        from IPython.display import display
        display(self._panel)
        return self


# ---------------------------------------------------------------------------
# File upload
# ---------------------------------------------------------------------------
def _in_colab() -> bool:
    try:
        import google.colab  # noqa: F401
        return True
    except ImportError:
        return False


def upload_files(dest: str = "uploads", accept: str = ".fsa") -> dict:
    """Show a file-picker button; save the chosen files and return {name: path}.

    On Google Colab this uses the native picker (``google.colab.files.upload``),
    which blocks until files are chosen. Elsewhere it falls back to an
    ``ipywidgets.FileUpload`` button -- in that case call :func:`collect_uploads`
    on the returned widget after picking files.

    Parameters
    ----------
    dest : str
        Directory to write the uploaded files into (created if needed).
    accept : str
        File-extension filter shown in the picker.
    """
    os.makedirs(dest, exist_ok=True)
    if _in_colab():                                     # pragma: no cover
        from google.colab import files as _files
        from IPython.display import display, HTML
        display(HTML(f"<div style='{_BOX}'><div style='{_TITLE}'>"
                     f"Upload your {accept} files</div>"
                     f"<div style='{_HINT}'>Select one or more files. "
                     f"They will be saved to <code>{dest}/</code>.</div></div>"))
        uploaded = _files.upload()
        paths = {}
        for name, data in uploaded.items():
            p = os.path.join(dest, name)
            with open(p, "wb") as fh:
                fh.write(data)
            paths[name] = p
        print(f"Saved {len(paths)} file(s) to {dest}/")
        return paths

    w = _w()
    fu = w.FileUpload(accept=accept, multiple=True,
                      description="Choose files",
                      layout=w.Layout(width="220px"))
    out = w.Output()
    fu._shapece_dest = dest          # noqa: SLF001 - stash for collect_uploads
    fu._shapece_paths = {}           # noqa: SLF001

    def _on(_change):
        with out:
            paths = collect_uploads(fu)
            print(f"Saved {len(paths)} file(s) to {dest}/")

    fu.observe(_on, names="value")
    from IPython.display import display
    display(_panel(f"Upload your {accept} files",
                   f"Select one or more files. They will be saved to {dest}/.",
                   fu, out))
    return fu


def collect_uploads(file_upload_widget) -> dict:
    """Write an ``ipywidgets.FileUpload`` payload to disk; return {name: path}.

    Handles both ipywidgets 7 (``value`` is a dict) and 8 (``value`` is a tuple).
    """
    dest = getattr(file_upload_widget, "_shapece_dest", "uploads")
    os.makedirs(dest, exist_ok=True)
    val = file_upload_widget.value
    items = []
    if isinstance(val, dict):                     # ipywidgets 7
        items = [(name, rec["content"]) for name, rec in val.items()]
    else:                                         # ipywidgets 8
        items = [(rec["name"], rec["content"]) for rec in val]
    paths = {}
    for name, content in items:
        p = os.path.join(dest, name)
        with open(p, "wb") as fh:
            fh.write(bytes(content))
        paths[name] = p
    file_upload_widget._shapece_paths = paths     # noqa: SLF001
    return paths


# ---------------------------------------------------------------------------
# Channel selection
# ---------------------------------------------------------------------------
class ChannelPicker:
    """Dropdowns for the reagent and size-standard channels of a `.fsa` file.

    The dropdown options are read from the file, so users pick from what actually
    exists instead of guessing a number. A short summary (length, max signal)
    accompanies each channel to make the choice obvious: the reagent channel has
    many sharp peaks; the size standard is an evenly spaced ladder.
    """

    def __init__(self, example_file: str,
                 title: str = "Which channel holds what?",
                 hint: str = "The reagent trace has many sharp peaks. "
                             "The size standard is an evenly spaced ladder. "
                             "Choose 'none' if your run has no size standard."):
        from .io import read_abif, list_data_channels, get_channel
        import numpy as np

        w = _w()
        entries = read_abif(example_file)
        chans = list_data_channels(entries)
        labels = []
        for c in chans:
            d = get_channel(entries, c)
            labels.append((f"DATA{c}   (n={len(d)}, max={int(np.max(d))})", c))

        default_data = 9 if 9 in chans else (2 if 2 in chans else chans[0])
        default_size = 205 if 205 in chans else None

        self.data = w.Dropdown(options=labels, value=default_data,
                               description="reagent:",
                               layout=w.Layout(width="360px"))
        self.size = w.Dropdown(options=[("none (no size standard)", None)] + labels,
                               value=default_size, description="size std:",
                               layout=w.Layout(width="360px"))
        self._panel = _panel(title, hint, self.data, self.size)

    def to_config(self) -> dict:
        return {"data_channel": self.data.value, "size_channel": self.size.value}

    def display(self):
        from IPython.display import display
        display(self._panel)
        return self


# ---------------------------------------------------------------------------
# Sample / pair assignment
# ---------------------------------------------------------------------------
class SamplePairer:
    """Assign each uploaded file a role, condition and replicate -- no typing.

    One row per file: a **role** dropdown ((+) reagent / (-) background / ladder /
    ignore), a **condition** box, and a **replicate** number. ``to_config()``
    matches each (+) with the (-) sharing its condition and replicate.
    """

    ROLES = [("(+) reagent", "plus"), ("(-) background", "minus"),
             ("ladder", "ladder"), ("ignore", "ignore")]

    def __init__(self, file_paths: dict, conditions=("condition 1", "condition 2"),
                 title: str = "Describe each file",
                 hint: str = "Set the role, condition and replicate for every file. "
                             "Each (+) is paired with the (-) that shares its "
                             "condition and replicate number."):
        w = _w()
        self.paths = dict(file_paths)
        self.rows = {}
        rows_ui = []
        for i, name in enumerate(self.paths):
            role = w.Dropdown(options=self.ROLES, value="plus" if i % 2 == 0 else "minus",
                              layout=w.Layout(width="150px"))
            cond = w.Combobox(options=list(conditions),
                              value=list(conditions)[0],
                              placeholder="condition",
                              layout=w.Layout(width="150px"))
            rep = w.BoundedIntText(value=(i // 2) + 1, min=1, max=99,
                                   layout=w.Layout(width="70px"))
            self.rows[name] = (role, cond, rep)
            rows_ui.append(w.HBox([
                w.HTML(f"<code style='font-size:12px'>{name}</code>",
                       layout=w.Layout(width="290px")),
                role, cond, rep]))
        header = w.HTML(
            "<div style='font-size:12px;color:#57606a'>"
            "<b style='display:inline-block;width:290px'>file</b>"
            "<b style='display:inline-block;width:150px'>role</b>"
            "<b style='display:inline-block;width:150px'>condition</b>"
            "<b>replicate</b></div>")
        self._panel = _panel(title, hint, header, w.VBox(rows_ui))

    def to_config(self) -> dict:
        """Return ``{"samples": {...}, "pairs": [...], "ladders": {...}}``."""
        samples, ladders = {}, {}
        plus, minus = {}, {}
        for name, (role, cond, rep) in self.rows.items():
            r = role.value
            if r == "ignore":
                continue
            if r == "ladder":
                ladders[name] = self.paths[name]
                continue
            key = (cond.value, int(rep.value))
            samples[name] = self.paths[name]
            (plus if r == "plus" else minus)[key] = name

        pairs, unmatched = [], []
        for key, pname in sorted(plus.items()):
            if key in minus:
                pairs.append((pname, minus[key]))
            else:
                unmatched.append(pname)
        return {"samples": samples, "pairs": pairs, "ladders": ladders,
                "unmatched_plus": unmatched}

    def display(self):
        from IPython.display import display
        display(self._panel)
        return self


# ---------------------------------------------------------------------------
# One form to rule them all
# ---------------------------------------------------------------------------
class ConfigPanel:
    """Combine sequence, channel and sample widgets into a single form.

    Examples
    --------
    >>> paths = upload_files()                    # doctest: +SKIP
    >>> panel = ConfigPanel(paths).display()      # doctest: +SKIP
    >>> CONFIG = panel.to_config()                # after filling the form
    """

    def __init__(self, file_paths: dict, conditions=("condition 1", "condition 2")):
        if not file_paths:
            raise ValueError("No files supplied -- run upload_files() first.")
        example = next(iter(file_paths.values()))
        self.sequence = SequenceInput()
        self.channels = ChannelPicker(example)
        self.samples = SamplePairer(file_paths, conditions)

    def display(self):
        from IPython.display import display
        display(self.sequence._panel, self.channels._panel, self.samples._panel)  # noqa: SLF001
        return self

    def to_config(self) -> dict:
        cfg = {"rna_sequence": self.sequence.value}
        cfg.update(self.channels.to_config())
        cfg.update(self.samples.to_config())
        cfg["reference"] = cfg["pairs"][0][0] if cfg["pairs"] else None
        cfg["reactivity_model"] = "area_difference"
        cfg["read_window"] = (None, None)
        return cfg

    def summary(self) -> str:
        c = self.to_config()
        lines = [
            f"sequence      : {len(c['rna_sequence'])} nt"
            if c["rna_sequence"] else "sequence      : (not set)",
            f"data channel  : DATA{c['data_channel']}",
            f"size standard : "
            + (f"DATA{c['size_channel']}" if c["size_channel"] else "none (DTW alignment)"),
            f"samples       : {len(c['samples'])}",
            f"pairs         : {len(c['pairs'])}",
        ]
        for p, m in c["pairs"]:
            lines.append(f"    {p}  (+)  ↔  {m}  (-)")
        if c["ladders"]:
            lines.append(f"ladders       : {', '.join(c['ladders'])}")
        if c["unmatched_plus"]:
            lines.append(f"⚠ unmatched (+) lanes: {', '.join(c['unmatched_plus'])}")
        return "\n".join(lines)
