#!/usr/bin/env python3
"""
make_braid.py -- draw open or closed braid diagrams, export SVG, and emit LaTeX.

A braid is given by its *word* in Artin generators (Tietze notation): a list of
non-zero integers, where ``i`` means the positive generator sigma_i (strand i
crosses over strand i+1) and ``-i`` its inverse.  Strands are 1-indexed in the
word, so ``[1, 1, -2, 1, -2]`` acts on at least 3 strands.

Two diagram styles are produced:

    * closed  -- the annular closure, strands drawn as concentric arcs.
    * open    -- the ordinary rectangular braid diagram.

A single ``--flow`` direction controls both: ``top-down``, ``bottom-up``,
``left-right`` or ``right-left``.  For an open braid it sets which way the word
runs across the page; for a closed braid it sets the winding sense and the
default angle of the start/end seam.  The seam (and the open braid's start/end
edges) can be marked with dashed lines via ``--mark-ends``.

Rendering uses matplotlib, so the script runs under a plain ``python3`` *or*
``sage``.  SageMath is only consulted (optionally) by ``--info``.

Examples
--------
    # Closed braid -> braid.svg
    python3 make_braid.py "1 1 -2 1 -2" -n 3 -o closed.svg

    # Open braid running left-to-right, mark the start/end edges
    python3 make_braid.py "1 -2 1 -2" --open --flow left-right --mark-ends

    # Closed braid, seam at 3 o'clock, dashed seam line
    python3 make_braid.py "1 1 1" -n 2 --start-angle 0 --mark-ends

    # Print the braid word's LaTeX (open form + closure) to stdout
    python3 make_braid.py "1 1 -2 1 -2" --latex

    # Launch the interactive GUI -- the LaTeX code is shown live in its log
    python3 make_braid.py --gui

Run ``python3 make_braid.py --help`` for the full parameter list.
"""

import argparse
import colorsys
import os
import re
import sys
from collections import namedtuple
from dataclasses import dataclass
from math import cos, pi, radians, sin

from matplotlib.colors import to_rgb
from matplotlib.figure import Figure


# --------------------------------------------------------------------------- #
# Braid combinatorics
# --------------------------------------------------------------------------- #
def parse_word(text):
    """Parse a braid word: "1 1 -2", "1,1,-2" or "[1, 1, -2]" -> [1, 1, -2]."""
    if isinstance(text, (list, tuple)):
        return [int(g) for g in text]
    text = text.strip().strip("[](){}").strip()
    if not text:
        return []
    word = []
    for p in re.split(r"[,\s]+", text):
        if not p:
            continue
        g = int(p)
        if g == 0:
            raise ValueError("Braid generators must be non-zero (got 0).")
        word.append(g)
    return word


def infer_strands(word):
    """Minimum number of strands implied by a braid word."""
    return max((abs(g) for g in word), default=0) + 1


def validate(word, n):
    if n < 1:
        raise ValueError("Number of strands must be >= 1.")
    for g in word:
        i = abs(g) - 1
        if i < 0 or i >= n - 1:
            raise ValueError(
                f"Generator {g} out of range for {n} strands "
                f"(need 1 <= |g| <= {n - 1})."
            )


def closed_braid_components(word, n):
    """Components of the standard closure: list of lists of 0-based labels."""
    labels = list(range(n))
    for g in word:
        i = abs(g) - 1
        labels[i], labels[i + 1] = labels[i + 1], labels[i]

    end_pos = [None] * n
    for pos, lab in enumerate(labels):
        end_pos[lab] = pos

    seen = [False] * n
    components = []
    for lab in range(n):
        if not seen[lab]:
            comp = []
            cur = lab
            while not seen[cur]:
                seen[cur] = True
                comp.append(cur)
                cur = end_pos[cur]
            components.append(comp)
    return components


# --------------------------------------------------------------------------- #
# Colours
# --------------------------------------------------------------------------- #
MARKER_COLOR = (0.45, 0.45, 0.45)   # dashed start/end indicators


def distinct_colors(k):
    if k <= 1:
        return [(0.0, 0.0, 0.0)]
    return [colorsys.hsv_to_rgb(j / k, 0.75, 0.85) for j in range(k)]


def resolve_color_mode(mode, closed):
    if mode == "auto":
        return "component" if closed else "strand"
    return mode


def label_colors(word, n, mode, closed, fg_color):
    """A colour per strand label (0..n-1)."""
    mode = resolve_color_mode(mode, closed)
    if mode == "single":
        return [to_rgb(fg_color)] * n
    if mode == "strand":
        return distinct_colors(n)
    comps = closed_braid_components(word, n)
    comp_colors = distinct_colors(len(comps))
    comp_id = [0] * n
    for j, comp in enumerate(comps):
        for lab in comp:
            comp_id[lab] = j
    return [comp_colors[comp_id[lab]] for lab in range(n)]


# --------------------------------------------------------------------------- #
# Style / parameters
# --------------------------------------------------------------------------- #
# flow -> (default closed seam angle in degrees, default winding: +1 CCW / -1 CW)
FLOW_DEFAULTS = {
    "top-down":   (90.0, -1),
    "bottom-up":  (270.0, +1),
    "left-right": (180.0, -1),
    "right-left": (0.0, +1),
}
FLOWS = tuple(FLOW_DEFAULTS)


@dataclass
class BraidStyle:
    # topology / direction
    closed: bool = True
    flow: str = "top-down"                  # top-down|bottom-up|left-right|right-left
    color_mode: str = "auto"                # auto|component|strand|single
    # closed seam / winding
    start_angle: float = None               # None -> flow default (degrees)
    winding: str = "auto"                   # auto|cw|ccw
    mark_ends: bool = False                  # dashed start/end indicators
    # line
    thickness: float = 4.0                   # matplotlib linewidth
    gap: float = 0.1                         # under-strand break, fraction of a cell
    samples: int = 40                        # samples per crossing arc
    # closed-braid geometry
    r0: float = 0.5                          # innermost radius
    dr: float = 0.25                         # radial spacing between strands
    # open-braid geometry
    spacing: float = 1.0                     # distance between strands
    level_height: float = 1.0                # size of one crossing
    # canvas
    figsize: float = 6.0
    margin: float = 0.08
    transparent: bool = False
    bg_color: str = "white"
    fg_color: str = "black"                  # used when color_mode == 'single'


def smoothstep(u):
    return u * u * (3.0 - 2.0 * u)


def resolve_start_dir(st):
    """Return (start angle in radians, winding direction +/-1) for closed braids."""
    default_angle, default_dir = FLOW_DEFAULTS[st.flow]
    angle = default_angle if st.start_angle is None else st.start_angle
    if st.winding == "cw":
        d = -1
    elif st.winding == "ccw":
        d = +1
    else:
        d = default_dir
    return radians(angle), d


# --------------------------------------------------------------------------- #
# Geometry -> paths (shared by the matplotlib and LaTeX back-ends)
# --------------------------------------------------------------------------- #
BraidPath = namedtuple("BraidPath", "pts color dashed width")


def _open_map(flow):
    """Map (pos across strands, prog along braid) -> (x, y) for a given flow."""
    if flow == "top-down":
        return lambda pos, prog: (pos, -prog)
    if flow == "bottom-up":
        return lambda pos, prog: (pos, prog)
    if flow == "left-right":
        return lambda pos, prog: (prog, -pos)
    if flow == "right-left":
        return lambda pos, prog: (-prog, -pos)
    raise ValueError(f"Unknown flow: {flow!r}")


def open_paths(word, n, st, colors):
    lw, sp, lh = st.thickness, st.spacing, st.level_height
    m = len(word)
    mp = _open_map(st.flow)
    paths = []

    def seg(pa, pb, k, u0=0.0, u1=1.0):
        out = []
        xa, xb = pa * sp, pb * sp
        for t in range(st.samples + 1):
            u = u0 + (u1 - u0) * t / st.samples
            pos = xa + (xb - xa) * smoothstep(u)
            out.append(mp(pos, (k + u) * lh))
        return out

    if m == 0:
        for p in range(n):
            paths.append(BraidPath([mp(p * sp, 0.0), mp(p * sp, lh)],
                                   colors[p], False, lw))
    else:
        labels = list(range(n))
        for k, g in enumerate(word):
            i = abs(g) - 1
            for p in range(n):
                if p != i and p != i + 1:
                    paths.append(BraidPath(seg(p, p, k), colors[labels[p]],
                                           False, lw))
            over_start, under_start = (i, i + 1) if g > 0 else (i + 1, i)

            def swapped(p):
                return i + 1 if p == i else (i if p == i + 1 else p)

            over_end, under_end = swapped(over_start), swapped(under_start)
            ol, ul = labels[over_start], labels[under_start]
            paths.append(BraidPath(seg(under_start, under_end, k, 0.0,
                                       0.5 - st.gap), colors[ul], False, lw))
            paths.append(BraidPath(seg(under_start, under_end, k, 0.5 + st.gap,
                                       1.0), colors[ul], False, lw))
            paths.append(BraidPath(seg(over_start, over_end, k), colors[ol],
                                   False, lw))
            labels[i], labels[i + 1] = labels[i + 1], labels[i]

    if st.mark_ends:
        total = (m if m > 0 else 1) * lh
        pos_lo, pos_hi = -0.5 * sp, (n - 1 + 0.5) * sp
        mw = max(1.0, lw * 0.55)
        paths.append(BraidPath([mp(pos_lo, 0.0), mp(pos_hi, 0.0)],
                               MARKER_COLOR, True, mw))
        paths.append(BraidPath([mp(pos_lo, total), mp(pos_hi, total)],
                               MARKER_COLOR, True, mw))
    return paths


def closed_paths(word, n, st, colors):
    lw = st.thickness
    radii = [st.r0 + j * st.dr for j in range(n)]
    m = len(word)
    theta0, d = resolve_start_dir(st)
    paths = []

    def polar(r, theta):
        return (r * cos(theta), r * sin(theta))

    if m == 0:
        steps = max(120, st.samples * 8)
        for p in range(n):
            r = radii[p]
            paths.append(BraidPath(
                [polar(r, theta0 + d * 2 * pi * t / steps)
                 for t in range(steps + 1)], colors[p], False, lw))
    else:
        def arc(p, q, k, u0=0.0, u1=1.0):
            out = []
            for t in range(st.samples + 1):
                u = u0 + (u1 - u0) * t / st.samples
                theta = theta0 + d * 2 * pi * (k + u) / m
                r = radii[p] + (radii[q] - radii[p]) * smoothstep(u)
                out.append(polar(r, theta))
            return out

        labels = list(range(n))
        for k, g in enumerate(word):
            i = abs(g) - 1
            for p in range(n):
                if p != i and p != i + 1:
                    paths.append(BraidPath(arc(p, p, k), colors[labels[p]],
                                           False, lw))
            over_start, under_start = (i, i + 1) if g > 0 else (i + 1, i)

            def swapped(p):
                return i + 1 if p == i else (i if p == i + 1 else p)

            over_end, under_end = swapped(over_start), swapped(under_start)
            ol, ul = labels[over_start], labels[under_start]
            paths.append(BraidPath(arc(under_start, under_end, k, 0.0,
                                       0.5 - st.gap), colors[ul], False, lw))
            paths.append(BraidPath(arc(under_start, under_end, k, 0.5 + st.gap,
                                       1.0), colors[ul], False, lw))
            paths.append(BraidPath(arc(over_start, over_end, k), colors[ol],
                                   False, lw))
            labels[i], labels[i + 1] = labels[i + 1], labels[i]

    if st.mark_ends:
        r_in = max(0.05, st.r0 - 0.6 * st.dr)
        r_out = radii[-1] + 0.6 * st.dr
        paths.append(BraidPath([polar(r_in, theta0), polar(r_out, theta0)],
                               MARKER_COLOR, True, max(1.0, lw * 0.55)))
    return paths


def braid_paths(word, n, st):
    validate(word, n)
    colors = label_colors(word, n, st.color_mode, st.closed, st.fg_color)
    if st.closed:
        return closed_paths(word, n, st, colors)
    return open_paths(word, n, st, colors)


# --------------------------------------------------------------------------- #
# matplotlib back-end
# --------------------------------------------------------------------------- #
def draw_into_figure(fig, word, n, st):
    fig.clear()
    ax = fig.add_subplot(111)
    for p in braid_paths(word, n, st):
        xs = [a[0] for a in p.pts]
        ys = [a[1] for a in p.pts]
        ax.plot(xs, ys, color=p.color, linewidth=p.width,
                linestyle="--" if p.dashed else "-",
                solid_capstyle="round", solid_joinstyle="round",
                dash_capstyle="round")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.margins(st.margin)
    face = "none" if st.transparent else st.bg_color
    fig.patch.set_alpha(0.0 if st.transparent else 1.0)
    if not st.transparent:
        fig.patch.set_facecolor(st.bg_color)
    ax.set_facecolor(face)
    return ax


def build_figure(word, n, st):
    fig = Figure(figsize=(st.figsize, st.figsize))
    draw_into_figure(fig, word, n, st)
    return fig


def save_figure(fig, path, st):
    fig.savefig(path, bbox_inches="tight", pad_inches=0.1,
                transparent=st.transparent,
                facecolor=("none" if st.transparent else st.bg_color))


# --------------------------------------------------------------------------- #
# LaTeX for the braid word (for display / copying, not a drawing)
# --------------------------------------------------------------------------- #
def word_to_tex(word):
    """Artin word as LaTeX, combining consecutive equal generators as powers.

    e.g. [1, 1, -2, 1, -2] -> \\sigma_{1}^{2}\\sigma_{2}^{-1}\\sigma_{1}\\sigma_{2}^{-1}
    """
    if not word:
        return "1"
    out = []
    i, m = 0, len(word)
    while i < m:
        g = word[i]
        j = i
        while j < m and word[j] == g:      # run of the same signed generator
            j += 1
        exp = (j - i) if g > 0 else -(j - i)
        if exp == 1:
            out.append(f"\\sigma_{{{abs(g)}}}")
        else:
            out.append(f"\\sigma_{{{abs(g)}}}^{{{exp}}}")
        i = j
    return "".join(out)


def latex_open(word):
    """LaTeX for the open braid word (the element of the braid group)."""
    return word_to_tex(word)


def latex_closed(word):
    """LaTeX for the closed braid (the closure of the word)."""
    return f"\\widehat{{{word_to_tex(word)}}}"


def latex_code(word):
    """Both forms as a small block of LaTeX, ready to paste."""
    return (f"% open braid word\n{latex_open(word)}\n"
            f"% closed braid (closure)\n{latex_closed(word)}")


# --------------------------------------------------------------------------- #
# Optional: link invariants and topology via SageMath / KnotInfo
# --------------------------------------------------------------------------- #
def _sage_link(word, n):
    """Build a Sage ``Link`` for the braid closure. Requires Sage; may raise."""
    from sage.all import BraidGroup, Link
    B = BraidGroup(n)
    b = B(word) if word else B.one()
    return Link(b)


def _identify_in_process(word, n):
    """KnotInfo lookup using an importable Sage (i.e. running under `sage`)."""
    try:
        L = _sage_link(word, n)
    except Exception as exc:
        return f"could not build link ({type(exc).__name__}: {exc})"
    try:
        return str(L.get_knotinfo())
    except Exception as exc:
        return f"not found in KnotInfo ({type(exc).__name__}: {exc})"


def _identify_via_sage_cli(word, n):
    """KnotInfo lookup by shelling out to the ``sage`` executable.

    Lets the GUI (which must run under the system Python for a working Tk)
    still identify topology, by delegating the Sage part to a subprocess.
    """
    import shutil
    import subprocess

    sage = shutil.which("sage")
    if not sage:
        return ("SageMath not found — install Sage (or run under `sage`) "
                "to identify topology.")
    script = (
        "from sage.all import BraidGroup, Link\n"
        f"w = {list(word)}\n"
        f"B = BraidGroup({n})\n"
        "b = B(w) if w else B.one()\n"
        "L = Link(b)\n"
        "try:\n"
        "    print('OK:' + str(L.get_knotinfo()))\n"
        "except Exception as e:\n"
        "    print('ERR:%s: %s' % (type(e).__name__, e))\n"
    )
    try:
        proc = subprocess.run([sage, "-python", "-c", script],
                              capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return "timed out waiting for `sage` KnotInfo lookup."
    except Exception as exc:
        return f"failed to run `sage` ({type(exc).__name__}: {exc})"

    for line in proc.stdout.splitlines():
        if line.startswith("OK:"):
            return line[3:]
        if line.startswith("ERR:"):
            return f"not found in KnotInfo ({line[4:].strip()})"
    tail = (proc.stderr.strip() or proc.stdout.strip()
            or f"exit code {proc.returncode}")
    return f"could not identify ({tail.splitlines()[-1] if tail else 'no output'})"


def identify_closure(word, n):
    """Name the topology of the braid closure via Sage's KnotInfo lookup.

    Returns a human-readable string.  If Sage is importable in this process
    (script run under ``sage``) it is used directly; otherwise the ``sage``
    executable is invoked as a subprocess so the plain-``python3`` GUI can
    still identify topology.  ``get_knotinfo()`` raises when the closure is
    not in the KnotInfo database or is not uniquely identifiable; that is
    caught and reported instead of propagating.
    """
    try:
        from sage.all import BraidGroup, Link  # noqa: F401
    except Exception:
        return _identify_via_sage_cli(word, n)
    return _identify_in_process(word, n)


def print_link_info(word, n):
    print(f"[info] strands           : {n}")
    print(f"[info] word              : {word}")
    print(f"[info] closure components: {len(closed_braid_components(word, n))}")
    try:
        from sage.all import BraidGroup, Link  # noqa: F401
    except Exception:
        print("[info] SageMath not available; run under `sage` for invariants.",
              file=sys.stderr)
        return
    try:
        L = _sage_link(word, n)
    except Exception as exc:
        print(f"[info] (could not build link: {exc})")
        return
    try:
        print(f"[info] Jones polynomial  : {L.jones_polynomial()}")
    except Exception as exc:
        print(f"[info] Jones polynomial  : (unavailable: {exc})")
    try:
        print(f"[info] KnotInfo          : {L.get_knotinfo()}")
    except Exception as exc:
        print(f"[info] KnotInfo          : not found "
              f"({type(exc).__name__}: {exc})")


# --------------------------------------------------------------------------- #
# GUI
# --------------------------------------------------------------------------- #
def run_gui(word_text="1 1 -2 1 -2", strands=None, style=None):
    import threading
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    try:
        from matplotlib.backends.backend_tkagg import (
            FigureCanvasTkAgg, NavigationToolbar2Tk)
    except Exception as exc:
        raise RuntimeError(
            "the Tk GUI backend failed to load under this Python "
            f"({type(exc).__name__}: {exc}).\n"
            "This is common with `sage -python` on macOS. Launch the GUI "
            "with your system Python instead:\n"
            "    python3 make_braid.py\n"
            "The 'Identify (KnotInfo)' button still calls `sage` in the "
            "background, so topology lookup keeps working.") from exc

    st0 = style or BraidStyle()

    root = tk.Tk()
    root.title("make_braid — braid diagram builder")
    root.geometry("1340x840")

    # --- preview on the LEFT ------------------------------------------- #
    canvas_frame = ttk.Frame(root, padding=6)
    canvas_frame.pack(side="left", fill="both", expand=True)
    fig = Figure(figsize=(st0.figsize, st0.figsize))
    canvas = FigureCanvasTkAgg(fig, master=canvas_frame)
    canvas.get_tk_widget().pack(fill="both", expand=True)
    NavigationToolbar2Tk(canvas, canvas_frame)

    # --- controls on the RIGHT (scrollable) ---------------------------- #
    right = ttk.Frame(root)
    right.pack(side="right", fill="y")
    scan = tk.Canvas(right, width=580, highlightthickness=0)
    vsb = ttk.Scrollbar(right, orient="vertical", command=scan.yview)
    controls = ttk.Frame(scan, padding=10)
    controls.bind("<Configure>",
                  lambda e: scan.configure(scrollregion=scan.bbox("all")))
    scan.create_window((0, 0), window=controls, anchor="nw")
    scan.configure(yscrollcommand=vsb.set)
    scan.pack(side="left", fill="y", expand=False)
    vsb.pack(side="right", fill="y")
    scan.bind_all("<MouseWheel>",
                  lambda e: scan.yview_scroll(-1 if e.delta > 0 else 1, "units"))

    v = {
        "word": tk.StringVar(value=word_text),
        "strands": tk.StringVar(value="" if strands is None else str(strands)),
        "closed": tk.BooleanVar(value=st0.closed),
        "flow": tk.StringVar(value=st0.flow),
        "color_mode": tk.StringVar(value=st0.color_mode),
        "start_angle": tk.StringVar(
            value="" if st0.start_angle is None else str(st0.start_angle)),
        "winding": tk.StringVar(value=st0.winding),
        "mark_ends": tk.BooleanVar(value=st0.mark_ends),
        "thickness": tk.DoubleVar(value=st0.thickness),
        "gap": tk.DoubleVar(value=st0.gap),
        "samples": tk.IntVar(value=st0.samples),
        "r0": tk.DoubleVar(value=st0.r0),
        "dr": tk.DoubleVar(value=st0.dr),
        "spacing": tk.DoubleVar(value=st0.spacing),
        "level_height": tk.DoubleVar(value=st0.level_height),
        "figsize": tk.DoubleVar(value=st0.figsize),
        "margin": tk.DoubleVar(value=st0.margin),
        "transparent": tk.BooleanVar(value=st0.transparent),
        "bg_color": tk.StringVar(value=st0.bg_color),
        "fg_color": tk.StringVar(value=st0.fg_color),
        "output": tk.StringVar(value="braid.svg"),
    }

    status = ttk.Label(controls, text="", foreground="#357", wraplength=540)
    topo_var = tk.StringVar(value="")
    save_dir = {"path": ""}          # remembers the last chosen save folder

    # LaTeX log: created early so render() can update it; gridded later.
    latex_box = tk.Text(controls, height=6, width=60, wrap="word",
                        relief="solid", borderwidth=1,
                        font=("Menlo", 10))
    latex_state = {"open": "", "closed": ""}

    def _set_latex(word):
        latex_state["open"] = latex_open(word)
        latex_state["closed"] = latex_closed(word)
        text = (f"open braid word:\n{latex_state['open']}\n\n"
                f"closed braid (closure):\n{latex_state['closed']}")
        latex_box.config(state="normal")
        latex_box.delete("1.0", "end")
        latex_box.insert("1.0", text)
        latex_box.config(state="disabled")

    def copy_latex():
        root.clipboard_clear()
        root.clipboard_append(f"{latex_state['open']}\n{latex_state['closed']}")
        status.config(text="LaTeX copied to clipboard.", foreground="#171")

    def identify():
        try:
            word, n, _ = collect()
        except Exception as exc:
            status.config(text=f"Error: {exc}", foreground="#a00")
            return
        status.config(text="Identifying topology via KnotInfo… "
                           "(Sage may take a while)", foreground="#357")
        topo_var.set("…")

        def work():                      # run Sage off the UI thread
            result = identify_closure(word, n)
            root.after(0, lambda: (
                topo_var.set(result),
                status.config(text="Topology lookup done.", foreground="#171")))

        threading.Thread(target=work, daemon=True).start()

    def collect():
        word = parse_word(v["word"].get())
        s = v["strands"].get().strip()
        n = int(s) if s else infer_strands(word)
        sa = v["start_angle"].get().strip()
        st = BraidStyle(
            closed=v["closed"].get(),
            flow=v["flow"].get(),
            color_mode=v["color_mode"].get(),
            start_angle=float(sa) if sa else None,
            winding=v["winding"].get(),
            mark_ends=v["mark_ends"].get(),
            thickness=float(v["thickness"].get()),
            gap=float(v["gap"].get()),
            samples=int(v["samples"].get()),
            r0=float(v["r0"].get()),
            dr=float(v["dr"].get()),
            spacing=float(v["spacing"].get()),
            level_height=float(v["level_height"].get()),
            figsize=float(v["figsize"].get()),
            margin=float(v["margin"].get()),
            transparent=v["transparent"].get(),
            bg_color=v["bg_color"].get().strip() or "white",
            fg_color=v["fg_color"].get().strip() or "black",
        )
        return word, n, st

    def render(*_):
        try:
            word, n, st = collect()
            draw_into_figure(fig, word, n, st)
            canvas.draw()
            _set_latex(word)
            kind = "closed" if st.closed else f"open/{st.flow}"
            status.config(text=f"OK — {kind}, {n} strands, "
                               f"{len(word)} crossings.", foreground="#357")
        except Exception as exc:
            status.config(text=f"Error: {exc}", foreground="#a00")

    def save():
        try:
            word, n, st = collect()
        except Exception as exc:
            messagebox.showerror("make_braid", f"Cannot render: {exc}")
            return
        # keep the dialog's directory and file name separate, so the folder
        # never ends up inside the file-name field.
        name = os.path.basename(v["output"].get().strip()) or "braid.svg"
        path = filedialog.asksaveasfilename(
            defaultextension=".svg",
            initialdir=save_dir["path"] or None,
            initialfile=name,
            filetypes=[("SVG", "*.svg"), ("PDF", "*.pdf"),
                       ("PNG", "*.png"), ("All files", "*.*")])
        if not path:
            return
        try:
            save_figure(build_figure(word, n, st), path, st)
            save_dir["path"] = os.path.dirname(path)
            v["output"].set(os.path.basename(path))   # store bare file name
            status.config(text=f"Saved: {path}", foreground="#171")
        except Exception as exc:
            messagebox.showerror("make_braid", f"Save failed: {exc}")

    # --- layout helpers: two columns of (label, widget) pairs ----------- #
    lay = {"row": 0, "side": 0}          # side 0 = left pair, 1 = right pair

    def _flush():                        # start the next block on a fresh row
        if lay["side"] == 1:
            lay["row"] += 1
            lay["side"] = 0

    def _cell():                         # (row, base col) for the next field
        col = 0 if lay["side"] == 0 else 2
        r = lay["row"]
        if lay["side"] == 0:
            lay["side"] = 1
        else:
            lay["side"] = 0
            lay["row"] += 1
        return r, col

    def section(title):
        _flush()
        ttk.Separator(controls).grid(row=lay["row"], column=0, columnspan=4,
                                     sticky="ew", pady=(9, 3))
        lay["row"] += 1
        ttk.Label(controls, text=title, font=("", 11, "bold")).grid(
            row=lay["row"], column=0, columnspan=4, sticky="w")
        lay["row"] += 1

    def _apply_entry(entry):
        # re-render when the user commits a text field (Enter or click away)
        entry.bind("<Return>", render)
        entry.bind("<FocusOut>", render)

    def add_entry(label, var, width=11):
        r, c = _cell()
        ttk.Label(controls, text=label).grid(row=r, column=c, sticky="w",
                                             pady=1, padx=(0, 4))
        e = ttk.Entry(controls, textvariable=var, width=width)
        e.grid(row=r, column=c + 1, sticky="w", pady=1, padx=(0, 14))
        _apply_entry(e)

    def add_combo(label, var, values, width=10):
        r, c = _cell()
        ttk.Label(controls, text=label).grid(row=r, column=c, sticky="w",
                                             pady=1, padx=(0, 4))
        ttk.Combobox(controls, textvariable=var, values=values, width=width,
                     state="readonly").grid(row=r, column=c + 1, sticky="w",
                                            pady=1, padx=(0, 14))

    def add_check(label, var):
        r, c = _cell()
        ttk.Checkbutton(controls, text=label, variable=var).grid(
            row=r, column=c, columnspan=2, sticky="w", pady=1)

    def add_wide(label, var, width=40):
        _flush()
        ttk.Label(controls, text=label).grid(row=lay["row"], column=0,
                                             sticky="w", pady=1, padx=(0, 4))
        e = ttk.Entry(controls, textvariable=var, width=width)
        e.grid(row=lay["row"], column=1, columnspan=3, sticky="ew", pady=1)
        _apply_entry(e)
        lay["row"] += 1

    section("Braid")
    add_wide("Word", v["word"])
    add_entry("Strands", v["strands"])
    add_check("Closed (annular)", v["closed"])
    add_combo("Flow", v["flow"], list(FLOWS))
    add_combo("Colour", v["color_mode"],
              ["auto", "component", "strand", "single"])

    section("Start / end")
    add_entry("Start°", v["start_angle"])
    add_combo("Winding", v["winding"], ["auto", "cw", "ccw"])
    add_check("Mark ends (dashed)", v["mark_ends"])

    section("Line")
    add_entry("Thickness", v["thickness"])
    add_entry("Under-gap", v["gap"])
    add_entry("Samples", v["samples"])

    section("Geometry")
    add_entry("r0 (closed)", v["r0"])
    add_entry("dr (closed)", v["dr"])
    add_entry("Spacing (open)", v["spacing"])
    add_entry("Level ht (open)", v["level_height"])

    section("Canvas")
    add_entry("Fig size", v["figsize"])
    add_entry("Margin", v["margin"])
    add_check("Transparent", v["transparent"])
    add_entry("BG colour", v["bg_color"])
    add_entry("FG colour", v["fg_color"])

    section("Output")
    add_wide("File", v["output"])

    _flush()
    btns = ttk.Frame(controls)
    btns.grid(row=lay["row"], column=0, columnspan=4, sticky="ew", pady=(12, 4))
    lay["row"] += 1
    ttk.Button(btns, text="Render", command=render).pack(side="left", padx=2)
    ttk.Button(btns, text="Save…", command=save).pack(side="left", padx=2)

    status.grid(row=lay["row"], column=0, columnspan=4, sticky="w", pady=(6, 0))
    lay["row"] += 1

    section("LaTeX (braid word)")
    latex_box.grid(row=lay["row"], column=0, columnspan=4, sticky="ew", pady=2)
    lay["row"] += 1
    ttk.Button(controls, text="Copy LaTeX", command=copy_latex).grid(
        row=lay["row"], column=0, columnspan=4, sticky="w", pady=(2, 6))
    lay["row"] += 1

    section("Topology (closure)")
    ttk.Button(controls, text="Identify (KnotInfo)", command=identify).grid(
        row=lay["row"], column=0, columnspan=2, sticky="w", pady=2)
    lay["row"] += 1
    ttk.Label(controls, textvariable=topo_var, wraplength=540,
              foreground="#334").grid(row=lay["row"], column=0, columnspan=4,
                                      sticky="w", pady=(0, 8))
    lay["row"] += 1

    for key in ("closed", "flow", "color_mode", "winding", "mark_ends",
                "transparent"):
        v[key].trace_add("write", render)

    render()
    root.mainloop()


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_arg_parser():
    p = argparse.ArgumentParser(
        prog="make_braid.py",
        description="Draw open or closed braid diagrams; export SVG and LaTeX.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("word", nargs="?",
                   help='Braid word in Tietze notation, e.g. "1 1 -2 1 -2". '
                        "If omitted, the GUI opens.")
    p.add_argument("-n", "--strands", type=int, default=None,
                   help="Number of strands (default: inferred from the word).")
    p.add_argument("-o", "--output", default="braid.svg",
                   help="Output image file (extension picks the format).")

    topo = p.add_argument_group("topology / direction")
    mode = topo.add_mutually_exclusive_group()
    mode.add_argument("--closed", dest="closed", action="store_true",
                      help="Draw the annular closure (default).")
    mode.add_argument("--open", dest="closed", action="store_false",
                      help="Draw the open (rectangular) braid.")
    p.set_defaults(closed=True)
    topo.add_argument("--flow", choices=list(FLOWS), default="top-down",
                      help="Direction the braid word runs (also sets closed "
                           "winding sense and default seam angle).")
    topo.add_argument("--start-angle", type=float, default=None,
                      help="Closed: seam angle in degrees "
                           "(0=right, 90=top). Default: from --flow.")
    topo.add_argument("--winding", choices=["auto", "cw", "ccw"], default="auto",
                      help="Closed: winding sense. Default: from --flow.")
    topo.add_argument("--mark-ends", action="store_true",
                      help="Mark the start/end seam (closed) or edges (open) "
                           "with dashed lines.")
    topo.add_argument("--color-mode",
                      choices=["auto", "component", "strand", "single"],
                      default="auto",
                      help="Colour by closure component, by strand, or single.")

    style = p.add_argument_group("style")
    style.add_argument("--thickness", type=float, default=4.0)
    style.add_argument("--gap", type=float, default=0.1,
                       help="Under-strand break size (fraction of a cell).")
    style.add_argument("--samples", type=int, default=40,
                       help="Samples per crossing arc (smoothness).")
    style.add_argument("--r0", type=float, default=0.5,
                       help="Closed: innermost radius.")
    style.add_argument("--dr", type=float, default=0.25,
                       help="Closed: radial spacing between strands.")
    style.add_argument("--spacing", type=float, default=1.0,
                       help="Open: distance between strands.")
    style.add_argument("--level-height", type=float, default=1.0,
                       help="Open: height of one crossing.")
    style.add_argument("--figsize", type=float, default=6.0)
    style.add_argument("--margin", type=float, default=0.08)
    style.add_argument("--transparent", action="store_true",
                       help="Transparent background.")
    style.add_argument("--bg-color", default="white")
    style.add_argument("--fg-color", default="black",
                       help="Colour used when --color-mode=single.")

    out = p.add_argument_group("other outputs")
    out.add_argument("--latex", "--tex", dest="latex", action="store_true",
                     help="Print LaTeX for the braid word (open form + "
                          "closure) to stdout, instead of drawing.")
    out.add_argument("--info", action="store_true",
                     help="Print closure invariants and KnotInfo topology "
                          "(needs Sage; run under `sage`).")
    p.add_argument("--gui", action="store_true", help="Launch the GUI.")
    return p


def style_from_args(args):
    return BraidStyle(
        closed=args.closed, flow=args.flow, color_mode=args.color_mode,
        start_angle=args.start_angle, winding=args.winding,
        mark_ends=args.mark_ends, thickness=args.thickness, gap=args.gap,
        samples=args.samples, r0=args.r0, dr=args.dr, spacing=args.spacing,
        level_height=args.level_height, figsize=args.figsize, margin=args.margin,
        transparent=args.transparent, bg_color=args.bg_color,
        fg_color=args.fg_color)


def main(argv=None):
    args = build_arg_parser().parse_args(argv)

    if args.gui or args.word is None:
        try:
            run_gui(word_text=args.word or "1 1 -2 1 -2", strands=args.strands,
                    style=style_from_args(args))
        except RuntimeError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1
        return 0

    try:
        word = parse_word(args.word)
        n = args.strands if args.strands is not None else infer_strands(word)
        validate(word, n)
        st = style_from_args(args)

        if args.latex:
            print(latex_code(word))
        else:
            save_figure(build_figure(word, n, st), args.output, st)
            kind = "closed" if st.closed else f"open ({st.flow})"
            print(f"Wrote {args.output}  ({kind} braid, {n} strands, "
                  f"{len(word)} crossings)")
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.info:
        print_link_info(word, n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
