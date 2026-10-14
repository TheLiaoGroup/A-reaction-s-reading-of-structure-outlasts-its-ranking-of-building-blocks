"""Extended Data Figs 3-8 and the shared figure style; plots the panel tables only."""

from __future__ import annotations

import argparse

import math

from pathlib import Path

import re

import numpy as np

import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

from matplotlib.lines import Line2D

from matplotlib.patches import Rectangle, PathPatch

from matplotlib.path import Path as MPath

from matplotlib.text import Annotation, Text

from rdkit import Chem, RDLogger

from rdkit.Chem import Draw, rdDepictor

from rdkit.Chem.Draw import rdMolDraw2D

from rdkit.Geometry import Point3D

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ROOT / "output/presentation/panel_data"
OUT = ROOT / "output/presentation/figures"

STRUCTURE_FEATURES = {"halogen next to ring N": ("Halogen beside a ring N", "Brc1ccccn1"),
                      "N–H": ("N–H in the aryl halide", "Brc1cccc2[nH]ncc12")}

STRUCTURE_SETS = {"first condition set": ("o", 4.6, True), "second condition set": ("D", 3.4, False)}

BLUE, ORANGE, TEAL, RED = "#28658A", "#B86C27", "#147A71", "#BD4939"

INK, GREY, LIGHT, MUTED = "#202932", "#929BA3", "#EEF1F3", "#59636C"

NAMES = {"Acylation": "Amide coupling", "Heck": "Heck", "Suzuki": "Suzuki–Miyaura",
         "Ullmann": "Ullmann ether", "Buchwald-Hartwig": "Buchwald–Hartwig",
         "Kumada": "Kumada", "Chan-Lam": "Chan–Lam", "SNAr": "S$_\\mathregular{N}$Ar",
         "Hiyama": "Hiyama", "Grignard": "Grignard addition", "C-O competitive reactions": "C–O coupling",
         "Reductive-amination": "Reductive amination", "Esterification": "Esterification",
         "Sonogashira": "Sonogashira", "Suzuki Bpin": "Suzuki–Miyaura (Bpin)"}

SHORT = {**NAMES, "Buchwald-Hartwig": "B–H", "Acylation": "Amide", "Suzuki": "Suzuki",
         "Grignard": "Grignard add.", "C-O competitive reactions": "C–O", "Reductive-amination": "Red. amin.",
         "Esterification": "Ester", "Sonogashira": "Sonog.", "Ullmann": "Ullmann", "Suzuki Bpin": "Suzuki (Bpin)"}

ROLE = {"Acylation": ("acid", "amine"), "Esterification": ("acid", "alcohol"),
        "Buchwald-Hartwig": ("aryl halide", "amine"), "SNAr": ("aryl halide", "amine"),
        "Chan-Lam": ("boronic acid", "amine"), "Reductive-amination": ("carbonyl", "amine"),
        "Suzuki": ("aryl halide", "boronic acid"), "Suzuki Bpin": ("aryl halide", "boronate"),
        "Heck": ("aryl bromide", "alkene"), "Sonogashira": ("aryl bromide", "alkyne"),
        "Hiyama": ("aryl halide", "silane"), "Kumada": ("aryl halide", "RMgX"),
        "Grignard": ("carbonyl", "RMgX"), "Ullmann": ("aryl bromide", "phenol"),
        "C-O competitive reactions": ("aryl bromide", "O-nucleophile")}

def dec(v):
    return f"{v:+.1f}".replace("-", "−")

def signed(v):
    return f"{v:+.0f}".replace("-", "−")

def table(name):
    return pd.read_csv(SOURCES / name)


def bool_series(series):
    return series.astype(str).str.lower().eq("true")

def style():
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.titlesize": 8.5,
                         "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                         "legend.fontsize": 7, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": INK, "ytick.color": INK, "axes.linewidth": .65,
                         "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
                         "mathtext.bf": "Arial:bold", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": SVG_TEXT,
                         "svg.hashsalt": "jacs-building-block-yield",
                         "savefig.facecolor": "white", "figure.facecolor": "white"})

def panel(ax, letter, title=None, x=-.1):
    """Panel letter and title share one baseline, 6 pt above the top of the axes; finalize() puts the letter on the left
    edge of the panel and the title after it."""
    a = ax.annotate(letter, xy=(x, 1), xycoords="axes fraction", xytext=(0, 6), textcoords="offset points",
                    fontsize=10, fontweight="bold", va="baseline", ha="left", gid="panel-letter")
    a.panel_ax = ax
    if title:
        b = ax.annotate(title, xy=(x, 1), xycoords="axes fraction", xytext=(GAP_PT, 6), textcoords="offset points",
                        fontsize=8.5, fontweight="bold", va="baseline", ha="left", gid="panel-title")
        b.panel_letter = a


# One type scale and one panel layout for every figure (applied by finalize() before saving).
LABEL_CASE = "lower"        # "upper" for journals that print capital panel letters
SVG_TEXT = "none"           # "path" writes SVG text as outlines, so any SVG reader draws it with the same spacing
GAP_PT = 13.0               # panel title starts this far after the left edge of its letter
CAP = .716                  # cap height of Arial as a fraction of the font size: titles are raised so capitals share one top


def tier(size):
    """Every text other than a panel letter (10 pt) or panel title (8.5 pt): 7 pt, or 6 pt for small annotations."""
    return 7.0 if size >= 6.5 else 6.0


def _shift(fig, t, dx, dy=0.0):
    """Move a text dx display pixels to the right and dy up."""
    if isinstance(t, Annotation):
        ox, oy = t.xyann
        t.xyann = (ox + dx * 72 / fig.dpi, oy + dy * 72 / fig.dpi)
    else:
        X, Y = t.get_transform().transform(t.get_position())
        t.set_position(tuple(t.get_transform().inverted().transform((X + dx, Y + dy))))


def finalize(fig):
    texts = [t for t in fig.findobj(Text) if t.get_text().strip() and t.get_visible()]
    for t in texts:
        gid = t.get_gid()
        if gid == "scheme":
            continue
        if gid == "panel-letter":
            t.set_text(t.get_text().upper() if LABEL_CASE == "upper" else t.get_text().lower())
            t.set_fontsize(10.0)
        elif gid == "panel-title":
            t.set_fontsize(8.5)
        else:
            t.set_fontsize(tier(t.get_fontsize()))
    letters = [t for t in texts if t.get_gid() == "panel-letter"]
    titles = [t for t in texts if t.get_gid() == "panel-title"]
    if not letters:
        return
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    for t in letters:   # a letter drawn by panel() sits on the left edge of its panel, tick labels included
        ax = getattr(t, "panel_ax", None)
        if ax is not None:
            _shift(fig, t, ax.get_tightbbox(r).x0 - t.get_window_extent(r).x0)
    fig.canvas.draw()
    for t in letters + titles:
        t.set_visible(False)
    left = fig.get_tightbbox(r).x0 * fig.dpi
    for t in letters + titles:
        t.set_visible(True)
    x0 = {t: t.get_window_extent(r).x0 for t in letters}
    first = min(x0.values())
    edge = min(left, first)
    for t in letters:   # the first column: every letter on the left edge of the drawing
        if x0[t] < first + 3 * fig.dpi / 72:
            _shift(fig, t, edge - x0[t])
    fig.canvas.draw()
    for t in titles:   # after its letter, on the letter's baseline raised so that the capitals of both share one top
        letter = t.panel_letter
        _shift(fig, t, letter.get_window_extent(r).x0 + GAP_PT * fig.dpi / 72 - t.get_window_extent(r).x0,
               CAP * (letter.get_fontsize() - t.get_fontsize()) * fig.dpi / 72 if not getattr(t, "raised", False) else 0.0)
        t.raised = True
    fig.canvas.draw()


def save_figure(fig, name):
    finalize(fig)
    for ext in ("pdf", "svg", "png"):
        metadata = {"Creator": "build_figures.py", "CreationDate": None, "ModDate": None} if ext == "pdf" else (
            {"Creator": "build_figures.py", "Date": None} if ext == "svg" else {"Software": "build_figures.py"})
        fig.savefig(OUT / f"{name}.{ext}", dpi=600 if ext == "png" else 300, bbox_inches="tight", pad_inches=.06, metadata=metadata)
    plt.close(fig)

def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=3, width=.65)

_NUM = r"-?\d+(?:\.\d+)?"

def _svg_path(d):
    tokens = re.findall(rf"[MLQCZ]|{_NUM}", d)
    verts, codes, i, cmd = [], [], 0, None
    while i < len(tokens):
        t = tokens[i]
        if t in "MLQCZ":
            cmd, i = t, i + 1
            if cmd == "Z":
                verts.append(verts[-1] if verts else (0, 0)); codes.append(MPath.CLOSEPOLY)
            continue
        n = {"M": 2, "L": 2, "Q": 4, "C": 6}[cmd]
        vals = [float(v) for v in tokens[i:i + n]]
        i += n
        if cmd in ("M", "L"):
            verts.append(tuple(vals)); codes.append(MPath.MOVETO if cmd == "M" else MPath.LINETO)
            cmd = "L" if cmd == "M" else cmd
        elif cmd == "Q":
            verts += [tuple(vals[:2]), tuple(vals[2:])]; codes += [MPath.CURVE3] * 2
        else:
            verts += [tuple(vals[:2]), tuple(vals[2:4]), tuple(vals[4:])]; codes += [MPath.CURVE4] * 3
    return MPath(np.array(verts), codes)

def upright(smiles, orient):
    """RDKit's 2D depiction turned rigidly so that its first six-membered ring stands with vertical side edges, then by
    orient x 60 degrees. Connectivity, bond lengths and atom labels are unchanged. Returns the molecule and that ring."""
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    ring = max(mol.GetRingInfo().AtomRings(), key=lambda r: (len(r) == 6, len(r)))
    p, q = conf.GetAtomPosition(ring[0]), conf.GetAtomPosition(ring[1])
    a = math.radians((30 - math.degrees(math.atan2(q.y - p.y, q.x - p.x))) % 60 + 60 * orient)
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(p.x * math.cos(a) - p.y * math.sin(a), p.x * math.sin(a) + p.y * math.cos(a), 0.0))
    return mol, ring

def depiction(smiles, orient=None):
    """The ACS 1996 SVG of a molecule, its width and height, and, with orient, the height of its ring in the SVG."""
    if orient is None:
        svg, ring_y = Draw.MolToACS1996SVG(Chem.MolFromSmiles(smiles)), None
    else:
        mol, ring = upright(smiles, orient)
        drawer = rdMolDraw2D.MolDraw2DSVG(-1, -1)
        rdMolDraw2D.SetACS1996Mode(drawer.drawOptions(), Draw.MeanBondLength(mol))   # as Draw.MolToACS1996SVG
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        ring_y = sum(drawer.GetDrawCoords(i).y for i in ring) / len(ring)
    w = float(re.search(r"width='([\d.]+)px'", svg).group(1))
    h = float(re.search(r"height='([\d.]+)px'", svg).group(1))
    return svg, w, h, ring_y

def structure(ax, smiles, cx, cy, scale, orient=None):
    """Draw a molecule in ACS 1996 style in the axes' point coordinates, centred at (cx, cy). With orient, the depiction
    of upright() is drawn, centred on its width and on the height of its ring, so rings drawn at one cy share a line."""
    svg, w, h, ring_y = depiction(smiles, orient)
    ox, oy = cx - w * scale / 2, cy + (h / 2 if ring_y is None else ring_y) * scale
    for el in re.finditer(r"<path ([^>]*?)/>", svg, re.S):
        attrs = el.group(1)
        d = re.search(r"d='([^']*)'", attrs, re.S).group(1)
        st = re.search(r"style='([^']*)'", attrs)
        st = st.group(1) if st else "fill:" + re.search(r"fill='([^']*)'", attrs).group(1)
        path = _svg_path(d)
        v = path.vertices.copy()
        v[:, 0], v[:, 1] = ox + v[:, 0] * scale, oy - v[:, 1] * scale
        path = MPath(v, path.codes)
        if re.search(r"fill:([^;]+)", st).group(1) != "none":
            ax.add_patch(PathPatch(path, facecolor="black", edgecolor="none", lw=0))
        else:
            sw = re.search(r"stroke-width:([\d.]+)px", st)
            ax.add_patch(PathPatch(path, facecolor="none", edgecolor="black", lw=(float(sw.group(1)) if sw else .6) * scale,
                                   capstyle="butt", joinstyle="miter"))
    return w * scale, h * scale

def structure_strip(fig, rect, rows, scale=.62, label_color=INK, centre=.58):
    """A strip of numbered structures under a scatter; rows carry smiles, number, a caption and optionally orient (see
    structure). Numbers and captions sit at fixed heights; the structures are drawn on the line at centre x height."""
    ax = fig.add_axes(rect)
    ax.set_axis_off()
    width_pt, height_pt = rect[2] * fig.get_figwidth() * 72, rect[3] * fig.get_figheight() * 72
    ax.set(xlim=(0, width_pt), ylim=(0, height_pt))
    step = width_pt / len(rows)
    line = height_pt * centre
    oriented = [depiction(r["smiles"], r["orient"])[1:] for r in rows if r.get("orient") is not None]
    if oriented:   # one ring line for the row, placed so that the row's drawings are centred at centre x height
        line -= (max(ry for _, _, ry in oriented) - max(h - ry for _, h, ry in oriented)) * scale / 2
    for k, r in enumerate(rows):
        cx = step * (k + .5)
        structure(ax, r["smiles"], cx, line, scale, r.get("orient"))
        ax.text(cx, height_pt * .1, r["caption"], ha="center", va="center", fontsize=6.3, color=MUTED)
        ax.text(cx - step * .42, height_pt * .95, str(r["number"]), ha="left", va="top", fontsize=7.5,
                fontweight="bold", color=r.get("color", label_color))
    return ax


PURPLE = "#766196"





def extended_data_fig2():
    defined = table("ExtendedDataFig2_defined_benchmark.csv"); pooled = table("ExtendedDataFig2_pooled_shortfall.csv")
    fig, ax = plt.subplots(figsize=(4.0, 3.7)); fig.subplots_adjust(left=.16, right=.97, top=.95, bottom=.14)
    ax.plot([0, 1], [0, 1], color=GREY, lw=.9, ls=(0, (4, 3)))
    ax.axhline(0, color="#CAD0D5", lw=.7)
    for side, color in (("elec", BLUE), ("nuc", ORANGE)):
        r = defined.loc[defined.side.eq(side)]
        ax.scatter(r.geometric_benchmark, r.rho_R_half, s=10 + 1.1 * r.molecules, c=color, alpha=.9, edgecolors="white", linewidths=.5,
                   label="Electrophiles" if side == "elec" else "Nucleophiles")
    ax.set(xlim=(0, 1), ylim=(-.6, 1), xlabel="Agreement within each reaction, new partners", ylabel="Agreement between the two reactions")
    for ra, rb, label, dx, dy in (("Heck", "Suzuki", "Heck / Suzuki", -.30, .22), ("Acylation", "Chan-Lam", "Amide / Chan–Lam", -.38, .20)):
        r = defined.loc[defined.rxn_a.eq(ra) & defined.rxn_b.eq(rb)].iloc[0]
        ax.annotate(label, (r.geometric_benchmark, r.rho_R_half), xytext=(r.geometric_benchmark + dx, r.rho_R_half + dy), fontsize=6.5,
                    arrowprops={"arrowstyle": "-", "lw": .5, "color": GREY})
    p = pooled.loc[pooled.minimum_partners.eq(2)].iloc[0]
    ax.text(.97, .04, f"{int(p.positive_benchmark_comparisons)} comparisons\nmean shortfall {p.D_G:.2f} ({p.low:.2f}–{p.high:.2f})",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.6)
    ax.legend(loc="upper left", frameon=False, fontsize=6.5, handletextpad=.2)
    clean(ax)
    save_figure(fig, "ExtendedDataFig3")

def figure_s1():
    inc = table("ExtendedDataFig7_prediction_all_reactions.csv")
    order = inc.groupby("reaction").scored_pairs.first().sort_values(ascending=False).index.tolist()
    fig, ax = plt.subplots(figsize=(6.5, 5.2)); fig.subplots_adjust(left=.27, right=.97, top=.93, bottom=.1)
    for y in range(len(order)):
        if y % 2 == 0:
            ax.axhspan(y - .48, y + .48, color=LIGHT, zorder=0)
    for i, rxn in enumerate(order):
        for side, shift, color in (("elec", -.15, BLUE), ("nuc", .15, ORANGE)):
            r = inc.loc[inc.reaction.eq(rxn) & inc.side.eq(side)].iloc[0]
            lo, hi = max(r.low, -2), min(r.high, 1)
            ax.plot([lo, hi], [i + shift] * 2, color=color, lw=1)
            ax.plot(np.clip(r.increment, -2, 1), i + shift, "o", ms=4, mfc=color if r.interval_above_zero else "white", mec=color)
    n = inc.groupby("reaction").scored_pairs.first()
    ax.set(yticks=range(len(order)), yticklabels=[f"{NAMES[r]} ({int(n[r])})" for r in order], ylim=(len(order) - .5, -.6),
           xlim=(-2.05, 1.05), xlabel="Share of prediction error removed by one partner's contributions")
    ax.axvline(0, color="#ABB3B9", lw=.8); ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False); clean(ax)
    ax.legend(handles=[Line2D([0], [0], marker="o", color=BLUE, label="Electrophile"), Line2D([0], [0], marker="o", color=ORANGE, label="Nucleophile")],
              loc="upper left", frameon=False)
    save_figure(fig, "ExtendedDataFig8")

FEATURE_SHORT = {"halogen next to ring N": "ring N beside X", "N–H": "N–H", "substituent beside the halogen": "ortho substituent",
                 "chloride leaving group": "X = Cl", "iodide leaving group": "X = I", "sulfonate leaving group": "X = OSO$_2$R",
                 "heteroaromatic ring at the halogen": "heteroaryl at X", "five-membered ring at the halogen": "five-membered ring at X",
                 "fused ring at the halogen": "fused ring at X", "electron-withdrawing group": "EWG", "electron-donating group": "EDG"}

def extended_data_profiles():
    camp = table("ExtendedDataFig6a_campaign_profiles.csv"); comp = table("ExtendedDataFig6bc_profile_comparisons.csv")
    fig = plt.figure(figsize=(7.0, 2.9))
    ax = fig.add_axes([.07, .17, .25, .66])
    lim = 1.15 * float(np.abs(camp[["first_campaign_pp", "second_campaign_pp"]].to_numpy()).max())
    ax.plot([-lim, lim], [-lim, lim], ls="--", color=GREY, lw=.7)
    ax.axhline(0, color="#CED4D8", lw=.6); ax.axvline(0, color="#CED4D8", lw=.6)
    ax.scatter(camp.first_campaign_pp, camp.second_campaign_pp, s=20, color=RED, ec="white", lw=.4, zorder=3)
    # Label placement only: the four features near the origin are spread to the four sides.
    place = {"electron-donating group": ((-4, 3), "right"), "electron-withdrawing group": ((5, -1), "left"),
             "heteroaromatic ring at the halogen": ((5, -3), "left"), "fused ring at the halogen": ((4, 2), "left")}
    for _, r in camp.iterrows():
        off, ha = place.get(r.feature, ((3, 2), "left"))
        ax.annotate(FEATURE_SHORT[r.feature], (r.first_campaign_pp, r.second_campaign_pp), xytext=off, textcoords="offset points",
                    fontsize=5.8, color=MUTED, ha=ha)
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel="First campaign (pp)", ylabel="Second campaign (pp)"); clean(ax)
    c = comp.loc[comp.comparison.eq("second Buchwald-Hartwig campaign against an analysed reaction") & comp.reaction_a.eq("Buchwald-Hartwig") & comp.reaction_b.eq("Buchwald-Hartwig(CMC)")].iloc[0]
    ax.text(.03, .97, f"structure: r = {c.pearson:.2f}\nbuilding blocks: {c.record_rank_agreement:.2f}".replace("-", "−"), transform=ax.transAxes, va="top", fontsize=6.4)
    panel(ax, "a", "Two B–H campaigns", x=-.2)
    ax = fig.add_axes([.48, .17, .17, .66])
    within = comp.loc[comp.comparison.str.startswith("same reaction")]
    order = sorted(set(within.reaction_a), key=lambda r: -within.loc[within.reaction_a.eq(r)].pearson.min())
    for i, r in enumerate(order):
        for kind, filled in (("same reaction, second condition set", True),):
            v = within.loc[within.reaction_a.eq(r) & within.comparison.eq(kind), "pearson"]
            if len(v):
                ax.plot(v.iloc[0], i, "o", ms=4, mfc=TEAL if filled else "white", mec=TEAL, mew=1)
    ax.set(yticks=range(len(order)), yticklabels=[SHORT[r] for r in order], ylim=(len(order) - .4, -.6), xlim=(-.5, 1.05), xlabel="Structural profile, r")
    ax.axvline(0, color="#ABB3B9", lw=.7); ax.tick_params(axis="y", length=0, labelsize=6.4); ax.spines["left"].set_visible(False); clean(ax)
    panel(ax, "b", "New conditions", x=-.62)
    ax = fig.add_axes([.77, .17, .21, .66])
    two = comp.loc[comp.comparison.eq("two analysed reactions")]
    ax.axhline(0, color="#CED4D8", lw=.6); ax.axvline(0, color="#CED4D8", lw=.6)
    ax.scatter(two.record_rank_agreement, two.pearson, s=12, color=GREY, ec="white", lw=.3, label="two reactions")
    ax.scatter([c.record_rank_agreement], [c.pearson], s=26, color=RED, ec="white", lw=.4, zorder=3, label="two B–H campaigns")
    ax.set(xlim=(-.6, 1), ylim=(-1, 1.05), xlabel="Building blocks: rank agreement", ylabel="Structural profile, r"); clean(ax)
    ax.legend(loc="lower right", frameon=False, fontsize=6.2, handletextpad=.1, borderaxespad=.1)
    panel(ax, "c", "Two reactions", x=-.32)
    save_figure(fig, "ExtendedDataFig7")

def figure_s3():
    und = table("ExtendedDataFig3_undefined_benchmark.csv").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(6.5, 3.6)); fig.subplots_adjust(left=.3, right=.97, top=.95, bottom=.24)
    for i, r in und.iterrows():
        if i % 2 == 0:
            ax.axhspan(i - .45, i + .45, color=LIGHT, zorder=0)
        ax.plot([r.rho_P_a, r.rho_P_b], [i, i], color="#BAC2C7", lw=.8)
        ax.plot(r.rho_P_a, i, marker="|", ms=8, mew=1.4, color=BLUE)
        ax.plot(r.rho_P_b, i, marker="|", ms=8, mew=1.4, color=ORANGE)
        ax.plot(r.rho_R_half, i, "o", color=RED, ms=4)
    labels = [f"{SHORT[r.rxn_a]} / {SHORT[r.rxn_b]} · {'E' if r.side == 'elec' else 'N'} ({int(r.molecules)})" for _, r in und.iterrows()]
    ax.axvline(0, color=GREY, lw=.7)
    ax.set(xlim=(-1, 1), ylim=(len(und) - .4, -.6), yticks=range(len(und)), yticklabels=labels, xlabel="Rank agreement on the same shared molecules")
    ax.tick_params(axis="y", length=0, labelsize=6.6); ax.spines["left"].set_visible(False)
    ax.legend(handles=[Line2D([0], [0], marker="|", ls="", ms=7, color=BLUE, label="First reaction, new partners"),
                       Line2D([0], [0], marker="|", ls="", ms=7, color=ORANGE, label="Second reaction, new partners"),
                       Line2D([0], [0], marker="o", ls="", ms=4, color=RED, label="Between reactions")],
              loc="upper center", bbox_to_anchor=(.4, -.17), ncol=3, frameon=False, fontsize=6.5)
    save_figure(fig, "ExtendedDataFig4")

def figure_s4():
    om = table("ExtendedDataFig4a_leave_one_reaction_out.csv"); sens = table("ExtendedDataFig4b_outcome_scale.csv")
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.6), gridspec_kw={"width_ratios": [1.1, 1]})
    fig.subplots_adjust(left=.17, right=.97, top=.9, bottom=.2, wspace=.9)
    ax = axes[0]; panel(ax, "a", "Leaving one reaction out")
    reactions = list(om.reaction_left_out.unique())
    for degree, shift, color in ((2, -.12, BLUE), (3, .12, ORANGE)):
        for _, r in om.loc[om.minimum_partners.eq(degree)].iterrows():
            ax.plot(r.D_G, reactions.index(r.reaction_left_out) + shift, "o", ms=3.6, color=color)
    ax.set(yticks=range(len(reactions)), yticklabels=[SHORT[x] for x in reactions], ylim=(len(reactions) - .4, -.6), xlim=(0, .5), xlabel="Mean shortfall")
    ax.tick_params(axis="y", labelsize=6.5)
    ax.legend(handles=[Line2D([0], [0], marker="o", ls="", color=BLUE, label="≥2 partners"), Line2D([0], [0], marker="o", ls="", color=ORANGE, label="≥3 partners")],
              loc="upper center", bbox_to_anchor=(.4, -.13), frameon=False, ncol=2)
    ax = axes[1]; panel(ax, "b", "Other yield scales")
    for i, r in sens.reset_index(drop=True).iterrows():
        ax.errorbar(r.D_G, i, xerr=[[r.D_G - r.low], [r.high - r.D_G]], fmt="o", ms=4, color=TEAL, capsize=2, lw=1)
        ax.text(.01, i + .28, f"{int(r.positive_benchmark_comparisons)}/{int(r.eligible_comparisons)} comparisons", fontsize=6.2, va="top", color=MUTED)
    ax.axvline(0, color=GREY, lw=.7)
    ax.set(yticks=range(len(sens)), yticklabels=[str(x).replace("Without SNAr and Chan–Lam", "Without S$_\\mathregular{N}$Ar\nand Chan–Lam") for x in sens.label],
           ylim=(len(sens) - .4, -.65), xlim=(-.02, .8), xlabel="Mean shortfall")
    save_figure(fig, "ExtendedDataFig5")

def figure_s2():
    data = table("ExtendedDataFig5_reaction_map.csv")
    order = sorted(set(data.rxn_a) | set(data.rxn_b)); arr = np.full((len(order), len(order)), np.nan); counts = {}
    for _, r in data.iterrows():
        a, b = order.index(r.rxn_a), order.index(r.rxn_b)
        i, j = (max(a, b), min(a, b)) if r.side == "elec" else (min(a, b), max(a, b))
        arr[i, j] = r.rho_obs; counts[i, j] = int(r.n_shared)
    fig, ax = plt.subplots(figsize=(6.5, 6.3)); fig.subplots_adjust(left=.2, right=.97, top=.97, bottom=.24)
    cmap = plt.get_cmap("RdBu").copy(); cmap.set_bad("#EDF0F2")
    im = ax.imshow(arr, vmin=-1, vmax=1, cmap=cmap)
    for (i, j), n in counts.items():
        ax.text(j, i, str(n), ha="center", va="center", fontsize=6.5, color="white" if abs(arr[i, j]) > .65 else INK)
    for i in range(len(order)):
        ax.add_patch(Rectangle((i - .5, i - .5), 1, 1, fc="#C8CDD1", ec="white", lw=.6))
    ax.set(xticks=range(len(order)), yticks=range(len(order)), xticklabels=[SHORT[x] for x in order], yticklabels=[SHORT[x] for x in order])
    plt.setp(ax.get_xticklabels(), rotation=90); ax.tick_params(length=0)
    ax.set_xticks(np.arange(-.5, len(order), 1), minor=True); ax.set_yticks(np.arange(-.5, len(order), 1), minor=True)
    ax.grid(which="minor", color="white", lw=.6); ax.tick_params(which="minor", length=0)
    cax = fig.add_axes([.3, .08, .5, .022]); cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[-1, 0, 1])
    cb.set_label("Rank agreement of shared molecules")
    save_figure(fig, "ExtendedDataFig6")

FIGURES = ["ExtendedDataFig3", "ExtendedDataFig4", "ExtendedDataFig5", "ExtendedDataFig6", "ExtendedDataFig7", "ExtendedDataFig8"]

def plot_figures():
    style()
    extended_data_fig2(); figure_s1(); figure_s2(); figure_s3(); figure_s4(); extended_data_profiles()
    print("Exported Extended Data Figs 3-8.")

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    plot_figures()
