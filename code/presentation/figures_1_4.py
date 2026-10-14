"""Figures 1 and 4 and Extended Data Fig. 2 of the Article, drawn from the Source Data
(article_sources.py) only. Figure 1a draws the fifteen reaction schemes of Extended Data Fig. 1 with their design counts.
Style, names, colours and structure drawing come from plot_figures.py; panel headers and the
text-overlap check from figures_2_3.py."""
import sys

sys.dont_write_bytecode = True

import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyBboxPatch, PathPatch, Rectangle  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC = ROOT / "output/presentation/source_data"
OUT = ROOT / "output/presentation/figures"
import plot_figures as BF  # noqa: E402  shared style, names, colours and structure drawing
import figures_2_3 as V1  # noqa: E402  panel headers and the text-overlap check

NAMES, ROLE = BF.NAMES, BF.ROLE
BLUE, ORANGE, TEAL, RED, PURPLE = BF.BLUE, BF.ORANGE, BF.TEAL, BF.RED, BF.PURPLE
INK, GREY, LIGHT, MUTED = BF.INK, BF.GREY, BF.LIGHT, BF.MUTED








def tab(name):
    return pd.read_csv(SRC / name)


def save(fig, name):
    BF.finalize(fig)
    for ext in ("pdf", "svg", "png"):
        meta = {"Creator": "figures_1_4.py", "CreationDate": None, "ModDate": None} if ext == "pdf" else (
            {"Creator": "figures_1_4.py", "Date": None} if ext == "svg" else {"Software": "figures_1_4.py"})
        fig.savefig(OUT / f"{name}.{ext}", dpi=600 if ext == "png" else 300, bbox_inches="tight", pad_inches=.06, metadata=meta)


CHANGE = {"Conditions": TEAL, "Partners": PURPLE, "Reaction": RED}   # one colour per change, as in Fig. 1c and Fig. 2a


SCHEMES = {s["reaction"]: s for s in json.loads((ROOT / "code" / "assets" / "Figure1a_schemes.json").read_text(encoding="utf-8"))["schemes"]}
SCHEME_SCALE, ARROW = .61, 28.0   # scale of the Extended Data Fig. 1 drawing in Fig. 1a; arrow length (source points) without labels


def scheme_shift(s):
    bar = next(p for p in s["paths"] if p.get("role") == "arrow bar")
    return bar["rect"][2] - bar["rect"][0] - ARROW, bar["rect"][2]


def scheme_size(s):
    """Width and height of a scheme in source points once its arrow is shortened."""
    shift, _ = scheme_shift(s)
    x0, y0, x1, y1 = s["bbox"]
    return x1 - x0 - shift, y1 - y0


def draw_scheme(ax, s, left, top, scale=SCHEME_SCALE):
    """Draw a reaction scheme of Extended Data Fig. 1 (code/assets/Figure1a_schemes.json) with the top-left corner of its box
    at (left, top) in point coordinates of ax (y up). Paths and glyph runs are the drawing's own; the arrow bar is drawn
    ARROW source points long, without its reagent labels, and what lies right of it moves left with it."""
    shift, end = scheme_shift(s)
    bx0, by0 = s["bbox"][0], s["bbox"][1]

    def to(x, y, right):
        return left + (x - bx0 - (shift if right else 0)) * scale, top - (y - by0) * scale
    for p in s["paths"]:
        right = p["rect"][0] >= end - 1.5
        bar = p.get("role") == "arrow bar"
        items = [["re", [p["rect"][0], p["rect"][1], p["rect"][0] + ARROW, p["rect"][3]]]] if bar else p["items"]
        verts, codes, cur = [], [], None
        for it in items:
            if it[0] == "re":
                x0, y0, x1, y1 = it[1]
                verts += [to(x, y, right) for x, y in ((x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0))]
                codes += [MPath.MOVETO, MPath.LINETO, MPath.LINETO, MPath.LINETO, MPath.CLOSEPOLY]
                cur = None
                continue
            if cur != it[1]:
                verts.append(to(*it[1], right))
                codes.append(MPath.MOVETO)
            if it[0] == "l":
                verts.append(to(*it[2], right))
                codes.append(MPath.LINETO)
                cur = it[2]
            else:
                verts += [to(*q, right) for q in it[2:5]]
                codes += [MPath.CURVE4] * 3
                cur = it[4]
        filled, stroked = "f" in p["type"], "s" in p["type"]
        if filled and codes[-1] != MPath.CLOSEPOLY:
            verts.append(verts[0])
            codes.append(MPath.CLOSEPOLY)
        ax.add_patch(PathPatch(MPath(verts, codes), facecolor=tuple(p["fill"]) if filled else "none",
                               edgecolor=tuple(p["color"]) if stroked else "none", lw=(p["width"] or 0) * scale if stroked else 0,
                               joinstyle="miter", capstyle="butt"))
    for t in s["texts"]:
        x, y = to(t["x"], t["y"], t["bbox"][0] >= end - 1.5)
        ax.text(x, y, t["s"], fontsize=t["size"] * scale, ha="left", va="baseline", family="Arial", color="black", gid="scheme")


def figure1():
    t = tab("Figure1a_reactions.csv")
    assert set(t.reaction) == set(SCHEMES), set(t.reaction) ^ set(SCHEMES)
    W = 7.0
    LOW = .49 * 5.7 * 72 + 14.0           # panels b and c keep the size they had in a 5.7-inch figure
    col = W * 72 / 3
    rows = [list(t.reaction.iloc[i:i + 3]) for i in range(0, len(t), 3)]
    band = [max(scheme_size(SCHEMES[r])[1] for r in row) * SCHEME_SCALE for row in rows]
    NAME, GAP, COUNT, PAD = 8.5, 4.5, 7.0, 4.0   # point heights inside a tile
    tile = [PAD + NAME + GAP + b + GAP + COUNT + PAD for b in band]
    a_h = sum(tile) + 5.0 * (len(rows) - 1)
    H = (LOW + a_h + 26.0) / 72
    fig = plt.figure(figsize=(W, H))
    k = 5.7 / H                            # rescales the lower panels' figure fractions
    ax = fig.add_axes([0, LOW / 72 / H, 1, a_h / 72 / H]); ax.set_axis_off(); ax.set(xlim=(0, W * 72), ylim=(0, a_h))
    y, rec = a_h, {r.reaction: r for r in t.itertuples()}
    for row, b, h in zip(rows, band, tile):
        for c, rx in enumerate(row):
            r, x0, w = rec[rx], c * col + 2.5, col - 5.0
            ax.add_patch(FancyBboxPatch((x0, y - h), w, h, boxstyle="round,pad=0,rounding_size=3.5", fc=LIGHT, ec="none"))
            ax.text(x0 + 6, y - PAD, NAMES[rx], fontsize=7.2, fontweight="bold", va="top")
            ax.text(x0 + w - 6, y - PAD, f"{r.shared_molecules} shared", fontsize=6.2, va="top", ha="right", color=MUTED)
            sw, sh = scheme_size(SCHEMES[rx])
            assert sw * SCHEME_SCALE <= w - 8, (rx, sw * SCHEME_SCALE, w)
            draw_scheme(ax, SCHEMES[rx], x0 + (w - sw * SCHEME_SCALE) / 2, y - PAD - NAME - GAP - (b - sh * SCHEME_SCALE) / 2)
            ax.text(x0 + 6, y - h + PAD, f"{r.electrophiles} × {r.nucleophiles} molecules, {r.pairs} pairs", fontsize=6.2,
                    va="bottom", color=MUTED)
        y -= h + 5.0
    V1.header(fig, .005, 1 - 18.5 / 72 / H, "a", "Fifteen reactions, crossed under twenty conditions each")

    rows = tab("Figure1b_three_changes.csv")
    ax = fig.add_axes([.02, .03 * k, .42, .43 * k]); ax.set_axis_off(); ax.set(xlim=(0, 1), ylim=(0, 1))

    def tiles(x, y, n, fc, ec, hatch=None):
        for j in range(n):
            ax.add_patch(Rectangle((x + j * .048, y), .04, .06, fc=fc, ec=ec, lw=.7, hatch=hatch))
    for i, r in rows.iterrows():
        y, c = .82 - i * .33, CHANGE[r.change]
        ax.text(0, y + .1, r.change, color=c, fontweight="bold", fontsize=8.0, va="center")
        ax.text(0, y + .02, r.question, fontsize=7.0, va="center")
        ax.text(0, y - .06, "kept: " + r.kept, fontsize=6.4, color=MUTED, va="center")
        x0 = .67
        if i == 0:
            ax.text(x0 - .02, y + .07, "conditions A", fontsize=5.8, ha="right", va="center", color=MUTED)
            ax.text(x0 - .02, y - .03, "conditions B", fontsize=5.8, ha="right", va="center", color=MUTED)
            tiles(x0, y + .04, 6, c, c); tiles(x0, y - .06, 6, "white", c, "////")
        elif i == 1:
            ax.text(x0 - .02, y + .02, "partners", fontsize=5.8, ha="right", va="center", color=MUTED)
            tiles(x0, y - .01, 3, c, c); tiles(x0 + 3 * .048, y - .01, 3, "white", c, "////")
            ax.text(x0 + .07, y - .06, "A", fontsize=6, ha="center", color=c); ax.text(x0 + .215, y - .06, "B", fontsize=6, ha="center", color=c)
        else:
            ax.text(x0 - .02, y + .07, "reaction 1", fontsize=5.8, ha="right", va="center", color=MUTED)
            ax.text(x0 - .02, y - .03, "reaction 2", fontsize=5.8, ha="right", va="center", color=MUTED)
            tiles(x0, y + .04, 6, c, c); tiles(x0, y - .06, 6, GREY, GREY)
    V1.header(fig, .005, .49 * k, "b", "Three ways to reuse a record")

    wb = tab("Figure1c_pooled_portability.csv")
    ax = fig.add_axes([.60, .165 * k, .385, .29 * k])
    color = {"conditions": TEAL, "partners": PURPLE, "reactions": RED}
    x, ticks, labels, centers = 0.0, [], [], {}
    for kind in ("conditions", "partners", "reactions"):
        start = x
        for r in wb.loc[wb.kind.eq(kind)].itertuples():
            filled = r.interval == "95% CI"
            ax.plot([x, x], [r.ci_low, r.ci_high], color=color[kind], lw=1.1)
            ax.plot(x, r.estimate, "o" if filled else "s", ms=4.6, mfc=color[kind] if filled else "white", mec=color[kind], mew=1.1)
            ticks.append(x); labels.append(r.label); x += 1
        centers[kind] = (start + x - 1) / 2
        x += .7
    ax.axhline(0, color="#ABB3B9", lw=.7)
    ax.set_xticks(ticks); ax.set_xticklabels(labels, fontsize=6.0, rotation=90)
    ax.set(ylim=(-.25, 1.02), xlim=(-.6, x - .1), ylabel="Rank agreement of contributions")
    for kind, text in (("conditions", "new\nconditions"), ("partners", "new\npartners"), ("reactions", "another\nreaction")):
        ax.text(centers[kind], 1.01, text, ha="center", va="bottom", fontsize=6.4, color=color[kind], fontweight="bold")
    BF.clean(ax)
    V1.header(fig, .5, .49 * k, "c", "How much of the record each change keeps")
    return fig


def figure4():
    pairs, metric = tab("Figure4a_heldout_amide_pairs.csv"), tab("Figure4a_metrics.csv").iloc[0]
    sel, pos = tab("Figure4b_selection_vs_agreement.csv"), tab("Figure4c_Heck_Suzuki_molecules.csv")
    fig = plt.figure(figsize=(7.0, 5.5))
    ax = fig.add_axes([.085, .6, .38, .32])
    ax.scatter(pairs.prediction_pp, pairs.target_pp, s=4, c=BLUE, alpha=.3, edgecolors="none", rasterized=True)
    lim = (-90, 60); ax.plot(lim, lim, color=GREY, lw=.8, ls=(0, (4, 3))); ax.set(xlim=lim, ylim=lim)
    ax.set_xlabel("Predicted (pp)"); ax.set_ylabel("Measured, condition-adjusted (pp)")
    ax.text(.04, .97, f"{int(metric.scored_pairs)} held-out amide pairs\nerror {metric.prediction_rmse_pp:.1f} pp\n(mean: {metric.reaction_mean_rmse_pp:.1f} pp)",
            transform=ax.transAxes, va="top", fontsize=6.4)
    BF.clean(ax)
    V1.header(fig, .005, .955, "a", "Held-out pairs")
    ax = fig.add_axes([.6, .6, .38, .32])
    ax.axhline(0, color=GREY, lw=.8)
    for side, color in (("elec", BLUE), ("nuc", ORANGE)):
        r = sel.loc[sel.side.eq(side)]
        ax.scatter(r.rho_obs, r.source_gain_pp, s=20, c=color, edgecolors="white", linewidths=.5, label="electrophiles" if side == "elec" else "nucleophiles")
    for ra, rb, side, label, dx, dy in (("Buchwald-Hartwig", "SNAr", "nuc", "B–H → S$_\\mathregular{N}$Ar", .08, 1.0),
                                        ("Heck", "Suzuki", "elec", "Heck → Suzuki", -.05, 2.2)):
        r = sel.loc[sel.rxn_a.eq(ra) & sel.rxn_b.eq(rb) & sel.side.eq(side)].iloc[0]
        ax.annotate(label, (r.rho_obs, r.source_gain_pp), xytext=(r.rho_obs + dx, r.source_gain_pp + dy), fontsize=6.3,
                    arrowprops={"arrowstyle": "-", "lw": .5, "color": GREY})
    ax.set(xlabel="Agreement of the two reactions", ylabel="Gain over a random ten (pp)", xlim=(-.35, .75))
    ax.legend(loc="lower right", frameon=False, fontsize=6.2, handletextpad=.1)
    BF.clean(ax)
    V1.header(fig, .515, .955, "b", "Borrowing a record")
    ax = fig.add_axes([.085, .08, .38, .36])
    ax.axhline(0, color="#CED4D8", lw=.7); ax.axvline(0, color="#CED4D8", lw=.7)
    ax.scatter(pos.x_contribution_pp, pos.y_contribution_pp, s=18, color=BLUE, alpha=.9, ec="white", lw=.4)
    for _, r in pos.loc[pos.drawn.gt(0)].iterrows():
        ax.annotate(str(int(r.drawn)), (r.x_contribution_pp, r.y_contribution_pp), xytext=(4, 3), textcoords="offset points", fontsize=7, fontweight="bold", color=INK)
    ax.set(xlabel="Heck (pp)", ylabel="Suzuki–Miyaura (pp)"); ax.margins(.12); BF.clean(ax)
    lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + .24 * (hi - lo))
    ax.text(.03, .97, f"{len(pos)} shared aryl bromides\nrank agreement {pos.rho_obs.iloc[0]:.2f}".replace("-", "−"), transform=ax.transAxes, va="top", fontsize=6.6)
    V1.header(fig, .005, .49, "c", "A pair of reactions whose records carry: Heck and Suzuki–Miyaura")
    strip = [{"smiles": r.mol, "orient": V1.ORIENT[r.mol], "number": int(r.drawn), "caption": f"{BF.signed(r.x_contribution_pp)} / {BF.signed(r.y_contribution_pp)}"}
             for _, r in pos.loc[pos.drawn.gt(0)].sort_values("drawn").iterrows()]
    assert len(strip) == 4
    BF.structure_strip(fig, [.53, .25, .42, .2], strip[:2], scale=.70, centre=.55)
    BF.structure_strip(fig, [.53, .03, .42, .2], strip[2:], scale=.70, centre=.55)
    return fig


def ed_fig2():
    m = tab("ExtendedDataFig2_amide_matrix.csv")
    fig = plt.figure(figsize=(4.8, 3.7))
    ax = fig.add_axes([.12, .14, .70, .80])
    n_acid, n_amine = m.acid_rank.max(), m.amine_rank.max()
    grid = np.full((n_acid, n_amine), np.nan)
    grid[m.acid_rank - 1, m.amine_rank - 1] = m.pair_effect_pp
    cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("white")
    lim = float(np.nanpercentile(np.abs(grid), 98))
    im = ax.imshow(grid, cmap=cmap, vmin=-lim, vmax=lim, aspect="auto", interpolation="nearest")
    ax.set(xticks=[], yticks=[])
    ax.set_xlabel(f"Amines and sulfonamides ({n_amine}), highest contribution at left")
    ax.set_ylabel(f"Carboxylic acids ({n_acid}), highest contribution at top")
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color("#C8CDD1"); sp.set_linewidth(.6)
    cax = fig.add_axes([.85, .14, .025, .80])
    cb = fig.colorbar(im, cax=cax); cb.set_label("Pair yield relative to reaction mean (pp)", fontsize=7)
    cb.ax.tick_params(labelsize=6.5, length=2)
    return fig




def main():
    BF.style()
    for name, make in (("Figure1", figure1), ("Figure4", figure4), ("ExtendedDataFig2", ed_fig2)):
        fig = make()
        BF.finalize(fig)   # text overlaps are checked on the finished layout
        bad = V1.overlaps(fig)
        save(fig, name)
        plt.close(fig)
        print(name, "text overlaps:", bad if bad else "none")


if __name__ == "__main__":
    main()
