"""Static figures for the README and write-up, in light and dark variants.

Figure 2 (benchmark_splits): balanced accuracy of each predictor under random, study-aware and
lineage-aware splits, one small panel per drug on a shared y-axis. A falling line is performance
that a random split attributed to the model but that does not survive grouping.

Figure 1 (colistin_mgrb): colistin sensitivity and specificity with 95% Wilson intervals for each
rule, with and without the targeted mgrB check (broth-microdilution labels, ADR-003). Only
mechanism-specified disruption is shown; the missense tier is excluded because leave-one-study-out
selection reduces its contribution from 19.1 to 0.9 points.

Figure 3 (colistin_transfer): sensitivity of the unfitted mgrB disruption rule in each held-out
study, showing that the mechanistic composition of colistin resistance differs between
collections.

Colors follow the project's validated categorical palette (3 slots, all-pairs CVD-safe in both
modes). Lines carry direct text labels in addition to the legend, so identity never depends on
color alone. Values: reports/benchmark/results.csv and reports/review/colistin_mgrb.md.

Outputs: reports/figures/<name>_light.png and <name>_dark.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

THEMES = {
    "light": {
        "surface": "#fcfcfb",
        "ink": "#0b0b0b",
        "ink2": "#52514e",
        "muted": "#898781",
        "grid": "#e1e0d9",
        "axis": "#c3c2b7",
        "series": ["#2a78d6", "#eb6834", "#1baf7a"],
    },
    "dark": {
        "surface": "#1a1a19",
        "ink": "#f0efec",
        "ink2": "#c3c2b7",
        "muted": "#898781",
        "grid": "#2c2c2a",
        "axis": "#383835",
        "series": ["#3987e5", "#d95926", "#199e70"],
    },
}
MODELS = [("rule", "Rule"), ("logreg", "Logistic regression"), ("gboost", "Gradient boosting")]
SCHEMES = [("random", "Random"), ("study", "Study"), ("lineage", "Lineage")]
DRUGS = ["meropenem", "ciprofloxacin", "gentamicin", "colistin", "tigecycline"]
FONT = ["system-ui", "-apple-system", "Segoe UI", "Helvetica Neue", "Arial", "DejaVu Sans"]


def style(ax, t: dict) -> None:
    ax.set_facecolor(t["surface"])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(t["axis"])
    ax.tick_params(colors=t["muted"], labelsize=9, length=0)
    ax.grid(axis="y", color=t["grid"], linewidth=0.8)
    ax.set_axisbelow(True)


def benchmark_figure(results: pd.DataFrame, theme: str, out: Path) -> None:
    t = THEMES[theme]
    r = results[(results["labels"] == "all")]
    fig, axes = plt.subplots(1, len(DRUGS), figsize=(13, 3.9), sharey=True)
    fig.patch.set_facecolor(t["surface"])
    xs = range(len(SCHEMES))
    for ax, drug in zip(axes, DRUGS, strict=True):
        style(ax, t)
        for color, (model, label) in zip(t["series"], MODELS, strict=True):
            d = r[(r.drug == drug) & (r.model == model)].set_index("scheme")
            ys = [d.loc[s, "balanced_accuracy"] for s, _ in SCHEMES]
            ax.plot(
                xs,
                ys,
                color=color,
                linewidth=2,
                marker="o",
                markersize=7,
                markeredgecolor=t["surface"],
                markeredgewidth=1.5,
                label=label,
                solid_capstyle="round",
            )
        ax.set_title(drug.capitalize(), color=t["ink"], fontsize=11, loc="left", pad=8)
        ax.set_xticks(list(xs), [s for _, s in SCHEMES])
        ax.set_xlim(-0.25, len(SCHEMES) - 0.25)
        ax.set_ylim(0.45, 1.0)
    # Direct labels on the panel with the widest spread, at the right end of each line.
    ax = axes[DRUGS.index("tigecycline")]
    d = r[r.drug == "tigecycline"]
    for model, label in MODELS:
        y = d[(d.model == model) & (d.scheme == "lineage")]["balanced_accuracy"].iloc[0]
        ax.annotate(
            label,
            (len(SCHEMES) - 1, y),
            xytext=(8, 0),
            textcoords="offset points",
            va="center",
            fontsize=8.5,
            color=t["ink2"],
            annotation_clip=False,
        )
    axes[0].set_ylabel("Balanced accuracy", color=t["ink2"], fontsize=10)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.06, 1.0),
        ncol=3,
        frameon=False,
        fontsize=9.5,
        labelcolor=t["ink2"],
    )
    fig.suptitle(
        "Random splits flatter models where rules fail",
        x=0.06,
        y=1.09,
        ha="left",
        color=t["ink"],
        fontsize=13,
        fontweight="bold",
    )
    fig.text(0.06, 0.985, "", color=t["ink2"])
    fig.text(
        0.06,
        -0.06,
        "Balanced accuracy on held-out genomes. Study = no BioProject in both "
        "training and test; lineage = no sequence type in both. 6,652 K. pneumoniae genomes, "
        "EUCAST v16.1 labels.",
        color=t["muted"],
        fontsize=8.5,
        ha="left",
    )
    fig.subplots_adjust(left=0.06, right=0.93, top=0.84, bottom=0.12, wspace=0.12)
    fig.savefig(out, dpi=200, facecolor=t["surface"], bbox_inches="tight")
    plt.close(fig)


def colistin_figure(table: pd.DataFrame, theme: str, out: Path) -> None:
    t = THEMES[theme]
    fig, axes = plt.subplots(
        1, 2, figsize=(11, 3.6), sharey=True, gridspec_kw={"width_ratios": [3, 2]}
    )
    fig.patch.set_facecolor(t["surface"])
    ys = list(range(len(table)))[::-1]
    color = t["series"][0]
    panels = [
        ("sensitivity", "sens_lo", "sens_hi", "Sensitivity (resistant called R)", (0, 1.08)),
        (
            "specificity",
            "spec_lo",
            "spec_hi",
            "Specificity (susceptible called S), zoomed",
            (0.85, 1.025),
        ),
    ]
    for ax, (metric, lo, hi, title, xlim) in zip(axes, panels, strict=True):
        style(ax, t)
        ax.grid(axis="y", visible=False)
        ax.grid(axis="x", color=t["grid"], linewidth=0.8)
        for y, (_, row) in zip(ys, table.iterrows(), strict=True):
            ax.plot([row[lo], row[hi]], [y, y], color=color, linewidth=2, solid_capstyle="round")
            ax.plot(
                row[metric],
                y,
                "o",
                color=color,
                markersize=8,
                markeredgecolor=t["surface"],
                markeredgewidth=1.5,
            )
            ax.annotate(
                f"{row[metric]:.0%}",
                (row[hi], y),
                xytext=(6, 0),
                textcoords="offset points",
                va="center",
                fontsize=8.5,
                color=t["ink2"],
            )
        ax.set_title(title, color=t["ink"], fontsize=10.5, loc="left", pad=8)
        ax.set_xlim(*xlim)
        lo_tick = 0.0 if xlim[0] == 0 else xlim[0]
        step = 0.2 if xlim[0] == 0 else 0.05
        ax.set_xticks(
            [round(lo_tick + i * step, 2) for i in range(int(round((1 - lo_tick) / step)) + 1)]
        )
        ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    axes[0].set_yticks(ys, table["rule"].tolist())
    axes[0].tick_params(axis="y", colors=t["ink2"], labelsize=9)
    fig.suptitle(
        "A targeted mgrB check recovers part of the colistin gap",
        x=0.0,
        y=1.06,
        ha="left",
        color=t["ink"],
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.0,
        -0.08,
        "Broth-microdilution colistin labels (n = 1,528; 367 R). Dots: estimate; "
        "bars: 95% Wilson interval. Specificity axis starts at 85%.",
        color=t["muted"],
        fontsize=8.5,
        ha="left",
    )
    fig.subplots_adjust(wspace=0.08)
    fig.savefig(out, dpi=200, facecolor=t["surface"], bbox_inches="tight")
    plt.close(fig)


def transfer_figure(studies: pd.DataFrame, pooled: float, theme: str, out: Path) -> None:
    """Per-held-out-study sensitivity of the unfitted mgrB disruption rule."""
    t = THEMES[theme]
    fig, ax = plt.subplots(figsize=(8.5, 3.6))
    fig.patch.set_facecolor(t["surface"])
    style(ax, t)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", color=t["grid"], linewidth=0.8)
    d = studies.sort_values("sens")
    ys = range(len(d))
    ax.axvline(pooled, color=t["muted"], linewidth=1.5, linestyle=(0, (4, 3)), zorder=1)
    ax.annotate(
        f"pooled {pooled:.0%}",
        (pooled, len(d) - 0.35),
        xytext=(5, 0),
        textcoords="offset points",
        fontsize=8.5,
        color=t["muted"],
    )
    for y, r in zip(ys, d.itertuples(), strict=True):
        ax.plot([0, r.sens], [y, y], color=t["grid"], linewidth=1.5, zorder=1)
        ax.plot(
            r.sens,
            y,
            "o",
            color=t["series"][0],
            markersize=9,
            markeredgecolor=t["surface"],
            markeredgewidth=1.5,
            zorder=3,
        )
        ax.annotate(
            f"{r.sens:.0%}",
            (r.sens, y),
            xytext=(9, 0),
            textcoords="offset points",
            va="center",
            fontsize=8.5,
            color=t["ink2"],
        )
    ax.set_yticks(list(ys), [f"{r.study}  (n={r.n_R} R)" for r in d.itertuples()])
    ax.tick_params(axis="y", colors=t["ink2"], labelsize=9)
    ax.set_xlim(0, 0.85)
    ax.set_ylim(-0.5, len(d) - 0.1)  # headroom above the top study for the pooled label
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    fig.suptitle(
        "The same unfitted rule transfers unevenly between collections",
        x=0.0,
        y=1.04,
        ha="left",
        color=t["ink"],
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.0,
        -0.09,
        "Sensitivity of mgrB disruption alone in each held-out study. Nothing is "
        "fitted, so the spread reflects\ndiffering mechanistic composition of colistin "
        "resistance, not overfitting.",
        color=t["muted"],
        fontsize=8.5,
        ha="left",
    )
    fig.savefig(out, dpi=200, facecolor=t["surface"], bbox_inches="tight")
    plt.close(fig)


def transfer_table(path: Path) -> tuple[pd.DataFrame, float]:
    """Parse per-study disruption-only sensitivity from reports/review/mgrb_robustness.md."""
    text = path.read_text().split("## 2. Leave-one-study-out", 1)[1].split("## 3.", 1)[0]
    rows = []
    for line in text.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 7 or not cells[0].startswith("PRJ"):
            continue
        rows.append(
            {"study": cells[0], "n_R": int(cells[1]), "sens": float(cells[3].rstrip("%")) / 100}
        )
    pooled = float(text.split("Disrupted only (nothing fitted): ")[1].split("%")[0]) / 100
    return pd.DataFrame(rows), pooled


def colistin_table(path: Path) -> pd.DataFrame:
    """Parse the broth-microdilution rules table from reports/review/colistin_mgrb.md."""
    text = path.read_text().split("### colistin: rules with and without the mgrB check", 1)[1]
    text = text.split("\n## ", 1)[0]
    rows = []
    for line in text.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 5 or cells[0] in ("Rule", "---"):
            continue

        def ci(cell: str) -> tuple[float, float, float]:
            est, rng = cell.split(" [")
            lo, hi = rng.rstrip("]").split("-")
            return tuple(float(v.rstrip("%")) / 100 for v in (est, lo, hi))

        se, sl, sh = ci(cells[1])
        sp, pl, ph = ci(cells[2])
        rows.append(
            {
                "rule": cells[0],
                "sensitivity": se,
                "sens_lo": sl,
                "sens_hi": sh,
                "specificity": sp,
                "spec_lo": pl,
                "spec_hi": ph,
            }
        )
    keep = [
        "AMRFinderPlus",
        "AMRFinderPlus + mgrB disrupted",
        "Kleborate",
        "Kleborate + mgrB disrupted",
        "mgrB disrupted alone",
    ]
    return pd.DataFrame(rows).set_index("rule").loc[keep].reset_index()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path("reports/benchmark/results.csv"))
    parser.add_argument("--colistin", type=Path, default=Path("reports/review/colistin_mgrb.md"))
    parser.add_argument(
        "--robustness", type=Path, default=Path("reports/review/mgrb_robustness.md")
    )
    parser.add_argument("--out", type=Path, default=Path("reports/figures"))
    args = parser.parse_args()

    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = FONT
    args.out.mkdir(parents=True, exist_ok=True)
    results = pd.read_csv(args.results)
    colistin = colistin_table(args.colistin)
    studies, pooled = transfer_table(args.robustness)
    for theme in THEMES:
        benchmark_figure(results, theme, args.out / f"benchmark_splits_{theme}.png")
        colistin_figure(colistin, theme, args.out / f"colistin_mgrb_{theme}.png")
        transfer_figure(studies, pooled, theme, args.out / f"colistin_transfer_{theme}.png")
    print(f"wrote {len(THEMES) * 3} figures to {args.out}/")


if __name__ == "__main__":
    main()
