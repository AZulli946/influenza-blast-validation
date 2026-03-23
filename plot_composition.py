#!/usr/bin/env python3
"""Plot influenza subtype composition by location and delivery date."""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from pathlib import Path

# --- Style ---
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

# Subtypes to show individually (rest grouped as "Other")
MAJOR_SUBTYPES = ["H3N2", "H1N1", "H3N8", "H5N1", "H13N6", "H16N3", "H1N2"]

SUBTYPE_COLORS = {
    "H3N2":  "#3574b0",
    "H1N1":  "#e8833a",
    "H3N8":  "#4ba34b",
    "H5N1":  "#d44040",
    "H13N6": "#8e6fbf",
    "H16N3": "#a0714f",
    "H1N2":  "#999999",
    "Other": "#cccccc",
}


def load_data():
    results_dir = Path(__file__).parent / "results"
    df = pd.read_csv(results_dir / "blast_validated.tsv", sep="\t", dtype=str)
    df["n_blast_hits"] = pd.to_numeric(df["n_blast_hits"], errors="coerce")
    df = df[df["n_blast_hits"] > 0].copy()

    # Group minor subtypes
    df["subtype_group"] = df["assigned_subtype"].where(
        df["assigned_subtype"].isin(MAJOR_SUBTYPES), "Other"
    )
    return df


def _short_delivery_label(d):
    """'2025-11-14' -> 'Nov 14'"""
    parts = d.split("-")
    months = ["", "Jan", "Feb", "Mar", "Apr", "May", "Jun",
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    return f"{months[int(parts[1])]} {int(parts[2])}"


def _delivery_sort_key(deliveries):
    """Return sorted deliveries and short labels."""
    deliveries = sorted(deliveries)
    labels = [_short_delivery_label(d) for d in deliveries]
    return deliveries, labels


def plot_per_site(df):
    """Faceted stacked bars: one panel per site, x = delivery date."""
    sites = sorted(df["City"].unique())
    deliveries, labels = _delivery_sort_key(df["delivery_date"].unique())
    all_subtypes = [s for s in MAJOR_SUBTYPES + ["Other"]
                    if s in df["subtype_group"].unique()]

    n_cols = 3
    n_rows = (len(sites) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 3.8 * n_rows))
    axes = axes.flatten()

    x = np.arange(len(deliveries))
    bar_w = 0.7

    for i, site in enumerate(sites):
        ax = axes[i]
        site_df = df[df["City"] == site]

        # Count reads per delivery x subtype
        counts = (site_df.groupby(["delivery_date", "subtype_group"])
                  .size().reset_index(name="n"))
        totals = (site_df.groupby("delivery_date").size()
                  .reindex(deliveries, fill_value=0))

        # Build proportion matrix (deliveries x subtypes)
        props = np.zeros((len(deliveries), len(all_subtypes)))
        for j, d in enumerate(deliveries):
            total = totals.get(d, 0)
            if total == 0:
                continue
            for k, s in enumerate(all_subtypes):
                row = counts[(counts["delivery_date"] == d) &
                             (counts["subtype_group"] == s)]
                if not row.empty:
                    props[j, k] = row["n"].values[0] / total

        # Stacked bars
        bottom = np.zeros(len(deliveries))
        for k, s in enumerate(all_subtypes):
            ax.bar(x, props[:, k], bar_w, bottom=bottom,
                   color=SUBTYPE_COLORS[s], label=s if i == 0 else "",
                   edgecolor="white", linewidth=0.4)
            bottom += props[:, k]

        # Read count annotation inside top of each bar
        for j, d in enumerate(deliveries):
            total = totals.get(d, 0)
            if total > 0:
                ax.text(x[j], 1.03, f"{int(total)}", ha="center", va="bottom",
                        fontsize=7, color="#555555")

        ax.set_title(site, fontsize=13, fontweight="bold", pad=12)
        ax.set_ylim(0, 1.15)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=9)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=50, ha="right", fontsize=8)
        ax.tick_params(axis="y", length=3)
        ax.tick_params(axis="x", length=0)

        # Light horizontal gridlines
        for yl in [0.25, 0.5, 0.75]:
            ax.axhline(yl, color="#e0e0e0", linewidth=0.5, zorder=0)

    # Hide empty panels
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    # Legend
    handles = [plt.Rectangle((0, 0), 1, 1, fc=SUBTYPE_COLORS[s], ec="white")
               for s in all_subtypes]
    fig.legend(handles, all_subtypes, loc="lower center",
               ncol=len(all_subtypes), fontsize=11, frameon=False,
               bbox_to_anchor=(0.5, -0.01), handlelength=1.2, handletextpad=0.5)

    fig.suptitle("Influenza A Subtype Composition by Site",
                 fontsize=18, fontweight="bold", y=1.005)
    fig.text(0.5, 0.993, "Bars show subtype proportion per delivery; numbers = total validated reads",
             ha="center", fontsize=10, color="#777777")

    fig.tight_layout(rect=[0, 0.03, 1, 0.99])
    return fig


def plot_aggregate(df):
    """Two-panel figure: stacked bars (proportion) + read count bars."""
    deliveries, labels = _delivery_sort_key(df["delivery_date"].unique())
    all_subtypes = [s for s in MAJOR_SUBTYPES + ["Other"]
                    if s in df["subtype_group"].unique()]

    counts = (df.groupby(["delivery_date", "subtype_group"])
              .size().reset_index(name="n"))
    totals_s = df.groupby("delivery_date").size().reindex(deliveries, fill_value=0)

    props = np.zeros((len(deliveries), len(all_subtypes)))
    totals = np.zeros(len(deliveries))
    for j, d in enumerate(deliveries):
        total = totals_s.get(d, 0)
        totals[j] = total
        if total == 0:
            continue
        for k, s in enumerate(all_subtypes):
            row = counts[(counts["delivery_date"] == d) &
                         (counts["subtype_group"] == s)]
            if not row.empty:
                props[j, k] = row["n"].values[0] / total

    x = np.arange(len(deliveries))
    bar_w = 0.65

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7),
                                    gridspec_kw={"height_ratios": [3, 1.2]},
                                    sharex=True)

    # Top: stacked proportion bars
    bottom = np.zeros(len(deliveries))
    for k, s in enumerate(all_subtypes):
        ax1.bar(x, props[:, k], bar_w, bottom=bottom,
                color=SUBTYPE_COLORS[s], label=s,
                edgecolor="white", linewidth=0.5)
        bottom += props[:, k]

    ax1.set_ylim(0, 1.0)
    ax1.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax1.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax1.set_ylabel("Subtype Proportion", fontsize=12)
    ax1.set_title("Influenza A Subtype Composition — All Sites Combined",
                  fontsize=16, fontweight="bold", pad=10)
    ax1.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=10,
               frameon=False, ncol=1, handlelength=1, handletextpad=0.4)
    for yl in [0.25, 0.5, 0.75]:
        ax1.axhline(yl, color="#e0e0e0", linewidth=0.5, zorder=0)
    ax1.tick_params(axis="x", length=0)

    # Bottom: total read counts
    ax2.bar(x, totals, bar_w, color="#555555", alpha=0.8)
    ax2.set_ylabel("Total Reads", fontsize=12)
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=40, ha="right", fontsize=10)
    ax2.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: f"{v / 1000:.0f}k" if v >= 1000 else f"{int(v)}"))

    # Label each bar with count
    for j in range(len(deliveries)):
        if totals[j] > 0:
            ax2.text(x[j], totals[j] + max(totals) * 0.02,
                     f"{int(totals[j]):,}", ha="center", va="bottom",
                     fontsize=8, color="#444444")
    ax2.set_ylim(0, max(totals) * 1.15)
    ax2.tick_params(axis="x", length=0)

    fig.tight_layout()
    return fig


def plot_heatmap(df):
    """Heatmap: rows = sites, columns = deliveries, color = dominant subtype,
    text = top 3 subtypes by percentage."""
    deliveries, labels = _delivery_sort_key(df["delivery_date"].unique())
    sites = sorted(df["City"].unique())

    # For each site x delivery: compute subtype breakdown
    grid_dominant = {}   # (site, delivery) -> dominant subtype name
    grid_top3 = {}       # (site, delivery) -> list of (subtype, pct) top 3
    grid_count = {}      # (site, delivery) -> total reads

    for site in sites:
        for d in deliveries:
            subset = df[(df["City"] == site) & (df["delivery_date"] == d)]
            if subset.empty:
                continue
            total = len(subset)
            sc = subset["subtype_group"].value_counts()
            grid_dominant[(site, d)] = sc.index[0]
            grid_count[(site, d)] = total
            # Top 3 subtypes with percentages
            top3 = []
            for s, n in sc.head(3).items():
                pct = n / total
                if pct >= 0.01:  # only show if >= 1%
                    top3.append((s, pct))
            grid_top3[(site, d)] = top3

    # Build color matrix
    all_subtypes = MAJOR_SUBTYPES + ["Other"]
    subtype_to_idx = {s: i for i, s in enumerate(all_subtypes)}
    cmap_colors = [SUBTYPE_COLORS.get(s, "#cccccc") for s in all_subtypes]

    from matplotlib.colors import ListedColormap, BoundaryNorm
    cmap = ListedColormap(cmap_colors)
    bounds = list(range(len(all_subtypes) + 1))
    norm = BoundaryNorm(bounds, cmap.N)

    color_matrix = np.full((len(sites), len(deliveries)), np.nan)
    for i, site in enumerate(sites):
        for j, d in enumerate(deliveries):
            s = grid_dominant.get((site, d))
            if s and s in subtype_to_idx:
                color_matrix[i, j] = subtype_to_idx[s] + 0.5

    fig, ax = plt.subplots(figsize=(18, 10))

    im = ax.imshow(color_matrix, cmap=cmap, norm=norm, aspect="auto",
                   interpolation="nearest")

    # Annotate cells with total reads + top 3 subtypes
    for i, site in enumerate(sites):
        for j, d in enumerate(deliveries):
            if (site, d) not in grid_count:
                ax.text(j, i, "-", ha="center", va="center",
                        fontsize=9, color="#bbbbbb")
                continue

            total = grid_count[(site, d)]
            top3 = grid_top3[(site, d)]

            # Build text: total on first line, then top 3
            lines = [f"n={int(total)}"]
            for s, pct in top3:
                lines.append(f"{s} {pct:.0%}")

            text = "\n".join(lines)
            ax.text(j, i, text, ha="center", va="center",
                    fontsize=6.5, color="white", linespacing=1.2)

    ax.set_xticks(range(len(deliveries)))
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=11)
    ax.set_yticks(range(len(sites)))
    ax.set_yticklabels(sites, fontsize=11)
    ax.tick_params(length=0)

    # Grid lines between cells
    for x in np.arange(-0.5, len(deliveries), 1):
        ax.axvline(x, color="white", linewidth=1.5)
    for y in np.arange(-0.5, len(sites), 1):
        ax.axhline(y, color="white", linewidth=1.5)

    # Legend
    handles = [plt.Rectangle((0, 0), 1, 1, fc=SUBTYPE_COLORS[s])
               for s in all_subtypes if s in df["subtype_group"].unique()]
    leg_labels = [s for s in all_subtypes if s in df["subtype_group"].unique()]
    ax.legend(handles, leg_labels, loc="upper left", bbox_to_anchor=(1.01, 1),
              fontsize=11, frameon=False, handlelength=1.2)

    ax.set_title("Influenza A Subtype Composition by Site and Delivery\n"
                 "Cell color = dominant subtype; text = top subtypes by read proportion",
                 fontsize=14, fontweight="bold", pad=14)

    fig.tight_layout()
    return fig


def main():
    df = load_data()
    out_dir = Path(__file__).parent / "results"

    fig1 = plot_per_site(df)
    p1 = out_dir / "subtype_composition_by_site.png"
    fig1.savefig(p1, dpi=150, bbox_inches="tight")
    print(f"Saved: {p1}")
    plt.close(fig1)

    fig2 = plot_aggregate(df)
    p2 = out_dir / "subtype_composition_aggregate.png"
    fig2.savefig(p2, dpi=150, bbox_inches="tight")
    print(f"Saved: {p2}")
    plt.close(fig2)

    fig3 = plot_heatmap(df)
    p3 = out_dir / "subtype_heatmap.png"
    fig3.savefig(p3, dpi=300, bbox_inches="tight")
    print(f"Saved: {p3}")
    plt.close(fig3)


if __name__ == "__main__":
    main()
