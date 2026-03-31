from interface import load_sections, get_section, require_columns
import matplotlib.pyplot as plt
import numpy as np

sections = load_sections("all_stats.tsv")

# scegli cosa plottare
plot_weighted = True      # usa count_per_bp
plot_normalized = True    # usa count_fraction

def plot_distance_curve(df, title, ycol, ylabel, outfile):
    plt.figure(figsize=(11, 4.8))
    plt.loglog(df["bin_center"], df[ycol], linewidth=2)

    plt.xlabel("Genomic distance")
    plt.ylabel(ylabel)
    plt.title(title)

    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(outfile, bbox_inches="tight")
    plt.show()


# BINNED_DISTANCE_GRAPH

# PS
df_ps = get_section(sections, "DIST_PS")
df_ps = require_columns(
    df_ps,
    ["bin_start", "bin_end", "count", "count_per_bp", "count_fraction"],
    "DIST_PS"
)
df_ps["bin_center"] = (df_ps["bin_start"] + df_ps["bin_end"]) / 2

if plot_weighted:
    plot_distance_curve(
        df_ps,
        "P(s) genomic distance weighted by bin width",
        "count_per_bp",
        "Count / bp",
        "P(s)_distance_plot_weighted.pdf"
    )

if plot_normalized:
    plot_distance_curve(
        df_ps,
        "P(s) genomic distance normalized",
        "count_fraction",
        "Fraction of total",
        "P(s)_distance_plot_normalized.pdf"
    )


# FF
df_ff = get_section(sections, "DIST_FF")
df_ff = require_columns(
    df_ff,
    ["bin_start", "bin_end", "count", "count_per_bp", "count_fraction"],
    "DIST_FF"
)
df_ff["bin_center"] = (df_ff["bin_start"] + df_ff["bin_end"]) / 2

if plot_weighted:
    plot_distance_curve(
        df_ff,
        "FF genomic distance weighted by bin width",
        "count_per_bp",
        "Count / bp",
        "FF_distance_plot_weighted.pdf"
    )

if plot_normalized:
    plot_distance_curve(
        df_ff,
        "FF genomic distance normalized",
        "count_fraction",
        "Fraction of total",
        "FF_distance_plot_normalized.pdf"
    )


# FR
df_fr = get_section(sections, "DIST_FR")
df_fr = require_columns(
    df_fr,
    ["bin_start", "bin_end", "count", "count_per_bp", "count_fraction"],
    "DIST_FR"
)
df_fr["bin_center"] = (df_fr["bin_start"] + df_fr["bin_end"]) / 2

if plot_weighted:
    plot_distance_curve(
        df_fr,
        "FR genomic distance weighted by bin width",
        "count_per_bp",
        "Count / bp",
        "FR_distance_plot_weighted.pdf"
    )

if plot_normalized:
    plot_distance_curve(
        df_fr,
        "FR genomic distance normalized",
        "count_fraction",
        "Fraction of total",
        "FR_distance_plot_normalized.pdf"
    )


# RF
df_rf = get_section(sections, "DIST_RF")
df_rf = require_columns(
    df_rf,
    ["bin_start", "bin_end", "count", "count_per_bp", "count_fraction"],
    "DIST_RF"
)
df_rf["bin_center"] = (df_rf["bin_start"] + df_rf["bin_end"]) / 2

if plot_weighted:
    plot_distance_curve(
        df_rf,
        "RF genomic distance weighted by bin width",
        "count_per_bp",
        "Count / bp",
        "RF_distance_plot_weighted.pdf"
    )

if plot_normalized:
    plot_distance_curve(
        df_rf,
        "RF genomic distance normalized",
        "count_fraction",
        "Fraction of total",
        "RF_distance_plot_normalized.pdf"
    )


# RR
df_rr = get_section(sections, "DIST_RR")
df_rr = require_columns(
    df_rr,
    ["bin_start", "bin_end", "count", "count_per_bp", "count_fraction"],
    "DIST_RR"
)
df_rr["bin_center"] = (df_rr["bin_start"] + df_rr["bin_end"]) / 2

if plot_weighted:
    plot_distance_curve(
        df_rr,
        "RR genomic distance weighted by bin width",
        "count_per_bp",
        "Count / bp",
        "RR_distance_plot_weighted.pdf"
    )

if plot_normalized:
    plot_distance_curve(
        df_rr,
        "RR genomic distance normalized",
        "count_fraction",
        "Fraction of total",
        "RR_distance_plot_normalized.pdf"
    )


print("PS sum:", df_ps["count"].sum())
print("FF sum:", df_ff["count"].sum())
print("FR sum:", df_fr["count"].sum())
print("RF sum:", df_rf["count"].sum())
print("RR sum:", df_rr["count"].sum())
#PAIRSTATS


categories = ["UU", "RU", "UR", "WW", "DD", "MU", "MR", "MM", "NM", "NU", "NR", "NN"]
pct_categories = ["pUU", "pRU", "pUR", "pWW", "pDD", "pMU", "pMR", "pMM", "pNM", "pNU", "pNR", "pNN"]

pair_df = get_section(sections, "PAIR_TYPES")
pair_df = require_columns(pair_df, ["run"] + categories + pct_categories, "PAIR_TYPES")

colors = {
    "UU": "#33a02c",
    "RU": "#b2df8a",
    "UR": "#a6cee3",
    "WW": "#1f78b4",
    "DD": "#e31a1c",
    "MU": "#ff7f00",
    "MR": "#fdbf6f",
    "MM": "#ffff33",
    "NM": "#b15928",
    "NU": "#6a3d9a",
    "NR": "#cab2d6",
    "NN": "#ff00ff",
}

y = np.arange(len(pair_df))
left = np.zeros(len(pair_df))

fig, ax = plt.subplots(figsize=(14, 2.8))

for cat in categories:
    vals = pair_df[cat].to_numpy()
    ax.barh(
        y,
        vals,
        left=left,
        color=colors[cat],
        label=cat,
        height=0.55
    )
    left += vals

ax.set_yticks(y)
ax.set_yticklabels(pair_df["run"])
ax.tick_params(axis="y", length=0)
ax.margins(x=0)
spazio_testo = 20 + (15 * len(pair_df))
ax.set_xlabel("Number of pairs", labelpad=spazio_testo)
ax.set_ylabel("")


# lascia spazio sotto
fig.subplots_adjust(bottom=0.32)

# una riga sotto al grafico per ogni run
for i in range(len(pair_df)):
    pct_text = "   ".join(
        f"{cat} {pair_df.iloc[i][pcat]:.1f}%"
        for cat, pcat in zip(categories, pct_categories)
        if pair_df.iloc[i][pcat] > 0
    )

    # coordinate assi: x in frazione asse, y sotto l'asse
    ax.text(
        0.0,
        -0.18,
        pct_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8
    )
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

plt.title("Pair types report")
plt.savefig("pair_types_stacked_with_pct.pdf", bbox_inches="tight")
plt.show()