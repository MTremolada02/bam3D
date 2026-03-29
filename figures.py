from interface import load_sections, get_section, require_columns
import matplotlib.pyplot as plt
import numpy as np

sections = load_sections("all_stats.tsv")

#BINNED_DISTANCE_GRAPH
#PS
df = get_section(sections, "DIST_PS")
df = require_columns(df, ["bin_start", "bin_end", "count"], "DIST_PS")

df["bin_center"] = (df["bin_start"] + df["bin_end"]) / 2

plt.figure(figsize=(11, 4.8))
plt.loglog(df["bin_center"], df["count"], linewidth=2)

plt.xlabel("Genomic distance")
plt.ylabel("Count")

plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("P(s)_distance_plot.pdf", bbox_inches="tight")
plt.show()

#FF
df = get_section(sections, "DIST_FF")
df = require_columns(df, ["bin_start", "bin_end", "count"], "DIST_FF")

df["bin_center"] = (df["bin_start"] + df["bin_end"]) / 2

plt.figure(figsize=(11, 4.8))
plt.loglog(df["bin_center"], df["count"], linewidth=2)

plt.xlabel("Genomic distance")
plt.ylabel("Count")

plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("FF_distance_plot.pdf", bbox_inches="tight")
plt.show()

#FR
df = get_section(sections, "DIST_FR")
df = require_columns(df, ["bin_start", "bin_end", "count"], "DIST_FR")

df["bin_center"] = (df["bin_start"] + df["bin_end"]) / 2

plt.figure(figsize=(11, 4.8))
plt.loglog(df["bin_center"], df["count"], linewidth=2)

plt.xlabel("Genomic distance")
plt.ylabel("Count")

plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("FR_distance_plot.pdf", bbox_inches="tight")
plt.show()

#RF
df = get_section(sections, "DIST_RF")
df = require_columns(df, ["bin_start", "bin_end", "count"], "DIST_RF")

df["bin_center"] = (df["bin_start"] + df["bin_end"]) / 2

plt.figure(figsize=(11, 4.8))
plt.loglog(df["bin_center"], df["count"], linewidth=2)

plt.xlabel("Genomic distance")
plt.ylabel("Count")

plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("RF_distance_plot.pdf", bbox_inches="tight")
plt.show()

#RR
df = get_section(sections, "DIST_RR")
df = require_columns(df, ["bin_start", "bin_end", "count"], "DIST_RR")

df["bin_center"] = (df["bin_start"] + df["bin_end"]) / 2

plt.figure(figsize=(11, 4.8))
plt.loglog(df["bin_center"], df["count"], linewidth=2)

plt.xlabel("Genomic distance")
plt.ylabel("Count")

plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("RR_distance_plot.pdf", bbox_inches="tight")
plt.show()

#PAIRSTATS
pair_df = get_section(sections, "PAIR_TYPES")
pair_df = require_columns(
    pair_df,
    ["run", "UU", "RU", "UR", "WW", "DD", "MU", "MR", "MM", "NM", "NU", "NR", "NN"],
    "PAIR_TYPES"
)


categories = ["UU", "RU", "UR", "WW", "DD", "MU", "MR", "MM", "NM", "NU", "NR", "NN"]
pct_categories = ["pUU", "pRU", "pUR", "pWW", "pDD", "pMU", "pMR", "pMM", "pNM", "pNU", "pNR", "pNN"]

df = get_section(sections, "PAIR_TYPES")
df = require_columns(df, ["run"] + categories + pct_categories, "PAIR_TYPES")

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

y = np.arange(len(df))
left = np.zeros(len(df))

fig, ax = plt.subplots(figsize=(14, 2.8))

for cat in categories:
    vals = df[cat].to_numpy()
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
ax.set_yticklabels(df["run"])
ax.tick_params(axis="y", length=0)
ax.margins(x=0)
spazio_testo = 20 + (15 * len(df))
ax.set_xlabel("Number of pairs", labelpad=spazio_testo)
ax.set_ylabel("")


# lascia spazio sotto
fig.subplots_adjust(bottom=0.32)

# una riga sotto al grafico per ogni run
for i in range(len(df)):
    pct_text = "   ".join(
        f"{cat} {df.iloc[i][pcat]:.1f}%"
        for cat, pcat in zip(categories, pct_categories)
        if df.iloc[i][pcat] > 0
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
plt.savefig("pair_types_stacked_with_pct.pdf", bbox_inches="tight")
plt.show()