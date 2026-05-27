from interface import load_sections, get_section, require_columns
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

sections = load_sections("all_stats.tsv")

#BINNED_DISTANCE_GRAPH
#PS
df_ps = get_section(sections, "DIST_PS")
df_ps = require_columns(df_ps, ["bin_start", "bin_end", "count"], "DIST_PS")

df_ps["bin_center"] = (df_ps["bin_start"] + df_ps["bin_end"]) / 2

plt.figure(figsize=(11, 4.8))
plt.loglog(df_ps["bin_center"], df_ps["count"], linewidth=2)
#plt.loglog(df_ps["bin_center"], df_ps["count_per_bp"], linewidth=2)

plt.xlabel("Genomic distance")
plt.ylabel("Count")
plt.title("P(s) genomic distance")

plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("P(s)_distance_plot.pdf", bbox_inches="tight")
plt.show()

#FF
df_ff = get_section(sections, "DIST_FF")
df_ff = require_columns(df_ff, ["bin_start", "bin_end", "count"], "DIST_FF")

df_ff["bin_center"] = (df_ff["bin_start"] + df_ff["bin_end"]) / 2

plt.figure(figsize=(11, 4.8))
plt.loglog(df_ff["bin_center"], df_ff["count"], linewidth=2)

plt.xlabel("Genomic distance")
plt.ylabel("Count")
plt.title("FF genomic distance")

plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("FF_distance_plot.pdf", bbox_inches="tight")
plt.show()

#FR
df_fr = get_section(sections, "DIST_FR")
df_fr = require_columns(df_fr, ["bin_start", "bin_end", "count"], "DIST_FR")
#df_fr = require_columns(df_fr, ["bin_start", "bin_end", "count", "count_per_bp", "count_fraction"], "DIST_FR")

df_fr["bin_center"] = (df_fr["bin_start"] + df_fr["bin_end"]) / 2

plt.figure(figsize=(11, 4.8))
plt.loglog(df_fr["bin_center"], df_fr["count"], linewidth=2)

plt.xlabel("Genomic distance")
plt.ylabel("Count")
plt.title("FR genomic distance")

plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("FR_distance_plot.pdf", bbox_inches="tight")
plt.show()

#RF
df_rf = get_section(sections, "DIST_RF")
df_rf = require_columns(df_rf, ["bin_start", "bin_end", "count"], "DIST_RF")

df_rf["bin_center"] = (df_rf["bin_start"] + df_rf["bin_end"]) / 2

plt.figure(figsize=(11, 4.8))
plt.loglog(df_rf["bin_center"], df_rf["count"], linewidth=2)

plt.xlabel("Genomic distance")
plt.ylabel("Count")
plt.title("RF genomic distance")

plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("RF_distance_plot.pdf", bbox_inches="tight")
plt.show()

#RR
df_rr = get_section(sections, "DIST_RR")
df_rr = require_columns(df_rr, ["bin_start", "bin_end", "count"], "DIST_RR")

df_rr["bin_center"] = (df_rr["bin_start"] + df_rr["bin_end"]) / 2

plt.figure(figsize=(11, 4.8))
plt.loglog(df_rr["bin_center"], df_rr["count"], linewidth=2)

plt.xlabel("Genomic distance")
plt.ylabel("Count")
plt.title("RR genomic distance")

plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("RR_distance_plot.pdf", bbox_inches="tight")
plt.show()



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

# ==============================================================================
# STRAND ORIENTATION STATS (SEPARATOR SUBPLOTS - NORMALIZED TO 100%)
# ==============================================================================

import matplotlib.pyplot as plt
import numpy as np

# Invertito l'ordine: RF viene prima di FR
orient_categories = ["FF", "RF", "FR", "RR"]

# Estraiamo la sezione dal file
orient_df = get_section(sections, "STRAND_ORIENTATION_BY_DISTANCE")
# Filtriamo per escludere la riga 'combined'
orient_df = orient_df[orient_df["distance_bin"] != "combined"].reset_index(drop=True)

# Palette di colori richiesta
orient_colors = {
    "FF": "#e31a1c",  # Rosso
    "RF": "#1f78b4",  # Blu
    "FR": "#33a02c",  # Verde
    "RR": "#6a3d9a",  # Viola
}

num_bins = len(orient_df)

# Creiamo una figura con N subplots verticali (uno sotto l'altro)
# Condividiamo l'asse X così la scala 0-100% è identica per tutti
fig, axes = plt.subplots(nrows=num_bins, ncols=1, figsize=(12, 1.2 * num_bins), sharex=True)

# Se c'è solo un bin, axes non è un array, lo forziamo ad esserlo
if num_bins == 1:
    axes = [axes]

for i, ax in enumerate(axes):
    row = orient_df.iloc[i]
    dist_label = row["distance_bin"]
    run_label = row["run"]
    
    # Recuperiamo i valori percentuali direttamente dalle colonne pFF, pRF... del tuo file
    p_vals = {cat: float(row[f"p{cat}"]) for cat in orient_categories}
    
    left = 0.0
    for cat in orient_categories:
        val = p_vals[cat]
        ax.barh(
            [0],  # Singola barra orizzontale in posizione 0
            val,
            left=left,
            color=orient_colors[cat],
            label=cat if i == 0 else "",  # Mettiamo la legenda solo sul primo grafico
            height=0.45
        )
        left += val
        
    # Configurazione di ogni singolo subplot
    ax.set_yticks([0])
    ax.set_yticklabels([f"{run_label}\n({dist_label})"], fontsize=9)
    ax.set_xlim(0, 100)  # Forza la scala da 0 a 100%
    
    # Rimuoviamo i bordi inutili per farlo pulito
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis="y", length=0)

# Configuriamo l'ultimo asse in basso per mostrare la percentuale
axes[-1].set_xlabel("Percentage (%)", fontsize=10, labelpad=10)
axes[-1].spines['bottom'].set_visible(True)

# Posizioniamo la legenda globale in alto a destra, orizzontale e pulita (segue l'ordine di orient_categories)
axes[0].legend(
    bbox_to_anchor=(1.0, 1.8), 
    loc="upper right", 
    ncol=4, 
    frameon=False,
    fontsize=10
)

plt.suptitle("Strand Orientation Composition by Distance Bin", fontsize=12, y=0.98, weight='bold')

# Aggiusta la disposizione per evitare sovrapposizioni
plt.tight_layout()
plt.savefig("strand_orientation_subplots_normalized.pdf", bbox_inches="tight")
plt.show()