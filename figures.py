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



#GRAFO NODI E ARCHI
nodes = get_section(sections, "GRAPH_NODES")
nodes = require_columns(nodes, ["id", "name", "length", "intra_count"], "GRAPH_NODES")

edges = get_section(sections, "GRAPH_EDGES")
edges = require_columns(edges, ["source", "target", "inter_count"], "GRAPH_EDGES")

# filtro opzionale: tieni solo nodi con intra_count > 0 oppure lunghi abbastanza
nodes = nodes[nodes["intra_count"] > 0].copy()

valid_ids = set(nodes["id"])
edges = edges[edges["source"].isin(valid_ids) & edges["target"].isin(valid_ids)].copy()

# filtro opzionale sugli archi
edges = edges[edges["inter_count"] >= 5].copy()

G = nx.Graph()

for _, row in nodes.iterrows():
    G.add_node(
        int(row["id"]),
        name=row["name"],
        length=int(row["length"]),
        intra_count=int(row["intra_count"])
    )

for _, row in edges.iterrows():
    G.add_edge(
        int(row["source"]),
        int(row["target"]),
        weight=int(row["inter_count"])
    )

# dimensione nodi: log dei contatti intra
node_sizes = []
for n in G.nodes():
    intra = G.nodes[n]["intra_count"]
    node_sizes.append(50 + 80 * np.log1p(intra))

# spessore archi: log dei contatti inter
edge_widths = []
for u, v in G.edges():
    w = G[u][v]["weight"]
    edge_widths.append(0.5 + 1.5 * np.log1p(w))

pos = nx.spring_layout(G, seed=42, k=0.8)

plt.figure(figsize=(12, 12))
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.4)
nx.draw_networkx_nodes(G, pos, node_size=node_sizes, alpha=0.8)

# etichette solo per i nodi più grandi
labels = {}
for n in G.nodes():
    if G.nodes[n]["intra_count"] >= 50:
        labels[n] = G.nodes[n]["name"]

nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)

plt.title("Reference contact graph")
plt.axis("off")
plt.tight_layout()
plt.savefig("reference_contact_graph.pdf", bbox_inches="tight")
plt.show()