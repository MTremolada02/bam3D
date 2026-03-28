from interface import load_binned_distance_file

bd = load_binned_distance_file("Pair_by_binned_distance.tsv", sep="\t")
df = bd.data