from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import pandas as pd


@dataclass
class BinnedDistanceData:
    input_file: Path
    data: pd.DataFrame

    @property
    def n_rows(self) -> int:
        return int(len(self.data))

    @property
    def total_counts(self) -> int:
        return int(self.data["count"].sum())

    @property
    def min_bin_start(self) -> int:
        return int(self.data["bin_start"].min())

    @property
    def max_bin_end(self) -> int:
        return int(self.data["bin_end"].max())


def load_binned_distance_file(file_path: str | Path, sep: str = "\t") -> BinnedDistanceData:
    path = Path(file_path)

    if not path.exists():
        raise FileNotFoundError(f"File non trovato: {path}")

    if not path.is_file():
        raise ValueError(f"Il path non è un file: {path}")

    df = pd.read_csv(
        path,
        sep=sep,
        comment="#"
    )

    expected_cols = {"bin_start", "bin_end", "count"}
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(
            f"Colonne mancanti nel file {path}. "
            f"Attese: {sorted(expected_cols)}. Trovate: {list(df.columns)}"
        )

    df = df[["bin_start", "bin_end", "count"]].copy()

    df["bin_start"] = pd.to_numeric(df["bin_start"], errors="coerce")
    df["bin_end"] = pd.to_numeric(df["bin_end"], errors="coerce")
    df["count"] = pd.to_numeric(df["count"], errors="coerce")

    if df[["bin_start", "bin_end", "count"]].isna().any().any():
        bad_rows = df[df[["bin_start", "bin_end", "count"]].isna().any(axis=1)]
        raise ValueError(
            "Ci sono valori non numerici nelle colonne "
            "'bin_start', 'bin_end' o 'count'. "
            f"Righe problematiche:\n{bad_rows}"
        )

    df["bin_start"] = df["bin_start"].astype("int64")
    df["bin_end"] = df["bin_end"].astype("int64")
    df["count"] = df["count"].astype("int64")

    if (df["bin_start"] < 0).any() or (df["bin_end"] < 0).any():
        raise ValueError("I bin contengono valori negativi.")

    if (df["count"] < 0).any():
        raise ValueError("La colonna 'count' contiene valori negativi.")

    if (df["bin_end"] < df["bin_start"]).any():
        raise ValueError("Ci sono righe con bin_end < bin_start.")

    df = df.sort_values(["bin_start", "bin_end"]).reset_index(drop=True)

    return BinnedDistanceData(input_file=path, data=df)


def print_binned_distance_summary(bd: BinnedDistanceData) -> None:
    print(f"File: {bd.input_file}")
    print(f"Numero di righe: {bd.n_rows}")
    print(f"Bin minimo: {bd.min_bin_start}")
    print(f"Bin massimo: {bd.max_bin_end}")
    print(f"Somma totale dei count: {bd.total_counts}")


if __name__ == "__main__":
    input_file = "Pair_by_binned_distance.tsv"

    bd = load_binned_distance_file(input_file, sep="\t")
    print_binned_distance_summary(bd)

    print("\nPrime righe:")
    print(bd.data.head())