from pathlib import Path
from io import StringIO
import pandas as pd


def load_sections(file_path: str | Path, sep: str = "\t") -> dict[str, pd.DataFrame]:
    path = Path(file_path)

    if not path.exists():
        raise FileNotFoundError(f"File non trovato: {path}")
    if not path.is_file():
        raise ValueError(f"Il path non è un file: {path}")

    sections = {}
    current_section = None
    current_lines = []

    with open(path, "r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")

            if not line.strip():
                continue

            if line.startswith("#"):
                if current_section is not None and current_lines:
                    block_text = "\n".join(current_lines)
                    sections[current_section] = pd.read_csv(StringIO(block_text), sep=sep)

                current_section = line[1:].strip()
                current_lines = []
            else:
                current_lines.append(line)

    if current_section is not None and current_lines:
        block_text = "\n".join(current_lines)
        sections[current_section] = pd.read_csv(StringIO(block_text), sep=sep)

    return sections


def get_section(sections: dict[str, pd.DataFrame], section_name: str) -> pd.DataFrame:
    if section_name not in sections:
        raise KeyError(
            f"Sezione {section_name!r} non trovata. "
            f"Sezioni disponibili: {list(sections.keys())}"
        )
    return sections[section_name].copy()


def require_columns(df: pd.DataFrame, required_cols: list[str], section_name: str) -> pd.DataFrame:
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"Nella sezione {section_name!r} mancano le colonne {missing}. "
            f"Colonne trovate: {list(df.columns)}"
        )
    return df[required_cols].copy()