"""Annotate a CSV file with ChEBI groundings using Gilda.

Usage
-----
Basic (CURIE only):
    python scripts/annotate_csv_chebi.py input.csv --column chemical_name

With score and matched name:
    python scripts/annotate_csv_chebi.py input.csv --column chemical_name --verbose

Custom output file:
    python scripts/annotate_csv_chebi.py input.csv --column chemical_name --output annotated.csv

With a context column for disambiguation:
    python scripts/annotate_csv_chebi.py input.csv --column chemical_name --context-column sentence
"""

import argparse

import pandas as pd

import gilda


def ground_to_chebi(text, context=None):
    """Ground a single text string to ChEBI, returning curie, name, and score."""
    if not isinstance(text, str):
        return None, None, None
    matches = gilda.ground(text, context=context, namespaces=["CHEBI"])
    if not matches:
        return None, None, None
    top = matches[0]
    return top.term.get_curie(), top.term.entry_name, round(top.score, 4)


def annotate_csv(
    input_path: str,
    column: str,
    output_path: str,
    context_column: str = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Load a CSV, annotate one column against ChEBI, and write the result.

    Parameters
    ----------
    input_path:
        Path to the input CSV file.
    column:
        Name of the column containing chemical entity text to ground.
    output_path:
        Path where the annotated CSV will be saved.
    context_column:
        Optional column containing sentence context for disambiguation.
    verbose:
        If True, also adds ``chebi_name`` and ``score`` columns in addition
        to ``chebi_curie``.
    """
    df = pd.read_csv(input_path)

    if column not in df.columns:
        raise ValueError(
            f"Column '{column}' not found in CSV. Available columns: {list(df.columns)}"
        )
    if context_column and context_column not in df.columns:
        raise ValueError(
            f"Context column '{context_column}' not found in CSV. "
            f"Available columns: {list(df.columns)}"
        )

    print(f"Grounding {len(df)} rows from column '{column}' against ChEBI...")

    if verbose or context_column:
        def _ground_row(row):
            context = row[context_column] if context_column else None
            return pd.Series(ground_to_chebi(row[column], context=context))

        df[["chebi_curie", "chebi_name", "score"]] = df.apply(_ground_row, axis=1)
    else:
        gilda.ground_df(df, source_column=column, target_column="chebi_curie", namespaces=["CHEBI"])

    df.to_csv(output_path, index=False)

    matched = df["chebi_curie"].notna().sum()
    print(f"Matched {matched}/{len(df)} rows ({100 * matched // len(df)}%)")
    print(f"Saved annotated CSV to: {output_path}")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Annotate a CSV column with ChEBI groundings using Gilda."
    )
    parser.add_argument("input", help="Path to the input CSV file.")
    parser.add_argument("--column", required=True, help="Column name containing chemical text.")
    parser.add_argument(
        "--output",
        default=None,
        help="Output CSV path. Defaults to <input>_chebi_annotated.csv.",
    )
    parser.add_argument(
        "--context-column",
        default=None,
        help="Optional column with sentence context for disambiguation.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Also output chebi_name and score columns.",
    )
    args = parser.parse_args()

    output_path = args.output
    if output_path is None:
        stem = args.input.removesuffix(".csv")
        output_path = f"{stem}_chebi_annotated.csv"

    annotate_csv(
        input_path=args.input,
        column=args.column,
        output_path=output_path,
        context_column=args.context_column,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
