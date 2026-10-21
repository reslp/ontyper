#!/usr/bin/env python3
"""Build machine-learning feature tables for genotype prediction."""

import os
import sys

import numpy as np
import pandas as pd


FEATURE_COLUMNS = [
    "allele_id",
    "allele_reads",
    "allele_fraction",
    "allele_rank",
    "distance_to_top_allele",
    "min_distance_to_other_alleles",
    "entropy",
    "higher_ranked_reads",
    "previous_allele_reads",
    "next_allele_reads",
    "higher_ranked_fraction",
    "previous_allele_fraction",
    "next_allele_fraction",
]

FEATURE_TABLE_COLUMNS = [
    "sample_id",
    "locus_id",
    "allele_id",
    "allele_reads",
    "total_reads",
    "other_reads",
    "allele_fraction",
    "allele_rank",
    "distance_to_top_allele",
    "min_distance_to_other_alleles",
    "entropy",
    "higher_ranked_reads",
    "previous_allele_reads",
    "next_allele_reads",
    "higher_ranked_fraction",
    "previous_allele_fraction",
    "next_allele_fraction",
]

TRAINING_TABLE_COLUMNS = [
    "sample_id",
    "locus_id",
    "allele_id",
    "label",
    "feature_disagreement",
    "allele_reads",
    "total_reads",
    "other_reads",
    "allele_fraction",
    "allele_rank",
    "distance_to_top_allele",
    "min_distance_to_other_alleles",
    "entropy",
    "n_alleles",
    "higher_ranked_reads",
    "previous_allele_reads",
    "next_allele_reads",
    "higher_ranked_fraction",
    "previous_allele_fraction",
    "next_allele_fraction",
]

METADATA_COLUMNS = {
    "sample_id",
    "locus_id",
    "locus-sample",
    "name",
    "nsequences",
    "total_reads",
    "other_reads",
    "label",
    "feature_disagreement",
    "n_alleles",
}


def sniff_delimiter(path, default=","):
    """Infer a simple CSV delimiter from the header line."""
    with open(path, "r", encoding="utf-8-sig") as handle:
        header = handle.readline()
    counts = {delimiter: header.count(delimiter) for delimiter in ["\t", ";", ","]}
    delimiter, count = max(counts.items(), key=lambda item: item[1])
    return delimiter if count > 0 else default


def read_table(path, delimiter=None):
    """Read a delimited table with delimiter inference when needed."""
    sep = delimiter if delimiter is not None else sniff_delimiter(path)
    return pd.read_csv(path, sep=sep)


def write_table(df, path):
    """Write a table as TSV for .tsv paths and CSV otherwise."""
    sep = "\t" if str(path).lower().endswith(".tsv") else ","
    df.to_csv(path, sep=sep, index=False)


def read_genotype_counts(path, delimiter=None):
    """Read a genotype count table produced by ontyper or the legacy scripts."""
    return read_table(path, delimiter=delimiter)


def read_mapping(path, delimiter=None):
    """Read and validate a name-to-sample/locus mapping table."""
    mapping = read_table(path, delimiter=delimiter)
    required = ["name", "sample_id", "locus_id"]
    missing = [column for column in required if column not in mapping.columns]
    if missing:
        raise ValueError("Mapping file is missing required columns: " + ", ".join(missing))
    if mapping["name"].duplicated().any():
        duplicates = sorted(mapping.loc[mapping["name"].duplicated(), "name"].astype(str).unique())
        raise ValueError("Mapping file contains duplicate names: " + ", ".join(duplicates))
    return mapping[required].copy()


def allele_id_from_column(column):
    """Convert an allele-count column name to an integer allele id, or return None."""
    text = str(column).strip()
    if text.startswith("X") and text[1:].isdigit():
        return int(text[1:])
    if text.endswith("x") and text[:-1].isdigit():
        return int(text[:-1])
    if text.isdigit():
        return int(text)
    return None


def allele_columns(df):
    """Return allele-count columns in numeric allele order."""
    columns = []
    for column in df.columns:
        if str(column) in METADATA_COLUMNS:
            continue
        allele_id = allele_id_from_column(column)
        if allele_id is not None:
            columns.append((allele_id, column))
    return [column for _, column in sorted(columns, key=lambda item: item[0])]


def add_mapping(counts_df, mapping_df):
    """Filter genotype rows to mapped names and add sample/locus identities."""
    if "sample_id" in counts_df.columns and "locus_id" in counts_df.columns:
        return counts_df.copy()
    name_column = counts_df.columns[0]
    counts_df = counts_df.rename(columns={name_column: "name"})
    genotype_names = set(counts_df["name"])
    mapping_names = set(mapping_df["name"])
    for names, message in [
        (genotype_names - mapping_names, "skipping genotype rows absent from mapping"),
        (mapping_names - genotype_names, "mapping entries absent from genotype table"),
    ]:
        if names:
            print("Warning: " + message + ": " + ", ".join(sorted(map(str, names))), file=sys.stderr)
    merged = counts_df.merge(mapping_df, on="name", how="inner")
    if merged.empty:
        raise ValueError("No matching sample names between genotype table and mapping file")
    missing = merged.loc[merged["sample_id"].isna() | merged["locus_id"].isna(), "name"].astype(str).tolist()
    if missing:
        raise ValueError("Mapping file has missing sample_id or locus_id for: " + ", ".join(missing))
    return merged


def wide_counts_to_long(counts_df, mapping_df=None):
    """Convert a wide allele-count table to one row per candidate allele."""
    if mapping_df is not None:
        counts_df = add_mapping(counts_df, mapping_df)
    elif "sample_id" not in counts_df.columns or "locus_id" not in counts_df.columns:
        raise ValueError("Genotype table needs sample_id/locus_id columns or a mapping file")

    count_columns = allele_columns(counts_df)
    if not count_columns:
        raise ValueError("No allele count columns found")

    id_columns = ["sample_id", "locus_id"]
    if "nsequences" in counts_df.columns:
        counts_df = counts_df.rename(columns={"nsequences": "total_reads"})
    if "total_reads" not in counts_df.columns:
        counts_df = counts_df.copy()
        counts_df["total_reads"] = counts_df[count_columns].sum(axis=1)
    id_columns.append("total_reads")

    long_df = counts_df.melt(
        id_vars=id_columns,
        value_vars=count_columns,
        var_name="allele_column",
        value_name="allele_reads",
    )
    long_df["allele_id"] = long_df["allele_column"].map(allele_id_from_column).astype(int)
    long_df["allele_reads"] = pd.to_numeric(long_df["allele_reads"], errors="coerce").fillna(0).astype(int)
    long_df["total_reads"] = pd.to_numeric(long_df["total_reads"], errors="coerce").fillna(0).astype(int)
    return long_df[["sample_id", "locus_id", "allele_id", "allele_reads", "total_reads"]].sort_values(
        ["sample_id", "locus_id", "allele_id"]
    ).reset_index(drop=True)


def safe_divide(numerator, denominator):
    """Divide numeric values and return zero where the denominator is zero."""
    numerator = pd.Series(numerator)
    denominator = pd.Series(denominator, index=numerator.index)
    result = numerator.div(denominator.replace(0, np.nan)).fillna(0)
    return result


def min_distances(values):
    """Return the nearest distance from each allele id to another allele id."""
    ids = list(values)
    distances = []
    for allele_id in ids:
        others = [other for other in ids if other != allele_id]
        distances.append(min((abs(allele_id - other) for other in others), default=0))
    return pd.Series(distances, index=values.index)


def group_entropy(group):
    """Compute Shannon entropy from positive allele fractions in one sample/locus group."""
    positive = group.loc[group["allele_fraction"] > 0, "allele_fraction"]
    if positive.empty:
        return 0.0
    return float(-(positive * np.log2(positive)).sum())


def add_group_context_features(df, group_cols):
    """Add neighbor-read and higher-ranked-read features within each sample/locus group."""
    df = df.copy()
    df["higher_ranked_reads"] = 0
    df["previous_allele_reads"] = 0
    df["next_allele_reads"] = 0
    for _, index in df.groupby(group_cols).groups.items():
        group = df.loc[index]
        rank_reads = group.groupby("allele_rank")["allele_reads"].sum().sort_index()
        cumulative_lower_rank_reads = rank_reads.cumsum().shift(fill_value=0)
        read_lookup = dict(zip(group["allele_id"], group["allele_reads"]))
        df.loc[index, "higher_ranked_reads"] = group["allele_rank"].map(cumulative_lower_rank_reads).fillna(0).astype(int)
        df.loc[index, "previous_allele_reads"] = group["allele_id"].map(lambda allele_id: read_lookup.get(allele_id - 1, 0))
        df.loc[index, "next_allele_reads"] = group["allele_id"].map(lambda allele_id: read_lookup.get(allele_id + 1, 0))
    return df


def build_allele_features(long_df):
    """Build model features for each candidate allele from long allele counts."""
    df = long_df.copy()
    group_cols = ["sample_id", "locus_id"]
    df["other_reads"] = df["total_reads"] - df["allele_reads"]
    df["allele_fraction"] = safe_divide(df["allele_reads"], df["total_reads"])
    df["allele_rank"] = df.groupby(group_cols)["allele_reads"].rank(method="min", ascending=False).astype(int)

    top_rows = df.loc[df.groupby(group_cols)["allele_reads"].idxmax(), group_cols + ["allele_id"]]
    top_rows = top_rows.rename(columns={"allele_id": "top_allele_id"})
    df = df.merge(top_rows, on=group_cols, how="left")
    df["distance_to_top_allele"] = (df["allele_id"] - df["top_allele_id"]).abs().astype(int)
    df = df.drop(columns=["top_allele_id"])

    df["min_distance_to_other_alleles"] = df.groupby(group_cols)["allele_id"].transform(min_distances).astype(int)
    entropy = df.groupby(group_cols, group_keys=False).apply(group_entropy).rename("entropy").reset_index()
    df = df.merge(entropy, on=group_cols, how="left")

    df = add_group_context_features(df, group_cols)
    df["higher_ranked_fraction"] = safe_divide(df["higher_ranked_reads"], df["total_reads"])
    df["previous_allele_fraction"] = safe_divide(df["previous_allele_reads"], df["total_reads"])
    df["next_allele_fraction"] = safe_divide(df["next_allele_reads"], df["total_reads"])

    numeric_int_columns = [
        "allele_id",
        "allele_reads",
        "total_reads",
        "other_reads",
        "allele_rank",
        "distance_to_top_allele",
        "min_distance_to_other_alleles",
        "higher_ranked_reads",
        "previous_allele_reads",
        "next_allele_reads",
    ]
    for column in numeric_int_columns:
        df[column] = df[column].fillna(0).astype(int)
    return df[FEATURE_TABLE_COLUMNS].sort_values(group_cols + ["allele_id"]).reset_index(drop=True)


def read_true_alleles(path, delimiter=None):
    """Read a true-allele table and return unique long-format true alleles."""
    df = read_table(path, delimiter=delimiter)
    required = ["sample_id", "locus_id"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError("True allele table is missing required columns: " + ", ".join(missing))
    allele_cols = [column for column in df.columns if str(column).lower().startswith("allele")]
    if "true_allele" in df.columns:
        allele_cols.append("true_allele")
    allele_cols = list(dict.fromkeys(allele_cols))
    if not allele_cols:
        raise ValueError("True allele table has no allele columns")
    long_df = df.melt(
        id_vars=["sample_id", "locus_id"],
        value_vars=allele_cols,
        value_name="true_allele",
    )
    long_df["true_allele"] = pd.to_numeric(long_df["true_allele"], errors="coerce")
    long_df = long_df.dropna(subset=["true_allele"])
    long_df["true_allele"] = long_df["true_allele"].astype(int)
    return long_df[["sample_id", "locus_id", "true_allele"]].drop_duplicates().reset_index(drop=True)


def add_training_labels(features_df, true_alleles_df):
    """Add binary true-allele labels to a feature table."""
    true_keys = set(
        zip(
            true_alleles_df["sample_id"].astype(str),
            true_alleles_df["locus_id"].astype(str),
            true_alleles_df["true_allele"].astype(int),
        )
    )
    df = features_df.copy()
    keys = zip(df["sample_id"].astype(str), df["locus_id"].astype(str), df["allele_id"].astype(int))
    df["label"] = [1 if key in true_keys else 0 for key in keys]
    df["n_alleles"] = df.groupby(["sample_id", "locus_id"])["label"].transform("sum").astype(int)
    df["feature_disagreement"] = False
    return df[TRAINING_TABLE_COLUMNS].reset_index(drop=True)


def build_feature_table(genotypes, mapping=None, delimiter=None, mapping_delimiter=None):
    """Build an unlabeled feature table for training or prediction from input paths."""
    counts = read_genotype_counts(genotypes, delimiter=delimiter)
    mapping_df = read_mapping(mapping, delimiter=mapping_delimiter) if mapping else None
    long_df = wide_counts_to_long(counts, mapping_df=mapping_df)
    return build_allele_features(long_df)


def main(args):
    """CLI entry point for feature table generation."""
    try:
        features = build_feature_table(
            args.genotypes,
            mapping=getattr(args, "mapping", None),
            delimiter=getattr(args, "delimiter", None),
            mapping_delimiter=getattr(args, "mapping_delimiter", None),
        )
    except ValueError as error:
        print(error, file=sys.stderr)
        return 1

    outdir = os.path.dirname(args.outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    write_table(features, args.outfile)
    print("Wrote feature table to", args.outfile, file=sys.stderr)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build genotype prediction feature tables")
    parser.add_argument("--genotypes", required=True, help="Genotype count table")
    parser.add_argument("--mapping", default=None, help="Optional name/sample/locus mapping table")
    parser.add_argument("--outfile", required=True, help="Output feature table")
    parser.add_argument("--delimiter", default=None, help="Genotype table delimiter")
    parser.add_argument("--mapping-delimiter", default=None, help="Mapping table delimiter")
    sys.exit(main(parser.parse_args()))
