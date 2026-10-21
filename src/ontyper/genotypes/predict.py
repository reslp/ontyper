#!/usr/bin/env python3
"""Predict true alleles from feature tables using trained scikit-learn models."""

import os
import sys

import joblib
import numpy as np
import pandas as pd

from ontyper.genotypes.features import read_table, write_table


OUTPUT_COLUMNS = ["sample_id", "locus_id", "allele_id", "probability", "pred_rank_within_locus"]


def load_model_bundle(path):
    """Load a joblib model bundle and validate its required fields."""
    bundle = joblib.load(path)
    required = ["model", "feature_cols"]
    missing = [field for field in required if field not in bundle]
    if missing:
        raise ValueError("Model bundle is missing required fields: " + ", ".join(missing))
    return bundle


def prepare_prediction_features(features_path, feature_cols, delimiter=None):
    """Read a feature table and return metadata plus ordered model features."""
    df = read_table(features_path, delimiter=delimiter)
    required = ["sample_id", "locus_id", "allele_id"] + list(feature_cols)
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError("Prediction feature table is missing required columns: " + ", ".join(missing))
    df = df.copy()
    for column in feature_cols:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    if df[feature_cols].isna().any().any():
        df[feature_cols] = df[feature_cols].fillna(df[feature_cols].median(numeric_only=True))
        df[feature_cols] = df[feature_cols].fillna(0)
    return df, df[feature_cols]


def positive_probabilities(model, X):
    """Return probabilities for the positive allele class."""
    if hasattr(model, "predict_proba"):
        probabilities = model.predict_proba(X)
        if probabilities.shape[1] == 1:
            return np.zeros(len(X)) if model.classes_[0] == 0 else np.ones(len(X))
        positive_index = list(model.classes_).index(1)
        return probabilities[:, positive_index]
    scores = model.decision_function(X)
    return 1.0 / (1.0 + np.exp(-scores))


def score_features(df, probabilities):
    """Build the full per-allele prediction table with within-locus ranks."""
    scored = df[["sample_id", "locus_id", "allele_id"]].copy()
    scored["probability"] = probabilities
    scored["pred_rank_within_locus"] = (
        scored.groupby(["sample_id", "locus_id"])["probability"].rank(method="first", ascending=False).astype(int)
    )
    return scored[OUTPUT_COLUMNS].sort_values(["sample_id", "locus_id", "pred_rank_within_locus"]).reset_index(drop=True)


def select_alleles(scored, top_k=None, prob_threshold=0.0):
    """Select predicted alleles by probability threshold and optional per-locus top-k cap."""
    selected = scored.copy()
    if prob_threshold and prob_threshold > 0:
        selected = selected[selected["probability"] >= prob_threshold]
    selected = selected.sort_values(["sample_id", "locus_id", "probability"], ascending=[True, True, False])
    if top_k is not None and top_k > 0:
        selected = selected.groupby(["sample_id", "locus_id"], group_keys=False).head(top_k)
    return selected[OUTPUT_COLUMNS].reset_index(drop=True)


def predict_features(features_path, modelfile, top_k=None, prob_threshold=0.0, delimiter=None):
    """Load a model and score a feature table."""
    bundle = load_model_bundle(modelfile)
    df, X = prepare_prediction_features(features_path, bundle["feature_cols"], delimiter=delimiter)
    probabilities = positive_probabilities(bundle["model"], X)
    scored = score_features(df, probabilities)
    selected = select_alleles(scored, top_k=top_k, prob_threshold=prob_threshold)
    return scored, selected


def main(args):
    """CLI entry point for allele prediction."""
    try:
        scored, selected = predict_features(
            args.features,
            args.modelfile,
            top_k=args.top_k,
            prob_threshold=args.prob_threshold,
            delimiter=getattr(args, "delimiter", None),
        )
    except ValueError as error:
        print(error, file=sys.stderr)
        return 1

    if args.predictions_out:
        outdir = os.path.dirname(args.predictions_out)
        if outdir:
            os.makedirs(outdir, exist_ok=True)
        write_table(scored, args.predictions_out)
        print("Wrote allele scores to", args.predictions_out, file=sys.stderr)

    outdir = os.path.dirname(args.outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    write_table(selected, args.outfile)
    print("Wrote selected allele calls to", args.outfile, file=sys.stderr)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Predict alleles with a trained model")
    parser.add_argument("--features", required=True, help="Prediction feature table")
    parser.add_argument("--modelfile", required=True, help="Trained .joblib model bundle")
    parser.add_argument("--predictions-out", default=None, help="Optional full scored allele table")
    parser.add_argument("--outfile", required=True, help="Selected allele call table")
    parser.add_argument("--top-k", type=int, default=6, help="Maximum alleles selected per sample/locus")
    parser.add_argument("--prob-threshold", type=float, default=0.0, help="Minimum probability for selected alleles")
    parser.add_argument("--delimiter", default=None, help="Feature table delimiter")
    sys.exit(main(parser.parse_args()))
