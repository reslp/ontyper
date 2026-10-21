#!/usr/bin/env python3
"""Train scikit-learn models for allele prediction."""

import json
import os
import sys

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier, GradientBoostingClassifier, HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss, precision_score, recall_score, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit

from ontyper import __version__
from ontyper.genotypes.features import FEATURE_COLUMNS, add_training_labels, read_table, read_true_alleles, write_table


MODEL_CHOICES = [
    "HistGradientBoostingClassifier",
    "RandomForestClassifier",
    "ExtraTreesClassifier",
    "GradientBoostingClassifier",
    "LogisticRegression",
]


def make_model(model_name, seed):
    """Create a supported scikit-learn classifier by name."""
    if model_name == "HistGradientBoostingClassifier":
        return HistGradientBoostingClassifier(random_state=seed, class_weight="balanced", max_iter=300)
    if model_name == "RandomForestClassifier":
        return RandomForestClassifier(n_estimators=300, random_state=seed, class_weight="balanced", n_jobs=-1)
    if model_name == "ExtraTreesClassifier":
        return ExtraTreesClassifier(n_estimators=300, random_state=seed, class_weight="balanced", n_jobs=-1)
    if model_name == "GradientBoostingClassifier":
        return GradientBoostingClassifier(random_state=seed)
    if model_name == "LogisticRegression":
        return LogisticRegression(random_state=seed, class_weight="balanced", max_iter=5000)
    raise ValueError("Unsupported model: " + model_name)


def prepare_training_table(path, delimiter=None, true_alleles=None, true_delimiter=None):
    """Read features and attach true-allele labels, or accept a legacy labeled table."""
    df = read_table(path, delimiter=delimiter)
    required = ["sample_id", "locus_id"] + FEATURE_COLUMNS
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError("Training feature table is missing required columns: " + ", ".join(missing))
    if true_alleles is not None:
        true_df = read_true_alleles(true_alleles, delimiter=true_delimiter)
        df = add_training_labels(df, true_df)
    elif "label" not in df.columns:
        raise ValueError("Training requires --true-alleles or a feature table containing a label column")
    df = df.copy()
    for column in FEATURE_COLUMNS:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    if df[FEATURE_COLUMNS].isna().any().any():
        df[FEATURE_COLUMNS] = df[FEATURE_COLUMNS].fillna(df[FEATURE_COLUMNS].median(numeric_only=True))
        df[FEATURE_COLUMNS] = df[FEATURE_COLUMNS].fillna(0)
    df["label"] = pd.to_numeric(df["label"], errors="coerce").fillna(0).astype(int)
    if not df["label"].isin([0, 1]).all():
        raise ValueError("Column label must contain only 0/1 values")
    return df


def split_by_group(df, group_column, train_size, seed):
    """Split rows into train and validation sets without splitting groups."""
    if group_column not in df.columns:
        raise ValueError("Group column not found in training table: " + group_column)
    groups = df[group_column]
    if groups.nunique() < 2:
        raise ValueError("At least two groups are required for group-wise train/validation splitting")
    splitter = GroupShuffleSplit(n_splits=1, train_size=train_size, random_state=seed)
    train_idx, val_idx = next(splitter.split(df[FEATURE_COLUMNS], df["label"], groups=groups))
    return train_idx, val_idx


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


def safe_auc(y_true, probabilities):
    """Compute ROC AUC, returning None when only one class is present."""
    if len(set(y_true)) < 2:
        return None
    return float(roc_auc_score(y_true, probabilities))


def safe_log_loss(y_true, probabilities):
    """Compute log loss with fixed binary labels."""
    clipped = np.clip(probabilities, 1e-7, 1 - 1e-7)
    return float(log_loss(y_true, clipped, labels=[0, 1]))


def add_prediction_ranks(df, probabilities):
    """Add probabilities and within-locus ranks to a prediction frame."""
    out = df[["sample_id", "locus_id", "allele_id", "label"]].copy()
    out["probability"] = probabilities
    out["pred_rank_within_locus"] = (
        out.groupby(["sample_id", "locus_id"])["probability"].rank(method="first", ascending=False).astype(int)
    )
    return out


def top1_accuracy(predictions):
    """Compute fraction of sample/locus groups whose top prediction is a true allele."""
    top = predictions.sort_values("probability", ascending=False).groupby(["sample_id", "locus_id"]).head(1)
    if top.empty:
        return None
    return float((top["label"] == 1).mean())


def topk_recall(predictions, top_k):
    """Compute recall of true alleles among top-k selected candidates per sample/locus."""
    positives = int((predictions["label"] == 1).sum())
    if positives == 0:
        return None
    selected = predictions.sort_values("probability", ascending=False).groupby(["sample_id", "locus_id"]).head(top_k)
    return float((selected["label"] == 1).sum() / positives)


def evaluate_predictions(y_true, probabilities, ranked_predictions, top_k, threshold):
    """Compute row-level and per-locus validation metrics."""
    predicted_labels = (probabilities >= threshold).astype(int)
    return {
        "roc_auc": safe_auc(list(y_true), probabilities),
        "log_loss": safe_log_loss(list(y_true), probabilities),
        "precision_at_threshold": float(precision_score(y_true, predicted_labels, zero_division=0)),
        "recall_at_threshold": float(recall_score(y_true, predicted_labels, zero_division=0)),
        "top1_accuracy": top1_accuracy(ranked_predictions),
        "topk_recall": topk_recall(ranked_predictions, top_k),
    }


def train_model(features_path, model_name="HistGradientBoostingClassifier", group="sample_id", train_size=0.8, seed=42, top_k=6, threshold=0.5, delimiter=None, true_alleles=None, true_delimiter=None):
    """Train a model and return a serializable model bundle plus predictions."""
    df = prepare_training_table(
        features_path, delimiter=delimiter, true_alleles=true_alleles, true_delimiter=true_delimiter
    )
    train_idx, val_idx = split_by_group(df, group, train_size, seed)

    X_train = df.iloc[train_idx][FEATURE_COLUMNS]
    y_train = df.iloc[train_idx]["label"]
    X_val = df.iloc[val_idx][FEATURE_COLUMNS]
    y_val = df.iloc[val_idx]["label"]

    if y_train.nunique() < 2:
        raise ValueError("Training split contains only one label class; choose a different group or train size")

    model = make_model(model_name, seed)
    model.fit(X_train, y_train)

    train_prob = positive_probabilities(model, X_train)
    val_prob = positive_probabilities(model, X_val)
    train_predictions = add_prediction_ranks(df.iloc[train_idx], train_prob)
    val_predictions = add_prediction_ranks(df.iloc[val_idx], val_prob)
    metrics = evaluate_predictions(y_val, val_prob, val_predictions, top_k, threshold)
    metrics.update({
        "model_type": model_name,
        "n_rows": int(len(df)),
        "n_train_rows": int(len(train_idx)),
        "n_validation_rows": int(len(val_idx)),
        "group_column": group,
        "train_size": float(train_size),
        "top_k": int(top_k),
        "threshold": float(threshold),
    })

    bundle = {
        "model": model,
        "feature_cols": FEATURE_COLUMNS,
        "model_type": model_name,
        "ontyper_version": __version__,
        "training_summary": metrics,
    }
    return bundle, train_predictions, val_predictions


def ensure_parent_dir(path):
    """Create the parent directory for an output path when needed."""
    outdir = os.path.dirname(path)
    if outdir:
        os.makedirs(outdir, exist_ok=True)


def write_json(data, path):
    """Write a JSON file with stable formatting."""
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main(args):
    """CLI entry point for model training."""
    try:
        bundle, train_predictions, val_predictions = train_model(
            args.features,
            model_name=args.model,
            group=args.group,
            train_size=args.train_size,
            seed=args.seed,
            top_k=args.top_k,
            threshold=args.prob_threshold,
            delimiter=getattr(args, "delimiter", None),
            true_alleles=getattr(args, "true_alleles", None),
            true_delimiter=getattr(args, "true_delimiter", None),
        )
    except ValueError as error:
        print(error, file=sys.stderr)
        return 1

    ensure_parent_dir(args.modelfile)
    joblib.dump(bundle, args.modelfile)
    print("Wrote model to", args.modelfile, file=sys.stderr)

    if args.metrics_out:
        ensure_parent_dir(args.metrics_out)
        write_json(bundle["training_summary"], args.metrics_out)
        print("Wrote metrics to", args.metrics_out, file=sys.stderr)
    else:
        print(json.dumps(bundle["training_summary"], indent=2, sort_keys=True), file=sys.stderr)

    if args.train_pred_out:
        ensure_parent_dir(args.train_pred_out)
        write_table(train_predictions, args.train_pred_out)
        print("Wrote training predictions to", args.train_pred_out, file=sys.stderr)
    if args.val_pred_out:
        ensure_parent_dir(args.val_pred_out)
        write_table(val_predictions, args.val_pred_out)
        print("Wrote validation predictions to", args.val_pred_out, file=sys.stderr)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Train an allele prediction model")
    parser.add_argument("--features", required=True, help="Feature table (legacy labeled tables also supported)")
    parser.add_argument("--true-alleles", default=None, help="Known allele table; required unless features already contain labels")
    parser.add_argument("--true-delimiter", default=None, help="True allele table delimiter (default: auto-detect)")
    parser.add_argument("--model", choices=MODEL_CHOICES, default="HistGradientBoostingClassifier", help="scikit-learn model type")
    parser.add_argument("--modelfile", required=True, help="Output .joblib model bundle")
    parser.add_argument("--group", choices=["sample_id", "locus_id"], default="sample_id", help="Group column for validation split")
    parser.add_argument("--train-size", type=float, default=0.8, help="Training fraction")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--top-k", type=int, default=6, help="Top-k value for validation metrics")
    parser.add_argument("--prob-threshold", type=float, default=0.5, help="Probability threshold for validation metrics")
    parser.add_argument("--delimiter", default=None, help="Feature table delimiter")
    parser.add_argument("--metrics-out", default=None, help="Optional metrics JSON output")
    parser.add_argument("--train-pred-out", default=None, help="Optional training predictions output")
    parser.add_argument("--val-pred-out", default=None, help="Optional validation predictions output")
    sys.exit(main(parser.parse_args()))
