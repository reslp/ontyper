#!/usr/bin/env python3
"""Run genotype identification and model-based allele calling in one workflow."""

import os
import sys
import tempfile
from types import SimpleNamespace

from ontyper.genotypes import identify
from ontyper.genotypes.features import build_feature_table, write_table
from ontyper.genotypes.predict import predict_features


def ensure_parent_dir(path):
    """Create the parent directory for an output path when needed."""
    outdir = os.path.dirname(path)
    if outdir:
        os.makedirs(outdir, exist_ok=True)


def run_identify_to_file(args, outfile):
    """Run the existing genotype identifier and write its count table to a file."""
    identify_args = SimpleNamespace(
        ntimes=args.ntimes,
        demfile=args.demfile,
        indir=args.indir,
        method=args.method,
        delimiter=args.demfile_delimiter,
        relative=False,
        detailed=False,
        outfile=outfile,
    )
    return identify.main(identify_args)


def main(args):
    """CLI entry point for integrated model-based allele calling."""
    with tempfile.TemporaryDirectory() as tmpdir:
        genotype_table = args.genotypes_out or os.path.join(tmpdir, "genotypes.tsv")
        feature_table = args.features_out or os.path.join(tmpdir, "features.csv")

        ensure_parent_dir(genotype_table)
        ensure_parent_dir(feature_table)
        ensure_parent_dir(args.outfile)
        if args.predictions_out:
            ensure_parent_dir(args.predictions_out)

        status = run_identify_to_file(args, genotype_table)
        if status != 0:
            return status

        try:
            features = build_feature_table(
                genotype_table,
                mapping=args.mapping,
                delimiter="\t",
                mapping_delimiter=args.mapping_delimiter,
            )
            write_table(features, feature_table)
            scored, selected = predict_features(
                feature_table,
                args.modelfile,
                top_k=args.top_k,
                prob_threshold=args.prob_threshold,
            )
        except ValueError as error:
            print(error, file=sys.stderr)
            return 1

        if args.predictions_out:
            write_table(scored, args.predictions_out)
            print("Wrote allele scores to", args.predictions_out, file=sys.stderr)
        write_table(selected, args.outfile)
        print("Wrote selected allele calls to", args.outfile, file=sys.stderr)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Identify repeat counts and call alleles with a trained model")
    parser.add_argument("--demfile", required=True, help="Demultiplex file with locus and motif info")
    parser.add_argument("--indir", required=True, help="Directory with demultiplexed FASTA files")
    parser.add_argument("--mapping", required=True, help="name/sample/locus mapping table")
    parser.add_argument("--modelfile", required=True, help="Trained .joblib model bundle")
    parser.add_argument("--outfile", required=True, help="Selected allele call table")
    parser.add_argument("--predictions-out", default=None, help="Optional full scored allele table")
    parser.add_argument("--genotypes-out", default=None, help="Optional genotype count table to keep")
    parser.add_argument("--features-out", default=None, help="Optional feature table to keep")
    parser.add_argument("--ntimes", default="50", help="Maximum repeat count for identify")
    parser.add_argument("--method", default="strict", help="Identify counting method")
    parser.add_argument("--demfile-delimiter", default=",", help="Demultiplex file delimiter for identify")
    parser.add_argument("--mapping-delimiter", default=None, help="Mapping table delimiter")
    parser.add_argument("--top-k", type=int, default=6, help="Maximum alleles selected per sample/locus")
    parser.add_argument("--prob-threshold", type=float, default=0.0, help="Minimum probability for selected alleles")
    sys.exit(main(parser.parse_args()))
