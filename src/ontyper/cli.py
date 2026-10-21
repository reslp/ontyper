#!/usr/bin/env python3
"""
ontyper CLI - Main entry point

Command structure:
    ontyper demultiplex [options]
    ontyper report [options]
    ontyper genotypes identify [options]
    ontyper genotypes extract [options]
    ontyper genotypes visualize [options]
"""

import sys
import argparse

from ontyper import __version__
from ontyper import demultiplex, report
from ontyper.genotypes import identify, extract, visualize, features, train, predict, call


def main(args=None):
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="ontyper",
        description="ONT microsatellite genotyping toolkit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  ontyper demultiplex -i reads.fasta -d demfile.txt -o output/
  ontyper report -d demultiplexed/ -g genotypes.tsv -o report.html
  ontyper genotypes identify -d demfile.txt -i demultiplexed/ --detailed
  ontyper genotypes extract -i reads.fasta -g detailed.tsv -l Locus1 -n 20x
  ontyper genotypes visualize -g genotypes.tsv -l Locus1
        """
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # -------------------------------------------------------------------------
    # ontyper demultiplex
    # -------------------------------------------------------------------------
    demux_parser = subparsers.add_parser(
        "demultiplex",
        help="Demultiplex ONT reads by primer and tag sequences",
        description="Demultiplex FASTA or FASTQ reads using primer/tag combinations from a demultiplex file."
    )
    demux_parser.add_argument("-i", "--infile", help="Input FASTA/FASTQ file with reads (.gz supported)", required=True)
    demux_parser.add_argument("-d", "--demfile", help="Demultiplex definition file", required=True)
    demux_parser.add_argument("-t", "--threads", help="Number of threads (default: 1)", type=int, default=1)
    demux_parser.add_argument("-m", "--mode", help="Demultiplex mode. Choices: single(default), dual. Single searches for tag+primer together, dual searches for primer first and then seperately for tag. Dual allows fine grained specification of nucleotide mismatches.", default="single")
    demux_parser.add_argument("-o", "--outdir", help="Output directory (default: ./)", default="./")
    demux_parser.add_argument("-c", "--cutoff", help="Length cutoff specification. Choices: global,range. Examples: -c global=100,700 (keep reads between 100 and 700bp), -r range=300,50 (keep reads between 250 and 350 bp). Default: keep all reads.", default="")
    demux_parser.add_argument("--primer-mismatches", type=int, default=5, help="Allowed primer edit distance in both modes (default: 5)")
    demux_parser.add_argument("--tag-mismatches", type=int, default=0, help="Allowed tag edit distance in dual mode (default: 0)")
    demux_parser.add_argument("--single-mismatches", type=int, default=5, help="Allowed combined-query edit distance in single mode (default: 5)")
    demux_parser.set_defaults(func=demultiplex.main)
    
    # -------------------------------------------------------------------------
    # ontyper report
    # -------------------------------------------------------------------------
    report_parser = subparsers.add_parser(
        "report",
        help="Generate a single-file HTML summary report",
        description="Create an interactive HTML report from demultiplexed reads and genotype results."
    )
    report_parser.add_argument("-d", "--demultiplexed", help="Directory with demultiplexed FASTA files", default=None)
    report_parser.add_argument("-g", "--genotypes", help="genotypes.tsv file", default=None)
    report_parser.add_argument("-p", "--predictions", help="Selected allele calls table", default=None)
    report_parser.add_argument("--true-alleles", help="Known allele table in the same format used for training", default=None)
    report_parser.add_argument("-o", "--outfile", help="Output HTML file (default: ontyper-report.html)", default="ontyper-report.html")
    report_parser.set_defaults(func=report.main)

    # -------------------------------------------------------------------------
    # ontyper genotypes (subcommand group)
    # -------------------------------------------------------------------------
    geno_parser = subparsers.add_parser(
        "genotypes",
        help="Genotype identification and extraction commands",
        description="Commands for microsatellite genotype identification and sequence extraction."
    )
    geno_subparsers = geno_parser.add_subparsers(dest="genotypes_command", help="Genotype subcommands")
    
    # ontyper genotypes identify
    identify_parser = geno_subparsers.add_parser(
        "identify",
        help="Identify microsatellite repeat genotypes from demultiplexed reads",
        description="Count repeat motifs in demultiplexed FASTA files to identify genotypes."
    )
    identify_parser.add_argument("-n", "--ntimes", help="Maximum repeat count to search for (default: 50)", default="50")
    identify_parser.add_argument("-d", "--demfile", help="Demultiplex file with locus and motif info", required=True)
    identify_parser.add_argument("-i", "--indir", help="Directory with demultiplexed FASTA files", required=True)
    identify_parser.add_argument("-m", "--method", help="Counting method. Options: strict (default),relaxed=X,Y. Examples: -m strict (allow no mismatches in repeat stretch). -m relaxed=1,2 (allow one mismatch per motif and two mismatches in whole sequence)", default="strict")
    identify_parser.add_argument("--delimiter", help="Delimiter in demfile (default: ,)", default=",")
    identify_parser.add_argument("--relative", action="store_true", help="Output relative frequencies", default=False)
    identify_parser.add_argument("--detailed", action="store_true", help="Write genotypes-detailed.tsv with sequence names", default=False)
    identify_parser.add_argument("-o", "--outfile", help="Write genotype table to this file instead of stdout", default=None)
    identify_parser.set_defaults(func=identify.main)
    
    # ontyper genotypes extract
    extract_parser = geno_subparsers.add_parser(
        "extract",
        help="Extract sequences belonging to a specific genotype bin",
        description="Extract all sequences for a given locus and repeat count from the original FASTA/FASTQ input."
    )
    extract_parser.add_argument("-i", "--infile", help="Original input FASTA/FASTQ file (.gz supported)", required=True)
    extract_parser.add_argument("-g", "--genotypes", help="genotypes-detailed.tsv file", required=True)
    extract_parser.add_argument("-l", "--locus", help="Locus name to extract", required=True)
    extract_parser.add_argument("-n", "--number", help="Repeat bin to extract (e.g., 20x)", required=True)
    extract_parser.add_argument("-o", "--outfile", help="Output FASTA file (default: <locus>_<number>_sequences.fas)", default=None)
    extract_parser.set_defaults(func=extract.main)

    # ontyper genotypes visualize
    visualize_parser = geno_subparsers.add_parser(
        "visualize",
        help="Visualize genotype count distributions",
        description="Create a publication-ready figure for one locus/sample from a genotype count table."
    )
    visualize_parser.add_argument("-g", "--genotypes", help="genotypes.tsv file", required=True)
    visualize_parser.add_argument("-l", "--locus", help="Locus/sample name to plot", required=True)
    visualize_parser.add_argument("-o", "--outfile", help="Output figure path (default: <locus>_genotypes.pdf)", default=None)
    visualize_parser.set_defaults(func=visualize.main)

    # ontyper genotypes features
    features_parser = geno_subparsers.add_parser(
        "features",
        help="Build allele prediction feature tables",
        description="Convert genotype count tables to machine-learning feature tables."
    )
    features_parser.add_argument("--genotypes", required=True, help="Genotype count table")
    features_parser.add_argument("--mapping", default=None, help="Optional name/sample/locus mapping table")
    features_parser.add_argument("-o", "--outfile", required=True, help="Output feature table")
    features_parser.add_argument("--delimiter", default=None, help="Genotype table delimiter")
    features_parser.add_argument("--mapping-delimiter", default=None, help="Mapping table delimiter")
    features_parser.set_defaults(func=features.main)

    # ontyper genotypes train
    train_parser = geno_subparsers.add_parser(
        "train",
        help="Train an allele prediction model",
        description="Train a scikit-learn model from allele features and known true alleles."
    )
    train_parser.add_argument("--features", required=True, help="Feature table (legacy labeled tables also supported)")
    train_parser.add_argument("--true-alleles", default=None, help="Known allele table; required unless features already contain labels")
    train_parser.add_argument("--true-delimiter", default=None, help="True allele table delimiter (default: auto-detect)")
    train_parser.add_argument("--model", choices=train.MODEL_CHOICES, default="HistGradientBoostingClassifier", help="scikit-learn model type")
    train_parser.add_argument("--modelfile", required=True, help="Output .joblib model bundle")
    train_parser.add_argument("--group", choices=["sample_id", "locus_id"], default="sample_id", help="Group column for validation split")
    train_parser.add_argument("--train-size", type=float, default=0.8, help="Training fraction")
    train_parser.add_argument("--seed", type=int, default=42, help="Random seed")
    train_parser.add_argument("--top-k", type=int, default=6, help="Top-k value for validation metrics")
    train_parser.add_argument("--prob-threshold", type=float, default=0.5, help="Probability threshold for validation metrics")
    train_parser.add_argument("--delimiter", default=None, help="Feature table delimiter")
    train_parser.add_argument("--metrics-out", default=None, help="Optional metrics JSON output")
    train_parser.add_argument("--train-pred-out", default=None, help="Optional training predictions output")
    train_parser.add_argument("--val-pred-out", default=None, help="Optional validation predictions output")
    train_parser.set_defaults(func=train.main)

    # ontyper genotypes predict
    predict_parser = geno_subparsers.add_parser(
        "predict",
        help="Predict true alleles from feature rows",
        description="Score allele feature rows with a trained model and select predicted alleles."
    )
    predict_parser.add_argument("--features", required=True, help="Prediction feature table")
    predict_parser.add_argument("--modelfile", required=True, help="Trained .joblib model bundle")
    predict_parser.add_argument("--predictions-out", default=None, help="Optional full scored allele table")
    predict_parser.add_argument("-o", "--outfile", required=True, help="Selected allele call table")
    predict_parser.add_argument("--top-k", type=int, default=6, help="Maximum alleles selected per sample/locus")
    predict_parser.add_argument("--prob-threshold", type=float, default=0.0, help="Minimum probability for selected alleles")
    predict_parser.add_argument("--delimiter", default=None, help="Feature table delimiter")
    predict_parser.set_defaults(func=predict.main)

    # ontyper genotypes call
    call_parser = geno_subparsers.add_parser(
        "call",
        help="Identify repeats and call alleles with a trained model",
        description="Run genotype identification, feature generation, and allele prediction in one workflow."
    )
    call_parser.add_argument("--demfile", required=True, help="Demultiplex file with locus and motif info")
    call_parser.add_argument("--indir", required=True, help="Directory with demultiplexed FASTA files")
    call_parser.add_argument("--mapping", required=True, help="name/sample/locus mapping table")
    call_parser.add_argument("--modelfile", required=True, help="Trained .joblib model bundle")
    call_parser.add_argument("-o", "--outfile", required=True, help="Selected allele call table")
    call_parser.add_argument("--predictions-out", default=None, help="Optional full scored allele table")
    call_parser.add_argument("--genotypes-out", default=None, help="Optional genotype count table to keep")
    call_parser.add_argument("--features-out", default=None, help="Optional feature table to keep")
    call_parser.add_argument("-n", "--ntimes", default="50", help="Maximum repeat count for identify")
    call_parser.add_argument("-m", "--method", default="strict", help="Identify counting method")
    call_parser.add_argument("--demfile-delimiter", default=",", help="Demultiplex file delimiter for identify")
    call_parser.add_argument("--mapping-delimiter", default=None, help="Mapping table delimiter")
    call_parser.add_argument("--top-k", type=int, default=6, help="Maximum alleles selected per sample/locus")
    call_parser.add_argument("--prob-threshold", type=float, default=0.0, help="Minimum probability for selected alleles")
    call_parser.set_defaults(func=call.main)
    
    # -------------------------------------------------------------------------
    # Parse and dispatch
    # -------------------------------------------------------------------------
    parsed = parser.parse_args(args)
    
    if parsed.command is None:
        parser.print_help()
        return 1
    
    if parsed.command == "genotypes" and parsed.genotypes_command is None:
        geno_parser.print_help()
        return 1
    
    if hasattr(parsed, "func"):
        return parsed.func(parsed)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
