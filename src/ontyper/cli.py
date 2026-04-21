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
from ontyper.genotypes import identify, extract, visualize


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
        description="Demultiplex FASTA reads using primer/tag combinations from a demultiplex file."
    )
    demux_parser.add_argument("-i", "--infile", help="Input FASTA file with reads", required=True)
    demux_parser.add_argument("-d", "--demfile", help="Demultiplex definition file", required=True)
    demux_parser.add_argument("-t", "--threads", help="Number of threads (default: 1)", type=int, default=1)
    demux_parser.add_argument("-m", "--mode", help="Demultiplexing mode: single or dual (default: single)", default="single")
    demux_parser.add_argument("-o", "--outdir", help="Output directory (default: ./)", default="./")
    demux_parser.add_argument("-c", "--cutoff", help="Length cutoff specification", default="")
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
    report_parser.add_argument("-o", "--output", help="Output HTML file (default: ontyper-report.html)", default="ontyper-report.html")
    report_parser.add_argument("--title", help="Report title", default="ontyper report")
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
    identify_parser.add_argument("-m", "--method", help="Counting method: strict or relaxed=X,Y (default: strict)", default="strict")
    identify_parser.add_argument("--delimiter", help="Delimiter in demfile (default: ,)", default=",")
    identify_parser.add_argument("--relative", action="store_true", help="Output relative frequencies", default=False)
    identify_parser.add_argument("--detailed", action="store_true", help="Write genotypes-detailed.tsv with sequence names", default=False)
    identify_parser.add_argument("-o", "--output", help="Write genotype table to this file instead of stdout", default=None)
    identify_parser.set_defaults(func=identify.main)
    
    # ontyper genotypes extract
    extract_parser = geno_subparsers.add_parser(
        "extract",
        help="Extract sequences belonging to a specific genotype bin",
        description="Extract all sequences for a given locus and repeat count from the original FASTA."
    )
    extract_parser.add_argument("-i", "--infile", help="Original input FASTA file", required=True)
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
