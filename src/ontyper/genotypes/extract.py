#!/usr/bin/env python3
"""
ontyper.genotypes.extract - Extract sequences for a specific genotype bin

This module extracts all sequences belonging to a specific locus and repeat count
from the original FASTA or FASTQ file using the genotypes-detailed.tsv file.
"""

import sys

from ontyper.sequence_io import iter_sequence_records


def read_sequence_names(genotype_file, locus_name, repeat_name):
    """
    Read sequence names for a specific locus and repeat bin from the detailed TSV.
    
    Args:
        genotype_file: Path to genotypes-detailed.tsv
        locus_name: Locus name to match
        repeat_name: Repeat bin to match (e.g., "20x")
    
    Returns:
        Set of sequence names
    """
    sequence_names = set()
    with open(genotype_file, "r", encoding="utf-8-sig") as infile:
        for line in infile:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("locus-name\t"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            if parts[0] != locus_name:
                continue
            if parts[1] != repeat_name:
                continue
            for name in parts[2].split(","):
                name = name.strip()
                if name:
                    sequence_names.add(name)
    return sequence_names


def write_matching_sequences(sequence_file, sequence_names, outfile):
    """
    Write matching sequences from FASTA/FASTQ input to output FASTA.
    
    Args:
        sequence_file: Input FASTA/FASTQ file path
        sequence_names: Set of sequence names to extract
        outfile: Output FASTA file path
    
    Returns:
        Number of sequences written
    """
    written = 0
    with open(outfile, "w") as out:
        for record in iter_sequence_records(sequence_file):
            if record.description not in sequence_names:
                continue
            out.write(">" + record.description + "\n")
            out.write(str(record.seq) + "\n")
            written += 1
    return written


def main(args):
    """
    Main entry point for sequence extraction.
    
    Args:
        args: Namespace object with infile, genotypes, locus, number, outfile
    """
    sequence_names = read_sequence_names(args.genotypes, args.locus, args.number)
    if not sequence_names:
        print("No sequence names found for locus", args.locus, "and genotype", args.number, file=sys.stderr)
        return 1

    outfile = args.outfile
    if not outfile:
        outfile = args.locus + "_" + args.number + "_sequences.fas"

    nwritten = write_matching_sequences(args.infile, sequence_names, outfile)
    print("Wrote", nwritten, "sequences to", outfile, file=sys.stderr)

    if nwritten == 0:
        print("No matching sequences were found in", args.infile, file=sys.stderr)
        return 1
    
    return 0


if __name__ == "__main__":
    # Allow running this module directly for testing
    import argparse
    parser = argparse.ArgumentParser(description="Extract sequences for a genotype bin")
    parser.add_argument("-i", "--infile", required=True, help="Input FASTA/FASTQ file")
    parser.add_argument("-g", "--genotypes", required=True)
    parser.add_argument("-l", "--locus", required=True)
    parser.add_argument("-n", "--number", required=True)
    parser.add_argument("-o", "--outfile", default=None)
    args = parser.parse_args()
    sys.exit(main(args))
