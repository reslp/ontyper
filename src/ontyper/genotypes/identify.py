#!/usr/bin/env python3
"""
ontyper.genotypes.identify - Identify microsatellite repeat genotypes

This module counts repeat motifs in demultiplexed FASTA files to identify genotypes.
"""

import sys
import glob
import os


def parse_method(method):
    """Parse the method string and return method name and parameters."""
    if method == "strict":
        return "strict", 0, 0
    if method.startswith("relaxed="):
        settings = method.split("=", 1)[1].split(",")
        if len(settings) != 2:
            print("Method relaxed requires two values: relaxed=X,Y", file=sys.stderr)
            sys.exit(1)
        try:
            per_motif = int(settings[0])
            total = int(settings[1])
        except ValueError:
            print("Method relaxed requires integer values: relaxed=X,Y", file=sys.stderr)
            sys.exit(1)
        return "relaxed", per_motif, total
    print("Unknown method:", method, file=sys.stderr)
    sys.exit(1)


def count_repeats_strict(sequence, motifes, motif):
    """Count repeats using exact string matching."""
    repeat_count = 0
    for m in motifes:
        if m in sequence:
            repeat_count = int(len(m) / len(motif))
            continue
        else:
            break
    return repeat_count


def count_mismatches(seq1, seq2):
    """Count character mismatches between two equal-length strings."""
    mismatches = 0
    for i in range(0, len(seq1)):
        if seq1[i] != seq2[i]:
            mismatches += 1
    return mismatches


def count_repeats_relaxed(sequence, motif, ntimes, allowed_mismatches_per_motif, allowed_total_mismatches):
    """Count repeats allowing mismatches."""
    motiflen = len(motif)
    best_repeat_count = 0
    sequence = sequence.strip()
    for start in range(0, len(sequence) - motiflen + 1):
        repeat_count = 0
        total_mismatches = 0
        pos = start
        while pos + motiflen <= len(sequence) and repeat_count < ntimes - 1:
            chunk = sequence[pos:pos+motiflen]
            mismatches = count_mismatches(chunk, motif)
            if mismatches > allowed_mismatches_per_motif:
                break
            if total_mismatches + mismatches > allowed_total_mismatches:
                break
            total_mismatches += mismatches
            repeat_count += 1
            pos += motiflen
        if repeat_count > best_repeat_count:
            best_repeat_count = repeat_count
        if best_repeat_count == ntimes - 1:
            break
    return best_repeat_count


def init_countdict(motifes):
    """Initialize count dictionary for all possible repeat counts."""
    countdict = {}
    countdict[0] = 0
    for m in motifes:
        countdict[float(len(m)) / 4] = 0
    return countdict


def write_detailed_output(detailed_bins, outfile):
    """Write detailed genotype output with sequence names per bin."""
    with open(outfile, "w") as out:
        out.write("locus-name\tnumber-of-repeats\tsequence-names\n")
        for locus_name, repeat_bins in detailed_bins.items():
            for repeat_count in sorted(repeat_bins.keys()):
                seqnames = repeat_bins[repeat_count]
                if not seqnames:
                    continue
                out.write(locus_name + "\t" + str(repeat_count) + "x\t" + ",".join(seqnames) + "\n")


def main(args):
    """
    Main entry point for genotype identification.
    
    Args:
        args: Namespace object with demfile, indir, ntimes, method, delimiter, relative, detailed
    """
    method_name, allowed_mismatches_per_motif, allowed_total_mismatches = parse_method(args.method)
    
    fastafiles = glob.glob(args.indir + "/*.fas")
    print("Directory", args.indir, "contains", len(fastafiles), "fasta files (*.fas)", file=sys.stderr)
    print("Using method:", args.method, file=sys.stderr)
    
    ntimes = int(args.ntimes) + 1
    detailed_bins = {}
    print("locus-sample\t", "\t".join([str(n)+"x" for n in range(0, ntimes)]), "\tnsequences")
    
    with open(args.demfile, "r", encoding='utf-8-sig') as demfile:
        for line in demfile:
            f = line.lstrip().split(args.delimiter)[0]
            motif = line.strip().split(args.delimiter)[-1]
            print("Sample:", f, "Motif:", motif, file=sys.stderr)
            nseqs = 0
            file = ""
            for fname in fastafiles:
                if f in fname:
                    file = fname
                    break
            if not os.path.isfile(file):
                print("File", file, "not found. No sequences for this?", file=sys.stderr)
                continue
            
            motifes = []
            for i in range(1, ntimes):
                motifes.append(motif * i)
            
            countdict = init_countdict(motifes)
            seqname = ""
            with open(file) as fh:
                for line in fh:
                    if ">" in line:
                        nseqs += 1
                        seqname = line.strip().lstrip(">")
                        continue
                    if method_name == "strict":
                        repeat_count = count_repeats_strict(line, motifes, motif)
                    else:
                        repeat_count = count_repeats_relaxed(line, motif, ntimes, allowed_mismatches_per_motif, allowed_total_mismatches)
                    countdict[repeat_count] += 1
                    if args.detailed:
                        locus_name = os.path.basename(file.split(".fa")[0])
                        if locus_name not in detailed_bins:
                            detailed_bins[locus_name] = {}
                        if repeat_count not in detailed_bins[locus_name]:
                            detailed_bins[locus_name][repeat_count] = []
                        detailed_bins[locus_name][repeat_count].append(seqname)
            
            outstr = os.path.basename(file.split(".fa")[0])
            for el in countdict.keys():
                outstr += "\t"
                if args.relative:
                    outstr += str(round((countdict[el] / nseqs), 2))
                else:
                    outstr += str(countdict[el])
            outstr += "\t" + str(nseqs)
            print(outstr)
    
    if args.detailed:
        detailed_outfile = "genotypes-detailed.tsv"
        write_detailed_output(detailed_bins, detailed_outfile)
        print("Wrote detailed genotype summary to", detailed_outfile, file=sys.stderr)
    
    return 0


if __name__ == "__main__":
    # Allow running this module directly for testing
    import argparse
    parser = argparse.ArgumentParser(description="Identify microsatellite genotypes")
    parser.add_argument("-n", "--ntimes", default="50")
    parser.add_argument("-d", "--demfile", required=True)
    parser.add_argument("-i", "--indir", required=True)
    parser.add_argument("-m", "--method", default="strict")
    parser.add_argument("--delimiter", default=",")
    parser.add_argument("--relative", action="store_true", default=False)
    parser.add_argument("--detailed", action="store_true", default=False)
    args = parser.parse_args()
    sys.exit(main(args))
