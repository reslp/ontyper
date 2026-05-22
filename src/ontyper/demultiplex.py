#!/usr/bin/env python3
"""
ontyper.demultiplex - Demultiplex ONT reads by primer and tag sequences

This module demultiplexes FASTA or FASTQ reads using primer/tag combinations
from a demultiplex file. Supports both single-pass and dual-pass demultiplexing
modes.
"""

import sys
import os
import csv
from datetime import datetime
import multiprocessing

from Bio.Seq import Seq
import edlib

from ontyper.sequence_io import detect_sequence_format, iter_sequence_records


# IUPAC ambiguity correspondences for edlib
CORRESPONDENCES = [
    ("R", "A"), ("R", "G"), ("Y", "C"), ("Y", "T"), ("M", "A"), ("M", "C"),
    ("K", "G"), ("K", "T"), ("S", "G"), ("S", "C"), ("W", "A"), ("W", "T"),
    ("B", "C"), ("B", "G"), ("B", "T"), ("D", "A"), ("D", "G"), ("D", "T"),
    ("H", "A"), ("H", "C"), ("H", "T"), ("V", "A"), ("V", "C"), ("V", "G"),
    ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T")
]

# Translation table for non-standard bases
TRANSLATION_TABLE = str.maketrans({"E": "A", "F": "G", "Q": "C", "P": "T"})

# Global worker data for multiprocessing
WORKER_DATA = {}

# Logging state
LOG_HANDLE = None
ORIGINAL_STDOUT = sys.stdout
ORIGINAL_STDERR = sys.stderr


class TeeOutput:
    """Tee output to multiple streams."""
    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for stream in self.streams:
            stream.write(text)
            stream.flush()
        return len(text)

    def flush(self):
        for stream in self.streams:
            stream.flush()


def start_logging(outdir, timestamp):
    """Start logging to file and stdout."""
    global LOG_HANDLE
    logfile = os.path.join(outdir, "demultiplex-" + timestamp + ".log")
    LOG_HANDLE = open(logfile, "w")
    sys.stdout = TeeOutput(ORIGINAL_STDOUT, LOG_HANDLE)
    sys.stderr = TeeOutput(ORIGINAL_STDERR, LOG_HANDLE)
    return logfile


def stop_logging():
    """Stop logging and restore stdout/stderr."""
    global LOG_HANDLE
    sys.stdout = ORIGINAL_STDOUT
    sys.stderr = ORIGINAL_STDERR
    if LOG_HANDLE:
        LOG_HANDLE.close()
        LOG_HANDLE = None


def revcomp(seq):
    """Return reverse complement of DNA sequence."""
    dna = Seq(seq)
    return str(dna.reverse_complement())


def clean_seq(seq):
    """Translate non-standard bases to standard DNA."""
    if "E" in seq or "F" in seq or "P" in seq or "Q" in seq:
        return seq.translate(TRANSLATION_TABLE)
    return seq


def parse_cutoff(text):
    """Parse length cutoff specification."""
    cutoff = {"type": "none"}
    if not text:
        return cutoff
    if text.startswith("global"):
        cutoff["type"] = "global"
        cutoff["min"] = int(text.split("=")[1].split(",")[0])
        cutoff["max"] = int(text.split(",")[-1])
        return cutoff
    if text.startswith("range"):
        center = int(text.split("=")[1].split(",")[0])
        tolerance = int(text.split(",")[-1])
        cutoff["type"] = "range"
        cutoff["min"] = center - tolerance
        cutoff["max"] = center + tolerance
        return cutoff
    if text.startswith("locus"):
        cutoff["type"] = "locus"
        return cutoff
    return cutoff


def passes_cutoff(length, cutoff):
    """Check if sequence length passes cutoff."""
    if cutoff["type"] == "none":
        return True
    if cutoff["type"] == "locus":
        return True
    return length >= cutoff["min"] and length <= cutoff["max"]


def detect_delimiter(lines, requested):
    """Auto-detect CSV delimiter."""
    if requested != "auto":
        return requested
    first = lines[0]
    if first.count(";") > first.count(","):
        return ";"
    return ","


def header_map(header):
    """Create column index map from header row."""
    mapped = {}
    for i, name in enumerate(header):
        mapped[name.strip().lower()] = i
    return mapped


def is_header_row(row):
    """Check if row looks like a header."""
    joined = " ".join(row).lower()
    if "primer" in joined or "demultiplex" in joined or "sample" in joined:
        return True
    return False


def get_field(row, indices, names, default=""):
    """Get field value by trying multiple column names."""
    for name in names:
        if name in indices and indices[name] < len(row):
            return row[indices[name]].strip()
    return default


def flag_demfile_line(line, delimiter, line_number):
    """Flag demultiplex rows that do not match the expected scheme."""
    fields = line.split(delimiter)
    problems = []
    if len(fields) != 7:
        problems.append("expected 7 fields, found " + str(len(fields)))
    if len(fields) > 0 and any(char.isspace() for char in fields[0]):
        problems.append("field 1 contains spaces")
    for field_number in range(2, min(len(fields), 6) + 1):
        value = fields[field_number - 1].strip().upper()
        if not value or any(base not in "ATCG" for base in value):
            problems.append("field " + str(field_number) + " contains non-ATCG characters")
    if problems:
        print(
            "Demultiplex file warning: line " + str(line_number) + " does not match scheme (" + "; ".join(problems) + ")",
            file=sys.stderr,
        )
    return problems


def read_demfile(path, delimiter):
    """Read and parse demultiplex file."""
    with open(path, "r", encoding="utf-8-sig") as infile:
        line_entries = [(line_number, line.strip()) for line_number, line in enumerate(infile, start=1) if line.strip()]
    if len(line_entries) == 0:
        print("No lines found")
        return []
    lines = [line for _, line in line_entries]
    delim = detect_delimiter(lines, delimiter)
    rows = [(line_number, line, next(csv.reader([line], delimiter=delim))) for line_number, line in line_entries]
    samples = []
    if is_header_row(rows[0][2]):
        indices = header_map(rows[0][2])
        data_rows = rows[1:]
    else:
        indices = None
        data_rows = rows
    for line_number, line, row in data_rows:
        if len(row) == 0:
            continue
        flag_demfile_line(line, delim, line_number)
        if indices is not None:
            name = get_field(row, indices, ["demultiplex_name", "name", "sample"], row[0].strip())
            tag_fwd = get_field(row, indices, ["tag_fwd", "tag1", "tag_forward"])
            tag_rev = get_field(row, indices, ["tag_rev", "tag2", "tag_reverse"])
            primer_fwd = get_field(row, indices, ["primer_fwd", "primer1", "forward_primer"])
            primer_rev = get_field(row, indices, ["primer_rev", "primer2", "reverse_primer"])
            extra = get_field(row, indices, ["tail", "extra", "motif"], "")
        else:
            if len(row) < 5:
                continue
            name = row[0].strip()
            tag_fwd = row[1].strip()
            tag_rev = row[2].strip()
            primer_fwd = row[3].strip()
            primer_rev = row[4].strip()
            extra = ""
            if len(row) > 5:
                extra = row[5].strip()
        if not name or not tag_fwd or not tag_rev or not primer_fwd or not primer_rev:
            continue
        samples.append([name, tag_fwd, tag_rev, primer_fwd, primer_rev, extra])
    return samples


def build_scheme(sample_rows):
    """Build demultiplexing scheme from sample rows."""
    samples = []
    groups = {}
    for row in sample_rows:
        name = row[0]
        tag_fwd = row[1]
        tag_rev = row[2]
        primer_fwd = row[3]
        primer_rev = row[4]
        extra = row[5]
        sample = {
            "name": name,
            "tag_fwd": tag_fwd,
            "tag_rev": tag_rev,
            "tag_rev_rc": revcomp(tag_rev),
            "primer_fwd": primer_fwd,
            "primer_rev": primer_rev,
            "single_left": tag_fwd + primer_fwd,
            "single_right": revcomp(tag_rev + primer_rev + extra),
        }
        samples.append(sample)
        key = primer_fwd + "\t" + primer_rev
        if key not in groups:
            groups[key] = {
                "primer_fwd": primer_fwd,
                "primer_rev": primer_rev,
                "primer_fwd_rc": revcomp(primer_fwd),
                "primer_rev_rc": revcomp(primer_rev),
                "samples": [],
            }
        groups[key]["samples"].append(sample)
    group_list = list(groups.values())
    return samples, group_list


def required_window(groups, mode):
    """Calculate minimum required search window."""
    minimum = 0
    for group in groups:
        minimum = max(minimum, len(group["primer_fwd"]))
        minimum = max(minimum, len(group["primer_rev"]))
        for sample in group["samples"]:
            if mode == "single":
                minimum = max(minimum, len(sample["single_left"]))
                minimum = max(minimum, len(sample["single_right"]))
            else:
                minimum = max(minimum, len(sample["tag_fwd"]))
                minimum = max(minimum, len(sample["tag_rev_rc"]))
    return minimum + 5


def load_records(seqfile, cutoff):
    """Load FASTA/FASTQ records, applying length cutoff."""
    records = []
    total = 0
    for record in iter_sequence_records(seqfile):
        total += 1
        seq = str(record.seq)
        if not passes_cutoff(len(seq), cutoff):
            continue
        seq = clean_seq(seq)
        records.append((record.description, seq))
    return records, total


def batch_records(records, batchsize):
    """Yield batches of records."""
    i = 0
    while i < len(records):
        yield records[i:i+batchsize]
        i += batchsize


def edlib_match(query, target, mismatches):
    """Check if query matches target within mismatch threshold."""
    if len(target) == 0:
        return False
    result = edlib.align(query, target, task="distance", mode="HW", k=mismatches, additionalEqualities=CORRESPONDENCES)
    if result["editDistance"] == -1:
        return False
    if result["editDistance"] <= mismatches:
        return True
    return False


def get_windows(seq, window):
    """Get left and right search windows from sequence."""
    if len(seq) <= window:
        return seq, seq
    return seq[:window], seq[-window:]


def init_worker(data):
    """Initialize worker process with shared data."""
    global WORKER_DATA
    WORKER_DATA = data


def find_primer_group(seq):
    """Find matching primer group for sequence."""
    left, right = get_windows(seq, WORKER_DATA["window"])
    candidates = []
    for group in WORKER_DATA["groups"]:
        if edlib_match(group["primer_fwd"], left, WORKER_DATA["primer_mismatches"]):
            if edlib_match(group["primer_rev_rc"], right, WORKER_DATA["primer_mismatches"]):
                candidates.append((group, "forward"))
                continue
        if edlib_match(group["primer_rev"], left, WORKER_DATA["primer_mismatches"]):
            if edlib_match(group["primer_fwd_rc"], right, WORKER_DATA["primer_mismatches"]):
                candidates.append((group, "reverse"))
    return candidates


def classify_record(record):
    """Classify a single record to a sample."""
    description = record[0]
    seq = record[1]
    primer_hits = find_primer_group(seq)
    if len(primer_hits) == 0:
        return None, None, "no_primer"
    if len(primer_hits) > 1:
        return None, None, "ambiguous_primer"
    group = primer_hits[0][0]
    orientation = primer_hits[0][1]
    if orientation == "reverse":
        seq = revcomp(seq)
    left, right = get_windows(seq, WORKER_DATA["window"])
    sample_hits = []
    if WORKER_DATA["mode"] == "dual":
        for sample in group["samples"]:
            if edlib_match(sample["tag_fwd"], left, WORKER_DATA["tag_mismatches"]):
                if edlib_match(sample["tag_rev_rc"], right, WORKER_DATA["tag_mismatches"]):
                    sample_hits.append(sample["name"])
    else:
        for sample in group["samples"]:
            if edlib_match(sample["single_left"], left, WORKER_DATA["single_mismatches"]):
                if edlib_match(sample["single_right"], right, WORKER_DATA["single_mismatches"]):
                    sample_hits.append(sample["name"])
    if len(sample_hits) == 0:
        return None, None, "no_sample"
    if len(sample_hits) > 1:
        return None, None, "ambiguous_sample"
    return sample_hits[0], ">" + description + "\n" + seq + "\n", "assigned"


def process_batch(records):
    """Process a batch of records."""
    results = {}
    stats = {"processed": 0, "assigned": 0, "no_primer": 0, "ambiguous_primer": 0, "no_sample": 0, "ambiguous_sample": 0}
    for record in records:
        stats["processed"] += 1
        sample_name, fasta_string, status = classify_record(record)
        stats[status] += 1
        if status == "assigned":
            if sample_name not in results:
                results[sample_name] = []
            results[sample_name].append(fasta_string)
    return results, stats


def combine_stats(total_stats, batch_stats):
    """Combine batch statistics into total."""
    for key in batch_stats.keys():
        total_stats[key] += batch_stats[key]


def run_demultiplex(records, samples, groups, args):
    """Run the demultiplexing process."""
    worker_data = {
        "mode": args.mode,
        "window": args.window,
        "groups": groups,
        "primer_mismatches": args.primer_mismatches,
        "tag_mismatches": args.tag_mismatches,
        "single_mismatches": args.single_mismatches,
    }
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    outfiles = {}
    sample_counts = {}
    for sample in samples:
        sample_counts[sample["name"]] = 0
        outfile = os.path.join(args.outdir, sample["name"] + ".fas")
        outfiles[sample["name"]] = outfile
        with open(outfile, "w") as f:
            pass
    total_stats = {"processed": 0, "assigned": 0, "no_primer": 0, "ambiguous_primer": 0, "no_sample": 0, "ambiguous_sample": 0}
    batches = batch_records(records, args.chunksize)
    if args.threads > 1:
        pool = multiprocessing.Pool(args.threads, initializer=init_worker, initargs=(worker_data,))
        try:
            for result_dict, batch_stats in pool.imap_unordered(process_batch, batches, chunksize=1):
                combine_stats(total_stats, batch_stats)
                for sample_name in result_dict.keys():
                    with open(outfiles[sample_name], "a") as f:
                        f.write("".join(result_dict[sample_name]))
                    sample_counts[sample_name] += len(result_dict[sample_name])
                print("\tProcessed reads:", total_stats["processed"], "/", len(records), end="\r")
        finally:
            pool.close()
            pool.join()
    else:
        init_worker(worker_data)
        for batch in batches:
            result_dict, batch_stats = process_batch(batch)
            combine_stats(total_stats, batch_stats)
            for sample_name in result_dict.keys():
                with open(outfiles[sample_name], "a") as f:
                    f.write("".join(result_dict[sample_name]))
                sample_counts[sample_name] += len(result_dict[sample_name])
            print("\tProcessed reads:", total_stats["processed"], "/", len(records), end="\r")
    print()
    return sample_counts, total_stats


def main(args):
    """
    Main entry point for demultiplexing.
    
    Args:
        args: Namespace object with infile, demfile, threads, mode, outdir, cutoff, etc.
    """
    # Set default values for optional args that CLI might not provide
    if not hasattr(args, 'window'):
        args.window = 150
    if not hasattr(args, 'chunksize'):
        args.chunksize = 1000
    if not hasattr(args, 'primer_mismatches'):
        args.primer_mismatches = 5
    if not hasattr(args, 'tag_mismatches'):
        args.tag_mismatches = 0
    if not hasattr(args, 'single_mismatches'):
        args.single_mismatches = 5
    if not hasattr(args, 'delimiter'):
        args.delimiter = "auto"
    
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    
    # Prepare output directory
    if not args.outdir:
        args.outdir = "demultiplex-" + timestamp
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    
    logfile = start_logging(args.outdir, timestamp)
    
    try:
        print("Log file:", logfile)
        print("Command: ontyper demultiplex")
        print("Parameters:")
        for key in sorted(vars(args).keys()):
            print("\t" + key + ":", getattr(args, key))
        print()
        
        print("Output directory:", args.outdir)
        print("Reading demultiplex file:", args.demfile)
        sample_rows = read_demfile(args.demfile, args.delimiter)
        if len(sample_rows) == 0:
            print("No demultiplex information found in", args.demfile, file=sys.stderr)
            return 1
        
        samples, groups = build_scheme(sample_rows)
        print("There are", len(samples), "samples in demfile.")
        print("There are", len(groups), "unique primer pairs.")
        
        minimum_window = required_window(groups, args.mode)
        if args.window < minimum_window:
            print("Requested window", args.window, "is too small for the longest query. Will use", minimum_window, "instead.")
            args.window = minimum_window
        
        cutoff = parse_cutoff(args.cutoff)
        print("Reading in sequences from file:", args.infile)
        print("Detected input format:", detect_sequence_format(args.infile).upper())
        
        if cutoff["type"] == "global":
            print("Using global length cutoff to filter sequences between", cutoff["min"], "and", cutoff["max"], ".")
        elif cutoff["type"] == "range":
            print("Using range as cutoff to filter sequences between", cutoff["min"], "and", cutoff["max"], ".")
        elif cutoff["type"] == "locus":
            print("Per locus length filtering was requested, but no per locus lengths are stored in the demultiplex file. All reads will be kept.")
        else:
            print("Length filtering not specified, will subject all sequences to demultiplexing")
        
        records, total_reads = load_records(args.infile, cutoff)
        print("There are", total_reads, "reads in the input file.")
        print("There are", len(records), "reads for demultiplexing after filtering.")
        print("Using", args.threads, "thread(s), window", args.window, "and batch size", args.chunksize)
        
        if args.mode == "dual":
            print("Starting two pass demultiplexing in one combined per-read workflow...")
            print("Number of allowed primer mismatches:", args.primer_mismatches)
            print("Number of allowed tag mismatches:", args.tag_mismatches)
        else:
            print("Starting single pass demultiplexing with primer pre-grouping...")
            print("Number of allowed primer mismatches for group detection:", args.primer_mismatches)
            print("Number of allowed mismatches for primer+tag matching:", args.single_mismatches)
        
        start = datetime.now()
        sample_counts, total_stats = run_demultiplex(records, samples, groups, args)
        end = datetime.now()
        duration = end - start
        
        print("Total assigned reads:", total_stats["assigned"])
        print("Reads with no matching primer pair:", total_stats["no_primer"])
        print("Reads with ambiguous primer pair:", total_stats["ambiguous_primer"])
        print("Reads with no matching sample inside primer group:", total_stats["no_sample"])
        print("Reads with ambiguous sample inside primer group:", total_stats["ambiguous_sample"])
        for sample in samples:
            print("\tTotal number of demultiplexed sequences for", sample["name"] + ":", sample_counts[sample["name"]])
        print("\tDuration:", duration.total_seconds(), "seconds.")
        print("Done demultiplexing")
        print("Log written to:", logfile)
        
        return 0
    
    finally:
        stop_logging()


if __name__ == "__main__":
    # Allow running this module directly for testing
    import argparse
    parser = argparse.ArgumentParser(description="Demultiplex ONT reads")
    parser.add_argument("-i", "--infile", required=True)
    parser.add_argument("-d", "--demfile", required=True)
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("-m", "--mode", default="single", choices=["single", "dual"])
    parser.add_argument("-o", "--outdir", default="")
    parser.add_argument("-c", "--cutoff", default="")
    parser.add_argument("--window", type=int, default=150)
    parser.add_argument("--chunksize", type=int, default=1000)
    parser.add_argument("--primer-mismatches", type=int, default=5)
    parser.add_argument("--tag-mismatches", type=int, default=0)
    parser.add_argument("--single-mismatches", type=int, default=5)
    parser.add_argument("--delimiter", default="auto")
    args = parser.parse_args()
    sys.exit(main(args))
