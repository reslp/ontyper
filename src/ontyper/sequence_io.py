"""Helpers for reading plain or gzipped FASTA/FASTQ files."""

import gzip
import os

from Bio import SeqIO

from ontyper import FASTA_SUFFIXES, FASTQ_SUFFIXES


def _strip_gzip_suffix(path):
    """Remove a trailing .gz suffix if present."""
    text = str(path)
    if text.lower().endswith(".gz"):
        return text[:-3]
    return text


def open_sequence_file(path):
    """Open a plain-text or gzipped sequence file for reading."""
    text = str(path)
    if text.lower().endswith(".gz"):
        return gzip.open(text, "rt", encoding="utf-8")
    return open(text, "r", encoding="utf-8")


def detect_sequence_format(path):
    """Detect whether a sequence file is FASTA or FASTQ."""
    name = _strip_gzip_suffix(path).lower()
    if name.endswith(FASTA_SUFFIXES):
        return "fasta"
    if name.endswith(FASTQ_SUFFIXES):
        return "fastq"

    with open_sequence_file(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                return "fasta"
            if line.startswith("@"):
                return "fastq"
            break

    raise ValueError(f"Could not determine sequence format for {path}")


def iter_sequence_records(path):
    """Yield records from a FASTA/FASTQ file, gzipped or plain."""
    file_format = detect_sequence_format(path)
    with open_sequence_file(path) as handle:
        yield from SeqIO.parse(handle, file_format)


def strip_fasta_suffix(path):
    """Remove a known FASTA suffix from a file name or path."""
    name = os.path.basename(str(path))
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name


def parse_genotype_header(line):
    """Return repeat-bin labels from a genotype TSV header row."""
    parts = [part.strip() for part in line.rstrip("\n").split("\t")]
    repeat_bins = []
    for part in parts[1:]:
        if part == "nsequences":
            break
        if part.endswith("x"):
            repeat_bins.append(part)
    return repeat_bins
