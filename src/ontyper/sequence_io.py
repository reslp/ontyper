"""Helpers for reading plain or gzipped FASTA/FASTQ files."""

import gzip

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
