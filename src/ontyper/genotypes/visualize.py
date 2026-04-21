#!/usr/bin/env python3
"""
ontyper.genotypes.visualize - Visualize genotype count distributions

This module reads a genotype count table and creates a publication-ready figure
for one locus/sample.
"""

import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def parse_header(line):
    """Return repeat-bin labels from the header row."""
    parts = [part.strip() for part in line.rstrip("\n").split("\t")]
    bins = []
    for part in parts[1:]:
        if part == "nsequences":
            break
        if part.endswith("x"):
            bins.append(part)
    return bins


def parse_row(line, bins):
    """Parse one genotype table row."""
    parts = [part.strip() for part in line.rstrip("\n").split("\t")]
    if len(parts) < len(bins) + 2:
        return None, None
    locus = parts[0]
    counts = []
    for value in parts[1:1 + len(bins)]:
        counts.append(int(value))
    return locus, counts


def read_counts(genotype_file, locus_name):
    """Read repeat-count distribution for one locus/sample."""
    with open(genotype_file, "r", encoding="utf-8-sig") as infile:
        header = infile.readline()
        bins = parse_header(header)
        for line in infile:
            if not line.strip():
                continue
            locus, counts = parse_row(line, bins)
            if locus == locus_name:
                return bins, counts
    return None, None


def trim_trailing_zeros(bins, counts):
    """Trim empty bins at the end of the distribution."""
    last_nonzero = 0
    for i, count in enumerate(counts):
        if count > 0:
            last_nonzero = i
    return bins[:last_nonzero + 1], counts[:last_nonzero + 1]


def tick_step(nbins):
    """Choose a readable x-axis tick interval."""
    if nbins <= 20:
        return 1
    if nbins <= 40:
        return 2
    return 5


def default_outfile(locus_name):
    """Build a default output file name."""
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", locus_name)
    return safe_name + "_genotypes.pdf"


def plot_distribution(locus_name, bins, counts, outfile):
    """Create a genotype distribution figure."""
    bins, counts = trim_trailing_zeros(bins, counts)
    xvalues = list(range(len(bins)))

    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
    })

    fig_width = max(6.0, min(10.0, 0.18 * len(bins) + 4.5))
    fig, ax = plt.subplots(figsize=(fig_width, 3.8))

    ax.bar(xvalues, counts, width=0.82, color="#4C78A8", edgecolor="black", linewidth=0.6)
    ax.set_title(locus_name)
    ax.set_xlabel("Repeat motif count")
    ax.set_ylabel("Number of sequences")
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    step = tick_step(len(bins))
    tick_positions = list(range(0, len(bins), step))
    if len(bins) - 1 not in tick_positions:
        tick_positions.append(len(bins) - 1)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([bins[i] for i in tick_positions])
    ax.set_xlim(-0.6, len(bins) - 0.4)

    total = sum(counts)
    ax.text(
        0.98,
        0.98,
        f"n = {total:,}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 2.5},
    )

    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)


def main(args):
    """Main entry point for genotype visualization."""
    bins, counts = read_counts(args.genotypes, args.locus)
    if bins is None:
        print("Locus not found in genotype table:", args.locus, file=sys.stderr)
        return 1

    outfile = args.outfile
    if not outfile:
        outfile = default_outfile(args.locus)

    outdir = os.path.dirname(outfile)
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir)

    plot_distribution(args.locus, bins, counts, outfile)
    print("Wrote genotype figure to", outfile, file=sys.stderr)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Visualize genotype count distributions")
    parser.add_argument("-g", "--genotypes", required=True, help="Genotype count table (TSV)")
    parser.add_argument("-l", "--locus", required=True, help="Locus/sample name to plot")
    parser.add_argument("-o", "--outfile", default=None, help="Output figure path (default: <locus>_genotypes.pdf)")
    sys.exit(main(parser.parse_args()))
