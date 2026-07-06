#!/usr/bin/env python3
"""
ontyper.report - Generate a single-file HTML summary report
"""

import html
import json
import os
import sys

from ontyper import FASTA_SUFFIXES
from ontyper.sequence_io import parse_genotype_header


REPORT_TITLE = "ontyper report"


def sample_name_from_path(path):
    """Return sample name from a FASTA file path."""
    name = os.path.basename(path)
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def count_fasta_records(path):
    """Count FASTA records in a file."""
    count = 0
    with open(path) as infile:
        for line in infile:
            if line.startswith(">"):
                count += 1
    return count


def read_demultiplexed_directory(indir):
    """Read per-sample read counts from a demultiplexed FASTA directory."""
    samples = {}
    if not indir:
        return samples
    if not os.path.isdir(indir):
        print("Demultiplexed directory not found:", indir, file=sys.stderr)
        return samples

    for name in sorted(os.listdir(indir)):
        path = os.path.join(indir, name)
        if not os.path.isfile(path):
            continue
        if not name.endswith(FASTA_SUFFIXES):
            continue
        sample_name = sample_name_from_path(path)
        samples[sample_name] = {"demultiplexed_reads": count_fasta_records(path)}
    return samples



def read_genotype_table(path):
    """Read genotype count data from a TSV table."""
    samples = {}
    repeat_bins = []
    if not path:
        return repeat_bins, samples
    if not os.path.isfile(path):
        print("Genotype file not found:", path, file=sys.stderr)
        return repeat_bins, samples

    with open(path, "r", encoding="utf-8-sig") as infile:
        header = infile.readline()
        if not header:
            return repeat_bins, samples
        repeat_bins = parse_genotype_header(header)
        for line in infile:
            if not line.strip():
                continue
            parts = [part.strip() for part in line.rstrip("\n").split("\t")]
            if len(parts) < len(repeat_bins) + 2:
                continue
            sample_name = parts[0]
            counts = []
            for value in parts[1:1 + len(repeat_bins)]:
                counts.append(int(value))
            samples[sample_name] = {
                "repeat_bins": repeat_bins,
                "repeat_counts": counts,
                "nsequences": int(parts[1 + len(repeat_bins)]),
            }
    return repeat_bins, samples


def build_samples(demultiplexed_dir, genotype_file):
    """Merge report data from the provided inputs."""
    merged = {}

    for sample_name, values in read_demultiplexed_directory(demultiplexed_dir).items():
        merged[sample_name] = dict(values)

    _, genotype_samples = read_genotype_table(genotype_file)
    for sample_name, values in genotype_samples.items():
        if sample_name not in merged:
            merged[sample_name] = {}
        merged[sample_name].update(values)

    sample_list = []
    for sample_name in sorted(merged):
        entry = {"name": sample_name}
        entry.update(merged[sample_name])
        if "demultiplexed_reads" not in entry and "nsequences" in entry:
            entry["demultiplexed_reads"] = entry["nsequences"]
        sample_list.append(entry)
    return sample_list


def build_summary(samples):
    """Build top-level report summary."""
    total_demultiplexed = 0
    total_genotyped = 0
    with_demultiplex = 0
    with_genotypes = 0

    for sample in samples:
        if "demultiplexed_reads" in sample:
            total_demultiplexed += sample["demultiplexed_reads"]
            with_demultiplex += 1
        if "nsequences" in sample:
            total_genotyped += sample["nsequences"]
            with_genotypes += 1

    return {
        "nsamples": len(samples),
        "samples_with_demultiplex": with_demultiplex,
        "samples_with_genotypes": with_genotypes,
        "total_demultiplexed_reads": total_demultiplexed,
        "total_genotyped_reads": total_genotyped,
    }


def render_html(samples, args):
    """Render the full single-file HTML report."""
    data = {
        "title": REPORT_TITLE,
        "summary": build_summary(samples),
        "samples": samples,
        "inputs": {
            "demultiplexed_dir": args.demultiplexed,
            "genotypes_file": args.genotypes,
        },
    }
    data_json = json.dumps(data).replace("</", "<\\/")
    title_html = html.escape(REPORT_TITLE)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title_html}</title>
<style>
:root {{
  --bg: #f5f7fb;
  --panel: #ffffff;
  --panel-2: #f8fafc;
  --text: #142033;
  --muted: #5e6b82;
  --line: #dbe2ea;
  --accent: #2563eb;
  --accent-soft: #dbeafe;
  --shadow: 0 10px 30px rgba(15, 23, 42, 0.08);
  --radius: 18px;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  background: linear-gradient(180deg, #f8fbff 0%, var(--bg) 100%);
  color: var(--text);
}}
.header {{
  padding: 28px 32px 20px;
  border-bottom: 1px solid var(--line);
  background: rgba(255,255,255,0.85);
  backdrop-filter: blur(10px);
  position: sticky;
  top: 0;
  z-index: 20;
}}
.header h1 {{ margin: 0 0 8px; font-size: 28px; }}
.header p {{ margin: 0; color: var(--muted); }}
.layout {{
  display: grid;
  grid-template-columns: 320px minmax(0, 1fr);
  gap: 24px;
  padding: 24px;
}}
.sidebar {{
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: var(--radius);
  box-shadow: var(--shadow);
  overflow: hidden;
  align-self: start;
  position: sticky;
  top: 120px;
}}
.sidebar-section {{ padding: 20px; border-bottom: 1px solid var(--line); }}
.sidebar h2, .content h2 {{ margin: 0 0 14px; font-size: 17px; }}
.search {{
  width: 100%;
  padding: 11px 13px;
  border: 1px solid var(--line);
  border-radius: 12px;
  background: var(--panel-2);
  font: inherit;
}}
.sample-list {{ max-height: calc(100vh - 280px); overflow: auto; padding: 10px; }}
.sample-button {{
  width: 100%;
  border: 0;
  background: transparent;
  text-align: left;
  padding: 12px 14px;
  border-radius: 14px;
  cursor: pointer;
  margin-bottom: 6px;
  transition: 0.15s ease;
}}
.sample-button:hover {{ background: var(--panel-2); }}
.sample-button.active {{ background: var(--accent-soft); color: var(--accent); }}
.sample-name {{ display: block; font-weight: 600; margin-bottom: 4px; }}
.sample-meta {{ color: var(--muted); font-size: 12px; }}
.content {{ display: grid; gap: 20px; }}
.cards {{ display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 16px; }}
.card {{
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: var(--radius);
  padding: 18px 20px;
  box-shadow: var(--shadow);
}}
.card-label {{ color: var(--muted); font-size: 13px; margin-bottom: 8px; }}
.card-value {{ font-size: 28px; font-weight: 700; }}
.panel {{
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: var(--radius);
  padding: 22px;
  box-shadow: var(--shadow);
}}
.panel-header {{ display: flex; align-items: baseline; justify-content: space-between; gap: 16px; margin-bottom: 16px; }}
.panel-subtitle {{ color: var(--muted); font-size: 14px; }}
.summary-grid {{ display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 14px; }}
.mini-stat {{ background: var(--panel-2); border-radius: 14px; padding: 16px; }}
.mini-stat strong {{ display: block; font-size: 22px; margin-bottom: 6px; }}
.kv {{ display: grid; grid-template-columns: 220px 1fr; gap: 10px 16px; }}
.kv-key {{ color: var(--muted); }}
.table-wrap {{ overflow: auto; max-height: 360px; border: 1px solid var(--line); border-radius: 14px; }}
table {{ width: 100%; border-collapse: collapse; font-size: 14px; }}
th, td {{ padding: 10px 12px; border-bottom: 1px solid var(--line); text-align: right; }}
th:first-child, td:first-child {{ text-align: left; }}
thead th {{ position: sticky; top: 0; background: var(--panel); }}
.empty {{ color: var(--muted); padding: 22px; text-align: center; border: 1px dashed var(--line); border-radius: 14px; }}
#chart {{ width: 100%; min-height: 340px; }}
svg text {{ fill: #334155; font-size: 11px; font-family: inherit; }}
.bar {{ fill: #4c78a8; }}
.bar:hover {{ fill: #2563eb; }}
.axis, .grid {{ stroke: #cbd5e1; stroke-width: 1; }}
.grid {{ stroke-dasharray: 3 4; }}
.footer {{ color: var(--muted); font-size: 13px; text-align: right; }}
@media (max-width: 1100px) {{
  .layout {{ grid-template-columns: 1fr; }}
  .sidebar {{ position: static; }}
  .cards, .summary-grid {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
}}
@media (max-width: 700px) {{
  .header {{ padding: 22px 18px 18px; }}
  .layout {{ padding: 18px; }}
  .cards, .summary-grid {{ grid-template-columns: 1fr; }}
  .kv {{ grid-template-columns: 1fr; }}
}}
</style>
</head>
<body>
<div class="header">
  <h1>{title_html}</h1>
  <p>Interactive summary of demultiplexing and genotype identification results.</p>
</div>
<div class="layout">
  <aside class="sidebar">
    <div class="sidebar-section">
      <h2>Samples</h2>
      <input id="search" class="search" type="search" placeholder="Search sample name...">
    </div>
    <div id="sampleList" class="sample-list"></div>
  </aside>
  <main class="content">
    <section class="panel">
      <div class="panel-header">
        <div>
          <h2>Report overview</h2>
          <div class="panel-subtitle">Input files and overall counts</div>
        </div>
      </div>
      <div class="summary-grid">
        <div class="mini-stat"><strong id="summarySamples"></strong><span>Samples in report</span></div>
        <div class="mini-stat"><strong id="summaryDemux"></strong><span>Samples with demultiplexed reads</span></div>
        <div class="mini-stat"><strong id="summaryGenotypes"></strong><span>Samples with genotype table entries</span></div>
        <div class="mini-stat"><strong id="summaryReads"></strong><span>Total genotyped reads</span></div>
      </div>
      <div style="height: 18px"></div>
      <div class="kv">
        <div class="kv-key">Demultiplexed directory</div><div id="inputDemultiplexed"></div>
        <div class="kv-key">Genotype table</div><div id="inputGenotypes"></div>
      </div>
    </section>

    <section class="cards">
      <div class="card"><div class="card-label">Selected sample</div><div id="cardSample" class="card-value">-</div></div>
      <div class="card"><div class="card-label">Demultiplexed reads</div><div id="cardDemultiplexed" class="card-value">-</div></div>
      <div class="card"><div class="card-label">Genotyped reads</div><div id="cardGenotyped" class="card-value">-</div></div>
      <div class="card"><div class="card-label">Peak repeat bin</div><div id="cardPeak" class="card-value">-</div></div>
    </section>

    <section class="panel">
      <div class="panel-header">
        <div>
          <h2>Sample details</h2>
          <div id="sampleSubtitle" class="panel-subtitle">Select a sample to inspect its repeat-count distribution.</div>
        </div>
      </div>
      <div id="sampleInfo" class="kv"></div>
    </section>

    <section class="panel">
      <div class="panel-header">
        <div>
          <h2>Repeat count histogram</h2>
          <div class="panel-subtitle">Distribution of reads across repeat motif counts</div>
        </div>
      </div>
      <div id="chart"></div>
    </section>

    <section class="panel">
      <div class="panel-header">
        <div>
          <h2>Repeat count table</h2>
          <div class="panel-subtitle">Counts by repeat bin for the selected sample</div>
        </div>
      </div>
      <div id="tableContainer"></div>
    </section>

    <div class="footer">Generated by ontyper report</div>
  </main>
</div>

<script>
const REPORT_DATA = {data_json};
const SAMPLE_LIST = document.getElementById('sampleList');
const SEARCH = document.getElementById('search');
const SAMPLE_INFO = document.getElementById('sampleInfo');
const TABLE_CONTAINER = document.getElementById('tableContainer');
const CHART = document.getElementById('chart');
let currentSample = null;

function formatNumber(value) {{
  if (value === null || value === undefined) return '-';
  return new Intl.NumberFormat().format(value);
}}

function fullDistribution(sample) {{
  if (!sample.repeat_bins || !sample.repeat_counts) return {{ bins: [], counts: [] }};
  return {{
    bins: sample.repeat_bins.slice(),
    counts: sample.repeat_counts.slice(),
  }};
}}

function nonZeroDistribution(sample) {{
  const distribution = fullDistribution(sample);
  let first = -1;
  let last = -1;
  for (let i = 0; i < distribution.counts.length; i += 1) {{
    if (distribution.counts[i] > 0) {{
      if (first === -1) first = i;
      last = i;
    }}
  }}
  if (first === -1) return distribution;
  return {{
    bins: distribution.bins.slice(first, last + 1),
    counts: distribution.counts.slice(first, last + 1),
  }};
}}

function peakBin(sample) {{
  if (!sample.repeat_bins || !sample.repeat_counts) return '-';
  let bestIndex = 0;
  for (let i = 1; i < sample.repeat_counts.length; i += 1) {{
    if (sample.repeat_counts[i] > sample.repeat_counts[bestIndex]) bestIndex = i;
  }}
  return sample.repeat_counts[bestIndex] > 0 ? sample.repeat_bins[bestIndex] : '-';
}}

function renderSummary() {{
  document.getElementById('summarySamples').textContent = formatNumber(REPORT_DATA.summary.nsamples);
  document.getElementById('summaryDemux').textContent = formatNumber(REPORT_DATA.summary.samples_with_demultiplex);
  document.getElementById('summaryGenotypes').textContent = formatNumber(REPORT_DATA.summary.samples_with_genotypes);
  document.getElementById('summaryReads').textContent = formatNumber(REPORT_DATA.summary.total_genotyped_reads);
  document.getElementById('inputDemultiplexed').textContent = REPORT_DATA.inputs.demultiplexed_dir || '-';
  document.getElementById('inputGenotypes').textContent = REPORT_DATA.inputs.genotypes_file || '-';
}}

function sampleMeta(sample) {{
  const bits = [];
  if (sample.demultiplexed_reads !== undefined) bits.push(`demultiplexed: ${{formatNumber(sample.demultiplexed_reads)}}`);
  if (sample.nsequences !== undefined) bits.push(`genotyped: ${{formatNumber(sample.nsequences)}}`);
  return bits.join(' · ');
}}

function renderSampleList(filterText = '') {{
  const query = filterText.trim().toLowerCase();
  SAMPLE_LIST.innerHTML = '';
  const filtered = REPORT_DATA.samples.filter(sample => sample.name.toLowerCase().includes(query));

  for (const sample of filtered) {{
    const button = document.createElement('button');
    button.className = 'sample-button' + (currentSample && currentSample.name === sample.name ? ' active' : '');
    button.innerHTML = `<span class="sample-name">${{sample.name}}</span><span class="sample-meta">${{sampleMeta(sample)}}</span>`;
    button.addEventListener('click', () => selectSample(sample.name));
    SAMPLE_LIST.appendChild(button);
  }}

  if (filtered.length === 0) {{
    SAMPLE_LIST.innerHTML = '<div class="empty">No samples match the current search.</div>';
  }}
}}

function renderSampleInfo(sample) {{
  document.getElementById('cardSample').textContent = sample.name;
  document.getElementById('cardDemultiplexed').textContent = formatNumber(sample.demultiplexed_reads);
  document.getElementById('cardGenotyped').textContent = formatNumber(sample.nsequences);
  document.getElementById('cardPeak').textContent = peakBin(sample);
  document.getElementById('sampleSubtitle').textContent = sampleMeta(sample) || 'Sample selected';

  const info = [
    ['Sample name', sample.name],
    ['Demultiplexed reads', formatNumber(sample.demultiplexed_reads)],
    ['Genotyped reads', formatNumber(sample.nsequences)],
    ['Peak repeat bin', peakBin(sample)],
  ];
  SAMPLE_INFO.innerHTML = info.map(([key, value]) => `<div class="kv-key">${{key}}</div><div>${{value}}</div>`).join('');
}}

function renderTable(sample) {{
  if (!sample.repeat_bins || !sample.repeat_counts) {{
    TABLE_CONTAINER.innerHTML = '<div class="empty">No genotype count table available for this sample.</div>';
    return;
  }}
  const distribution = nonZeroDistribution(sample);
  const rows = distribution.bins.map((bin, i) => `<tr><td>${{bin}}</td><td>${{formatNumber(distribution.counts[i])}}</td></tr>`).join('');
  TABLE_CONTAINER.innerHTML = `
    <div class="table-wrap">
      <table>
        <thead><tr><th>Repeat bin</th><th>Reads</th></tr></thead>
        <tbody>${{rows}}</tbody>
      </table>
    </div>`;
}}

function renderChart(sample) {{
  if (!sample.repeat_bins || !sample.repeat_counts) {{
    CHART.innerHTML = '<div class="empty">No genotype count table available for this sample.</div>';
    return;
  }}

  const distribution = fullDistribution(sample);
  const bins = distribution.bins;
  const counts = distribution.counts;
  const width = CHART.clientWidth || 800;
  const height = 340;
  const margin = {{ top: 18, right: 18, bottom: 58, left: 56 }};
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;
  const maxCount = Math.max(...counts, 1);
  const barWidth = innerWidth / Math.max(bins.length, 1);
  const yTicks = 5;
  const step = bins.length <= 20 ? 1 : (bins.length <= 40 ? 2 : 5);

  let svg = `<svg viewBox="0 0 ${{width}} ${{height}}" width="100%" height="${{height}}" aria-label="Repeat count histogram">`;

  for (let i = 0; i <= yTicks; i += 1) {{
    const y = margin.top + (innerHeight * i / yTicks);
    const value = Math.round(maxCount - (maxCount * i / yTicks));
    svg += `<line class="grid" x1="${{margin.left}}" y1="${{y}}" x2="${{width - margin.right}}" y2="${{y}}"></line>`;
    svg += `<text x="${{margin.left - 8}}" y="${{y + 4}}" text-anchor="end">${{value}}</text>`;
  }}

  svg += `<line class="axis" x1="${{margin.left}}" y1="${{margin.top + innerHeight}}" x2="${{width - margin.right}}" y2="${{margin.top + innerHeight}}"></line>`;
  svg += `<line class="axis" x1="${{margin.left}}" y1="${{margin.top}}" x2="${{margin.left}}" y2="${{margin.top + innerHeight}}"></line>`;

  for (let i = 0; i < bins.length; i += 1) {{
    const x = margin.left + i * barWidth + 1;
    const barHeight = innerHeight * counts[i] / maxCount;
    const y = margin.top + innerHeight - barHeight;
    svg += `<rect class="bar" x="${{x}}" y="${{y}}" width="${{Math.max(barWidth - 2, 1)}}" height="${{barHeight}}"><title>${{bins[i]}}: ${{counts[i]}} reads</title></rect>`;
    if (i % step === 0 || i === bins.length - 1) {{
      const tx = margin.left + i * barWidth + barWidth / 2;
      svg += `<text x="${{tx}}" y="${{height - 18}}" text-anchor="middle">${{bins[i]}}</text>`;
    }}
  }}

  svg += `<text x="${{width / 2}}" y="${{height - 2}}" text-anchor="middle">Repeat motif count</text>`;
  svg += `<text transform="translate(16, ${{height / 2}}) rotate(-90)" text-anchor="middle">Number of sequences</text>`;
  svg += `</svg>`;
  CHART.innerHTML = svg;
}}

function selectSample(name) {{
  currentSample = REPORT_DATA.samples.find(sample => sample.name === name) || null;
  if (!currentSample) return;
  renderSampleList(SEARCH.value);
  renderSampleInfo(currentSample);
  renderTable(currentSample);
  renderChart(currentSample);
}}

SEARCH.addEventListener('input', () => renderSampleList(SEARCH.value));
window.addEventListener('resize', () => {{ if (currentSample) renderChart(currentSample); }});

renderSummary();
renderSampleList();
if (REPORT_DATA.samples.length > 0) {{
  selectSample(REPORT_DATA.samples[0].name);
}}
</script>
</body>
</html>
"""


def main(args):
    """Main entry point for report generation."""
    if not args.demultiplexed and not args.genotypes:
        print("Please provide at least one input with -d/--demultiplexed or -g/--genotypes", file=sys.stderr)
        return 1

    samples = build_samples(args.demultiplexed, args.genotypes)
    if not samples:
        print("No report data found in the provided inputs", file=sys.stderr)
        return 1

    outfile = args.outfile or "ontyper-report.html"
    html = render_html(samples, args)
    with open(outfile, "w", encoding="utf-8") as out:
        out.write(html)
    print("Wrote report to", outfile, file=sys.stderr)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate a single-file HTML report")
    parser.add_argument("-d", "--demultiplexed", default=None, help="Directory with demultiplexed FASTA files")
    parser.add_argument("-g", "--genotypes", default=None, help="genotypes.tsv file")
    parser.add_argument("-o", "--outfile", default="ontyper-report.html", help="Output HTML file")
    sys.exit(main(parser.parse_args()))
