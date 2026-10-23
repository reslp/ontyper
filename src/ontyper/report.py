#!/usr/bin/env python3
"""
ontyper.report - Generate a single-file HTML summary report
"""

import base64
import csv
import html
import json
import os
import sys
from pathlib import Path

from ontyper import FASTA_SUFFIXES
from ontyper.genotypes.features import read_true_alleles
from ontyper.sequence_io import parse_genotype_header, strip_fasta_suffix


REPORT_TITLE = "ontyper report"



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
        sample_name = strip_fasta_suffix(path)
        samples[sample_name] = {"demultiplexed_reads": count_fasta_records(path)}
    return samples



def detect_table_delimiter(header):
    """Detect a simple delimiter from a table header line."""
    counts = {delimiter: header.count(delimiter) for delimiter in ["\t", ",", ";"]}
    delimiter, count = max(counts.items(), key=lambda item: item[1])
    return delimiter if count > 0 else "\t"


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


def prediction_sample_name(row):
    """Return the report sample name for one prediction row."""
    name = (row.get("name") or "").strip()
    if name:
        return name
    sample_id = (row.get("sample_id") or "").strip()
    locus_id = (row.get("locus_id") or "").strip()
    if sample_id and locus_id:
        return sample_id + "_" + locus_id
    return sample_id or locus_id


def parse_prediction_row(row):
    """Normalize one prediction row for embedding in the report JSON."""
    prediction = {
        "sample_id": (row.get("sample_id") or "").strip(),
        "locus_id": (row.get("locus_id") or "").strip(),
        "allele_id": (row.get("allele_id") or "").strip(),
    }
    probability = (row.get("probability") or row.get("pred") or "").strip()
    if probability:
        prediction["probability"] = float(probability)
    rank = (row.get("pred_rank_within_locus") or row.get("rank") or "").strip()
    if rank:
        prediction["rank"] = int(float(rank))
    return prediction


def read_prediction_table(path):
    """Read model prediction results grouped by report sample name."""
    samples = {}
    if not path:
        return samples
    if not os.path.isfile(path):
        print("Prediction file not found:", path, file=sys.stderr)
        return samples

    with open(path, "r", encoding="utf-8-sig", newline="") as infile:
        header = infile.readline()
        if not header:
            return samples
        delimiter = detect_table_delimiter(header)
        infile.seek(0)
        reader = csv.DictReader(infile, delimiter=delimiter)
        for row in reader:
            sample_name = prediction_sample_name(row)
            if not sample_name:
                continue
            if sample_name not in samples:
                samples[sample_name] = {"predictions": []}
            samples[sample_name]["predictions"].append(parse_prediction_row(row))

    for values in samples.values():
        values["predictions"].sort(key=lambda item: (item.get("locus_id", ""), item.get("rank", 999999), -item.get("probability", -1)))
        values["npredictions"] = len(values["predictions"])
    return samples


def read_true_allele_table(path):
    """Read the training truth format and group distinct alleles by sample/locus."""
    if not path:
        return {}
    truth = read_true_alleles(path)
    samples = {}
    for (sample_id, locus_id), group in truth.groupby(["sample_id", "locus_id"]):
        sample_id, locus_id = str(sample_id), str(locus_id)
        name = prediction_sample_name({"sample_id": sample_id, "locus_id": locus_id})
        samples[name] = {
            "sample_id": sample_id,
            "locus_id": locus_id,
            "true_alleles": sorted(group["true_allele"].astype(int).unique().tolist()),
        }
    return samples


def compare_alleles(sample, predictions_supplied):
    """Compare unique calls with known alleles; absent truth is not a negative label."""
    calls = {}
    for prediction in sample.get("predictions", []):
        allele = int(prediction["allele_id"])
        if allele not in calls or prediction.get("probability", -1) > calls[allele].get("probability", -1):
            calls[allele] = prediction
    truth = set(sample.get("true_alleles", []))
    has_truth = "true_alleles" in sample
    evaluated = has_truth and predictions_supplied
    comparison = []
    for allele in sorted(set(calls) | truth):
        called = allele in calls
        known_true = allele in truth if has_truth else None
        status = "called" if called else "known"
        if evaluated:
            status = "correct" if called and known_true else "incorrect" if called else "missed"
        row = dict(calls.get(allele, {}))
        row.update({
            "allele_id": allele,
            "locus_id": row.get("locus_id", sample.get("locus_id", "")),
            "called": called,
            "known_true": known_true,
            "status": status,
        })
        comparison.append(row)
    overlap = None
    if evaluated:
        true_positives = len(set(calls) & truth)
        overlap = {
            "correct": true_positives,
            "incorrect": len(set(calls) - truth),
            "missed": len(truth - set(calls)),
            "precision": true_positives / len(calls) if calls else None,
            "recall": true_positives / len(truth) if truth else None,
        }
    return comparison, overlap


def build_samples(demultiplexed_dir, genotype_file, prediction_file=None, true_allele_file=None):
    """Merge report data from the provided inputs."""
    merged = {}

    for sample_name, values in read_demultiplexed_directory(demultiplexed_dir).items():
        merged[sample_name] = dict(values)

    _, genotype_samples = read_genotype_table(genotype_file)
    for sample_name, values in genotype_samples.items():
        if sample_name not in merged:
            merged[sample_name] = {}
        merged[sample_name].update(values)

    prediction_samples = read_prediction_table(prediction_file)
    for sample_name, values in prediction_samples.items():
        if sample_name not in merged:
            merged[sample_name] = {}
        merged[sample_name].update(values)

    for sample_name, values in read_true_allele_table(true_allele_file).items():
        merged.setdefault(sample_name, {}).update(values)

    sample_list = []
    for sample_name in sorted(merged):
        entry = {"name": sample_name}
        entry.update(merged[sample_name])
        entry["allele_comparison"], entry["allele_overlap"] = compare_alleles(entry, bool(prediction_file))
        if "demultiplexed_reads" not in entry and "nsequences" in entry:
            entry["demultiplexed_reads"] = entry["nsequences"]
        sample_list.append(entry)
    return sample_list


def build_truth_statistics(samples):
    """Summarize call/truth overlap across evaluated sample/locus entries only."""
    evaluated = [sample["allele_overlap"] for sample in samples if sample.get("allele_overlap") is not None]
    total = len(evaluated)
    counts = {
        "exact_matches": sum(values["correct"] > 0 and values["incorrect"] == 0 and values["missed"] == 0 for values in evaluated),
        "with_missed_true_alleles": sum(values["missed"] > 0 for values in evaluated),
        "with_no_true_allele_called": sum(values["correct"] == 0 for values in evaluated),
    }
    return {
        "evaluated": total,
        **{key: {"count": count, "proportion": count / total if total else None} for key, count in counts.items()},
    }


def build_summary(samples):
    """Build top-level report summary."""
    total_demultiplexed = 0
    total_genotyped = 0
    with_demultiplex = 0
    with_genotypes = 0
    with_predictions = 0
    total_predictions = 0

    for sample in samples:
        if "demultiplexed_reads" in sample:
            total_demultiplexed += sample["demultiplexed_reads"]
            with_demultiplex += 1
        if "nsequences" in sample:
            total_genotyped += sample["nsequences"]
            with_genotypes += 1
        if "predictions" in sample:
            with_predictions += 1
            total_predictions += len(sample["predictions"])

    return {
        "nsamples": len(samples),
        "samples_with_demultiplex": with_demultiplex,
        "samples_with_genotypes": with_genotypes,
        "samples_with_predictions": with_predictions,
        "total_demultiplexed_reads": total_demultiplexed,
        "total_genotyped_reads": total_genotyped,
        "total_predictions": total_predictions,
        "truth_statistics": build_truth_statistics(samples),
    }


def report_logo_data_uri():
    """Embed the packaged logo, with a source-tree fallback for editable installs."""
    logo = Path(__file__).parent / "assets" / "ontyper-icon-ont.png"
    if not logo.is_file():
        logo = Path(__file__).resolve().parents[2] / "assets" / "ontyper-icon-ont.png"
    return "data:image/png;base64," + base64.b64encode(logo.read_bytes()).decode("ascii")


def render_html(samples, args):
    """Render the full single-file HTML report."""
    data = {
        "title": REPORT_TITLE,
        "summary": build_summary(samples),
        "samples": samples,
        "inputs": {
            "demultiplexed_dir": getattr(args, "demultiplexed", None),
            "genotypes_file": getattr(args, "genotypes", None),
            "predictions_file": getattr(args, "predictions", None),
            "true_alleles_file": getattr(args, "true_alleles", None),
        },
    }
    data_json = json.dumps(data).replace("</", "<\\/")
    title_html = html.escape(REPORT_TITLE)
    logo_uri = report_logo_data_uri()

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
  display: flex;
  align-items: center;
  gap: 14px;
  padding: 16px 24px;
  background: var(--panel);
}}
.report-logo {{ width: 48px; height: 48px; flex-shrink: 0; }}
.header h1 {{ margin: 0 0 3px; font-size: 24px; }}
.header p {{ margin: 0; color: var(--muted); font-size: 12px; }}
.report-overview {{ padding: 0 24px 12px; background: var(--panel); border-bottom: 1px solid var(--line); }}
.report-overview h2 {{ margin: 0 0 8px; color: var(--muted); font-size: 12px; font-weight: 600; }}
.overview-details {{ margin-top: 10px; font-size: 12px; }}
.overview-details summary, .sample-details summary {{ cursor: pointer; color: var(--muted); }}
.overview-details .kv {{ margin-top: 10px; gap: 5px 16px; overflow-wrap: anywhere; }}
.sample-details .kv {{ margin-top: 14px; }}
#histogramPanel {{ scroll-margin-top: 16px; }}
#histogramPanel h2 {{ margin: 3px 0 5px; font-size: 20px; overflow-wrap: anywhere; }}
.eyebrow {{ font-size: 12px; color: var(--muted); }}
.sample-metrics {{ display: flex; flex-wrap: wrap; gap: 6px 20px; margin-bottom: 14px; color: var(--muted); font-size: 12px; }}
.sample-metrics strong {{ color: var(--text); font-size: 14px; }}
.layout {{
  display: grid;
  grid-template-columns: 280px minmax(0, 1fr);
  gap: 18px;
  padding: 18px 24px;
}}
.sidebar {{
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: var(--radius);
  box-shadow: var(--shadow);
  overflow: hidden;
  align-self: start;
  position: sticky;
  top: 16px;
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
.sample-list {{ max-height: calc(100vh - 180px); overflow: auto; padding: 10px; }}
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
.content {{ display: grid; align-content: start; gap: 16px; min-width: 0; }}
.panel {{
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: var(--radius);
  padding: 22px;
  box-shadow: var(--shadow);
}}
.panel-header {{ display: flex; align-items: baseline; justify-content: space-between; gap: 16px; margin-bottom: 16px; }}
.panel-subtitle {{ color: var(--muted); font-size: 14px; }}
.summary-grid {{ display: flex; flex-wrap: wrap; gap: 8px 24px; }}
.mini-stat {{ display: flex; align-items: baseline; gap: 6px; font-size: 12px; }}
.mini-stat strong {{ font-size: 16px; }}
.truth-grid .mini-stat {{ flex: 1; min-width: 180px; flex-direction: column; gap: 2px; }}
#truthStatistics {{ margin-top: 12px; }}
#truthStatisticsNote {{ color: var(--muted); font-size: 12px; }}
.stat-green strong {{ color: #16a34a; }}
.stat-red strong {{ color: #dc2626; }}
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
.bar.called {{ fill: #f97316; }}
.bar.called:hover {{ fill: #ea580c; }}
.bar.correct {{ fill: #16a34a; }}
.bar.correct:hover {{ fill: #15803d; }}
.bar.missed, .bar.known {{ fill: #dc2626; }}
.bar.missed:hover, .bar.known:hover {{ fill: #b91c1c; }}
.allele-legend {{ display: flex; flex-wrap: wrap; gap: 10px 20px; margin-bottom: 14px; font-size: 13px; }}
.allele-legend[hidden] {{ display: none; }}
.legend-item {{ display: inline-flex; align-items: center; gap: 7px; }}
.swatch {{ display: inline-block; width: 12px; height: 12px; border-radius: 3px; background: #4c78a8; }}
.swatch.called, .swatch.incorrect {{ background: #f97316; }}
.swatch.correct {{ background: #16a34a; }}
.swatch.missed, .swatch.known {{ background: #dc2626; }}
#alleleOverlap {{ margin: 0 0 14px; color: var(--muted); font-size: 14px; }}
.axis, .grid {{ stroke: #cbd5e1; stroke-width: 1; }}
.grid {{ stroke-dasharray: 3 4; }}
.footer {{ color: var(--muted); font-size: 13px; text-align: right; }}
@media (max-width: 1000px) {{
  .layout {{ grid-template-columns: 220px minmax(0, 1fr); }}
}}
@media (max-width: 700px) {{
  .header {{ padding: 12px; }}
  .report-overview {{ padding: 0 12px 12px; }}
  .layout {{ grid-template-columns: 1fr; padding: 12px; }}
  .sidebar {{ position: static; }}
  .sample-list {{ max-height: 200px; }}
  .panel {{ padding: 16px; }}
  .kv {{ grid-template-columns: 1fr; }}
}}
</style>
</head>
<body>
<header class="header">
  <img class="report-logo" src="{logo_uri}" alt="ontyper logo" width="48" height="48">
  <div>
    <h1>{title_html}</h1>
    <p>Interactive summary of demultiplexing and genotype identification results. Allele predictions are included when provided.</p>
  </div>
</header>
<section class="report-overview" aria-labelledby="overviewTitle">
  <h2 id="overviewTitle">Report overview</h2>
  <div class="summary-grid">
    <div class="mini-stat"><strong id="summarySamples"></strong><span>Samples in report</span></div>
    <div class="mini-stat"><strong id="summaryDemux"></strong><span>Samples with demultiplexed reads</span></div>
    <div class="mini-stat"><strong id="summaryGenotypes"></strong><span>Samples with genotype table entries</span></div>
    <div class="mini-stat"><strong id="summaryReads"></strong><span>Total genotyped reads</span></div>
    <div class="mini-stat"><strong id="summaryPredictions"></strong><span>Predicted allele rows</span></div>
  </div>
  <div id="truthStatistics" hidden>
    <h2>Allele call accuracy</h2>
    <div class="summary-grid truth-grid">
      <div class="mini-stat stat-green"><strong id="summaryExactMatches"></strong><span>Exact matches: all true alleles called, no extra calls (all green)</span></div>
      <div class="mini-stat stat-red"><strong id="summaryMissedTrue"></strong><span>Entries with uncalled true alleles</span></div>
      <div class="mini-stat stat-red"><strong id="summaryNoTrueCalled"></strong><span>Entries with no known true allele called</span></div>
    </div>
  </div>
  <details class="overview-details">
    <summary>Input files and evaluation details</summary>
    <div class="kv">
      <div class="kv-key">Demultiplexed directory</div><div id="inputDemultiplexed"></div>
      <div class="kv-key">Genotype table</div><div id="inputGenotypes"></div>
      <div class="kv-key">Prediction table</div><div id="inputPredictions"></div>
      <div class="kv-key">True allele table</div><div id="inputTrueAlleles"></div>
    </div>
    <p id="truthStatisticsNote" hidden></p>
  </details>
</section>
<div class="layout">
  <aside class="sidebar">
    <div class="sidebar-section">
      <h2>Samples</h2>
      <input id="search" class="search" type="search" placeholder="Search sample name...">
    </div>
    <div id="sampleList" class="sample-list"></div>
  </aside>
  <main class="content">
    <section id="histogramPanel" class="panel">
      <div class="panel-header">
        <div>
          <div class="eyebrow">Repeat count histogram</div>
          <h2 id="cardSample">Selected sample</h2>
          <div id="sampleSubtitle" class="panel-subtitle">Select a sample to inspect its repeat-count distribution.</div>
        </div>
      </div>
      <div class="sample-metrics">
        <span>Demultiplexed reads <strong id="cardDemultiplexed">-</strong></span>
        <span>Genotyped reads <strong id="cardGenotyped">-</strong></span>
        <span>Peak repeat bin <strong id="cardPeak">-</strong></span>
        <span>Predicted alleles <strong id="cardPredictions">-</strong></span>
      </div>
      <div id="alleleLegend" class="allele-legend" aria-label="Allele color legend" hidden></div>
      <div id="chart"></div>
    </section>

    <section class="panel">
      <div class="panel-header">
        <div>
          <h2>Predicted alleles</h2>
          <div class="panel-subtitle">Model-selected calls and known true alleles, when supplied</div>
        </div>
      </div>
      <p id="alleleOverlap" hidden></p>
      <div id="predictionContainer"></div>
    </section>

    <details class="panel sample-details">
      <summary>Sample details</summary>
      <div id="sampleInfo" class="kv"></div>
    </details>

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
const PREDICTION_CONTAINER = document.getElementById('predictionContainer');
const CHART = document.getElementById('chart');
let currentSample = null;

function formatNumber(value) {{
  if (value === null || value === undefined) return '-';
  return new Intl.NumberFormat().format(value);
}}

function fullDistribution(sample) {{
  if (!sample.repeat_bins || !sample.repeat_counts) return {{ bins: [], counts: [] }};
  const counts = new Map(sample.repeat_bins.map((bin, i) => [Number(bin.replace(/x$/, '')), sample.repeat_counts[i]]));
  for (const allele of sample.allele_comparison || []) {{
    if (!counts.has(allele.allele_id)) counts.set(allele.allele_id, 0);
  }}
  const alleles = Array.from(counts.keys()).sort((a, b) => a - b);
  return {{ bins: alleles.map(allele => `${{allele}}x`), counts: alleles.map(allele => counts.get(allele)) }};
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
  document.getElementById('summaryPredictions').textContent = formatNumber(REPORT_DATA.summary.total_predictions);
  document.getElementById('inputDemultiplexed').textContent = REPORT_DATA.inputs.demultiplexed_dir || '-';
  document.getElementById('inputGenotypes').textContent = REPORT_DATA.inputs.genotypes_file || '-';
  document.getElementById('inputPredictions').textContent = REPORT_DATA.inputs.predictions_file || '-';
  document.getElementById('inputTrueAlleles').textContent = REPORT_DATA.inputs.true_alleles_file || '-';
  document.getElementById('truthStatistics').hidden = !REPORT_DATA.inputs.true_alleles_file;
  document.getElementById('truthStatisticsNote').hidden = !REPORT_DATA.inputs.true_alleles_file;
  const stats = REPORT_DATA.summary.truth_statistics;
  for (const [id, key] of [['summaryExactMatches', 'exact_matches'], ['summaryMissedTrue', 'with_missed_true_alleles'], ['summaryNoTrueCalled', 'with_no_true_allele_called']]) {{
    const value = stats[key];
    document.getElementById(id).textContent = value.proportion === null
      ? 'Not evaluated'
      : `${{(100 * value.proportion).toFixed(1)}}% (${{value.count}}/${{stats.evaluated}})`;
  }}
  document.getElementById('truthStatisticsNote').textContent = stats.evaluated
    ? `Based on ${{stats.evaluated}} evaluated sample/locus entries; entries without truth are excluded. Entries with no true allele called also count as having uncalled true alleles, so percentages overlap.`
    : 'No evaluated sample/locus entries. Both prediction results and matching truth data are needed.';
}}

function sampleMeta(sample) {{
  const bits = [];
  if (sample.demultiplexed_reads !== undefined) bits.push(`demultiplexed: ${{formatNumber(sample.demultiplexed_reads)}}`);
  if (sample.nsequences !== undefined) bits.push(`genotyped: ${{formatNumber(sample.nsequences)}}`);
  if (sample.predictions !== undefined) bits.push(`predictions: ${{formatNumber(sample.predictions.length)}}`);
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
    button.addEventListener('click', () => selectSample(sample.name, true));
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
  document.getElementById('cardPredictions').textContent = sample.predictions ? formatNumber(sample.predictions.length) : '-';
  document.getElementById('sampleSubtitle').textContent = sampleMeta(sample) || 'Sample selected';

  const info = [
    ['Sample name', sample.name],
    ['Demultiplexed reads', formatNumber(sample.demultiplexed_reads)],
    ['Genotyped reads', formatNumber(sample.nsequences)],
    ['Peak repeat bin', peakBin(sample)],
    ['Predicted alleles', sample.predictions ? formatNumber(sample.predictions.length) : '-'],
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

function formatProbability(value) {{
  if (value === null || value === undefined) return '-';
  return Number(value).toFixed(3);
}}

function escapeHtml(value) {{
  const element = document.createElement('span');
  element.textContent = String(value);
  return element.innerHTML;
}}

function renderAlleleLegend(sample) {{
  const legend = document.getElementById('alleleLegend');
  legend.hidden = !REPORT_DATA.inputs.predictions_file && !REPORT_DATA.inputs.true_alleles_file;
  const items = [['other', 'Other bins (blue)']];
  if (sample.allele_overlap) {{
    items.push(['correct', 'Called and known true (green)'],
               ['incorrect', 'Called but not known true (orange)'],
               ['missed', 'Known true but not called (red)']);
  }} else {{
    if (REPORT_DATA.inputs.predictions_file) items.push(['called', 'Called; truth unavailable (orange)']);
    if (sample.true_alleles !== undefined) items.push(['known', 'Known true; calls not supplied (red)']);
  }}
  legend.innerHTML = items.map(([status, label]) => `<span class="legend-item"><span class="swatch ${{status}}"></span>${{label}}</span>`).join('');
}}

function renderPredictions(sample) {{
  const overlap = document.getElementById('alleleOverlap');
  overlap.hidden = !REPORT_DATA.inputs.true_alleles_file;
  if (sample.allele_overlap) {{
    const values = sample.allele_overlap;
    const percent = value => value === null ? '-' : `${{(100 * value).toFixed(1)}}%`;
    overlap.textContent = `Correct: ${{values.correct}} · Incorrect: ${{values.incorrect}} · Missed: ${{values.missed}} · Precision: ${{percent(values.precision)}} · Recall: ${{percent(values.recall)}}`;
  }} else {{
    overlap.textContent = sample.true_alleles === undefined
      ? 'Not evaluated: no truth data for this sample/locus.'
      : 'Not evaluated: no prediction file supplied.';
  }}
  const alleles = sample.allele_comparison || [];
  if (alleles.length === 0) {{
    PREDICTION_CONTAINER.innerHTML = '<div class="empty">No allele calls or known true alleles available for this sample.</div>';
    return;
  }}
  const showTruth = Boolean(REPORT_DATA.inputs.true_alleles_file);
  const rows = alleles.map(allele => `
    <tr>
      <td>${{escapeHtml(allele.locus_id || '-')}}</td>
      <td><span class="swatch ${{allele.status}}"></span> ${{allele.allele_id}}</td>
      ${{showTruth ? `<td>${{allele.called ? 'Yes' : 'No'}}</td><td>${{allele.known_true === null ? 'Not evaluated' : allele.known_true ? 'Yes' : 'No'}}</td>` : ''}}
      <td>${{formatProbability(allele.probability)}}</td>
      <td>${{allele.rank ?? '-'}}</td>
    </tr>`).join('');
  PREDICTION_CONTAINER.innerHTML = `
    <div class="table-wrap">
      <table>
        <thead><tr><th>Locus</th><th>Allele</th>${{showTruth ? '<th>Called</th><th>Known true</th>' : ''}}<th>Probability</th><th>Rank</th></tr></thead>
        <tbody>${{rows}}</tbody>
      </table>
    </div>`;
}}

function calledAlleleMap(sample) {{
  const called = new Map();
  if (!sample.predictions) return called;
  for (const prediction of sample.predictions) {{
    const allele = String(prediction.allele_id ?? '').replace(/x$/, '');
    if (!allele) continue;
    const probability = prediction.probability;
    if (!called.has(allele) || probability > called.get(allele)) {{
      called.set(allele, probability);
    }}
  }}
  return called;
}}

function renderChart(sample) {{
  if (!sample.repeat_bins || !sample.repeat_counts) {{
    CHART.innerHTML = '<div class="empty">No genotype count table available for this sample.</div>';
    return;
  }}

  const distribution = fullDistribution(sample);
  const bins = distribution.bins;
  const counts = distribution.counts;
  const called = calledAlleleMap(sample);
  const comparison = new Map((sample.allele_comparison || []).map(allele => [String(allele.allele_id), allele]));
  const statusLabels = {{ correct: 'called and known true', incorrect: 'called but not known true', missed: 'known true but not called', known: 'known true; calls not supplied', called: 'called; truth not evaluated' }};
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
    const allele = String(bins[i]).replace(/x$/, '');
    const entry = comparison.get(allele);
    // Keep a small marker visible for known/called alleles with zero supporting reads.
    const barHeight = Math.max(innerHeight * counts[i] / maxCount, entry ? 3 : 0);
    const y = margin.top + innerHeight - barHeight;
    const isCalled = called.has(allele);
    const probability = called.get(allele);
    const className = (isCalled ? 'bar called' : 'bar') + (entry && entry.status !== 'called' ? ` ${{entry.status}}` : '');
    let title = isCalled
      ? `${{bins[i]}}: ${{counts[i]}} reads · called allele probability ${{formatProbability(probability)}}`
      : `${{bins[i]}}: ${{counts[i]}} reads`;
    if (entry) title += ` · ${{statusLabels[entry.status]}}`;
    svg += `<rect class="${{className}}" x="${{x}}" y="${{y}}" width="${{Math.max(barWidth - 2, 1)}}" height="${{barHeight}}"><title>${{escapeHtml(title)}}</title></rect>`;
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

function selectSample(name, revealChart = false) {{
  currentSample = REPORT_DATA.samples.find(sample => sample.name === name) || null;
  if (!currentSample) return;
  renderSampleList(SEARCH.value);
  renderSampleInfo(currentSample);
  renderTable(currentSample);
  renderPredictions(currentSample);
  renderAlleleLegend(currentSample);
  renderChart(currentSample);
  if (revealChart) document.getElementById('histogramPanel').scrollIntoView({{ block: 'start', behavior: 'auto' }});
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
    predictions = getattr(args, "predictions", None)
    if not args.demultiplexed and not args.genotypes and not predictions:
        print("Please provide at least one input with -d/--demultiplexed, -g/--genotypes, or -p/--predictions", file=sys.stderr)
        return 1

    try:
        samples = build_samples(args.demultiplexed, args.genotypes, predictions, getattr(args, "true_alleles", None))
    except (ValueError, OSError) as error:
        print(error, file=sys.stderr)
        return 1
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
    parser.add_argument("-p", "--predictions", default=None, help="Selected allele calls table")
    parser.add_argument("--true-alleles", default=None, help="Known allele table in the same format used for training")
    parser.add_argument("-o", "--outfile", default="ontyper-report.html", help="Output HTML file")
    sys.exit(main(parser.parse_args()))
